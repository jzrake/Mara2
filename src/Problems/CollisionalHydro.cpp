#include <cmath>
#include "../Mara.hpp"
#include "../BoundaryConditions.hpp"
#include "../CartesianMeshGeometry.hpp"
#include "../CellCenteredFieldCT.hpp"
#include "../Checkpoint.hpp"
#include "../ConservationLaws.hpp"
#include "../FieldOperator.hpp"
#include "../IntercellFluxSchemes.hpp"
#include "../MeshData.hpp"
#include "../MeshOperator.hpp"
#include "../ParticleData.hpp"
#include "../SolutionSchemes.hpp"
#include "../TaskScheduler.hpp"
#include "../TimeSeriesManager.hpp"
#include "CollisionalHydro.hpp"
#include "MPI.hpp"
#include "HDF5.hpp"
#include "Timer.hpp"
#include "Logger.hpp"

using namespace Cow;




// ============================================================================
class CollisionalHydroScheme : public GenericSolutionScheme
{
public:
    CollisionalHydroScheme();
    void advance (MeshData& solution, double dt) const override;
    void advance (MeshData& solution, ParticleData& particles, double dt) const override;
};




// ============================================================================
CollisionalHydroScheme::CollisionalHydroScheme()
{
    fluxScheme = std::make_shared<MethodOfLines>();
}

void CollisionalHydroScheme::advance (MeshData& solution, double dt) const
{
    throw std::logic_error ("CollisionalHydroScheme does not implement advance() without particle data");
}

void CollisionalHydroScheme::advance (MeshData& solution, ParticleData& particleData, double dt) const
{
    if (! fieldOperator)     throw std::logic_error ("No FieldOperator instance");
    if (! meshOperator)      throw std::logic_error ("No MeshOperator instance");
    if (! boundaryCondition) throw std::logic_error ("No BoundaryCondition instance");
    if (! fluxScheme)        throw std::logic_error ("No IntercellFluxScheme instance");


    // Figure out the scheme footprint, and if we have enough guard zones
    // ------------------------------------------------------------------------
    auto footprint = Shape3D();
    auto startIndex = Index();
    auto interior = Region();

    SchemeHelpers::makeFootprint (
        fluxScheme->getStencilSize(),
        solution.P.shape(),
        solution.getBoundaryShape(),
        footprint, startIndex, interior);


    // Setup callback to compute Godunov fluxes
    // ------------------------------------------------------------------------
    auto cl = fieldOperator->getConservationLaw();
    auto nq = cl->getNumConserved();

    auto D = IntercellFluxScheme::FaceData();
    D.conservationLaw = cl;

    auto Fhat = [&] (GodunovStencil& stencil)
    {
        D.areaElement = stencil.faceNormal.cartesian();
        D.stencilData = stencil.cellData;

        auto S = fluxScheme->intercellFlux (D);

        for (int q = 0; q < nq; ++q)
        {
            stencil.godunovFlux[q] = S.F[q];
        }
    };


    // Update
    // ------------------------------------------------------------------------
    auto U0 = fieldOperator->generateConserved (solution.P);
    auto U = U0;

    auto F = meshOperator->godunov (Fhat, solution.P, solution.B, footprint, startIndex);
    auto L = meshOperator->divergence (F, -1.0);


    // Update particles
    // ------------------------------------------------------------------------
    const double meanFreePath = 0.02;
    const double particleMass = 1.0 / particleData.particles.size();
    std::vector<Index> indexes;
    std::vector<double> weights;

    int numCollisions = 0;
    const auto& geometry = meshOperator->getMeshGeometry();

    for (auto& p : particleData.particles)
    {
        auto x = Coordinate {{ p.position, 0.0, 0.0 }};
        auto xind = geometry.indexAtCoordinate(x)[0] - startIndex[0];
        double v = solution.P (xind, 0, 0, 1);

        p.position += dt * p.velocity;
        p.opticalDepth += dt * std::fabs (p.velocity - v) / meanFreePath;


        // Collisions
        // --------------------------------------------------------------------
        if (p.opticalDepth > 1)
        {
            numCollisions += 1;
            meshOperator->weight (x, indexes, weights);

            for (unsigned int i = 0; i < indexes.size(); ++i)
            {
                auto I = indexes[i];
                I[0] -= startIndex[0];
                I[1] -= startIndex[1];
                I[2] -= startIndex[2];

                L(I[0], I[1], I[2], 1) += 2 * p.velocity * particleMass * weights[i];
            }

            p.velocity *= -1.0;
            p.opticalDepth -= 1.0;
        }


        // Boundary conditions
        // --------------------------------------------------------------------
        if (p.position < 0)
        {
            p.position += 1;
        }
        else if (p.position >= 1)
        {
            p.position -= 1;
        }
    }

    for (int n = 0; n < L.size(); ++n)
    {
        U[n] += dt * L[n];
    }

    fieldOperator->recoverPrimitive (U[interior], solution.P[interior]);
    solution.applyBoundaryCondition (*boundaryCondition);
}




// ============================================================================
int CollisionalHydroProgram::run (int argc, const char* argv[])
{
    auto cs = Shape {{256, 1, 1 }}; // cells shape
    auto bs = Shape {{  2, 0, 0 }};
    auto ss = std::make_shared<CollisionalHydroScheme>();
    auto bc = std::make_shared<PeriodicBoundaryCondition>();
    auto mg = std::make_shared<CartesianMeshGeometry>();
    auto mo = std::make_shared<MeshOperator>();
    auto cl = std::make_shared<NewtonianHydro>();
    auto fo = std::make_shared<FieldOperator>();
    auto md = std::make_shared<MeshData> (cs, bs, 5);
    auto pd = std::make_shared<ParticleData>();

    mg->setCellsShape (cs);
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setFieldOperator (fo);
    ss->setMeshOperator (mo);
    ss->setBoundaryCondition (bc);

    auto status = SimulationStatus();
    auto cflParameter = 0.5;
    auto finalTime = 0.04;
    auto L = mo->linearCellDimension();
    auto P = md->getPrimitive();

    auto timestep  = [&] () { return cflParameter * fo->getCourantTimestep (P, L); };
    auto condition = [&] () { return status.simulationTime < finalTime; };
    auto advance   = [&] (double dt) { return ss->advance (*md, *pd, dt); };
    auto scheduler = std::make_shared<TaskScheduler>();
    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();

    writer->setFilenamePrefix ("chkpt");
    writer->setMeshDecomposition (nullptr);
    writer->setTimeSeriesManager (nullptr);

    scheduler->schedule ([&] (SimulationStatus, int rep)
    {
        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);
    }, TaskScheduler::Recurrence (finalTime));

    status = SimulationStatus();
    status.totalCellsInMesh = mg->totalCellsInMesh();


    int N = 1000;

    for (int i = 0; i < cs[0]; ++i)
    {
        double x = mg->coordinateAtIndex (i, 0, 0)[0];

        for (int n = 0; n < N; ++n)
        {
            auto p = ParticleData::Particle();
            double pressure = 1.0 + 0.0 * std::sin (4 * M_PI * x);
            p.position = x;
            p.velocity = 1.0 + std::sqrt (pressure) * (-1 + 2 * double (rand()) / RAND_MAX);
            p.opticalDepth = 0;
            pd->particles.push_back(p);
        }
    }

    auto initialData = [&] (double x, double, double) -> std::vector<double>
    {
        return std::vector<double> {1, 0, 0, 0, 1};
    };

    md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    md->applyBoundaryCondition (*bc);

    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);
}
