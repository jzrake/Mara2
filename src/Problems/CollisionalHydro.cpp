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
    void advance (MeshData& solution, ParticleData& particleData, double dt) const override;
    double getCollisionTime (MeshData& solution, ParticleData& particleData) const;
private:
    double crossSection;
};




// ============================================================================
CollisionalHydroScheme::CollisionalHydroScheme()
{
    fluxScheme = std::make_shared<MethodOfLines>();
    crossSection = 500.0; // per unit mass
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
    const double particleMass = 5.0 / particleData.particles.size();
    std::vector<Index> indexes;
    std::vector<double> weights;

    int numCollisions = 0;
    const auto& geometry = meshOperator->getMeshGeometry();

    for (auto& p : particleData.particles)
    {
        auto x = Coordinate {{ p.position, 0.0, 0.0 }};
        auto xind = geometry.indexAtCoordinate(x)[0] - startIndex[0];
        double df = solution.P (xind, 0, 0, 0);
        double vf = solution.P (xind, 0, 0, 1);
        double dl = 1.0 / (crossSection * df);

        p.position += dt * p.velocity;
        p.opticalDepth += std::fabs (p.velocity - vf) / (dl / dt);


        // Collisions
        // --------------------------------------------------------------------
        if (p.opticalDepth > 1.0)
        {
            numCollisions += 1;
            meshOperator->weight (x, indexes, weights);

            const double v0 = p.velocity;
            const double u0 = v0 - vf; // initial velocity in fluid frame
            const double u1 = -u0;     // final velocity in fluid frame
            const double v1 = u1 + vf; // final velocity in system frame
            const double dv = v1 - v0; // change in particle velocity (both frames)

            p.velocity = v1;
            p.opticalDepth -= 1.0;

            for (unsigned int i = 0; i < indexes.size(); ++i)
            {
                auto I = indexes[i];
                I[0] -= startIndex[0];
                I[1] -= startIndex[1];
                I[2] -= startIndex[2];

                L(I[0], I[1], I[2], 1) -= dv * particleMass * weights[i];
            }
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

double CollisionalHydroScheme::getCollisionTime (MeshData& solution, ParticleData& particleData) const
{
    const auto& geometry = meshOperator->getMeshGeometry();
    double collisionTime = 1.0;
    auto footprint = Shape3D();
    auto startIndex = Index();
    auto interior = Region();

    SchemeHelpers::makeFootprint (
        fluxScheme->getStencilSize(),
        solution.P.shape(),
        solution.getBoundaryShape(),
        footprint, startIndex, interior);

    for (auto& p : particleData.particles)
    {
        auto x = Coordinate {{ p.position, 0.0, 0.0 }};
        auto xind = geometry.indexAtCoordinate(x)[0] - startIndex[0];
        double df = solution.P (xind, 0, 0, 0);
        double vf = solution.P (xind, 0, 0, 1);
        double dl = 1.0 / (crossSection * df);
        double dt = dl / std::fabs (p.velocity - vf);

        collisionTime = std::min (collisionTime, dt);
    }
    return collisionTime;
}




// ============================================================================
int CollisionalHydroProgram::run (int argc, const char* argv[])
{
    auto cs = Shape {{128, 1, 1 }}; // cells shape
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
    auto cflParameter = 0.25;
    auto finalTime = 12.;
    auto L = mo->linearCellDimension();
    auto timestepSize = 0.0;

    auto timestep  = [&] () { return timestepSize; };
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
    }, TaskScheduler::Recurrence (0.1));

    scheduler->schedule ([&] (SimulationStatus, int rep)
    {
        auto L = mo->linearCellDimension();
        auto P = md->getPrimitive();
        timestepSize = ss->getCollisionTime (*md, *pd);
        timestepSize = std::min (timestepSize, fo->getCourantTimestep (P, L));
        timestepSize *= cflParameter;
    }, TaskScheduler::Recurrence (0.0, 0.0, 10));

    status = SimulationStatus();
    status.totalCellsInMesh = mg->totalCellsInMesh();

    int N = 10000;

    for (int i = 0; i < cs[0]; ++i)
    {
        double x = mg->coordinateAtIndex (i, 0, 0)[0];

        for (int n = 0; n < N; ++n)
        {
            auto p = ParticleData::Particle();
            double pressure = 1.0 + 0.1 * std::sin (4 * M_PI * x);
            p.position = x;
            p.velocity = std::sqrt (pressure) * (-1 + 2 * double (rand()) / RAND_MAX);
            p.opticalDepth = 1.0;
            pd->particles.push_back(p);
        }
    }

    auto initialData = [&] (double x, double, double) -> std::vector<double>
    {
        return std::vector<double> {1, 0, 0, 0, 0.01};
        // return std::vector<double> {1, 0, 0, 0, 0.10 + 0.001 * std::sin (4 * M_PI * x)};
    };

    md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    md->applyBoundaryCondition (*bc);

    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);
}
