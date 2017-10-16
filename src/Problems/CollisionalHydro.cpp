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

    cl->addSourceTerms (solution.P, L);


    const double meanFreePath = 1.0;

    for (auto& p : particleData.particles)
    {
        p.position += dt * p.velocity;
        p.opticalDepth += dt * p.velocity / meanFreePath;

        if (p.opticalDepth > 1)
        {
            int meshIndexOfParticle = 0; // TODO
            double cellVolume = 1.0; // TODO

            L[meshIndexOfParticle] += 2 * p.momentum / cellVolume;

            p.velocity *= -1.0;
            p.momentum *= -1.0;
            p.opticalDepth -= 1.0;
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
    auto finalTime = 1.0;
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

    auto initialData = [&] (double x, double, double) -> std::vector<double>
    {
        return std::vector<double> {1, 0, 0, 0, 1};
    };
    md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    md->applyBoundaryCondition (*bc);

    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);
}
