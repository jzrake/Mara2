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
    virtual ~CollisionalHydroScheme() {}
    void advance (MeshData& solution, double dt) const override;
    void advance (MeshData& solution, ParticleData& pd, double dt) const override;
    double getCollisionTime (const MeshData& solution, const ParticleData& pd) const;
private:

    /**
    Move particles forward in time by the given time step dt. This updates
    their position and optical depth. Also, update the time remaining in their
    scattering events, and remove scattering events whose durations are
    expired.
    */
    void advanceParticles (ParticleData& pd, double dt) const;

    /**
    Add scattering events to any particles whose optical depth now exceeds 1,
    and subtract 1 from their optical depth.
    */
    void doCollisions (ParticleData& pd) const;

    /**
    Sample the velocity and density field in the given MeshData instance at
    particle positions, and record those values in the particle data
    structure.
    */
    void sampleGridToParticles (MeshData& md, ParticleData& pd) const;

    /**
    Generate an array of source terms (energy and momentum per unit volume)
    from all scattering events in the particle data. Source terms are already
    multiplied by time. The array size does not include guard zones.
    */
    Array getSourceTermsFromParticles (const ParticleData& pd, double dt) const;

    double scatteringEventDuration;
    double crossSection;
    double particleMass;
};




// ============================================================================
void CollisionalHydroScheme::advanceParticles (ParticleData& pd, double dt) const
{
    const auto& geometry = meshOperator->getMeshGeometry();

    for (auto& p : pd.particles)
    {
        double df = p.fluidDensity;
        double vf = p.fluidVelocity;
        double dl = 1.0 / (crossSection * df);
        p.position     += dt * p.velocity;
        p.opticalDepth += dt / dl * std::fabs (p.velocity - vf);


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

        p.meshIndex = geometry.indexAtCoordinate (Coordinate ({{p.position, 0.0, 0.0}}));


        // Update scattering events
        // --------------------------------------------------------------------
        auto it = p.events.begin();

        while (it != p.events.end())
        {
            it->remaining -= dt;

            if (it->remaining < 0.0)
            {
                it = p.events.erase (it);
            }
            else
            {
                ++it;
            }
        }
    }
}

void CollisionalHydroScheme::doCollisions (ParticleData& pd) const
{
    for (auto& p : pd.particles)
    {
        if (p.opticalDepth > 1.0)
        {
            const double vf = p.fluidVelocity;
            const double v0 = p.velocity;
            const double u0 = v0 - vf; // initial velocity in fluid frame
            const double u1 = -u0;     // final velocity in fluid frame
            const double v1 = u1 + vf; // final velocity in system frame
            const double dv = v1 - v0; // change in particle velocity (both frames)
            const double dp = dv * particleMass; // impulse given to particle
            const double dW = vf * dp; // work done on particle

            auto event = ParticleData::ScatteringEvent();
            event.energy = dW * p.weight;
            event.momentum = dp * p.weight;
            event.duration = scatteringEventDuration;
            event.remaining = event.duration;

            p.opticalDepth -= 1.0;
            p.events.push_back (event);
        }
    }
}

void CollisionalHydroScheme::sampleGridToParticles (MeshData& md, ParticleData& pd) const
{
    for (auto& p : pd.particles)
    {
        auto i = p.meshIndex;
        double df = md.P (i[0], 0, 0, 0);
        double vf = md.P (i[0], 0, 0, 1);

        p.fluidDensity = df;
        p.fluidVelocity = vf;
    }
}

Array CollisionalHydroScheme::getSourceTermsFromParticles (const ParticleData& pd, double dt) const
{
    const auto& geometry = meshOperator->getMeshGeometry();
    const auto cellsShape = Shape3D (geometry.cellsShape());

    auto S = Array (cellsShape.withComponents(5));

    int numEvents = 0;

    for (const auto& p : pd.particles)
    {
        const auto i = p.meshIndex;
        const auto V = geometry.cellVolume (i[0], i[1], i[2]);

        for (const auto& event : p.events)
        {
            // Scattering events indicate the momentum taken by the particle.

            const double ds = std::min (dt, event.remaining) / event.duration;

            S (i[0], i[1], i[2], 1) -= ds / V * event.momentum;
            S (i[0], i[1], i[2], 4) -= ds / V * event.energy;

            ++numEvents;
        }
    }

    // std::cout << "There are " << double (numEvents) / pd.particles.size() << " scattering events per particle\n";

    return S;
}

double CollisionalHydroScheme::getCollisionTime (const MeshData& md, const ParticleData& pd) const
{
    double collisionTime = 1.0;

    for (const auto& p : pd.particles)
    {
        double df = p.fluidDensity;
        double vf = p.fluidVelocity;
        double dl = 1.0 / (crossSection * df);
        double dt = dl / std::fabs (p.velocity - vf);

        collisionTime = std::min (collisionTime, dt);
    }
    return collisionTime;
}




// ============================================================================
CollisionalHydroScheme::CollisionalHydroScheme()
{
    fluxScheme = std::make_shared<MethodOfLines>();
    crossSection = 40.0; // per unit mass
    particleMass = 1e-5;
    scatteringEventDuration = 1e-1; // time over which scattering events are communicated to the grid
}

void CollisionalHydroScheme::advance (MeshData& solution, double dt) const
{
    throw std::logic_error ("CollisionalHydroScheme does not implement advance() without particle data");
}

void CollisionalHydroScheme::advance (MeshData& md, ParticleData& pd, double dt) const
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
        md.P.shape(),
        md.getBoundaryShape(),
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
    auto U0 = fieldOperator->generateConserved (md.P);
    auto U = U0;

    auto F = meshOperator->godunov (Fhat, md.P, md.B, footprint, startIndex, nullptr);
    auto L = meshOperator->divergence (F, -1.0);
    auto S = Array (L.shape());

    sampleGridToParticles (md, pd);
    doCollisions (pd);
    S[interior] = getSourceTermsFromParticles (pd, dt);
    advanceParticles (pd, dt);

    for (int n = 0; n < L.size(); ++n)
    {
        U[n] += dt * L[n] + S[n];
    }

    fieldOperator->recoverPrimitive (U[interior], md.P[interior]);
    md.applyBoundaryCondition (*boundaryCondition);
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
    auto finalTime = 16.0;
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
    }, TaskScheduler::Recurrence (0.05));

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

    int N = 1000;

    for (int i = 0; i < cs[0]; ++i)
    {
        double x = mg->coordinateAtIndex (i, 0, 0)[0];
        double dx = mg->cellLength (i, 0, 0, 0);

        for (int n = 0; n < N; ++n)
        {
            auto p = ParticleData::Particle();
            double pressure = 1.0 + 0.1 * std::sin (4 * M_PI * x);
            double b = -0.5 + 1.0 * n / (N - 1);
            p.position = x - b * dx;
            p.velocity = std::sqrt (pressure) * (-1 + 2 * double (rand()) / RAND_MAX);
            p.opticalDepth = double (rand()) / RAND_MAX;
            p.weight = 1.0 / N;
            pd->particles.push_back(p);
        }
    }

    auto initialData = [&] (double x, double, double) -> std::vector<double>
    {
        return std::vector<double> {1, 0, 0, 0, 0.01};
    };

    md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    md->applyBoundaryCondition (*bc);

    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);
}
