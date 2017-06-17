#include <iomanip>
#include "Problems.hpp"
#include "Mara.hpp"
#include "TaskScheduler.hpp"

using namespace Cow;

#include "Timer.hpp"
#include "Logger.cpp"




// ============================================================================
static void maraMainLoop (
    SimulationStatus& status,
    std::function<double ()> timestep,
    std::function<bool ()> condition,
    std::function<void (double)> advance,
    TaskScheduler& scheduler,
    Logger& logger)
{
    auto simulationTimer = Timer();

    while (condition())
    {
        scheduler.dispatch (status);

        auto stepTimer = Timer();
        double dt = timestep();
        advance (dt);

        status.simulationTime += dt;
        status.simulationIter += 1;
        status.wallMinutes = simulationTimer.minutes();

        double kzps = 1e-3 * status.totalCellsInMesh / stepTimer.age();
        logger.log() << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        logger.log() << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        logger.log() << "dt=" << std::setprecision (4) << std::scientific << dt << " ";
        logger.log() << "kzps=" << std::setprecision (2) << std::fixed << kzps << std::endl;
    }
}




// ============================================================================
#include "BoundaryConditions.hpp"
#include "CartesianMeshGeometry.hpp"
#include "Checkpoint.hpp"
#include "ConservationLaws.hpp"
#include "MeshData.hpp"
#include "MeshOperator.hpp"
#include "FieldOperator.hpp"
#include "SolutionSchemes.hpp"
#include "TimeSeriesManager.hpp"

int SimpleTestProgram::run (int argc, const char* argv[]) const
{
    auto cs = Shape {{128, 1, 1 }}; // cells shape
    auto bs = Shape {{  1, 0, 0 }}; // boundary shape
    auto mg = std::make_shared<CartesianMeshGeometry>(cs);
    auto mo = std::make_shared<MeshOperator>(mg);
    auto cl = std::make_shared<NewtonianHydro>();
    auto fo = std::make_shared<FieldOperator>(cl);
    auto ss = std::make_shared<MethodOfLinesTVD>();
    auto bc = std::make_shared<PeriodicBoundaryCondition>();
    auto md = std::make_shared<MeshData>(cs, bs, 5);
    auto id = [] (double x, double, double) { return std::vector<double> {1., 0., 0., 0., 1.}; };

    ss->setFieldOperator (fo);
    ss->setMeshOperator (mo);
    ss->setBoundaryCondition (bc);

    md->assignPrimitive (mo->generate (id, MeshLocation::cell));
    ss->applyBoundaryCondition (*md);

    auto status = SimulationStatus();
    auto finalTime = 0.5;
    auto cflParameter = 0.5;

    auto L = mo->linearCellDimension();
    auto P = md->getPrimitive();

    auto timestep  = [&] () { return cflParameter * fo->getCourantTimestep (P, L); };
    auto condition = [&] () { return status.simulationTime < finalTime; };
    auto advance   = [&] (double dt) { return ss->advance (dt, *md); };
    auto scheduler = std::make_shared<TaskScheduler>();
    auto logger    = std::make_shared<Logger>();

    status.totalCellsInMesh = mg->totalCellsInMesh();

    scheduler->schedule ([&] (SimulationStatus status, int rep)
    {
    	// place-holder
    }, TaskScheduler::Recurrence (0.2));

    maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);

    auto tseries = std::make_shared<TimeSeriesManager>();
    auto writer = CheckpointWriter();

    writer.writeCheckpoint (0, status, *cl, *md, *mg, *tseries, *logger);

    return 0;
}
