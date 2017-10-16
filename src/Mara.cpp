#include <iomanip>

// Mara includes
#include "Mara.hpp"
#include "Problems.hpp"
#include "BackwardsCompat.hpp"
#include "Checkpoint.hpp"
#include "TaskScheduler.hpp"
#include "Timer.hpp"

// Cow includes
#include "MPI.hpp"
#include "DebugHelper.hpp"

// Problems
#include "Problems/CollisionalHydro.hpp"
#include "Problems/MagneticBraiding.hpp"
#include "Problems/UnstablePinch.hpp"

using namespace Cow;




// ============================================================================
int maraMainLoop (
    SimulationStatus& status,
    std::function<double ()> timestep,
    std::function<bool ()> condition,
    std::function<void (double)> advance,
    TaskScheduler& scheduler,
    Logger& logger)
{
    auto simulationTimer = Timer();
    auto sessionKzps = 0.0;
    auto sessionIter = 0;

    logger.log ("Mara") << "Beginning main loop" << std::endl;

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
        logger.log() << "["     << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        logger.log() << "t="    << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        logger.log() << "dt="   << std::setprecision (4) << std::scientific << dt << " ";
        logger.log() << "kzps=" << std::setprecision (2) << std::fixed << kzps << std::endl;

        sessionKzps += kzps;
        sessionIter +=1;
    }

    scheduler.dispatch (status);

    if (sessionIter > 0)
    {
        logger.log ("Mara") << "Completed main loop" << std::endl;
        logger.log ("Mara") << "Mean kzps was " << sessionKzps / sessionIter << std::endl;
    }
    return 0;
}




// ============================================================================
SimulationSetup::SimulationSetup()
{
    finalTime = 1.0;
    checkpointInterval = 1.0;
    vtkOutputInterval = 1.0;
    timeSeriesInterval = 0.1;
    vtkUseBinary = true;
    cflParameter = 0.25;
    rungeKuttaOrder = 2;
    outputDirectory = ".";
    runName = "test";
    disableCT = false;
}




// ============================================================================
SimulationStatus::SimulationStatus()
{
    simulationIter = 0;
    numCheckpoints = 0;
    numVtkOutputs = 0;
    numTimeSeriesEntries = 0;
    simulationTime = 0.0;
    lastCheckpoint = 0.0;
    lastVtkOutput = 0.0;
    lastTimeSeriesEntry = 0.0;
}

void SimulationStatus::update (const Variant::NamedValues& values)
{
    simulationIter = values.at ("simulationIter");
    numCheckpoints = values.at ("numCheckpoints");
    numVtkOutputs = values.at ("numVtkOutputs");
    numTimeSeriesEntries = values.at ("numTimeSeriesEntries");
    simulationTime = values.at ("simulationTime");
    lastCheckpoint = values.at ("lastCheckpoint");
    lastVtkOutput = values.at ("lastVtkOutput");
    lastTimeSeriesEntry = values.at ("lastTimeSeriesEntry");
}

Variant::NamedValues SimulationStatus::pack() const
{
    auto values = Variant::NamedValues();
    values["simulationIter"] = simulationIter;
    values["numCheckpoints"] = numCheckpoints;
    values["numVtkOutputs"] = numVtkOutputs;
    values["numTimeSeriesEntries"] = numTimeSeriesEntries;
    values["simulationTime"] = simulationTime;
    values["lastCheckpoint"] = lastCheckpoint;
    values["lastVtkOutput"] = lastVtkOutput;
    values["lastTimeSeriesEntry"] = lastTimeSeriesEntry;
    return values;
}

void SimulationStatus::print (std::ostream& stream)
{
    stream << "run status:\n" << pack();
}




// ============================================================================
class HelpProgram : public SubProgram
{
public:
    int run (int argc, const char* argv[]) override
    {
        std::cout <<
        "Mara is an astrophysics code for multi-dimensional "
        "gas and magnetofluid dynamics.\n";
        return 0;
    }
};




// ============================================================================
#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_FAST_COMPILE
#include "Catch.hpp"

class UnitTestProgram : public SubProgram
{
public:
    int run (int argc, const char* argv[]) override
    {
        int result = Catch::Session().run (argc, argv);
        return (result < 0xff ? result : 0xff);
    }
};




// ============================================================================
int main (int argc, const char* argv[])
{
    MpiSession mpiSession;
    std::set_terminate (Cow::terminateWithBacktrace);

    auto programs = std::map<std::string, std::unique_ptr<SubProgram>>();

    programs["help"]        = std::make_unique<HelpProgram>();
    programs["test"]        = std::make_unique<UnitTestProgram>();
    programs["regress-1d"]  = std::make_unique<Hydro1DTestProgram>();
    programs["regress-2d"]  = std::make_unique<Hydro2DTestProgram>();
    programs["regress-rel"] = std::make_unique<Relativistic1DTestProgram>();
    programs["regress-mhd"] = std::make_unique<NewtonianMHD2DTestProgram>();
    programs["braid"]       = std::make_unique<MagneticBraidingProgram>();
    programs["pinch"]       = std::make_unique<UnstablePinchProgram>();
    programs["kinetic"]     = std::make_unique<CollisionalHydroProgram>();
    programs["stitch"]      = std::make_unique<CheckpointStitcherProgram>();
    programs["tovtk"]       = std::make_unique<CheckpointToVtkProgram>();

    if (argc == 1)
    {
        std::cout << "usages: \n";

        for (auto& prog : programs)
        {
            std::cout << "    mara " << prog.first << std::endl;
        }
        return 0;
    }
    else if (programs.find (argv[1]) != programs.end())
    {
        return programs[argv[1]]->run (argc - 1, argv + 1);
    }
    else
    {
        return BackwardsCompatProgram().run (argc - 1, argv + 1);
    }
}
