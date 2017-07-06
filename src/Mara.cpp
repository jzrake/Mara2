// Mara includes
#include "Mara.hpp"
#include "Problems.hpp"
#include "BackwardsCompat.hpp"
#include "Checkpoint.hpp"

// Cow includes
#include "MPI.hpp"
#include "DebugHelper.hpp"

using namespace Cow;




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

    programs["help"]        .reset (new HelpProgram);
    programs["test"]        .reset (new UnitTestProgram);
    programs["regress-1d"]  .reset (new Hydro1DTestProgram);
    programs["regress-2d"]  .reset (new Hydro2DTestProgram);
    programs["regress-mhd"] .reset (new NewtonianMHD2DTestProgram);
    programs["braid"]       .reset (new MagneticBraidingProgram);
    programs["stitch"]      .reset (new CheckpointStitcherProgram);
    programs["tovtk"]       .reset (new CheckpointToVtkProgram);

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
