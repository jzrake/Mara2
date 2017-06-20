#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

// Mara includes
#include "Mara.hpp"
#include "Configuration.hpp"
#include "FluxConservativeSystem.hpp"
#include "BlockDecomposition.hpp"
#include "CartesianMeshGeometry.hpp"
#include "TimeSeriesManager.hpp"
#include "TestSuite.hpp"

// Cow includes
#include "HDF5.hpp"
#include "Timer.hpp"
#include "VTK.hpp"
#include "MPI.hpp"
#include "FileSystem.hpp"
#include "DebugHelper.hpp"
#include "Logger.hpp"

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
void checkpointToVtk (std::string filename, Logger& logger)
{
    using namespace Cow;

    auto h5Filename = filename;
    auto doth5 = filename.find (".h5");
    auto vtkFilename = filename.replace (doth5, 5, ".vtk");

    auto file = H5::File (h5Filename, "r");
    auto points = file.getGroup ("mesh").getGroup ("points");
    auto X = points.readArray ("x");
    auto Y = points.readArray ("y");
    auto Z = points.readArray ("z");
    auto cellsShape = Shape {{X.size() - 1, Y.size() - 1, Z.size() - 1, 1, 1}};
    auto primitive = file.getGroup ("primitive");

    auto vtkStream = std::ofstream (vtkFilename);
    auto vtkMesh = VTK::RectilinearGrid (cellsShape);
    vtkMesh.setTitle (file.readString ("run_name"));
    vtkMesh.setUseBinaryFormat (true);
    vtkMesh.setPointCoordinates (X, 0);
    vtkMesh.setPointCoordinates (Y, 1);
    vtkMesh.setPointCoordinates (Z, 2);

    auto velocityNames = std::vector<std::string> {"velocity1", "velocity2", "velocity3"};
    auto magneticNames = std::vector<std::string> {"magnetic1", "magnetic2", "magnetic3"};

    if (primitive.hasDataSet ("density"))
    {
        vtkMesh.addScalarField ("density", primitive.readArray ("density"));
    }
    if (primitive.hasDataSet ("pressure"))
    {
        vtkMesh.addScalarField ("pressure", primitive.readArray ("pressure"));
    }
    if (primitive.hasDataSets (velocityNames))
    {
        vtkMesh.addVectorField ("velocity", primitive.readArrays (velocityNames, 3));
    }
    if (primitive.hasDataSets (magneticNames))
    {
        vtkMesh.addVectorField ("magnetic", primitive.readArrays (magneticNames, 3));
    }

    logger.log ("Mara") << h5Filename << " -> " << vtkFilename << std::endl;
    vtkMesh.write (vtkStream);
}




// ============================================================================
void writeVtkOutput (
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system,
    BlockDecomposition& block,
    Logger& logger)
{
    using namespace Cow;

    if (block.getCommunicator().size() != 1)
    {
        block.getCommunicator().onMasterOnly ([&] ()
        {
            logger.log ("Mara") << "Warning: no support for VTK output in parallel runs\n";
        });

        status.lastVtkOutput += setup.vtkOutputInterval;
        status.numVtkOutputs += 1;
        return;
    }

    auto geometry = dynamic_cast<const CartesianMeshGeometry*> (setup.meshGeometry.get());

    if (geometry == nullptr)
    {
        logger.log ("Mara") << "Warning: only CartesianMeshGeometry supports VTK output\n";
        return;
    }

    auto dir = setup.outputDirectory;
    FileSystem::ensureDirectoryExists (dir);

    auto vtkFilename = FileSystem::makeFilename (dir, "mesh", ".vtk", status.numVtkOutputs);
    auto vtkStream = std::ofstream (vtkFilename);
    auto vtkMesh = VTK::RectilinearGrid (setup.meshGeometry->cellsShape());
    vtkMesh.setTitle (setup.runName);
    vtkMesh.setUseBinaryFormat (setup.vtkUseBinary);
    vtkMesh.setPointCoordinates (geometry->getPointCoordinates(0), 0);
    vtkMesh.setPointCoordinates (geometry->getPointCoordinates(1), 1);
    vtkMesh.setPointCoordinates (geometry->getPointCoordinates(2), 2);

    using VT = ConservationLaw::VariableType;
    auto claw = setup.conservationLaw;

    auto indexD = claw->getIndexFor (VT::density);
    auto indexP = claw->getIndexFor (VT::pressure);
    auto indexV = claw->getIndexFor (VT::velocity);
    auto indexB = claw->getIndexFor (VT::magnetic);

    if (indexD != -1)
    {
        auto D = system.getPrimitive (indexD);
        vtkMesh.addScalarField ("density", D);
    }
    if (indexP != -1)
    {
        auto P = system.getPrimitive (indexP);
        vtkMesh.addScalarField ("pressure", P);
    }
    if (indexV != -1)
    {
        auto V = system.getPrimitiveVector (indexV);
        vtkMesh.addVectorField ("velocity", V);
    }
    if (indexB != -1)
    {
        auto B = system.getPrimitiveVector (indexB);
        vtkMesh.addVectorField ("magnetic", B);

        // Write divergence of magnetic field (at mesh vertices)
        auto M = setup.constrainedTransport->computeMonopole (MeshLocation::vert);
        vtkMesh.addScalarField ("monopole", M, VTK::RectilinearGrid::MeshLocation::vert);
    }

    vtkMesh.addScalarField ("health", system.getZoneHealth());
    logger.log ("Mara") << "writing VTK file " << vtkFilename << std::endl;
    vtkMesh.write (vtkStream);

    status.lastVtkOutput += setup.vtkOutputInterval;
    status.numVtkOutputs += 1;
}




// ============================================================================
void writeCheckpoint (
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system,
    BlockDecomposition& block,
    TimeSeriesManager& tseries,
    Logger& logger)
{
    using VT = ConservationLaw::VariableType; // For VT::magnetic
    using ML = MeshLocation; // For ML::cell

    auto timer = Timer();
    auto comm = block.getCommunicator();
    auto outdir = setup.outputDirectory;
    auto filename = FileSystem::makeFilename (outdir, "chkpt", ".h5", status.numCheckpoints);


    // Update the status before writing it
    // ------------------------------------------------------------------------
    status.lastCheckpoint += setup.checkpointInterval;
    status.numCheckpoints += 1;


    comm.onMasterOnly ([&] ()
    {
        logger.log ("Mara") << "writing checkpoint file " << filename << std::endl;


        // Create data directory, HDF5 checkpoint file, and groups
        // --------------------------------------------------------------------
        FileSystem::ensureDirectoryExists (outdir);
        auto file            = H5::File (filename, "w");
        auto statusGroup     = file.createGroup ("status");
        auto primitiveGroup  = file.createGroup ("primitive");
        auto diagnosticGroup = file.createGroup ("diagnostic");
        auto meshGroup       = file.createGroup ("mesh");
        auto tseriesGroup    = file.createGroup ("time_series");


        // Write misc data like date and run script
        // --------------------------------------------------------------------
        auto t = std::time (nullptr);
        auto tm = *std::localtime (&t);
        auto oss = std::ostringstream();
        oss << std::put_time (&tm, "%m-%d-%Y %H:%M:%S");

        file.writeString ("run_name", setup.runName);
        file.writeString ("date", oss.str());
        file.writeString ("lua_script", setup.luaScript);
        file.writeString ("lua_command_line", setup.luaCommandLine);


        // Create data sets for the primitive and diagnostic data
        // --------------------------------------------------------------------
        auto globalShape = Array::vectorFromShape (block.getGlobalShape());
        auto localsShape = Array::vectorFromShape (setup.meshGeometry->cellsShape());
        auto dtype = H5::DataType::nativeDouble();
        auto plist = H5::PropertyList::DataSetCreate().setChunk (localsShape);

        for (int q = 0; q < setup.conservationLaw->getNumConserved(); ++q)
        {
            auto field = setup.conservationLaw->getPrimitiveName(q);
            primitiveGroup.createDataSet (field, globalShape, dtype, plist);
        }

        if (setup.conservationLaw->getIndexFor (VT::magnetic) != -1)
        {
            diagnosticGroup.createDataSet ("monopole", globalShape, dtype, plist);
        }


        // Write the simulation status data into the checkpoint
        // --------------------------------------------------------------------
        for (auto member : status.pack())
        {
            statusGroup.writeVariant (member.first, member.second);
        }


        // Write the mesh information (needs to be generalized to non-rectilinear meshes)
        // --------------------------------------------------------------------
        meshGroup.writeString ("type", "cartesian");

        auto geometry = dynamic_cast<const CartesianMeshGeometry*> (block.getGlobalGeometry().get());
        auto pointGroup = meshGroup.createGroup ("points");
        pointGroup.writeArray ("x", geometry->getPointCoordinates (0));
        pointGroup.writeArray ("y", geometry->getPointCoordinates (1));
        pointGroup.writeArray ("z", geometry->getPointCoordinates (2));


        // Write the time series data
        // --------------------------------------------------------------------
        tseries.write (tseriesGroup);
    });

    comm.inSequence ([&] (int rank)
    {
        // Open the file and groups
        // --------------------------------------------------------------------
        auto file            = H5::File (filename, "a");
        auto primitiveGroup  = file.getGroup ("primitive");
        auto diagnosticGroup = file.getGroup ("diagnostic");


        // Create data sets for the primitive and diagnostic data
        // --------------------------------------------------------------------
        auto targetRegion = block.getPatchRegion();

        for (int q = 0; q < setup.conservationLaw->getNumConserved(); ++q)
        {
            auto field = setup.conservationLaw->getPrimitiveName(q);
            auto dset = primitiveGroup.getDataSet (field);
            dset[targetRegion] = system.getPrimitive(q);
        }

        if (setup.conservationLaw->getIndexFor (VT::magnetic) != -1)
        {
            auto dset = diagnosticGroup.getDataSet ("monopole");
            auto M = setup.constrainedTransport->computeMonopole (ML::cell);
            dset[targetRegion] = M;
        }        
    });

    logger.log ("Mara") << "required " << timer.ageInSeconds() << std::endl;
}




// ============================================================================
void readCheckpoint (
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system,
    BlockDecomposition& block,
    TimeSeriesManager& tseries,
    Logger& logger)
{
    logger.log ("Mara") << "reading checkpoint file " << setup.restartFile << std::endl;

    // block.getCommunicator().inSequence ([&] (int rank)
    // Note: we read all procs at once here, because BC's will hang otherwise

    auto file            = H5::File (setup.restartFile, "r");
    auto statusGroup     = file.getGroup ("status");
    auto primitiveGroup  = file.getGroup ("primitive");
    auto tseriesGroup    = file.getGroup ("time_series");

    system.assignPrimitive (primitiveGroup.readArrays (
        setup.conservationLaw->getPrimitiveNames(), 3, block.getPatchRegion()));

    auto packedStatus = status.pack();

    for (auto& member : packedStatus)
    {
        member.second = statusGroup.readVariant (member.first);
    }
    status.update (packedStatus);
    // status.update (statusGroup.readNamedValues()); <-- add Cow feature

    status.print (logger.log ("Mara"));
    tseries.load (tseriesGroup);
}




// ============================================================================
void doOutputStage(
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system,
    BlockDecomposition& block,
    TimeSeriesManager& tseries,
    Logger& logger)
{
    auto isTime = [=] (double interval, double next) -> bool
    {
        return interval > 0 && status.simulationTime >= next;
    };

    double nextVtk   = status.lastVtkOutput       + setup.vtkOutputInterval;
    double nextChkpt = status.lastCheckpoint      + setup.checkpointInterval;
    double nextEntry = status.lastTimeSeriesEntry + setup.timeSeriesInterval;

    if (isTime (setup.timeSeriesInterval, nextEntry))
    {
        auto volumeIntegrated = system.volumeIntegratedDiagnostics();
        auto volumeAveraged = block.volumeAverageOverPatches (volumeIntegrated);
        auto fieldNames = setup.conservationLaw->getDiagnosticNames();
        auto entry = Variant::NamedValues();

        for (unsigned int n = 0; n < fieldNames.size(); ++n)
        {
            entry[fieldNames[n]] = volumeAveraged[n];
        }
        tseries.append (status, entry);

        status.lastTimeSeriesEntry += setup.timeSeriesInterval;
        status.numTimeSeriesEntries += 1;
    }

    if (isTime (setup.vtkOutputInterval, nextVtk))
    {
        writeVtkOutput (setup, status, system, block, logger);
    }

    if (isTime (setup.checkpointInterval, nextChkpt))
    {
        writeCheckpoint (setup, status, system, block, tseries, logger);
    }
}




// ============================================================================
MaraSession::MaraSession()
{
    if (MpiCommunicator::world().isThisMaster())
    {
        logger->setLogToStdout();
    }
    else
    {
        logger->setLogToNull();
    }
}

SimulationStatus MaraSession::launch (SimulationSetup& setup)
{
    // "Dependency injection"
    // ------------------------------------------------------------------------
    auto block = BlockDecomposition (setup.meshGeometry, *logger);
    auto tseries = TimeSeriesManager();


    setup.meshGeometry = block.decompose();
    setup.boundaryCondition->setMeshGeometry (setup.meshGeometry);
    setup.boundaryCondition->setConservationLaw (setup.conservationLaw);
    setup.boundaryCondition->setBoundaryValueFunction (setup.boundaryValueFunction);


    setup.boundaryCondition = block.createBoundaryCondition (setup.boundaryCondition);
    setup.constrainedTransport->setMeshGeometry (setup.meshGeometry);
    setup.constrainedTransport->setBoundaryCondition (setup.boundaryCondition);


    tseries.setLogger (logger);


    auto status = SimulationStatus();
    auto system = FluxConservativeSystem (setup);


    // Initial solution data
    // ------------------------------------------------------------------------
    if (setup.initialDataFunction != nullptr)
    {
        system.setInitialData (setup.initialDataFunction, setup.vectorPotentialFunction);
    }
    else if (! setup.restartFile.empty())
    {
        readCheckpoint (setup, status, system, block, tseries, *logger);
    }
    else
    {
        throw std::runtime_error ("No initial data function or restart file was given");
    }


    // Improve on this thing here:
    bool writeInitial = status.simulationIter == 0;

    if (writeInitial)
    {
        status.lastCheckpoint      -= setup.checkpointInterval;
        status.lastVtkOutput       -= setup.vtkOutputInterval;
        status.lastTimeSeriesEntry -= setup.timeSeriesInterval;
    }
    doOutputStage (setup, status, system, block, tseries, *logger);


    while (status.simulationTime < setup.finalTime)
    {
        // Invoke the solver to advance the solution
        // --------------------------------------------------------------------
        auto timer = Timer();
        double dt = block.getCommunicator().minimum (setup.cflParameter * system.getCourantTimestep());
        system.advance (dt);
        status.simulationTime += dt;
        status.simulationIter += 1;


        // Generate iteration output message
        // --------------------------------------------------------------------
        double kzps = 1e-3 * setup.meshGeometry->totalCellsInMesh() / timer.age();
        logger->log() << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        logger->log() << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        logger->log() << "dt=" << std::setprecision (4) << std::scientific << dt << " ";
        logger->log() << "kzps=" << std::setprecision (2) << std::fixed << kzps << std::endl;


        doOutputStage (setup, status, system, block, tseries, *logger);
    }

    return status;
}




// ============================================================================
#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_FAST_COMPILE
#include "Catch.hpp"
#include "TestSuite.hpp"

int TestSuite::runAllTests (int argc, const char* argv[])
{
    int result = Catch::Session().run (argc, argv);
    return (result < 0xff ? result : 0xff);
}




#include "Problems.hpp"




// ============================================================================
int main (int argc, const char* argv[])
{
    using namespace Cow;
    MpiSession mpiSession;

    std::set_terminate (Cow::terminateWithBacktrace);
    //std::set_terminate (Cow::terminateWithPrintException);

    if (argc == 1)
    {
        std::cout << "usages: \n";
        std::cout << "\tmara config.lua\n";
        std::cout << "\tmara run script.lua\n";
        std::cout << "\tmara tovtk chkpt.*.h5\n";
        std::cout << "\tmara chkpt.0000.h5\n";
        std::cout << "\tmara test\n";
        std::cout << "\tmara help\n";
        return 0;
    }

    auto session = MaraSession();
    auto configuration = Configuration (argc, argv);
    auto command = std::string (argv[1]);
    auto extension = FileSystem::fileExtension (command);

    if (command == "help")
    {
        std::cout <<
        "Mara is an astrophysics code for multi-dimensional "
        "gas and magnetofluid dynamics.\n";
        return 0;
    }
    if (command == "test")
    {
        return TestSuite().runAllTests (argc - 1, argv + 1);
    }
    if (command == "regress")
    {
        auto program = SimpleTestProgram();
        return program.run (argc - 1, argv + 1);
    }
    if (command == "run")
    {
        if (argc < 3)
        {
            std::cout << "'run': no script provided" << std::endl;
            return 0;
        }
        return configuration.launchFromScript (session, argv[2]);
    }
    if (command == "tovtk")
    {
        Logger logger;

        for (int n = 2; n < argc; ++n)
        {
            checkpointToVtk (argv[n], logger);
        }
        return 0;
    }
    if (extension == ".lua")
    {
        auto setup = configuration.fromLuaFile (argv[1]);
        session.launch (setup);
        return 0;
    }
    if (extension == ".h5")
    {
        auto setup = configuration.fromCheckpoint (argv[1]);
        session.launch (setup);
        return 0;
    }

    std::cout << "unrecognized command " << command << std::endl;
    return 0;
}
