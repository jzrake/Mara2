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

// Cow includes
#include "HDF5.hpp"
#include "Timer.hpp"
#include "VTK.hpp"
#include "MPI.hpp"
#include "FileSystem.hpp"
#include "DebugHelper.hpp"




// ============================================================================
SimulationSetup::SimulationSetup()
{
    finalTime = 1.0;
    checkpointInterval = 1.0;
    vtkOutputInterval = 1.0;
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
    simulationTime = 0.0;
    simulationIter = 0;
    checkpointsWrittenSoFar = 0;
    vtkOutputsWrittenSoFar = 0;
}




// ============================================================================
void checkpointToVtk (std::string filename)
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

    std::cout << "[Mara] " << h5Filename << " -> " << vtkFilename << std::endl;
    vtkMesh.write (vtkStream);
}




// ============================================================================
void writeVtkOutput (
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system,
    BlockDecomposition& block)
{
    using namespace Cow;

    if (block.getCommunicator().size() != 1)
    {
        block.getCommunicator().onMasterOnly ([&] ()
        {
            // NOTE: possible future syntax for logging
            // Log::warning ("Mara") << "no support for VTK output\n";
            std::cout << "[Mara] Warning: no support for VTK output in parallel runs\n";
        });

        ++status.vtkOutputsWrittenSoFar;
        return;
    }

    auto geometry = dynamic_cast<const CartesianMeshGeometry*> (setup.meshGeometry.get());

    if (geometry == nullptr)
    {
        std::cout << "[Mara] Warning: only CartesianMeshGeometry supports VTK output\n";
        return;
    }

    auto dir = setup.outputDirectory;
    FileSystem::ensureDirectoryExists (dir);

    auto vtkFilename = FileSystem::makeFilename (dir, "mesh", ".vtk", status.vtkOutputsWrittenSoFar);
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
        auto M = setup.constrainedTransport->computeMonopole (ConstrainedTransport::MeshLocation::vert);
        vtkMesh.addScalarField ("monopole", M, VTK::RectilinearGrid::MeshLocation::vert);
    }

    vtkMesh.addScalarField ("health", system.getZoneHealth());

    std::cout << "[Mara] writing VTK file " << vtkFilename << std::endl;
    vtkMesh.write (vtkStream);
    ++status.vtkOutputsWrittenSoFar;
}




// ============================================================================
void writeCheckpoint (
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system,
    BlockDecomposition& block)
{
    using namespace Cow;
    using VT = ConservationLaw::VariableType; // For VT::magnetic
    using ML = ConstrainedTransport::MeshLocation; // For ML::cell

    auto comm = block.getCommunicator();
    auto outdir = setup.outputDirectory;
    auto filename = FileSystem::makeFilename (outdir, "chkpt", ".h5", status.checkpointsWrittenSoFar);

    comm.onMasterOnly ([&] ()
    {
        std::cout << "[Mara] writing checkpoint file " << filename << std::endl;


        // Create data directory, HDF5 checkpoint file, and groups
        // --------------------------------------------------------------------
        FileSystem::ensureDirectoryExists (outdir);
        auto file            = H5::File (filename, "w");
        auto statusGroup     = file.createGroup ("status");
        auto primitiveGroup  = file.createGroup ("primitive");
        auto diagnosticGroup = file.createGroup ("diagnostic");
        auto meshGroup       = file.createGroup ("mesh");


        // Write misc data like date and run script
        // --------------------------------------------------------------------
        auto t = std::time (nullptr);
        auto tm = *std::localtime (&t);
        std::ostringstream oss;
        oss << std::put_time (&tm, "%m-%d-%Y %H:%M:%S");

        file.writeString ("run_name", setup.runName);
        file.writeString ("date", oss.str());
        file.writeString ("script", setup.luaScript);


        // Create data sets for the primitive and diagnostic data
        // --------------------------------------------------------------------
        auto globalShape = Array::vectorFromShape (block.getGlobalShape());

        for (int q = 0; q < setup.conservationLaw->getNumConserved(); ++q)
        {
            auto field = setup.conservationLaw->getPrimitiveName(q);
            primitiveGroup.createDataSet (field, globalShape);
        }

        if (setup.conservationLaw->getIndexFor (VT::magnetic) != -1)
        {
            diagnosticGroup.createDataSet ("monopole", globalShape);
        }


        // Write the simulation status data into the checkpoint
        // --------------------------------------------------------------------
        statusGroup.writeInt ("vtkOutputsWrittenSoFar", status.vtkOutputsWrittenSoFar);
        statusGroup.writeInt ("checkpointsWrittenSoFar", status.checkpointsWrittenSoFar);
        statusGroup.writeInt ("simulationIter", status.simulationIter);
        statusGroup.writeDouble ("simulationTime", status.simulationTime);


        // Write the mesh information (needs to be generalized to non-rectilinear meshes)
        // --------------------------------------------------------------------
        meshGroup.writeString ("type", "CartesianMeshGeometry");

        auto geometry = dynamic_cast<const CartesianMeshGeometry*> (block.getGlobalGeometry().get());
        auto pointGroup = meshGroup.createGroup ("points");
        pointGroup.writeArray ("x", geometry->getPointCoordinates (0));
        pointGroup.writeArray ("y", geometry->getPointCoordinates (1));
        pointGroup.writeArray ("z", geometry->getPointCoordinates (2));
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

    ++status.checkpointsWrittenSoFar;
}




// ============================================================================
void readCheckpoint (
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system,
    BlockDecomposition& block)
{
    block.getCommunicator().onMasterOnly ([&] ()
    {
        std::cout << "[Mara] reading checkpoint file " << setup.restartFile << std::endl;
    });

    // block.getCommunicator().inSequence ([&] (int rank)
    // Note: we read all procs at once here, because BC's will hang otherwise
    {
        auto file            = Cow::H5::File (setup.restartFile, "r");
        auto statusGroup     = file.getGroup ("status");
        auto primitiveGroup  = file.getGroup ("primitive");

        system.assignPrimitive (primitiveGroup.readArrays (
            setup.conservationLaw->getPrimitiveNames(), 3, block.getPatchRegion()));

        status.vtkOutputsWrittenSoFar  = statusGroup.readInt ("vtkOutputsWrittenSoFar");
        status.checkpointsWrittenSoFar = statusGroup.readInt ("checkpointsWrittenSoFar");
        status.simulationIter          = statusGroup.readInt ("simulationIter");
        status.simulationTime          = statusGroup.readDouble ("simulationTime");
    }//);
}




// ============================================================================
int MaraSession::launch (SimulationSetup& setup)
{
    // "Dependency injection"
    // ------------------------------------------------------------------------
    auto blockDecomposition = BlockDecomposition (setup.meshGeometry);
    auto blockCommunicator = blockDecomposition.getCommunicator();

    setup.meshGeometry = blockDecomposition.decompose();

    setup.boundaryCondition->setMeshGeometry (setup.meshGeometry);
    setup.boundaryCondition->setConservationLaw (setup.conservationLaw);
    setup.boundaryCondition->setBoundaryValueFunction (setup.boundaryValueFunction);
    setup.boundaryCondition = blockDecomposition.createBoundaryCondition (setup.boundaryCondition);

    setup.constrainedTransport->setMeshGeometry (setup.meshGeometry);
    setup.constrainedTransport->setBoundaryCondition (setup.boundaryCondition);

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
        readCheckpoint (setup, status, system, blockDecomposition);
    }
    else
    {
        throw std::runtime_error ("No initial data function or restart file was given");
    }


    // Initial output stage
    // --------------------------------------------------------------------
    if (setup.vtkOutputInterval > 0)
    {
        writeVtkOutput (setup, status, system, blockDecomposition);
    }
    if (setup.checkpointInterval > 0)
    {
        writeCheckpoint (setup, status, system, blockDecomposition);
    }


    while (status.simulationTime < setup.finalTime)
    {
        // Perform output of different formats if necessary
        // --------------------------------------------------------------------
        double nextVtk = status.vtkOutputsWrittenSoFar * setup.vtkOutputInterval;
        double nextChkpt = status.checkpointsWrittenSoFar * setup.checkpointInterval;

        if (setup.vtkOutputInterval > 0 && status.simulationTime >= nextVtk)
        {
            writeVtkOutput (setup, status, system, blockDecomposition);
        }

        if (setup.checkpointInterval > 0 && status.simulationTime >= nextChkpt)
        {
            writeCheckpoint (setup, status, system, blockDecomposition);
        }


        // Invoke the solver to advance the solution
        // --------------------------------------------------------------------
        auto timer = Cow::Timer();
        double dt = blockCommunicator.minimum (setup.cflParameter * system.getCourantTimestep());
        system.advance (dt);
        status.simulationTime += dt;
        status.simulationIter += 1;


        // Generate iteration output message
        // --------------------------------------------------------------------
        blockDecomposition.getCommunicator().onMasterOnly ([&] ()
        {
            double kzps = 1e-3 * setup.meshGeometry->totalCellsInMesh() / timer.age();
            std::cout << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
            std::cout << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
            std::cout << "dt=" << std::setprecision (4) << std::scientific << dt << " ";
            std::cout << "kzps=" << std::setprecision (2) << std::fixed << kzps << std::endl;
        });
    }

    return 0;
}




// ============================================================================
int main (int argc, const char* argv[])
{
    using namespace Cow;
    MpiSession mpiSession;

    // std::set_terminate (Cow::terminateWithBacktrace);
    std::set_terminate (Cow::terminateWithPrintException);

    if (argc == 1)
    {
        std::cout << "usages: \n";
        std::cout << "\tmara config.lua\n";
        std::cout << "\tmara run script.lua\n";
        std::cout << "\tmara tovtk chkpt.*.h5\n";
        std::cout << "\tmara chkpt.0000.h5\n";
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
    else if (command == "run")
    {
        if (argc < 3)
        {
            std::cout << "'run': no script provided" << std::endl;
            return 0;
        }
        return configuration.launchFromScript (session, argv[2]);
    }
    else if (command == "tovtk")
    {
        for (int n = 2; n < argc; ++n)
        {
            checkpointToVtk (argv[n]);
        }
    }
    else if (extension == ".lua")
    {
        auto setup = configuration.fromLuaFile (argv[1]);
        return session.launch (setup);
    }
    else if (extension == ".h5")
    {
        auto setup = configuration.fromCheckpoint (argv[1]);
        return session.launch (setup);
    }
    else
    {
        std::cout << "unrecognized command " << command << std::endl;
    }
    return 0;
}
