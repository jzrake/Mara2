#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// Mara includes
#include "Mara.hpp"
#include "Configuration.hpp"
#include "FluxConservativeSystem.hpp"
#include "BlockDecomposition.hpp"

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






void writeVtkOutput (
    SimulationSetup& setup,
    SimulationStatus& status,
    FluxConservativeSystem& system)
{
    using namespace Cow;
    auto mpiWorld = MpiCommunicator::world();
    auto dir = setup.outputDirectory;
    mpiWorld.onMasterOnly ([&] () { FileSystem::ensureDirectoryExists (dir); });

    auto vtkFilename = FileSystem::makeFilename (dir, "mesh", ".vtk", status.vtkOutputsWrittenSoFar);
    auto vtkStream = std::ofstream (vtkFilename);
    auto vtkDataSet = VTK::DataSet (setup.meshGeometry->cellsShape());
    vtkDataSet.setTitle (setup.runName);
    vtkDataSet.setUseBinaryFormat (setup.vtkUseBinary);

    using VT = ConservationLaw::VariableType;
    auto claw = setup.conservationLaw;

    auto indexD = claw->getIndexFor (VT::density);
    auto indexP = claw->getIndexFor (VT::pressure);
    auto indexV = claw->getIndexFor (VT::velocity);
    auto indexB = claw->getIndexFor (VT::magnetic);

    if (indexD != -1)
    {
        auto D = system.getPrimitive (indexD);
        vtkDataSet.addScalarField ("density", D);
    }
    if (indexP != -1)
    {
        auto P = system.getPrimitive (indexP);
        vtkDataSet.addScalarField ("pressure", P);
    }
    if (indexV != -1)
    {
        auto V = system.getPrimitiveVector (indexV);
        vtkDataSet.addVectorField ("velocity", V);
    }
    if (indexB != -1)
    {
        auto B = system.getPrimitiveVector (indexB);
        vtkDataSet.addVectorField ("magnetic", B);

        // Write divergence of magnetic field (at mesh vertices)
        auto M = setup.constrainedTransport->computeMonopole (ConstrainedTransport::MeshLocation::vert);
        vtkDataSet.addScalarField ("monopole", M, VTK::DataSet::MeshLocation::vert);
    }

    vtkDataSet.addScalarField ("health", system.getZoneHealth());

    std::cout << "writing VTK file " << vtkFilename << std::endl;
    vtkDataSet.write (vtkStream);
    ++status.vtkOutputsWrittenSoFar;
}




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
    auto indexB = setup.conservationLaw->getIndexFor (VT::magnetic);

    comm.onMasterOnly ([&] ()
    {
        std::cout << "[Mara] writing checkpoint file " << filename << std::endl;


        // Create data directory, HDF5 checkpoint file, and groups
        // --------------------------------------------------------------------
        FileSystem::ensureDirectoryExists (outdir);
        auto file = H5::File (filename, "w");
        auto statusGroup = file.createGroup ("status");
        auto primitiveGroup = file.createGroup ("primitive");
        auto diagnosticGroup = file.createGroup ("diagnostic");
        auto globalShape = Array::vectorFromShape (block.getGlobalShape());


        // Create data sets for the primitive and diagnostic data
        // --------------------------------------------------------------------
        for (int q = 0; q < setup.conservationLaw->getNumConserved(); ++q)
        {
            auto field = setup.conservationLaw->getPrimitiveName(q);
            primitiveGroup.createDataSet (field, globalShape);
        }

        if (indexB != -1)
        {
            diagnosticGroup.createDataSet ("monopole", globalShape);
        }


        // Write the simulation status data into the checkpoint
        // --------------------------------------------------------------------
        statusGroup.write ("vtkOutputsWrittenSoFar", status.vtkOutputsWrittenSoFar);
        statusGroup.write ("checkpointsWrittenSoFar", status.checkpointsWrittenSoFar);
        statusGroup.write ("simulationIter", status.simulationIter);
        statusGroup.write ("simulationTime", status.simulationTime);
    });

    comm.inSequence ([&] (int rank)
    {
        // Open the file and groups
        // --------------------------------------------------------------------
        auto file = H5::File (filename, "a");
        auto primitiveGroup = file.getGroup ("primitive");
        auto diagnosticGroup = file.getGroup ("diagnostic");
        auto targetRegion = block.getPatchRegion();


        // Create data sets for the primitive and diagnostic data
        // --------------------------------------------------------------------
        for (int q = 0; q < setup.conservationLaw->getNumConserved(); ++q)
        {
            auto field = setup.conservationLaw->getPrimitiveName(q);
            auto dset = primitiveGroup.getDataSet (field);
            dset[targetRegion] = system.getPrimitive(q);
        }

        if (indexB != -1)
        {
            auto dset = diagnosticGroup.getDataSet ("monopole");
            auto M = setup.constrainedTransport->computeMonopole (ML::cell);
            dset[targetRegion] = M;
        }        
    });

    ++status.checkpointsWrittenSoFar;
}




// ============================================================================
int MaraSession::launch (SimulationSetup& setup)
{
    using namespace Cow;


    // More general setup validation code should go here
    // ------------------------------------------------------------------------
    if (setup.initialDataFunction == nullptr)
    {
        throw std::runtime_error ("No initial data function was provided");
    }


    // Mesh decomposition steps
    // ------------------------------------------------------------------------
    auto blockDecomposition = BlockDecomposition (setup.meshGeometry);
    setup.meshGeometry = blockDecomposition.decompose();
    setup.boundaryCondition = blockDecomposition.createBoundaryCondition (setup.boundaryCondition);

    auto status = SimulationStatus();
    auto system = FluxConservativeSystem (setup); // This also initializes CT.

    system.setInitialData (setup.initialDataFunction, setup.vectorPotentialFunction);

    if (setup.vtkOutputInterval > 0)
    {
        writeVtkOutput (setup, status, system);
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
            writeVtkOutput (setup, status, system);
        }

        if (setup.checkpointInterval > 0 && status.simulationTime >= nextChkpt)
        {
            writeCheckpoint (setup, status, system, blockDecomposition);
        }


        // Invoke the solver to advance the solution
        // --------------------------------------------------------------------
        auto timer = Cow::Timer();
        double dt = setup.cflParameter * system.getCourantTimestep();
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

    std::set_terminate (Cow::terminateWithBacktrace);

    if (argc == 1)
    {
        std::cout << "usages: \n";
        std::cout << "\tmara config.lua\n";
        std::cout << "\tmara run script.lua\n";
        std::cout << "\tmara help\n";
        return 0;
    }

    auto session = MaraSession();
    auto configuration = Configuration();
    auto command = std::string (argv[1]);

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
            std::cout << "'run': no script provided\n";
            return 0;
        }
        return configuration.launchFromScript (session, argv[2]);
    }
    else
    {
        auto setup = configuration.fromLuaFile (argv[1]);
        return session.launch (setup);
    }

    return 0;
}
