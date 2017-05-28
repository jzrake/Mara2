#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// Mara includes
#include "Mara.hpp"
#include "Configuration.hpp"
#include "FluxConservativeSystem.hpp"
#include "RiemannSolver.hpp"
#include "Matrix.hpp"
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







// ============================================================================
ConservationLaw::State ScalarUpwind::intercellFlux (const FaceData& faceData) const
{
    ConservationLaw::Request request;
    request.areaElement = faceData.areaElement;

    const double* PL = &faceData.stencilData (0);
    const double* PR = &faceData.stencilData (1);

    auto L = faceData.conservationLaw->fromPrimitive (request, PL);
    auto R = faceData.conservationLaw->fromPrimitive (request, PR);

    UpwindRiemannSolver riemannSolver;
    return riemannSolver.solve (L, R, faceData.areaElement);
}

int ScalarUpwind::getStencilSize() const
{
    return 1;
}




// ============================================================================
MethodOfLines::MethodOfLines()
{
    riemannSolver.reset (new HlleRiemannSolver);
}

void MethodOfLines::setRiemannSolver (std::shared_ptr<RiemannSolver> solverToUse)
{
    riemannSolver = solverToUse;
}

ConservationLaw::State MethodOfLines::intercellFlux (const FaceData& faceData) const
{
    ConservationLaw::Request request;
    request.areaElement = faceData.areaElement;

    const double* P0 = &faceData.stencilData(0);
    const double* P1 = &faceData.stencilData(1);

    auto L = faceData.conservationLaw->fromPrimitive (request, P0);
    auto R = faceData.conservationLaw->fromPrimitive (request, P1);

    return riemannSolver->solve (L, R, faceData.areaElement);
}

int MethodOfLines::getStencilSize() const
{
    return 1;
}




// ============================================================================
MethodOfLinesPlm::MethodOfLinesPlm (double plmTheta)
{
    plm.setPlmTheta (plmTheta);
    riemannSolver.reset (new HlleRiemannSolver);
}

void MethodOfLinesPlm::setRiemannSolver (std::shared_ptr<RiemannSolver> solverToUse)
{
    riemannSolver = solverToUse;
}

ConservationLaw::State MethodOfLinesPlm::intercellFlux (const FaceData& faceData) const
{
    ConservationLaw::Request request;
    request.areaElement = faceData.areaElement;

    const double* P0 = &faceData.stencilData(0);
    const double* P1 = &faceData.stencilData(1);
    const double* P2 = &faceData.stencilData(2);
    const double* P3 = &faceData.stencilData(3);
    const int nq = faceData.conservationLaw->getNumConserved();

    std::vector<double> PL(nq);
    std::vector<double> PR(nq);

    for (int q = 0; q < nq; ++q)
    {
        const double ps[4] = {P0[q], P1[q], P2[q], P3[q]};
        PL[q] = plm.reconstruct (&ps[1], Reconstruction::PLM_C2R);
        PR[q] = plm.reconstruct (&ps[2], Reconstruction::PLM_C2L);
    }
    auto L = faceData.conservationLaw->fromPrimitive (request, &PL[0]);
    auto R = faceData.conservationLaw->fromPrimitive (request, &PR[0]);

    return riemannSolver->solve (L, R, faceData.areaElement);
}

int MethodOfLinesPlm::getStencilSize() const
{
    return 2;
}




// ============================================================================
MethodOfLinesWeno::MethodOfLinesWeno (double shenZhaA)
{
    weno.setSmoothnessIndicator (Reconstruction::ImprovedShenZha10);
    weno.setShenZha10A (shenZhaA);
}

ConservationLaw::State MethodOfLinesWeno::intercellFlux (const FaceData& faceData) const
{
    auto request = ConservationLaw::Request();
    request.areaElement = faceData.areaElement;

    const auto claw = faceData.conservationLaw;
    const auto states = claw->fromPrimitive (request, faceData.stencilData);
    const auto maxLam = claw->maxEigenvalueMagnitude (states);
    const auto Shat = claw->averageStates (request, states[2], states[3]);
    const int nq = claw->getNumConserved();

    Cow::Matrix Fp (nq, 6);
    Cow::Matrix Fm (nq, 6);

    for (int n = 0; n < 6; ++n)
    {
        for (int q = 0; q < nq; ++q)
        {
            const auto& S = states[n];
            Fp (q, n) = S.F[q] + maxLam * S.U[q]; // Lax-Friedrichs flux splitting
            Fm (q, n) = S.F[q] - maxLam * S.U[q];
        }
    }

    auto fp = Shat.L * Fp;
    auto fm = Shat.L * Fm;
    auto fhat = Cow::Matrix (nq, 1); // column vector (nq rows)

    for (int q = 0; q < nq; ++q)
    {
        double fhatp = weno.reconstruct (&fp (q, 2), Reconstruction::WENO5_FD_C2R);
        double fhatm = weno.reconstruct (&fm (q, 3), Reconstruction::WENO5_FD_C2L);
        fhat (q, 0) = 0.5 * (fhatp + fhatm);
    }

    auto Fhat = Shat.R * fhat;
    auto S = ConservationLaw::State();

    for (int q = 0; q < nq; ++q)
    {
        S.F[q] = Fhat (q, 0);
    }
    return S;
}

int MethodOfLinesWeno::getStencilSize() const
{
    return 3;
}




void writeVtkOutput (SimulationSetup& setup, SimulationStatus& status, FluxConservativeSystem& system)
{
    using namespace Cow;
    auto mpiWorld = MpiCommunicator::world();
    auto dir = setup.outputDirectory;
    mpiWorld.onMasterOnly ([&] () { FileSystem::ensureDirectoryExists (dir); });

    auto vtkFilename = FileSystem::makeFilename (dir, "mesh", ".vtk", status.vtkOutputsWrittenSoFar);
    auto vtkStream = std::ofstream (vtkFilename);
    auto vtkDataSet = VTK::DataSet (setup.meshGeometry->domainShape());
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




void writeCheckpoint (SimulationSetup& setup, SimulationStatus& status, FluxConservativeSystem& system)
{
    using namespace Cow;
    using VT = ConservationLaw::VariableType; // For VT::magnetic

    auto mpiWorld = MpiCommunicator::world();
    auto dir = setup.outputDirectory;

    mpiWorld.onMasterOnly ([&] () { FileSystem::ensureDirectoryExists (dir); });

    auto h5Filename = FileSystem::makeFilename (dir, "chkpt", ".h5", status.checkpointsWrittenSoFar);
    auto file = H5::File (h5Filename, "w");

    std::cout << "writing checkpoint file " << h5Filename << std::endl;

    for (int q = 0; q < setup.conservationLaw->getNumConserved(); ++q)
    {
        auto field = setup.conservationLaw->getPrimitiveName(q);
        file.write (field, system.getPrimitive(q));
    }


    auto indexB = setup.conservationLaw->getIndexFor (VT::magnetic);

    if (indexB != -1)
    {
        auto M = setup.constrainedTransport->computeMonopole (ConstrainedTransport::MeshLocation::vert);
        file.write ("monopole", M);
    }

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
    auto localGeometry = blockDecomposition.decompose();
    setup.meshGeometry = localGeometry; // over-write global mesh geometry


    auto status = SimulationStatus();
    auto system = FluxConservativeSystem (setup); // This also initializes CT.

    system.setInitialData (setup.initialDataFunction, setup.vectorPotentialFunction);
    writeVtkOutput (setup, status, system);

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
            writeCheckpoint (setup, status, system);
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
        double kzps = 1e-3 * setup.meshGeometry->totalCellsInMesh() / timer.age();
        std::cout << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        std::cout << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        std::cout << "dt=" << std::setprecision (4) << std::scientific << dt << " ";
        std::cout << "kzps=" << std::setprecision (2) << std::fixed << kzps << std::endl;
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
        "Mara is an astrophysics code for gas and magnetofluid "
        "dynamics simulations.\n";
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
