#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "Mara.hpp"
#include "Configuration.hpp"
#include "FileSystem.hpp"
#include "FluxConservativeSystem.hpp"
#include "RiemannSolver.hpp"
#include "HDF5.hpp"
#include "VTK.hpp"
#include "Matrix.hpp"
#include "DebugHelper.hpp"




// ============================================================================
SimulationSetup::SimulationSetup()
{
    finalTime = 1.0;
    checkpointInterval = 1.0;
    vtkOutputInterval = 1.0;
    cflParameter = 0.25;
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
void PeriodicBoundaryCondition::apply (Cow::Array& P, int numGuard) const
{
    int ng = numGuard;
    int ni = P.size(0) - 2 * ng;
    int nj = P.size(1) - 2 * ng;
    int nk = P.size(2) - 2 * ng;
    int nq = P.size(3);

    if (P.size(0) > 1)
    {
        for (int j = 0; j < P.size(1); ++j)
        {
            for (int k = 0; k < P.size(2); ++k)
            {
                for (int i = 0; i < ng; ++i)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i,           j, k, q) = P (i + ni, j, k, q);
                        P (i + ng + ni, j, k, q) = P (i + ng, j, k, q);
                    }
                }
            }
        }
    }

    if (P.size(1) > 1)
    {
        for (int k = 0; k < P.size(2); ++k)
        {
            for (int i = 0; i < P.size(0); ++i)
            {
                for (int j = 0; j < ng; ++j)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i, j,           k, q) = P (i, j + nj, k, q);
                        P (i, j + ng + nj, k, q) = P (i, j + ng, k, q);
                    }
                }
            }
        }
    }

    if (P.size(2) > 1)
    {
        for (int i = 0; i < P.size(0); ++i)
        {
            for (int j = 0; j < P.size(1); ++j)
            {
                for (int k = 0; k < ng; ++k)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i, j, k,           q) = P (i, j, k + nk, q);
                        P (i, j, k + ng + nk, q) = P (i, j, k + ng, q);
                    }
                }
            }
        }
    }
}




// ============================================================================
ConservationLaw::State ConservationLaw::averageStates (const Request& request,
    const State& L, const State& R) const
{
    int nq = getNumConserved();
    auto Paverage = std::vector<double> (nq);

    for (int q = 0; q < nq; ++q)
    {
        Paverage[q] = 0.5 * (L.P[q] + R.P[q]);
    }
    return fromPrimitive (request, &Paverage[0]);
}

std::vector<ConservationLaw::State> ConservationLaw::fromPrimitive
(const Request& request, const Cow::Array& P) const
{
    assert (P.size (1) == getNumConserved());
    auto states = std::vector<ConservationLaw::State>(P.size (0));

    for (int n = 0; n < P.size (0); ++n)
    {
        states[n] = fromPrimitive (request, &P(n));
    }
    return states;
}

double ConservationLaw::maxEigenvalueMagnitude (const State& state) const
{
    double maxLambda = 0.0;
    int nq = getNumConserved();

    for (int n = 0; n < nq; ++n)
    {
        if (maxLambda < std::fabs (state.A[n]))
        {
            maxLambda = std::fabs (state.A[n]);
        }
    }
    return maxLambda;
}

double ConservationLaw::maxEigenvalueMagnitude (const StateVector& states) const
{
    double maxLambda = 0.0;

    for (int n = 0; n < states.size(); ++n)
    {
        double A = maxEigenvalueMagnitude (states[n]);

        if (maxLambda < A) maxLambda = A;
    }
    return maxLambda;
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
    S.F.resize (nq);

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

void writeOutput (SimulationSetup& setup, SimulationStatus& status, FluxConservativeSystem& system)
{
    using namespace Cow;

    auto dir = setup.outputDirectory;

    FileSystem::ensureDirectoryExists (dir);
    auto h5fname = FileSystem::makeFilename (dir, "chkpt", ".h5", status.checkpointsWrittenSoFar);
    auto vtkname = FileSystem::makeFilename (dir, "mesh", ".vtk", status.vtkOutputsWrittenSoFar);

    auto file = H5::File (h5fname, "w");
    auto vtkStream = std::ofstream (vtkname);
    auto vtkDataSet = VTK::DataSet (setup.meshGeometry->domainShape());
    vtkDataSet.setTitle (setup.runName);

    for (int q = 0; q < setup.conservationLaw->getNumConserved(); ++q)
    {
        std::string field = setup.conservationLaw->getPrimitiveName(q);

        auto P = system.getPrimitive (q);
        file.write (field, P);

        vtkDataSet.addField (field, P);
    }
    vtkDataSet.write (vtkStream);

    ++status.checkpointsWrittenSoFar;
    ++status.vtkOutputsWrittenSoFar;
}




// ============================================================================
int main(int argc, const char* argv[])
{
    using namespace Cow;
    std::set_terminate (Cow::terminateWithBacktrace);


    if (argc == 1)
    {
        std::cout << "usage: mara config.lua\n";
        return 0;
    }
    auto configuration = Configuration();
    auto setup = configuration.fromLuaFile (argv[1]);


    // Setup lines specific to conservation law problems
    auto system = FluxConservativeSystem (setup);
    system.setInitialData (setup.initialDataFunction);


    auto status = SimulationStatus();
    writeOutput (setup, status, system);

    while (status.simulationTime < setup.finalTime)
    {
        if (status.simulationTime >= status.checkpointsWrittenSoFar * setup.checkpointInterval)
        {
            writeOutput (setup, status, system);
        }

        double dt = setup.cflParameter * system.getCourantTimestep();
        system.advance (dt);
        status.simulationTime += dt;
        status.simulationIter += 1;

        std::cout << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        std::cout << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        std::cout << "dt=" << std::setprecision (2) << std::scientific << dt << "\n";
    }

    return 0;
}
