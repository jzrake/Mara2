#include <iostream>
#include <iomanip>
#include <cmath>
#include "Mara.hpp"
#include "Configuration.hpp"
#include "FluxConservativeSystem.hpp"
#include "HDF5.hpp"




// ============================================================================
SimulationSetup::SimulationSetup()
{
    finalTime = 1.0;
    checkpointInterval = 1.0;
    cflParameter = 0.25;
}




// ============================================================================
SimulationStatus::SimulationStatus()
{
    simulationTime = 0.0;
    simulationIter = 0;
    outputsWrittenSoFar = 0;
}




// ============================================================================
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
MethodOfLines::MethodOfLines (double plmTheta)
{
    plm.setPlmTheta (plmTheta);
}

ConservationLaw::State MethodOfLines::intercellFlux (const FaceData& faceData) const
{
    ConservationLaw::Request request;
    request.areaElement = faceData.areaElement;

    const double* P0 = &faceData.stencilData (0);
    const double* P1 = &faceData.stencilData (1);
    const double* P2 = &faceData.stencilData (2);
    const double* P3 = &faceData.stencilData (3);
    const int nq = faceData.conservationLaw->getNumConserved();

    std::vector<double> PL (nq);
    std::vector<double> PR (nq);

    for (int q = 0; q < nq; ++q)
    {
        const double ps[4] = {P0[q], P1[q], P2[q], P3[q]};
        PL[q] = plm.reconstruct (&ps[1], Reconstruction::PLM_C2R);
        PR[q] = plm.reconstruct (&ps[2], Reconstruction::PLM_C2L);
    }
    auto L = faceData.conservationLaw->fromPrimitive (request, &PL[0]);
    auto R = faceData.conservationLaw->fromPrimitive (request, &PR[0]);

    UpwindRiemannSolver riemannSolver;
    return riemannSolver.solve (L, R, faceData.areaElement);
}

int MethodOfLines::getStencilSize() const
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
    auto claw = faceData.conservationLaw;
    auto request = ConservationLaw::Request();
    request.areaElement = faceData.areaElement;

    const auto states = claw->fromPrimitive (request, faceData.stencilData);
    const auto maxLam = claw->maxEigenvalueMagnitude (states);

    claw->maxEigenvalueMagnitude (states);

    double Fp[6];
    double Fm[6];

    for (int n = 0; n < 6; ++n)
    {
        const auto& S = states[n];
        Fp[n] = S.F[0] + maxLam * S.U[0]; // Lax-Friedrichs flux splitting
        Fm[n] = S.F[0] - maxLam * S.U[0];
    }

    double Fhatp = weno.reconstruct (Fp + 2, Reconstruction::WENO5_FD_C2R);
    double Fhatm = weno.reconstruct (Fm + 3, Reconstruction::WENO5_FD_C2L);
    auto S = ConservationLaw::State();
    S.F = {0.5 * (Fhatp + Fhatm)};
    return S;
}

int MethodOfLinesWeno::getStencilSize() const
{
    return 3;
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
    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0000.h5", "w");
        file.write ("primitive", P);
    }


    while (status.simulationTime < setup.finalTime)
    {
        double dt = setup.cflParameter * system.getCourantTimestep();
        system.advance (dt);
        status.simulationTime += dt;
        status.simulationIter += 1;

        std::cout << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        std::cout << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        std::cout << "dt=" << std::setprecision (2) << std::scientific << dt << "\n";
    }

    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0001.h5", "w");
        file.write ("primitive", P);
        file.write ("t", status.simulationTime);
    }

    return 0;
}
