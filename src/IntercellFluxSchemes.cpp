#include <cassert>
#include "IntercellFluxSchemes.hpp"
#include "RiemannSolvers.hpp"




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
    riemannSolver = std::make_shared<HlleRiemannSolver>();
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
MethodOfLinesPlm::MethodOfLinesPlm()
{
    plm.setPlmTheta (1.5);
    riemannSolver = std::make_shared<HlleRiemannSolver>();
}

void MethodOfLinesPlm::setPlmTheta (double plmTheta)
{
    plm.setPlmTheta (plmTheta);    
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

    assert (nq <= 8);

    std::array<double, 8> PL;
    std::array<double, 8> PR;

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
MethodOfLinesWeno::MethodOfLinesWeno()
{
    weno.setSmoothnessIndicator (Reconstruction::ImprovedShenZha10);
    weno.setShenZha10A (50.0);
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

