#include <cmath>
#include <cassert>
#include "CellCenteredFieldCT.hpp"
#include "FieldOperator.hpp"
#include "IntercellFluxSchemes.hpp"
#include "MeshData.hpp"
#include "MeshOperator.hpp"
#include "SolutionSchemes.hpp"
#include <cstdio>

using namespace Cow;




// ============================================================================
void SchemeHelpers::makeFootprint (
    int stencilSize,
    Shape3D arrayShape,
    Shape3D boundaryShape,
    Shape3D& footprint,
    Index& startIndex,
    Region& interior)
{
    int ng = stencilSize;

    for (int axis = 0; axis < 3; ++axis)
    {
        if (arrayShape[axis] > 1)
        {
            footprint[axis] = 2 * ng;
            startIndex[axis] = -boundaryShape[axis];
            interior.lower[axis] =  boundaryShape[axis];
            interior.upper[axis] = -boundaryShape[axis];
        }
        else
        {
            footprint[axis] = 0;
            startIndex[axis] = 0;
            interior.lower[axis] = 0;
            interior.upper[axis] = 0;
        }
    }

    if (! boundaryShape.contains (footprint / 2))
    {
        throw std::logic_error ("Boundary region of mesh data is smaller than the scheme's stencil");
    }
}




// ============================================================================
int GenericSolutionScheme::getStencilSize() const
{
    if (! fluxScheme) throw std::logic_error ("No IntercellFluxScheme instance");
    return fluxScheme->getStencilSize();
}

void GenericSolutionScheme::setIntercellFluxScheme (std::shared_ptr<IntercellFluxScheme> fs)
{
    fluxScheme = fs;
}




// ============================================================================
MethodOfLinesTVD::MethodOfLinesTVD()
{
    fluxScheme = std::make_shared<MethodOfLines>();
    rungeKuttaOrder = 1;
    disableFieldCT = false;
}

void MethodOfLinesTVD::setRungeKuttaOrder (int rungeKuttaOrderToUse)
{
    if (0 >= rungeKuttaOrderToUse || rungeKuttaOrderToUse > 3)
    {
        throw std::logic_error ("Runge-Kutta order not 1, 2, or 3");
    }
    rungeKuttaOrder = rungeKuttaOrderToUse;
}

void MethodOfLinesTVD::setDisableFieldCT (bool shouldDisableFieldCT)
{
    disableFieldCT = shouldDisableFieldCT;
}

void MethodOfLinesTVD::setViscousFluxFunction (std::function<void(const Cow::Array&, Cow::Array&, double)> viscousFluxToUse)
{
    viscousFlux = viscousFluxToUse;
}

void MethodOfLinesTVD::setStarParticleDerivatives (
    std::function<std::vector<double>(const Cow::Array&, const std::vector<double>&)>
    starParticleDerivativesToUse)
{
    starParticleDerivatives = starParticleDerivativesToUse;
}

void MethodOfLinesTVD::setSourceTermsWithParticles (SourceTermsWithParticles sourceTermsWithParticlesToUse)
{
    sourceTermsWithParticles = sourceTermsWithParticlesToUse;
}

void MethodOfLinesTVD::advance (MeshData& solution, double t0, double dt) const
{
    check_valid();

    // Figure out the scheme footprint, and if we have enough guard zones
    // ------------------------------------------------------------------------
    auto footprint = Shape3D();
    auto startIndex = Index();
    auto interior = Region();

    SchemeHelpers::makeFootprint (
        fluxScheme->getStencilSize(),
        solution.P.shape(),
        solution.getBoundaryShape(),
        footprint, startIndex, interior);


    // Setup callback to compute Godunov fluxes
    // ------------------------------------------------------------------------
    auto cl = fieldOperator->getConservationLaw();
    auto nq = cl->getNumConserved();

    auto D = IntercellFluxScheme::FaceData();
    D.conservationLaw = cl;


    // RK parameters
    // ------------------------------------------------------------------------
    double b1[1] = {1.0};
    double b2[2] = {1.0, 0.5};
    double b3[3] = {1.0, 1./4, 2./3};
    double *b = nullptr;

    switch (rungeKuttaOrder)
    {
        case 1: b = b1; break;
        case 2: b = b2; break;
        case 3: b = b3; break;
    }


    // RK updates
    // ------------------------------------------------------------------------
    auto U0 = fieldOperator->generateConserved (solution.P, t0);
    auto H0 = solution.starParticles;
    auto U = U0;
    auto t = t0;
    auto H = H0; // H is for 'holes'

    auto Fhat = [&] (GodunovStencil& stencil)
    {
        D.areaElement = stencil.faceNormal.cartesian();
        D.stencilData = stencil.cellData;

        auto S = fluxScheme->intercellFlux (D, t);

        for (int q = 0; q < nq; ++q)
        {
            stencil.godunovFlux[q] = S.F[q];
        }
    };


    auto fluxCorrection = [&] (const Array& P, Array& fluxes)
    {
        auto ib = cl->getIndexFor (ConservationLaw::VariableType::magnetic);

        if (ib != -1 && ! disableFieldCT)
        {
            auto ct = CellCenteredFieldCT();
            ct.correctGodunovFluxes (fluxes, ib);
        }

        if (viscousFlux)
        {
            viscousFlux (P, fluxes, t);
        }
    };



    for (int rk = 0; rk < rungeKuttaOrder; ++rk)
    {
        auto F = meshOperator->godunov (Fhat, solution.P, solution.B, footprint, startIndex, fluxCorrection);
        auto L = meshOperator->divergence (F, -1.0, startIndex);

        cl->addSourceTerms (solution.P, L);

        if (sourceTermsFunction)
        {
            auto S = meshOperator->generateSourceTerms (sourceTermsFunction, solution.P, t, startIndex);

            for (int n = 0; n < L.size(); ++n)
            {
                L[n] += S[n];
            }
        }
        else if (sourceTermsWithParticles)
        {
            auto S = meshOperator->generateSourceTerms (sourceTermsWithParticles, solution.P, t, H, startIndex);

            for (int n = 0; n < L.size(); ++n)
            {
                L[n] += S[n];
            }
        }

        for (int n = 0; n < L.size(); ++n)
        {
            U[n] = U0[n] * (1 - b[rk]) + (U[n] + dt * L[n]) * b[rk];
        }

        // Update the inter-timestep time
        t = t0 * (1 - b[rk]) + (t + dt) * b[rk];

        // Update the star particle positions
        if (starParticleDerivatives)
        {
            auto Hdot = starParticleDerivatives (solution.P, H);

            assert(Hdot.size() == H.size());

            for (std::size_t n = 0; n < H.size(); ++n)
            {
                H[n] = H0[n] * (1 - b[rk]) + (H[n] + dt * Hdot[n]) * b[rk];
            }
        }

        fieldOperator->recoverPrimitive (U[interior], solution.P[interior], t);
        solution.starParticles = H;
        solution.applyBoundaryCondition (*boundaryCondition);
    }
}

Cow::Array MethodOfLinesTVD::computeAdvectiveFluxes (MeshData& solution, double t0) const
{
    check_valid();

    // Figure out the scheme footprint, and if we have enough guard zones
    // ------------------------------------------------------------------------
    auto footprint = Shape3D();
    auto startIndex = Index();
    auto interior = Region();

    SchemeHelpers::makeFootprint (
        fluxScheme->getStencilSize(),
        solution.P.shape(),
        solution.getBoundaryShape(),
        footprint, startIndex, interior);


    // Setup callback to compute Godunov fluxes
    // ------------------------------------------------------------------------
    auto cl = fieldOperator->getConservationLaw();
    auto nq = cl->getNumConserved();

    auto D = IntercellFluxScheme::FaceData();
    D.conservationLaw = cl;

    auto Fhat = [&] (GodunovStencil& stencil)
    {
        D.areaElement = stencil.faceNormal.cartesian();
        D.stencilData = stencil.cellData;

        auto S = fluxScheme->intercellFlux (D, t0);

        for (int q = 0; q < nq; ++q)
        {
            stencil.godunovFlux[q] = S.F[q];
        }
    };
    return meshOperator->godunov (Fhat, solution.P, solution.B, footprint, startIndex, nullptr);
}

Cow::Array MethodOfLinesTVD::computeViscousFluxes (MeshData& solution, double t0) const
{
    check_valid();

    auto footprint = Shape3D();
    auto startIndex = Index();
    auto interior = Region();

    SchemeHelpers::makeFootprint (
        fluxScheme->getStencilSize(),
        solution.P.shape(),
        solution.getBoundaryShape(),
        footprint, startIndex, interior);

    auto Fhat = [] (GodunovStencil&) {};
    auto visc = [this, t0] (const Array& P, Array& F) { if (viscousFlux) viscousFlux (P, F, t0); };

    return meshOperator->godunov (Fhat, solution.P, solution.B, footprint, startIndex, visc);
}

void MethodOfLinesTVD::check_valid() const
{
    if (! fieldOperator)     throw std::logic_error ("No FieldOperator instance");
    if (! meshOperator)      throw std::logic_error ("No MeshOperator instance");
    if (! boundaryCondition) throw std::logic_error ("No BoundaryCondition instance");
    if (! fluxScheme)        throw std::logic_error ("No IntercellFluxScheme instance");
}
