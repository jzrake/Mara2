#include <iostream> //DEBUG
#include <cassert>
#include <cmath>
#include "FluxConservativeSystem.hpp"
#define MIN2(a, b) ((a) < (b) ? a : b)
#define MIN3(a, b, c) ((a) < (b) ? MIN2(a, c) : MIN2(b, c))
using namespace Cow;




// ============================================================================
ConservationLaw::Request::Request()
{
    getPrimitive = false;
    getConserved = false;
    getFluxes = false;
    getEigenvalues = false;
    areaElement[0] = 1.0;
    areaElement[1] = 0.0;
    areaElement[2] = 0.0;
}




// ============================================================================
FluxConservativeSystem::FluxConservativeSystem (SimulationSetup setup)
{
    meshGeometry           = setup.meshGeometry;
    conservationLaw        = setup.conservationLaw;
    intercellFluxScheme    = setup.intercellFluxScheme;
    boundaryCondition      = setup.boundaryCondition;
    rungeKuttaOrder        = setup.rungeKuttaOrder;

    domainShape   = meshGeometry->domainShape();
    numConserved  = conservationLaw->getNumConserved();
    stencilSize   = intercellFluxScheme->getStencilSize();

    Shape shapeU  = domainShape;
    Shape shapeP  = domainShape;
    Shape shapeL  = domainShape;
    Shape shapeF1 = domainShape;
    Shape shapeF2 = domainShape;
    Shape shapeF3 = domainShape;

    shapeU [3] = numConserved;
    shapeP [3] = numConserved;
    shapeL [3] = numConserved;
    shapeF1[3] = numConserved;
    shapeF2[3] = numConserved;
    shapeF3[3] = numConserved;

    for (int n = 0; n < 3; ++n)
    {
        if (domainShape[n] > 1) // We accommodate domain shapes such as (1, 1, 128)
        {
            shapeU[n] += 2 * stencilSize;
            shapeP[n] += 2 * stencilSize;
            shapeL[n] += 2 * stencilSize;
            updateableRegion.lower[n] =  stencilSize;
            updateableRegion.upper[n] = -stencilSize;
        }
    }

    // Each array of fluxes has one more element along its own coordinate
    // axis.
    shapeF1[0] += 1;
    shapeF2[1] += 1;
    shapeF3[2] += 1;

    U  = Array (shapeU);
    P  = Array (shapeP);
    L  = Array (shapeL);
    U0 = Array (shapeU);
    F1 = Array (shapeF1);
    F2 = Array (shapeF2);
    F3 = Array (shapeF3);
}

Cow::Array::Reference FluxConservativeSystem::getPrimitive()
{
    return P[updateableRegion];
}

void FluxConservativeSystem::setInitialData (InitialDataFunction F)
{
    auto request = ConservationLaw::Request();
    auto Preg = P[updateableRegion];
    auto Ureg = U[updateableRegion];
    auto pit = Preg.begin();
    auto uit = Ureg.begin();

    for ( ; uit != Ureg.end(); ++uit, ++pit)
    {
        auto index = pit.relativeIndex();
        auto coord = meshGeometry->coordinateAtIndex (index[0], index[1], index[2]);
        auto P = F (coord[0], coord[1], coord[2]);
        auto S = conservationLaw->fromPrimitive (request, &P[0]);

        for (int q = 0; q < numConserved; ++q)
        {
            pit[q] = S.P[q];
            uit[q] = S.U[q];
        }
    }
}

double FluxConservativeSystem::getCourantTimestep()
{
    auto request = ConservationLaw::Request();
    auto Preg = P[updateableRegion];

    double courantTimestep = std::numeric_limits<double>::max();

    for (auto pit = Preg.begin(); pit != Preg.end(); ++pit)
    {
        const auto index = pit.relativeIndex();
        const int i = index[0];
        const int j = index[1];
        const int k = index[2];

        const double dx1 = meshGeometry->cellLength (i, j, k, 0);
        const double dx2 = meshGeometry->cellLength (i, j, k, 1);
        const double dx3 = meshGeometry->cellLength (i, j, k, 2);

        auto S = conservationLaw->fromPrimitive (request, pit);
        double maxWaveSpeed = conservationLaw->maxEigenvalueMagnitude (S);
        double minLength = MIN3(dx1, dx2, dx3);

        courantTimestep = MIN2(courantTimestep, minLength / maxWaveSpeed);
    }
    return courantTimestep;
}

void FluxConservativeSystem::advance (double dt)
{
    cacheConserved();

    switch (rungeKuttaOrder)
    {
        case 1:
        {
            takeRungeKuttaSubstep (dt, 1.);
            break;
        }
        case 2:
        {
            takeRungeKuttaSubstep (dt, 1.);
            takeRungeKuttaSubstep (dt, 1./2);
            break;
        }
        case 3:
        {
            takeRungeKuttaSubstep (dt, 1.);
            takeRungeKuttaSubstep (dt, 1./4);
            takeRungeKuttaSubstep (dt, 2./3);
            break;
        }
    }
}

void FluxConservativeSystem::computeIntercellFluxes()
{
    int N1 = domainShape[0];
    int N2 = domainShape[1];
    int N3 = domainShape[2];

    // Flux sweep along axis 1
    // ------------------------------------------------------------------------
    auto faceData = IntercellFluxScheme::FaceData();
    faceData.stencilData = Cow::Array (stencilSize * 2, numConserved);
    faceData.areaElement[0] = 1.0;
    faceData.areaElement[1] = 0.0;
    faceData.areaElement[2] = 0.0;
    faceData.conservationLaw = conservationLaw;

    for (int i = 0; i < N1 + 1; ++i)
    {
        for (int j = 0; j < N2; ++j)
        {
            for (int k = 0; k < N3; ++k)
            {
                for (int q = 0; q < numConserved; ++q)
                {
                    for (int n = 0; n < stencilSize * 2; ++n)
                    {
                        faceData.stencilData (n, q) = P (i + n, j, k, q);
                    }
                }

                auto flux = intercellFluxScheme->intercellFlux (faceData);

                for (int q = 0; q < numConserved; ++q)
                {
                    F1 (i, j, k, q) = flux.F[q];
                }
            }
        }
    }
}

void FluxConservativeSystem::computeTimeDerivative()
{
    auto reg = L[updateableRegion];

    for (auto it = reg.begin(); it != reg.end(); ++it)
    {
        const auto index = it.relativeIndex();
        const int i = index[0];
        const int j = index[1];
        const int k = index[2];

        const double A1L = meshGeometry->faceArea (i - 0, j, k, 1);
        const double A1R = meshGeometry->faceArea (i + 1, j, k, 1);
        const double A2L = meshGeometry->faceArea (i, j - 0, k, 2);
        const double A2R = meshGeometry->faceArea (i, j + 1, k, 2);
        const double A3L = meshGeometry->faceArea (i, j, k - 0, 3);
        const double A3R = meshGeometry->faceArea (i, j, k + 1, 3);
        const double Vol = meshGeometry->cellVolume (i, j, k);

        for (int q = 0; q < numConserved; ++q)
        {
            const double dAF1 = F1 (i + 1, j, k, q) * A1R - F1 (i, j, k, q) * A1L;
            const double dAF2 = F2 (i, j + 1, k, q) * A2R - F2 (i, j, k, q) * A2L;
            const double dAF3 = F3 (i, j, k + 1, q) * A3R - F3 (i, j, k, q) * A3L;

            it[q] = -(dAF1 + dAF2 + dAF3) / Vol;
        }
    }
}

void FluxConservativeSystem::applyBoundaryCondition()
{
    boundaryCondition->apply (P, stencilSize);
}

void FluxConservativeSystem::cacheConserved()
{
    for (int n = 0; n < U.size(); ++n)
    {
        U0[n] = U[n];
    }
}

void FluxConservativeSystem::updateConserved (double dt)
{
    for (int n = 0; n < U.size(); ++n)
    {
        U[n] += dt * L[n];
    }
}

void FluxConservativeSystem::averageRungeKutta (double b)
{
    for (int n = 0; n < U.size(); ++n)
    {
        U[n] = U0[n] * (1 - b) + U[n] * b;
    }
}

void FluxConservativeSystem::recoverPrimitive()
{
    auto request = ConservationLaw::Request();
    auto Ureg = U[updateableRegion];
    auto Preg = P[updateableRegion];
    auto pit = Preg.begin();
    auto uit = Ureg.begin();

    for ( ; uit != Ureg.end(); ++uit, ++pit)
    {
        auto S = conservationLaw->fromConserved (request, uit);

        for (int q = 0; q < numConserved; ++q)
        {
            pit[q] = S.P[q];
        }
    }
}

void FluxConservativeSystem::takeRungeKuttaSubstep (double dt, double b)
{
    computeIntercellFluxes();
    computeTimeDerivative();
    updateConserved (dt);
    averageRungeKutta (b);
    recoverPrimitive();
    applyBoundaryCondition();
}
