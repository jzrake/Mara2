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

    domainShape = meshGeometry->domainShape();
    numDimensions = Array::vectorFromShape (domainShape).size();
    numConserved = conservationLaw->getNumConserved();
    schemeOrder = intercellFluxScheme->getSchemeOrder();

    // This is the shape of the updateableRegion, with space for conserved
    // quantities on the final axis.
    domainShapeAndState = domainShape;
    domainShapeAndState[numDimensions] = numConserved;

    Shape shapeU = domainShape;
    Shape shapeP = domainShape;
    Shape shapeL = domainShape;
    Shape shapeF1 = domainShape;
    Shape shapeF2 = domainShape;
    Shape shapeF3 = domainShape;

    for (int n = 0; n < numDimensions; ++n)
    {
        if (domainShape[n] > 1) // We tolerate domain shapes such as (1, 1, 128)
        {
            shapeU[n] += 2 * schemeOrder;
            shapeP[n] += 2 * schemeOrder;
            shapeL[n] += 2 * schemeOrder;
            updateableRegion.lower[n] =  schemeOrder;
            updateableRegion.upper[n] = -schemeOrder;
        }
    }

    shapeU[numDimensions] = numConserved;
    shapeP[numDimensions] = numConserved;
    shapeL[numDimensions] = numConserved;
    shapeF1[numDimensions] = numConserved;
    shapeF2[numDimensions] = numConserved;
    shapeF3[numDimensions] = numConserved;

    // Each array of fluxes has one more element along its own coordinate
    // axis.
    shapeF1[0] += 1;
    shapeF2[1] += 1;
    shapeF3[2] += 1;

    U = Array (shapeU);
    P = Array (shapeP);
    L = Array (shapeL);
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

void FluxConservativeSystem::computeIntercellFluxes()
{
    int N1 = domainShape[0];
    int N2 = domainShape[1];
    int N3 = domainShape[2];

    auto stateVector = IntercellFluxScheme::StateVector (schemeOrder * 2);
    auto faceData = IntercellFluxScheme::FaceData();
    auto request = ConservationLaw::Request();



    request.getPrimitive = true;
    request.getConserved = false;
    request.getFluxes = false;
    request.getEigenvalues = false;


    // Flux sweep along axis 1
    // ------------------------------------------------------------------------
    request.areaElement[0] = 1.0;
    request.areaElement[1] = 0.0;
    request.areaElement[2] = 0.0;

    faceData.stencilData = Cow::Array (schemeOrder * 2, numConserved);
    faceData.areaElement = request.areaElement;
    faceData.conservationLaw = conservationLaw;

    for (int i = 0; i < N1 + 1; ++i)
    {
        for (int j = 0; j < N2; ++j)
        {
            for (int k = 0; k < N3; ++k)
            {
                for (int n = 0; n < schemeOrder * 2; ++n)
                {
                    stateVector[n] = conservationLaw->fromPrimitive (request, &P (i + n, j, k));
                    faceData.stencilData (n, 0) = P (i + n, j, k, 0);
                }


                auto flux = intercellFluxScheme->intercellFlux (stateVector, conservationLaw.get());
                auto flux2 = intercellFluxScheme->intercellFlux (faceData);

                for (int q = 0; q < numConserved; ++q)
                {
                    F1 (i, j, k, q) = flux2.F[q];
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

void FluxConservativeSystem::applyBoundaryConditions()
{
    boundaryCondition->apply (P, schemeOrder);
}

void FluxConservativeSystem::updateConserved (double dt, double rungeKuttaParameter)
{
    double b = rungeKuttaParameter;

    for (int n = 0; n < U.size(); ++n)
    {
        U[n] += U[n] * b + L[n] * dt * (1 - b);
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
        double maxWaveSpeed = std::fabs (S.A[0]);
        double minLength = MIN3(dx1, dx2, dx3);

        courantTimestep = MIN2(courantTimestep, minLength / maxWaveSpeed);
    }
    return courantTimestep;
}
