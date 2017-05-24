#include <iostream> //DEBUG
#include <cassert>
#include <sstream>
#include <cmath>
#include "FluxConservativeSystem.hpp"
#define MIN2(a, b) ((a) < (b) ? a : b)
#define MIN3(a, b, c) ((a) < (b) ? MIN2(a, c) : MIN2(b, c))
using namespace Cow;




// ============================================================================
const char* FluxConservativeSystem::SolverFailure::what() const noexcept
{
    std::stringstream stream;
    static std::string message;

    for (auto& e : failedStates)
    {
        stream << e.what() << std::endl;
    }
    message = stream.str();

    return message.c_str();
}




// ============================================================================
FluxConservativeSystem::FluxConservativeSystem (SimulationSetup setup)
{
    meshGeometry           = setup.meshGeometry;
    conservationLaw        = setup.conservationLaw;
    intercellFluxScheme    = setup.intercellFluxScheme;
    constrainedTransport   = setup.constrainedTransport;
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

    updateableRegion.stride[3] = numConserved;

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

    // Initialize CT.
    auto ct = getCT();
    ct->setMeshGeometry (meshGeometry);
    ct->setBoundaryCondition (boundaryCondition);

    int imag = conservationLaw->getIndexFor (ConservationLaw::VariableType::magnetic);

    if (imag != -1)
    {
        magneticIndices.lower[3] = imag;
        magneticIndices.upper[3] = imag + 3;
    }
}

Cow::Array::Reference FluxConservativeSystem::getPrimitive (int fieldIndex)
{
    Region R = updateableRegion;

    R.stride[3] = 1;

    if (fieldIndex != -1)
    {
        assert (0 <= fieldIndex && fieldIndex < numConserved);
        R.lower[3] = fieldIndex;
        R.upper[3] = fieldIndex + 1;
    }
    return P[R];
}

Cow::Array::Reference FluxConservativeSystem::getPrimitiveVector (int fieldIndex)
{
    Region R = updateableRegion;
    R.stride[3] = 1;
    R.lower[3] = fieldIndex;
    R.upper[3] = fieldIndex + 3;
    return P[R];
}

void FluxConservativeSystem::setInitialData (InitialDataFunction F, InitialDataFunction A)
{
    auto request = ConservationLaw::Request();
    auto Preg = P[updateableRegion];
    auto Ureg = U[updateableRegion];
    auto pit = Preg.begin();
    auto uit = Ureg.begin();
    int imag = conservationLaw->getIndexFor (ConservationLaw::VariableType::magnetic);
    bool useVectorPotential = A != nullptr && imag != -1;

    // If there is a vector potential function, then cache the resulting
    // magnetic field in the U array (cell centers).
    if (useVectorPotential)
    {
        auto ct = getCT();

        ct->assignVectorPotential (A, ConstrainedTransport::MeshLocation::face);

        auto ctFluxes = ct->computeGodunovFluxesFieldCT();
        //auto ctFluxes = ct->getGodunovFluxes();
        F1[magneticIndices] = ctFluxes.F1;
        F2[magneticIndices] = ctFluxes.F2;
        F3[magneticIndices] = ctFluxes.F3;

        computeTimeDerivative();
        updateConserved (1.0);
    }

    for ( ; uit != Ureg.end(); ++uit, ++pit)
    {
        auto index = pit.relativeIndex();
        auto coord = meshGeometry->coordinateAtIndex (index[0], index[1], index[2]);
        auto P = F (coord[0], coord[1], coord[2]);

        if (P.size() != numConserved)
        {
            throw std::runtime_error ("initial data function returned vector of length != nq");
        }

        if (useVectorPotential) // then fields are already computed and stored in U
        {
            P[imag + 0] = uit[imag + 0];
            P[imag + 1] = uit[imag + 1];
            P[imag + 2] = uit[imag + 2];
        }

        auto S = conservationLaw->fromPrimitive (request, &P[0]);

        for (int q = 0; q < numConserved; ++q)
        {
            pit[q] = S.P[q];
            uit[q] = S.U[q];
        }
    }

    applyBoundaryCondition();
    uploadFieldsToCT();
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

        // !!! Need to compute wave speeds in other directions
        auto S = conservationLaw->fromPrimitive (request, pit);
        double maxWaveSpeed = conservationLaw->maxEigenvalueMagnitude(S);
        double minLength = MIN3(dx1, dx2, dx3);

        courantTimestep = MIN2(courantTimestep, minLength / maxWaveSpeed);
    }
    return courantTimestep;
}

void FluxConservativeSystem::advance (double dt)
{
    U0 = U; // cache conserved

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

    uploadFieldsToCT();
}

void FluxConservativeSystem::computeIntercellFluxes()
{
    if (domainShape[0] > 1) intercellFluxSweep(0);
    if (domainShape[1] > 1) intercellFluxSweep(1);
    if (domainShape[2] > 1) intercellFluxSweep(2);

    int imag = conservationLaw->getIndexFor (ConservationLaw::VariableType::magnetic);

    if (imag != -1)
    {
        auto ct = getCT();
        ct->assignGodunovFluxes (F1[magneticIndices], F2[magneticIndices], F3[magneticIndices]);
        auto ctFluxes = ct->computeGodunovFluxesFieldCT();

        F1[magneticIndices] = ctFluxes.F1;
        F2[magneticIndices] = ctFluxes.F2;
        F3[magneticIndices] = ctFluxes.F3;
    }
}

void FluxConservativeSystem::intercellFluxSweep (int axis)
{
    int N1 = domainShape[0];
    int N2 = domainShape[1];
    int N3 = domainShape[2];
    int i0 = domainShape[0] > 1 ? stencilSize : 0; // account for guard zones in P
    int j0 = domainShape[1] > 1 ? stencilSize : 0;
    int k0 = domainShape[2] > 1 ? stencilSize : 0;

    auto faceData = IntercellFluxScheme::FaceData();
    faceData.stencilData = Cow::Array (stencilSize * 2, numConserved);
    faceData.areaElement[0] = axis == 0 ? 1.0 : 0.0;
    faceData.areaElement[1] = axis == 1 ? 1.0 : 0.0;
    faceData.areaElement[2] = axis == 2 ? 1.0 : 0.0;
    faceData.conservationLaw = conservationLaw;

    for (int i = 0; i < N1 + (axis == 0); ++i)
    {
        for (int j = 0; j < N2 + (axis == 1); ++j)
        {
            for (int k = 0; k < N3 + (axis == 2); ++k)
            {
                for (int q = 0; q < numConserved; ++q)
                {
                    for (int n = 0; n < stencilSize * 2; ++n)
                    {
                        switch (axis)
                        {
                            case 0: faceData.stencilData (n, q) = P (i + n,  j + j0, k + k0, q); break;
                            case 1: faceData.stencilData (n, q) = P (i + i0, j + n,  k + k0, q); break;
                            case 2: faceData.stencilData (n, q) = P (i + i0, j + j0, k + n,  q); break;
                        }
                    }
                }

                auto flux = intercellFluxScheme->intercellFlux (faceData);

                for (int q = 0; q < numConserved; ++q)
                {
                    switch (axis)
                    {
                        case 0: F1 (i, j, k, q) = flux.F[q]; break;
                        case 1: F2 (i, j, k, q) = flux.F[q]; break;
                        case 2: F3 (i, j, k, q) = flux.F[q]; break;
                    }
                }
            }
        }
    }



    // This is a possible revision to the method of flux sweeps, but requires
    // some new Array methods:
    // ------------------------------------------------------------------------

    // auto fluxableRegion = Region::strided (3, numConserved);
    // auto Freg = F1[updateableRegion];
    // auto Fit = Freg.begin();

    // for ( ; Fit != Freg.end(); ++Fit)
    // {
    //     auto P0 = P.getIterator (Fit);

    //     for (int q = 0; q < numConserved; ++q)
    //     {
    //         for (int n = 0; n < stencilSize * 2; ++n)
    //         {
    //             faceData.stencilData (n, q) = P0 (n, 0, 0, q);
    //         }
    //     }
    //     auto flux = intercellFluxScheme->intercellFlux (faceData);

    //     for (int q = 0; q < numConserved; ++q)
    //     {
    //         Fit[q] = flux.F[q];
    //     }
    // }
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

        const double A1L = meshGeometry->faceArea (i - 0, j, k, 0);
        const double A1R = meshGeometry->faceArea (i + 1, j, k, 0);
        const double A2L = meshGeometry->faceArea (i, j - 0, k, 1);
        const double A2R = meshGeometry->faceArea (i, j + 1, k, 1);
        const double A3L = meshGeometry->faceArea (i, j, k - 0, 2);
        const double A3R = meshGeometry->faceArea (i, j, k + 1, 2);
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
    boundaryCondition->apply (P, *conservationLaw, stencilSize);
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
    auto solverFailure = SolverFailure();

    for ( ; uit != Ureg.end(); ++uit, ++pit)
    {
        try
        {
            auto S = conservationLaw->fromConserved (request, uit);

            for (int q = 0; q < numConserved; ++q)
            {
                pit[q] = S.P[q];
            }
        }
        catch (ConservationLaw::StateFailure& stateFailure)
        {
            stateFailure.zoneIndex = uit.index();
            solverFailure.failedStates.push_back (stateFailure);
        }
    }

    if (solverFailure.failedStates.size() > 0)
    {
        throw solverFailure;
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

void FluxConservativeSystem::uploadFieldsToCT()
{
    int imag = conservationLaw->getIndexFor (ConservationLaw::VariableType::magnetic);

    if (imag != -1)
    {
        auto ct = getCT();
        ct->assignCellCenteredB (getPrimitiveVector (imag));
    }
}

UniformCartesianCT* FluxConservativeSystem::getCT()
{
    if (auto ct = dynamic_cast<UniformCartesianCT*> (constrainedTransport.get()))
    {
        return ct;
    }

    throw std::logic_error ("constrainedTransport must be UniformCartesianCT");
}
