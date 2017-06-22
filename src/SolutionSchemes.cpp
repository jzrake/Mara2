#include "BoundaryConditions.hpp"
#include "CartesianMeshGeometry.hpp"
#include "CellCenteredFieldCT.hpp"
#include "ConservationLaws.hpp"
#include "FieldOperator.hpp"
#include "IntercellFluxSchemes.hpp"
#include "MeshData.hpp"
#include "MeshOperator.hpp"
#include "SolutionSchemes.hpp"
#include "RiemannSolvers.hpp"

using namespace Cow;




// ============================================================================
MethodOfLinesTVD::MethodOfLinesTVD()
{
    fluxScheme = std::make_shared<MethodOfLines>();
    rungeKuttaOrder = 1;
    disableFieldCT = false;
}

void MethodOfLinesTVD::setDisableFieldCT (bool shouldDisableFieldCT)
{
    disableFieldCT = shouldDisableFieldCT;
}

int MethodOfLinesTVD::getStencilSize() const
{
    if (! fluxScheme) throw std::logic_error ("No IntercellFluxScheme instance");
    return fluxScheme->getStencilSize();
}

void MethodOfLinesTVD::setRungeKuttaOrder (int rungeKuttaOrderToUse)
{
    if (0 >= rungeKuttaOrderToUse || rungeKuttaOrderToUse > 3)
    {
        throw std::logic_error ("Runge-Kutta order not 1, 2, or 3");
    }
    rungeKuttaOrder = rungeKuttaOrderToUse;
}

void MethodOfLinesTVD::setIntercellFluxScheme (std::shared_ptr<IntercellFluxScheme> fs)
{
    fluxScheme = fs;
}

void MethodOfLinesTVD::advance (MeshData& solution, double dt) const
{
    if (! fieldOperator)     throw std::logic_error ("No FieldOperator instance");
    if (! meshOperator)      throw std::logic_error ("No MeshOperator instance");
    if (! boundaryCondition) throw std::logic_error ("No BoundaryCondition instance");
    if (! fluxScheme)        throw std::logic_error ("No IntercellFluxScheme instance");


    // Figure out the scheme footprint, and if we have enough guard zones
    // ------------------------------------------------------------------------
    auto footprint = Shape3D();
    auto startIndex = Index();
    auto interior = Region();
    makeFootprint (solution, *fluxScheme, footprint, startIndex, interior);


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

        auto S = fluxScheme->intercellFlux (D);

        for (int q = 0; q < nq; ++q)
        {
            stencil.godunovFlux[q] = S.F[q];
        }
    };

    auto fieldCT = [&] (Array& fluxes)
    {
        auto ib = cl->getIndexFor (ConservationLaw::VariableType::magnetic);

        if (ib != -1 && ! disableFieldCT)
        {
            auto ct = CellCenteredFieldCT();
            ct.correctGodunovFluxes (fluxes, ib);
        }
    };


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
    auto U0 = fieldOperator->generateConserved (solution.P);
    auto U = U0;

    for (int rk = 0; rk < rungeKuttaOrder; ++rk)
    {
        auto F = meshOperator->godunov (Fhat, solution.P, solution.B, footprint, startIndex, fieldCT);
        auto L = meshOperator->divergence (F);

        for (int n = 0; n < L.size(); ++n)
        {
            const double uold = U0[n];
            const double unew = U[n] - dt * L[n];
            U[n] = uold * (1 - b[rk]) + unew * b[rk];
        }

        solution.Z = fieldOperator->recoverPrimitive (U[interior], solution.P[interior]);
        solution.applyBoundaryCondition (*boundaryCondition);
    }
}

void MethodOfLinesTVD::makeFootprint (
    const MeshData& solution,
    const IntercellFluxScheme& fs,
    Shape3D& footprint,
    Index& startIndex,
    Region& interior) const
{
    int ng = fs.getStencilSize();

    for (int axis = 0; axis < 3; ++axis)
    {
        if (solution.P.size (axis) > 1)
        {
            footprint[axis] = 2 * ng;
            startIndex[axis] = -ng;
            interior.lower[axis] =  ng;
            interior.upper[axis] = -ng;
        }
        else
        {
            footprint[axis] = 0;
            startIndex[axis] = 0;
            interior.lower[axis] = 0;
            interior.upper[axis] = 0;
        }
    }

    if (! solution.getBoundaryShape().contains (footprint / 2))
    {
        throw std::logic_error ("Boundary region of mesh data is smaller than the scheme's stencil");
    }
}
