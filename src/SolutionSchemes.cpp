#include "BoundaryConditions.hpp"
#include "CartesianMeshGeometry.hpp"
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
    auto fs = std::make_shared<MethodOfLines>();
    setIntercellFluxScheme (fs);

    rungeKuttaOrder = 1;
}

int MethodOfLinesTVD::getStencilSize() const
{
    if (! fluxScheme) throw std::logic_error ("No IntercellFluxScheme instance");
    return fluxScheme->getStencilSize();
}

void MethodOfLinesTVD::setIntercellFluxScheme (std::shared_ptr<IntercellFluxScheme> fs)
{
    auto ng = fs->getStencilSize();
    footprint = Shape {{ 2 * ng, 0, 0, }};
    startIndex = Index {{ -ng, 0, 0 }};
    fluxScheme = fs;
}

void MethodOfLinesTVD::advance (MeshData& solution, double dt) const
{
    if (! fieldOperator)     throw std::logic_error ("No FieldOperator instance");
    if (! meshOperator)      throw std::logic_error ("No MeshOperator instance");
    if (! boundaryCondition) throw std::logic_error ("No BoundaryCondition instance");
    if (! fluxScheme)        throw std::logic_error ("No IntercellFluxScheme instance");
    if (! solution.getBoundaryShape().contains (footprint / 2))
    {
        throw std::logic_error ("Boundary region of mesh data is smaller than the scheme's stencil");
    }


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
        auto F = meshOperator->godunov (Fhat, solution.P, solution.B, footprint, startIndex);
        auto L = meshOperator->divergence (F);

        for (int n = 0; n < L.size(); ++n)
        {
            const double uold = U0[n];
            const double unew = U[n] - dt * L[n];
            U[n] = uold * (1 - b[rk]) + unew * b[rk];
        }

        solution.P = std::move (fieldOperator->recoverPrimitive (U));
        solution.applyBoundaryCondition (*boundaryCondition);
    }
}
