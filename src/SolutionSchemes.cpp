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
    auto rs = std::make_shared<HlleRiemannSolver>();
    auto fs = std::make_shared<MethodOfLinesPlm>();
    auto ng = fs->getStencilSize();

    fs->setRiemannSolver (rs);
    footprint = Shape {{ 2 * ng, 0, 0, }};
    startIndex = Index {{ -ng, 0, 0 }};
    fluxScheme = fs;
}

int MethodOfLinesTVD::getStencilSize() const
{
    return fluxScheme->getStencilSize();
}

void MethodOfLinesTVD::advance (MeshData& solution, double dt) const
{
    if (! fieldOperator)     throw std::logic_error ("No FieldOperator instance");
    if (! meshOperator)      throw std::logic_error ("No MeshOperator instance");
    if (! boundaryCondition) throw std::logic_error ("No BoundaryCondition instance");
    if (! solution.getBoundaryShape().contains (footprint / 2))
    {
        throw std::logic_error ("Boundary region of mesh data is smaller than the scheme's stencil");
    }

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

    auto U = fieldOperator->generateConserved (solution.P);
    auto F = meshOperator->godunov (Fhat, solution.P, solution.B, footprint, startIndex);
    auto L = meshOperator->divergence (F);

    for (int n = 0; n < L.size(); ++n)
    {
        U[n] -= dt * L[n];
    }

    auto newP = fieldOperator->recoverPrimitive (U);
    solution.P = std::move (newP);
    solution.applyBoundaryCondition (*boundaryCondition);
}
