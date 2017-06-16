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

}

void MethodOfLinesTVD::advance (double dt, MeshData& solution) const
{
    if (! fieldOperator)     throw std::logic_error ("No FieldOperator instance");
    if (! meshOperator)      throw std::logic_error ("No MeshOperator instance");
    if (! boundaryCondition) throw std::logic_error ("No BoundaryCondition instance");

    auto cl = fieldOperator->getConservationLaw();
    auto rs = std::make_shared<HlleRiemannSolver>();
    auto fs = std::make_shared<MethodOfLines>();

    fs->setRiemannSolver (rs);

    int nq = cl->getNumConserved();
    auto footprint = Shape {{ 2, 0, 0, }};
    auto start = Shape {{ -1, 0, 0 }};

    auto Fhat = [&] (GodunovStencil& stencil)
    {
        IntercellFluxScheme::FaceData D;
        D.areaElement = stencil.faceNormal.cartesian();
        D.stencilData = stencil.cellData;
        D.conservationLaw = cl;

        auto S = fs->intercellFlux (D);

        for (int q = 0; q < nq; ++q)
        {
            stencil.godunovFlux[q] = S.F[q];
        }
    };

    auto U = fieldOperator->generateConserved (solution.P);
    auto F = meshOperator->godunov (Fhat, solution.P, solution.B, footprint, start);
    auto L = meshOperator->divergence (F);

    for (int n = 0; n < L.size(); ++n)
    {
        U[n] -= dt * L[n];
    }

    auto newP = fieldOperator->recoverPrimitive (U);

    for (int axis = 0; axis < 3; ++axis)
    {
        if (footprint[axis] > 1)
        {
            boundaryCondition->apply (newP,
                MeshLocation::cell,
                MeshBoundary::left,
                axis,
                footprint[axis] / 2);

            boundaryCondition->apply (newP,
                MeshLocation::cell,
                MeshBoundary::right,
                axis,
                footprint[axis] / 2);
        }
    }

    solution.P = std::move (newP);
}
