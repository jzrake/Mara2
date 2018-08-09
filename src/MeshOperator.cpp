#include <cmath>
#include <cassert>
#include <algorithm>
#include "MeshOperator.hpp"

#define ENSURE_GEOMETRY_IS_VALID do {\
if (! geometry)\
{\
    throw std::logic_error ("MeshOperator has no MeshGeometry instance");\
} } while (0)\

using namespace Cow;




// ============================================================================
GodunovStencil::GodunovStencil (int footprint, int numCellQ, int numFaceQ) :
cellCoords (footprint),
cellData (footprint, numCellQ),
faceData (numFaceQ),
godunovFlux (numCellQ)
{

}




// ============================================================================
MeshOperator::MeshOperator (std::shared_ptr<MeshGeometry> geometry) : geometry (geometry)
{

}

void MeshOperator::setMeshGeometry (std::shared_ptr<MeshGeometry> g)
{
    geometry = g;
}

const MeshGeometry& MeshOperator::getMeshGeometry() const
{
    ENSURE_GEOMETRY_IS_VALID;
    return *geometry;
}

Array MeshOperator::generate (InitialDataFunction F, MeshLocation location, VectorMode vectorMode, Shape3D boundaryShape) const
{
    ENSURE_GEOMETRY_IS_VALID;

    int nq = F (0.0, 0.0, 0.0).size();
    auto totalCellsShape = Shape3D (geometry->cellsShape()).increased (boundaryShape * 2);

    switch (vectorMode)
    {
        case VectorMode::scalars:
            break;
        case VectorMode::fluxish:
            if (nq != 3) throw std::logic_error ("Flux-ish data function returns size(4)=3");
            nq = 1;
            break;
        case VectorMode::emflike:
            if (nq != 3) throw std::logic_error ("EMF-like data function returns size(4)=3");
            nq = 1;
            break;
    }

    switch (location)
    {
        case MeshLocation::edge:
        {
            auto E = Array (totalCellsShape.increased(1).withComponents(nq).withRank(3));

            E.shape3D().deploy ([&] (int i, int j, int k)
            {
                int ic = i - boundaryShape[0];
                int jc = j - boundaryShape[1];
                int kc = k - boundaryShape[2];

                auto coord0 = geometry->coordinateAtIndex (ic, jc - 0.5, kc - 0.5);
                auto coord1 = geometry->coordinateAtIndex (ic - 0.5, jc, kc - 0.5);
                auto coord2 = geometry->coordinateAtIndex (ic - 0.5, jc - 0.5, kc);

                auto P0 = F (coord0[0], coord0[1], coord0[2]);
                auto P1 = F (coord1[0], coord1[1], coord1[2]);
                auto P2 = F (coord2[0], coord2[1], coord2[2]);

                switch (vectorMode)
                {
                    case VectorMode::scalars:
                    {
                        for (int q = 0; q < nq; ++q)
                        {
                            E (i, j, k, q, 0) = P0[q];
                            E (i, j, k, q, 1) = P1[q];
                            E (i, j, k, q, 2) = P2[q];
                        }
                        break;
                    }
                    case VectorMode::emflike:
                    {
                        auto nhat0 = geometry->edgeVector (ic, jc, kc, 0);
                        auto nhat1 = geometry->edgeVector (ic, jc, kc, 1);
                        auto nhat2 = geometry->edgeVector (ic, jc, kc, 2);

                        E (i, j, k, 0, 0) = nhat0.project (P0[0], P0[1], P0[2]);
                        E (i, j, k, 0, 1) = nhat1.project (P1[0], P1[1], P1[2]);
                        E (i, j, k, 0, 2) = nhat2.project (P2[0], P2[1], P2[2]);
                        break;
                    }
                    default:
                    {
                        throw std::logic_error ("Invalid vector mode for generating data on edges");
                    }
                }
            });
            return E;
        }
        case MeshLocation::face:
        {
            auto B = Array (totalCellsShape.increased(1).withComponents(nq).withRank(3));

            B.shape3D().deploy ([&] (int i, int j, int k)
            {
                int ic = i - boundaryShape[0];
                int jc = j - boundaryShape[1];
                int kc = k - boundaryShape[2];

                auto coord0 = geometry->coordinateAtIndex (ic - 0.5, jc, kc);
                auto coord1 = geometry->coordinateAtIndex (ic, jc - 0.5, kc);
                auto coord2 = geometry->coordinateAtIndex (ic, jc, kc - 0.5);

                auto P0 = F (coord0[0], coord0[1], coord0[2]);
                auto P1 = F (coord1[0], coord1[1], coord1[2]);
                auto P2 = F (coord2[0], coord2[1], coord2[2]);

                switch (vectorMode)
                {
                    case VectorMode::scalars:
                    {
                        for (int q = 0; q < nq; ++q)
                        {
                            B (i, j, k, q, 0) = P0[q];
                            B (i, j, k, q, 1) = P1[q];
                            B (i, j, k, q, 2) = P2[q];
                        }
                        break;
                    }
                    case VectorMode::fluxish:
                    {
                        auto nhat0 = geometry->faceNormal (ic, jc, kc, 0);
                        auto nhat1 = geometry->faceNormal (ic, jc, kc, 1);
                        auto nhat2 = geometry->faceNormal (ic, jc, kc, 2);

                        B (i, j, k, 0, 0) = nhat0.project (P0[0], P0[1], P0[2]);
                        B (i, j, k, 0, 1) = nhat1.project (P1[0], P1[1], P1[2]);
                        B (i, j, k, 0, 2) = nhat2.project (P2[0], P2[1], P2[2]);
                        break;
                    }
                    default:
                    {
                        throw std::logic_error ("Invalid vector mode for generating data on faces");
                    }
                }
            });
            return B;
        }
        case MeshLocation::cell:
        {
            int nq = F (0.0, 0.0, 0.0).size();
            auto P = Array (totalCellsShape.withComponents (nq));
            auto R = Region().withStride (3, nq);

            for (auto it = P[R].begin(); it != P.end(); ++it)
            {
                int ic = it.index()[0] - boundaryShape[0];
                int jc = it.index()[1] - boundaryShape[1];
                int kc = it.index()[2] - boundaryShape[2];

                auto coord = geometry->coordinateAtIndex (ic, jc, kc);
                auto P0 = F (coord[0], coord[1], coord[2]);

                for (int q = 0; q < nq; ++q)
                {
                    it[q] = P0[q];
                }
            }
            return P;
        }
        default:
        {
            throw std::logic_error ("Unsupported mesh location for MeshOperator::generate");
        }
    }
}

Array MeshOperator::generateSourceTerms (SourceTermsFunction S, const Array& P, Index start) const
{
    ENSURE_GEOMETRY_IS_VALID;
    auto sourceTermArray = Array (P.shape());

    P.shape3D().deploy ([&] (int i, int j, int k)
    {
        auto p = StateArray();
        auto x = geometry->coordinateAtIndex (i + start[0], j + start[1], k + start[2]);

        for (int q = 0; q < P.size(3); ++q)
        {
            p[q] = P(i, j, k, q);
        }
        auto s = S (x[0], x[1], x[2], p);

        for (int q = 0; q < P.size(3); ++q)
        {
            sourceTermArray (i, j, k, q) = s[q];
        }
    });
    return sourceTermArray;
}

Array MeshOperator::measure (MeshLocation location) const
{
    switch (location)
    {
        case MeshLocation::vert:
        {
            throw std::logic_error ("No measure associated with vertices");
        }
        case MeshLocation::edge:
        {
            auto lengths = Array (Shape3D (geometry->cellsShape()).increased().withRank(3));

            Array::deploy (lengths.shape(), [&] (int i, int j, int k)
            {
                lengths (i, j, k, 0, 0) = geometry->edgeLength (i, j, k, 0);
                lengths (i, j, k, 0, 1) = geometry->edgeLength (i, j, k, 1);
                lengths (i, j, k, 0, 2) = geometry->edgeLength (i, j, k, 2);
            });
            return lengths;
        }
        case MeshLocation::face:
        {
            auto areas = Array (Shape3D (geometry->cellsShape()).increased().withRank(3));

            Array::deploy (areas.shape(), [&] (int i, int j, int k)
            {
                areas (i, j, k, 0, 0) = geometry->faceArea (i, j, k, 0);
                areas (i, j, k, 0, 1) = geometry->faceArea (i, j, k, 1);
                areas (i, j, k, 0, 2) = geometry->faceArea (i, j, k, 2);
            });
            return areas;
        }
        case MeshLocation::cell:
        {
            auto volumes = Array (geometry->cellsShape());

            Array::deploy (volumes.shape(), [&] (int i, int j, int k)
            {
                volumes (i, j, k) = geometry->cellVolume (i, j, k);
            });
            return volumes;
        }
    }
    assert (false);
}

Array MeshOperator::linearCellDimension() const
{
    auto L = Array (geometry->cellsShape());

    Array::deploy (L.shape(), [&] (int i, int j, int k)
    {
        auto L3 = std::array<double, 3>
        {{
            geometry->cellLength (i, j, k, 0),
            geometry->cellLength (i, j, k, 1),
            geometry->cellLength (i, j, k, 2)
        }};
        L (i, j, k) = *std::min_element (L3.begin(), L3.end());
    });
    return L;
}

Array MeshOperator::cellCentroidCoordinates() const
{
    auto X = Array (Shape3D (geometry->cellsShape()).withComponents (3));

    Array::deploy (X.shape(), [&] (int i, int j, int k)
    {
        auto C = geometry->coordinateAtIndex (i, j, k);
        X (i, j, k, 0) = C[0];
        X (i, j, k, 1) = C[1];
        X (i, j, k, 2) = C[2];
    });
    return X;
}

Array MeshOperator::divergence (const Array& flux, double factor, Index start) const
{
    ENSURE_GEOMETRY_IS_VALID;

    if (flux.size(4) != 3)
    {
        throw std::logic_error ("Attempt to compute divergence from data whose size(4) != 3");
    }
    /*
                        F (i, j + 1, k) : 1
                                ^
                                |
                        +-----------------+
                        |                 |
                        |                 |
    F (i, j, k) : 0  <--|      div F      |--> F (i + 1, j, k) : 0
                        |                 |
                        |                 |
                        +-----------------+
                                |
                                v
                        F (i, j + 0, k) : 1
    */
#define A(di, dj, dk, a   ) geometry->faceArea   (i + di + start[0], j + dj + start[1], k + dk + start[2], a)
#define V(di, dj, dk      ) geometry->cellVolume (i + di + start[0], j + dj + start[1], k + dk + start[2])
#define F(di, dj, dk, q, a) flux                 (i + di, j + dj, k + dk, q, a)
#define D(di, dj, dk, q   ) divergence           (i + di, j + dj, k + dk, q, 0)

    const int nq = flux.size(3);
    auto divergence = Array (Shape3D (flux).reduced(1).withComponents (nq).withRank(1));

    Array::deploy (divergence.shape(), [&] (int i, int j, int k)
    {
        const double A1000 = A (1,0,0,0);
        const double A0000 = A (0,0,0,0);
        const double A0101 = A (0,1,0,1);
        const double A0001 = A (0,0,0,1);
        const double A0012 = A (0,0,1,2);
        const double A0002 = A (0,0,0,2);
        const double vol = V (0,0,0);

        for (int q = 0; q < nq; ++q)
        {
            D (0,0,0,q) = factor / vol * (
                +F (1,0,0,q,0) * A1000
                -F (0,0,0,q,0) * A0000
                +F (0,1,0,q,1) * A0101
                -F (0,0,0,q,1) * A0001
                +F (0,0,1,q,2) * A0012
                -F (0,0,0,q,2) * A0002);
        }
    });
    return divergence;

#undef A
#undef V
#undef F
#undef D
}

Array MeshOperator::curl (const Array& potential, Index start) const
{
    ENSURE_GEOMETRY_IS_VALID;

    if (potential.size(3) != 1)
    {
        throw std::logic_error ("Attempt to compute curl data with more than one component");
    }
    if (potential.size(4) != 3)
    {
        throw std::logic_error ("Attempt to compute curl of EMF-like data whose size(4) != 3");
    }
    /*
                    E (i, j + 1, k) : 0

                    +--------<--------+
                    |                 |
                    |                 |
    E (i, j, k) : 1 v      curl E     ^ E (i + 1, j, k) : 1
                    |                 |
                    |                 |
                    +-------->--------+

                    E (i, j + 0, k) : 0
    */
#define A(di, dj, dk, a) geometry->faceArea   (i + di + start[0], j + dj + start[1], k + dk + start[2], a)
#define L(di, dj, dk, a) geometry->edgeLength (i + di + start[0], j + dj + start[1], k + dk + start[2], a)
#define E(di, dj, dk, a) potential (i + di, j + dj, k + dk, 0, a)
#define B(di, dj, dk, a) curlfield (i + di, j + dj, k + dk, 0, a)

    auto curlfield = Array (potential.shape());

    Array::deploy (Shape3D (curlfield).reduced(1), [&] (int i, int j, int k)
    {
        B (0,0,0,0) = 1. / A (0,0,0,0) * (
            +E (0,1,0,2) * L (0,1,0,2)
            -E (0,0,0,2) * L (0,0,0,2)
            +E (0,0,1,1) * L (0,0,1,1)
            -E (0,0,0,1) * L (0,0,0,1));

        B (0,0,0,1) = 1. / A (0,0,0,1) * (
            +E (0,0,1,0) * L (0,0,1,0)
            -E (0,0,0,0) * L (0,0,0,0)
            +E (1,0,0,2) * L (1,0,0,2)
            -E (0,0,0,2) * L (0,0,0,2));

        B (0,0,0,2) = 1. / A (0,0,0,2) * (
            +E (1,0,0,1) * L (1,0,0,1)
            -E (0,0,0,1) * L (0,0,0,1)
            +E (0,1,0,0) * L (0,1,0,0)
            -E (0,0,0,0) * L (0,0,0,0));
    });
    return curlfield;

#undef A
#undef L
#undef E
#undef B
}

Array MeshOperator::godunov (
    GodunovStencil::IntercellFluxFunction Fhat,
    const Array& cellData,
    const Array& faceData,
    Shape footprint,
    Index start,
    FluxCorrection fluxCorrection) const
{
    ENSURE_GEOMETRY_IS_VALID;

    const int nq = cellData.size(3);
    const int np = faceData.size(3);

    auto& fp = footprint; // Footprint is always even (total cells in stencil)
    auto& mg = geometry;
    auto fluxShape = Shape3D (cellData).increased(1).withRank(3);
    auto result = Array (fluxShape);

    for (int sweep = 0; sweep < 3; ++sweep)
    {
        if (fp[sweep] == 0)
        {
            continue; // It's a zero-stencil in this direction
        }

        auto shape = cellData.shape3D().reduced (sweep, fp[sweep] - 1);
        auto gs = GodunovStencil (fp[sweep], nq, np);

        shape.deploy ([&] (int i, int j, int k)
        {
            const int ic = i + start[0]; // Indexes for querying geometry
            const int jc = j + start[1];
            const int kc = k + start[2];

            for (int n = 0; n < fp[sweep]; ++n)
            {
                for (int q = 0; q < nq; ++q)
                {
                    switch (sweep)
                    {
                        case 0: gs.cellData (n, q) = cellData (i + n, j, k, q); break;
                        case 1: gs.cellData (n, q) = cellData (i, j + n, k, q); break;
                        case 2: gs.cellData (n, q) = cellData (i, j, k + n, q); break;
                    }
                }
                switch (sweep)
                {
                    case 0: gs.cellCoords[n] = mg->coordinateAtIndex (ic + n, jc, kc); break;
                    case 1: gs.cellCoords[n] = mg->coordinateAtIndex (ic, jc + n, kc); break;
                    case 2: gs.cellCoords[n] = mg->coordinateAtIndex (ic, jc, kc + n); break;
                }
            }
            for (int p = 0; p < np; ++p)
            {
                switch (sweep)
                {
                    case 0: gs.faceData(p) = faceData (i + fp[0] / 2, j, k, p); break;
                    case 1: gs.faceData(p) = faceData (i, j + fp[1] / 2, k, p); break;
                    case 2: gs.faceData(p) = faceData (i, j, k + fp[2] / 2, p); break;
                }
            }

            gs.faceNormal = mg->faceNormal (ic, jc, kc, sweep);
            Fhat (gs);

            for (int q = 0; q < nq; ++q)
            {
                switch (sweep)
                {
                    case 0: result (i + fp[0] / 2, j, k, q, sweep) = gs.godunovFlux[q]; break;
                    case 1: result (i, j + fp[1] / 2, k, q, sweep) = gs.godunovFlux[q]; break;
                    case 2: result (i, j, k + fp[2] / 2, q, sweep) = gs.godunovFlux[q]; break;
                }
            }
        });
    }

    if (fluxCorrection)
    {
        fluxCorrection (result);
    }
    return result;
}

void MeshOperator::weight (Coordinate point, std::vector<Index>& indexes, std::vector<double>& weights) const
{
    ENSURE_GEOMETRY_IS_VALID;
    const int particleShape = 0;


    if (particleShape == 0)
    {
        const auto i = geometry->indexAtCoordinate (point);
        const auto V = geometry->cellVolume (i[0], i[1], i[2]);

        indexes.resize(1);
        weights.resize(1);
        indexes[0] = i;
        weights[0] = 1 / V;
    }

    if (particleShape == 1)
    {
        const auto i0 = geometry->indexAtCoordinate (point);
        const auto im = Index {{i0[0] - 1, i0[1], i0[2], 0, 0}};
        const auto ip = Index {{i0[0] + 1, i0[1], i0[2], 0, 0}};

        const double x = point[0];
        const double dx = geometry->cellLength (i0[0], i0[1], i0[2], 0);
        const double xm = geometry->coordinateAtIndex (im)[0];
        const double x0 = geometry->coordinateAtIndex (i0)[0];
        const double xp = geometry->coordinateAtIndex (ip)[0];

        const double wm = std::max (0.0, 1 - std::fabs (xm - x) / dx);
        const double w0 = std::max (0.0, 1 - std::fabs (x0 - x) / dx);
        const double wp = std::max (0.0, 1 - std::fabs (xp - x) / dx);
        const double wt = wm + w0 + wp;

        indexes.resize(3);
        weights.resize(3);

        const auto Vm = geometry->cellVolume (im[0], im[1], im[2]);
        const auto V0 = geometry->cellVolume (i0[0], i0[1], i0[2]);
        const auto Vp = geometry->cellVolume (ip[0], ip[1], ip[2]);

        indexes[0] = im;
        indexes[1] = i0;
        indexes[2] = ip;
        weights[0] = wm / Vm / wt;
        weights[1] = w0 / V0 / wt;
        weights[2] = wp / Vp / wt;
    }

    if (particleShape == 2)
    {
        const auto i0 = geometry->indexAtCoordinate (point);
        const auto im1 = Index {{i0[0] - 1, i0[1], i0[2], 0, 0}};
        const auto ip1 = Index {{i0[0] + 1, i0[1], i0[2], 0, 0}};
        const auto im2 = Index {{i0[0] - 2, i0[1], i0[2], 0, 0}};
        const auto ip2 = Index {{i0[0] + 2, i0[1], i0[2], 0, 0}};

        const double x = point[0];
        const double dx = geometry->cellLength (i0[0], i0[1], i0[2], 0);
        const double x0 = geometry->coordinateAtIndex (i0)[0];
        const double xm1 = geometry->coordinateAtIndex (im1)[0];
        const double xp1 = geometry->coordinateAtIndex (ip1)[0];
        const double xm2 = geometry->coordinateAtIndex (im2)[0];
        const double xp2 = geometry->coordinateAtIndex (ip2)[0];

        const double w0 = std::max (0.0, 1 - std::fabs (x0 - x) / dx / 2);
        const double wm1 = std::max (0.0, 1 - std::fabs (xm1 - x) / dx / 2);
        const double wp1 = std::max (0.0, 1 - std::fabs (xp1 - x) / dx / 2);
        const double wm2 = std::max (0.0, 1 - std::fabs (xm2 - x) / dx / 2);
        const double wp2 = std::max (0.0, 1 - std::fabs (xp2 - x) / dx / 2);
        const double wt = w0 + wm1 + wp1 + wm2 + wp2;

        indexes.resize(5);
        weights.resize(5);

        const auto V0 = geometry->cellVolume (i0[0], i0[1], i0[2]);

        indexes[0] = i0;
        indexes[1] = im1;
        indexes[2] = ip1;
        indexes[3] = im2;
        indexes[4] = ip2;
        weights[0] = w0 / V0 / wt;
        weights[1] = wm1 / V0 / wt;
        weights[2] = wp1 / V0 / wt;
        weights[3] = wm2 / V0 / wt;
        weights[4] = wp2 / V0 / wt;
    }
}
