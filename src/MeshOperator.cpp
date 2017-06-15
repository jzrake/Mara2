#include "MeshOperator.hpp"

#define ENSURE_GEOMETRY_IS_VALID do {\
if (! geometry)\
{\
    throw std::logic_error ("MeshOperator has no MeshGeometry instance");\
} } while (0)\

using namespace Cow;




// ============================================================================
class MeshOperator::RichShape
{
public:
    RichShape (Shape S) : S (S) {}
    RichShape (const Array& A) : S (A.shape()) {}
    operator Shape() { return S; }
    RichShape reduced (int delta=1) const
    {
        return Shape {{ S[0] - delta, S[1] - delta, S[2] - delta, S[3], S[4] }};
    }
    RichShape increased (int delta=1) const
    {
        return Shape {{ S[0] + delta, S[1] + delta, S[2] + delta, S[3], S[4] }};
    }
    RichShape reduced (Shape delta) const
    {
        return Shape {{ S[0] - delta[0], S[1] - delta[1], S[2] - delta[2], S[3], S[4] }};
    }
    RichShape increased (Shape delta) const
    {
        return Shape {{ S[0] + delta[0], S[1] + delta[1], S[2] + delta[2], S[3], S[4] }};
    }
    RichShape reduced (int axis, int delta) const
    {
        auto s = S;
        s[axis] -= delta;
        return s;
    }
    RichShape increased (int axis, int delta) const
    {
        auto s = S;
        s[axis] += delta;
        return s;
    }
    RichShape withComponents (int numComponents)
    {
        return Shape {{ S[0], S[1], S[2], numComponents, S[4] }};
    }
    RichShape withRank (int rank)
    {
        return Shape {{ S[0], S[1], S[2], S[3], rank }};
    }
private:
    Shape S;
};




// ============================================================================
GodunovStencil::GodunovStencil (int footprint, int numCellQ, int numFaceQ) :
cellCoords (footprint),
cellData (footprint, numCellQ),
faceData (numFaceQ),
faceFlux (numCellQ)
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

Array MeshOperator::generate (InitialDataFunction F, MeshLocation location, VectorMode vectorMode) const
{
    ENSURE_GEOMETRY_IS_VALID;

    int nq = F (0.0, 0.0, 0.0).size();

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
            auto E = Array (RichShape (geometry->cellsShape()).increased(1).withComponents(nq).withRank(3));

            Array::deploy (E.shape(), [&] (int i, int j, int k)
            {
                auto coord0 = geometry->coordinateAtIndex (i, j - 0.5, k - 0.5);
                auto coord1 = geometry->coordinateAtIndex (i - 0.5, j, k - 0.5);
                auto coord2 = geometry->coordinateAtIndex (i - 0.5, j - 0.5, k);

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
                        auto nhat0 = geometry->edgeVector (i, j, k, 0);
                        auto nhat1 = geometry->edgeVector (i, j, k, 1);
                        auto nhat2 = geometry->edgeVector (i, j, k, 2);

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
            auto B = Array (RichShape (geometry->cellsShape()).increased(1).withComponents(nq).withRank(3));

            Array::deploy (B.shape(), [&] (int i, int j, int k)
            {
                auto coord0 = geometry->coordinateAtIndex (i - 0.5, j, k);
                auto coord1 = geometry->coordinateAtIndex (i, j - 0.5, k);
                auto coord2 = geometry->coordinateAtIndex (i, j, k - 0.5);

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
                        auto nhat0 = geometry->faceNormal (i, j, k, 0);
                        auto nhat1 = geometry->faceNormal (i, j, k, 1);
                        auto nhat2 = geometry->faceNormal (i, j, k, 2);

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
            auto P = Array (RichShape (geometry->cellsShape()).withComponents (nq));
            auto R = Region().withStride (3, nq);

            for (auto it = P[R].begin(); it != P.end(); ++it)
            {
                auto coord = geometry->coordinateAtIndex (it.index());
                auto P0 = F (coord[0], coord[1], coord[2]);

                for (unsigned int q = 0; q < nq; ++q)
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

Array MeshOperator::divergence (const Array& flux, Index start) const
{
    ENSURE_GEOMETRY_IS_VALID;

    if (flux.size(4) != 3)
    {
        throw std::logic_error ("Attempt to compute divergence from data whose size(4) != 3");
    }
    /*
                    F (i, j + 1, k) : 1

                    +-----------------+
                    |                 |
                    |                 |
    F (i, j, k) : 0 |      div F      | F (i + 1, j, k) : 0
                    |                 |
                    |                 |
                    +-----------------+

                    F (i, j + 0, k) : 1
    */
#define A(di, dj, dk, a   ) geometry->faceArea   (i + di + start[0], j + dj + start[1], k + dk + start[2], a)
#define V(di, dj, dk      ) geometry->cellVolume (i + di + start[0], j + dj + start[1], k + dk + start[2])
#define F(di, dj, dk, q, a) flux                 (i + di, j + dj, k + dk, q, a)
#define D(di, dj, dk, q   ) divergence           (i + di, j + dj, k + dk, q, 0)

    const int nq = flux.size(3);
    auto divergence = Array (RichShape (flux).reduced(1).withComponents (nq).withRank(1));

    Array::deploy (divergence.shape(), [&] (int i, int j, int k)
    {
        for (int q = 0; q < nq; ++q)
        {
            D (0,0,0,q) = 1. / V (0,0,0) * (
                +F (1,0,0,q,0) * A (0,0,0,0)
                -F (0,0,0,q,0) * A (1,0,0,0)
                +F (0,1,0,q,1) * A (0,1,0,1)
                -F (0,0,0,q,1) * A (0,0,0,1)
                +F (0,0,1,q,2) * A (0,0,1,2)
                -F (0,0,0,q,2) * A (0,0,0,2));
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

                    +-----------------+
                    |                 |
                    |                 |
    E (i, j, k) : 1 |      curl E     | E (i + 1, j, k) : 1
                    |                 |
                    |                 |
                    +-----------------+

                    E (i, j + 0, k) : 0
    */
#define A(di, dj, dk, a) geometry->faceArea   (i + di + start[0], j + dj + start[1], k + dk + start[2], a)
#define L(di, dj, dk, a) geometry->edgeLength (i + di + start[0], j + dj + start[1], k + dk + start[2], a)
#define E(di, dj, dk, a) potential (i + di, j + dj, k + dk, 0, a)
#define B(di, dj, dk, a) curlfield (i + di, j + dj, k + dk, 0, a)

    auto curlfield = Array (potential.shape());

    Array::deploy (RichShape (curlfield).reduced(1), [&] (int i, int j, int k)
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
    Index start) const
{
    const int nq = cellData.size(3);
    const int np = faceData.size(3);

    auto& fp = footprint; // Footprint is always even (total cells in stencil)
    auto& mg = geometry;
    auto fluxShape = RichShape (cellData).increased(1).withRank(3);
    auto result = Array (fluxShape);

    for (int sweep = 0; sweep < 3; ++sweep)
    {
        auto shape = RichShape (mg->cellsShape()).reduced (sweep, fp[sweep]);

        Array::deploy (shape, [&] (int i, int j, int k)
        {
            const int ic = i + start[0]; // Indexes for querying geometry
            const int jc = j + start[1];
            const int kc = k + start[2];
            
            auto gs = GodunovStencil (fp[sweep], nq, np);

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
                    case 0: gs.faceData (p) = faceData (i + fp[0] / 2, j, k, p); break;
                    case 1: gs.faceData (p) = faceData (i, j + fp[1] / 2, k, p); break;
                    case 2: gs.faceData (p) = faceData (i, j, k + fp[2] / 2, p); break;
                }
            }

            gs.faceNormal = mg->faceNormal (ic, jc, kc, sweep);
            Fhat (gs);

            for (int q = 0; q < nq; ++q)
            {
                switch (sweep)
                {
                    case 0: result (i + fp[0] / 2, j, k, q, sweep) = gs.faceFlux[q]; break;
                    case 1: result (i, j + fp[1] / 2, k, q, sweep) = gs.faceFlux[q]; break;
                    case 2: result (i, j, k + fp[2] / 2, q, sweep) = gs.faceFlux[q]; break;
                }
            }
        });
    }



    // auto fs0 = GodunovStencil (fp[0], nq, np);
    // auto fs1 = GodunovStencil (fp[1], nq, np);
    // auto fs2 = GodunovStencil (fp[2], nq, np);


    // Array::deploy (RichShape (mg->cellsShape()).reduced (fp), [&] (int i, int j, int k)
    // {
    //     const int ic = i + start[0]; // Indexes for querying geometry
    //     const int jc = j + start[1];
    //     const int kc = k + start[2];

    //     const int ir = i + fp[0] / 2; // Indexes into resulting flux array
    //     const int jr = j + fp[1] / 2;
    //     const int kr = k + fp[2] / 2;

    //     for (int q = 0; q < nq; ++q)
    //     {
    //         for (int n = 0; n < fp[0]; ++n) fs0.cellData (n, q) = cellData (i + n, j, k, q);
    //         for (int n = 0; n < fp[1]; ++n) fs1.cellData (n, q) = cellData (i, j + n, k, q);
    //         for (int n = 0; n < fp[2]; ++n) fs2.cellData (n, q) = cellData (i, j, k + n, q);
    //     }

    //     for (int n = 0; n < fp[0]; ++n) fs0.cellCoords[n] = mg->coordinateAtIndex (ic + n, jc, kc);
    //     for (int n = 0; n < fp[1]; ++n) fs1.cellCoords[n] = mg->coordinateAtIndex (ic, jc + n, kc);
    //     for (int n = 0; n < fp[2]; ++n) fs2.cellCoords[n] = mg->coordinateAtIndex (ic, jc, kc + n);

    //     fs0.faceNormal = mg->faceNormal (ic, jc, kc, 0);
    //     fs1.faceNormal = mg->faceNormal (ic, jc, kc, 1);
    //     fs2.faceNormal = mg->faceNormal (ic, jc, kc, 2);

    //     Fhat (fs0);
    //     Fhat (fs1);
    //     Fhat (fs2);

    //     for (int q = 0; q < nq; ++q)
    //     {
    //         result (ir, jr, kr, q, 0) = fs0.faceFlux[q];
    //         result (ir, jr, kr, q, 1) = fs1.faceFlux[q];
    //         result (ir, jr, kr, q, 2) = fs2.faceFlux[q];
    //     }
    // });

    return result;
}
