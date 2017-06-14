#include "MeshOperator.hpp"
#include "Stencil.hpp"

#define ENSURE_GEOMETRY_IS_VALID do {\
if (! geometry)\
{\
    throw std::logic_error ("MeshOperator has no MeshGeometry instance");\
} } while (0)\

using namespace Cow;





MeshOperator::MeshOperator()
{

}

void MeshOperator::setMeshGeometry (std::shared_ptr<MeshGeometry> g)
{
    geometry = g;

    auto edgesShape = makeEdgesShape (1, VectorMode::scalars);
    auto facesShape = makeFacesShape (1, VectorMode::scalars);
    auto cellsShape = makeCellsShape (1);

    edgeLengths = Array (edgesShape);
    faceAreas   = Array (facesShape);
    cellVolumes = Array (cellsShape);

    for (int i = 0; i < edgesShape[0]; ++i)
    {
        for (int j = 0; j < edgesShape[1]; ++j)
        {
            for (int k = 0; k < edgesShape[2]; ++k)
            {
                edgeLengths (i, j, k, 0, 0) = geometry->edgeLength (i, j, k, 0);
                edgeLengths (i, j, k, 0, 1) = geometry->edgeLength (i, j, k, 1);
                edgeLengths (i, j, k, 0, 2) = geometry->edgeLength (i, j, k, 2);
            }
        }
    }

    for (int i = 0; i < facesShape[0]; ++i)
    {
        for (int j = 0; j < facesShape[1]; ++j)
        {
            for (int k = 0; k < facesShape[2]; ++k)
            {
                faceAreas (i, j, k, 0, 0) = geometry->faceArea (i, j, k, 0);
                faceAreas (i, j, k, 0, 1) = geometry->faceArea (i, j, k, 1);
                faceAreas (i, j, k, 0, 2) = geometry->faceArea (i, j, k, 2);
            }
        }
    }

    for (int i = 0; i < cellsShape[0]; ++i)
    {
        for (int j = 0; j < cellsShape[1]; ++j)
        {
            for (int k = 0; k < cellsShape[2]; ++k)
            {
                cellVolumes (i, j, k) = geometry->cellVolume (i, j, k);
            }
        }
    }

}

Array MeshOperator::generate (InitialDataFunction F, MeshLocation location, VectorMode vectorMode) const
{
    ENSURE_GEOMETRY_IS_VALID;

    switch (location)
    {
        case MeshLocation::edge:
        {
            auto shape = makeEdgesShape (F (0.0, 0.0, 0.0).size(), vectorMode);
            auto E = Array (shape);

            for (int i = 0; i < shape[0]; ++i)
            {
                for (int j = 0; j < shape[1]; ++j)
                {
                    for (int k = 0; k < shape[2]; ++k)
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
                                for (int q = 0; q < shape[3]; ++q)
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
                    }
                }
            }
            return E;
        }
        case MeshLocation::face:
        {
            auto shape = makeFacesShape (F (0.0, 0.0, 0.0).size(), vectorMode);
            auto B = Array (shape);

            for (int i = 0; i < shape[0]; ++i)
            {
                for (int j = 0; j < shape[1]; ++j)
                {
                    for (int k = 0; k < shape[2]; ++k)
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
                                for (int q = 0; q < shape[3]; ++q)
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
                    }
                }
            }
            return B;
        }
        case MeshLocation::cell:
        {
            auto shape = makeCellsShape (F (0.0, 0.0, 0.0).size());
            auto P = Array (shape);
            auto R = Region().withStride (3, shape[3]);

            for (auto it = P[R].begin(); it != P.end(); ++it)
            {
                auto coord = geometry->coordinateAtIndex (it.index());
                auto P0 = F (coord[0], coord[1], coord[2]);

                for (unsigned int q = 0; q < shape[3]; ++q)
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

Array MeshOperator::divergence (const Array& flux, MeshLocation location) const
{
    ENSURE_GEOMETRY_IS_VALID;

    const int numComponents = flux.size(3);

    auto stencil = Stencil();
    stencil.setFootprintLower (0, 0, 0);
    stencil.setFootprintUpper (1, 1, 1);
    stencil.setCodomainRank (numComponents, 1);

    auto div = [&] (const Array& F, const Array &A, Array& d)
    {
        for (int q = 0; q < numComponents; ++q)
        {
            d[q] += F (1,0,0,q,0) * A (1,0,0,0,0) - F (0,0,0,q,0) * A (0,0,0,0,0);
            d[q] += F (0,1,0,q,1) * A (0,1,0,0,1) - F (0,0,0,q,1) * A (0,0,0,0,1);
            d[q] += F (0,0,1,q,2) * A (0,0,1,0,2) - F (0,0,0,q,2) * A (0,0,0,0,2);
        }
    };
    auto D = stencil.evaluate (div, flux, faceAreas);

    for (int i = 0; i < D.size(); ++i)
    {
        D[i] /= cellVolumes[i / numComponents];
    }
    return D;
}

Array MeshOperator::curl (const Array& potential, MeshLocation location) const
{
    ENSURE_GEOMETRY_IS_VALID;

    if (potential.size(4) != 3)
    {
        throw std::logic_error ("Attempt to compute curl of EMF-like data whose size(4) != 3");
    }

    auto stencil0 = Stencil();
    auto stencil1 = Stencil();
    auto stencil2 = Stencil();
    stencil0.setFootprintUpper (0, 1, 1);
    stencil1.setFootprintUpper (1, 0, 1);
    stencil2.setFootprintUpper (1, 1, 0);

    auto rot0 = [&] (const Array& E, const Array &L, const Array &A, Array& r)
    {
        r[0] += E (0,1,0,0,2) * L (0,1,0,0,2) - E (0,0,0,0,2) * L (0,0,0,0,2);
        r[0] += E (0,0,1,0,1) * L (0,0,1,0,1) - E (0,0,0,0,1) * L (0,0,0,0,1);
        r[0] /= A (0,0,0,0,0);
    };
    auto rot1 = [&] (const Array& E, const Array &L, const Array &A, Array& r)
    {
        r[0] += E (0,0,1,0,0) * L (0,0,1,0,0) - E (0,0,0,0,0) * L (0,0,0,0,0);
        r[0] += E (1,0,0,0,2) * L (1,0,0,0,2) - E (0,0,0,0,2) * L (0,0,0,0,2);
        r[0] /= A (0,0,0,0,1);
    };
    auto rot2 = [&] (const Array& E, const Array &L, const Array &A, Array& r)
    {
        r[0] += E (1,0,0,0,1) * L (1,0,0,0,1) - E (0,0,0,0,1) * L (0,0,0,0,1);
        r[0] += E (0,1,0,0,0) * L (0,1,0,0,0) - E (0,0,0,0,0) * L (0,0,0,0,0);
        r[0] /= A (0,0,0,0,2);
    };

    auto B0 = stencil0.evaluate (rot0, potential, edgeLengths, faceAreas);
    auto B1 = stencil1.evaluate (rot1, potential, edgeLengths, faceAreas);
    auto B2 = stencil2.evaluate (rot2, potential, edgeLengths, faceAreas);

    auto B = Array (makeFacesShape(3, VectorMode::fluxish));

    B.insert (B0, Region().withLower (4, 0).withUpper (4, 1).withUpper (1, -1).withUpper (2, -1));
    B.insert (B1, Region().withLower (4, 1).withUpper (4, 2).withUpper (2, -1).withUpper (0, -1));
    B.insert (B2, Region().withLower (4, 2).withUpper (4, 3).withUpper (0, -1).withUpper (1, -1));

    return B;

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
    // auto stencil = Stencil();
    // stencil.setFootprintLower (0, 0, 0);
    // stencil.setFootprintUpper (1, 1, 1);
    // stencil.setCodomainRank (1, 3);

    // auto rot = [&] (const Array& E, const Array &L, const Array &A, Array& r)
    // {
    //     r[0] += E (0,1,0,0,2) * L (0,1,0,0,2) - E (0,0,0,0,2) * L (0,0,0,0,2);
    //     r[0] += E (0,0,1,0,1) * L (0,0,1,0,1) - E (0,0,0,0,1) * L (0,0,0,0,1);

    //     r[1] += E (0,0,1,0,0) * L (0,0,1,0,0) - E (0,0,0,0,0) * L (0,0,0,0,0);
    //     r[1] += E (1,0,0,0,2) * L (1,0,0,0,2) - E (0,0,0,0,2) * L (0,0,0,0,2);

    //     r[2] += E (1,0,0,0,1) * L (1,0,0,0,1) - E (0,0,0,0,1) * L (0,0,0,0,1);
    //     r[2] += E (0,1,0,0,0) * L (0,1,0,0,0) - E (0,0,0,0,0) * L (0,0,0,0,0);

    //     r[0] /= A (0,0,0,0,0);
    //     r[1] /= A (0,0,0,0,1);
    //     r[2] /= A (0,0,0,0,2);
    // };
    // return stencil.evaluate (rot, potential, edgeLengths, faceAreas);
}

Shape MeshOperator::makeEdgesShape (int numComponents,VectorMode mode) const
{
    if (mode == VectorMode::emflike)
    {
        if (numComponents != 3)
        {
            throw std::logic_error ("Attempt to assign EMF-like data whose shape is not 3");
        }
        numComponents = 1;
    }

    auto S = geometry->cellsShape();
    S[0] += 1;
    S[1] += 1;
    S[2] += 1;
    S[3] = numComponents;
    S[4] = 3;
    return S;
}

Shape MeshOperator::makeFacesShape (int numComponents, VectorMode mode) const
{
    if (mode == VectorMode::fluxish)
    {
        if (numComponents != 3)
        {
            throw std::logic_error ("Attempt to assign flux-like data whose shape is not 3");
        }
        numComponents = 1;
    }

    auto S = geometry->cellsShape();
    S[0] += 1;
    S[1] += 1;
    S[2] += 1;
    S[3] = numComponents;
    S[4] = 3;
    return S;
}

Shape MeshOperator::makeCellsShape (int numComponents) const
{
    auto S = geometry->cellsShape();
    S[3] = numComponents;
    return S;
}
