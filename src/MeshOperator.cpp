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

    auto cellsShape = geometry->cellsShape();
    auto facesShape = makeFacesShape(1);

    cellVolumes = Array (cellsShape);
    faceAreas   = Array (facesShape);

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

Array MeshOperator::generate (InitialDataFunction F, MeshLocation location) const
{
    ENSURE_GEOMETRY_IS_VALID;

    switch (location)
    {
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
        case MeshLocation::face:
        {
            auto shape = makeFacesShape (F (0.0, 0.0, 0.0).size());
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

                        for (int q = 0; q < shape[3]; ++q)
                        {
                            B(i, j, k, q, 0) = P0[q];
                            B(i, j, k, q, 1) = P1[q];
                            B(i, j, k, q, 2) = P2[q];
                        }
                    }
                }
            }
            return B;
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
            d[q] += F(1, 0, 0, q, 0) * A (1, 0, 0, 0, 0) - F(0, 0, 0, q, 0) * A (0, 0, 0, 0, 0);
            d[q] += F(0, 1, 0, q, 1) * A (0, 1, 0, 0, 1) - F(0, 0, 0, q, 1) * A (0, 0, 0, 0, 1);
            d[q] += F(0, 0, 1, q, 2) * A (0, 0, 1, 0, 2) - F(0, 0, 0, q, 2) * A (0, 0, 0, 0, 2);
        }
    };
    auto D = stencil.evaluate (div, flux, faceAreas);

    for (int i = 0; i < D.size(); ++i)
    {
        D[i] /= cellVolumes[i / numComponents];
    }
    return D;
}

Shape MeshOperator::makeFacesShape (int numComponents) const
{
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
