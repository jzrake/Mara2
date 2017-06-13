#include "InitialDataGenerator.hpp"

using namespace Cow;



// ============================================================================
Array InitialDataGenerator::generatePrimitive (InitialDataFunction F, const MeshGeometry& geometry) const
{
    auto shape = geometry.cellsShape();
    shape[3] = F (0.0, 0.0, 0.0).size();

    auto P = Array (shape);
    auto R = Region().withStride (3, shape[3]);

    for (auto it = P[R].begin(); it != P.end(); ++it)
    {
        auto coord = geometry.coordinateAtIndex (it.index());
        auto P0 = F (coord[0], coord[1], coord[2]);

        for (unsigned int q = 0; q < P0.size(); ++q)
        {
            it[q] = P0[q];
        }
    }
    return P;
}
