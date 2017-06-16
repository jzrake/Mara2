#include "MeshData.hpp"




MeshData::MeshData (Shape baseShape, Shape guardZone, int numComponents)
{
	auto cellsShape = baseShape;
	auto facesShape = baseShape;

    for (int n = 0; n < 3; ++n)
    {
        cellsShape[n] += 2 * guardZone[n];
        facesShape[n] += 2 * guardZone[n] + 1;
        updateableRegion.lower[n] =  guardZone[n];
        updateableRegion.upper[n] = -guardZone[n];
    }

    P = Array (cellsShape);
    B = Array (facesShape);
    Z = Array (cellsShape);
}

Cow::Array::Reference MeshData::getPrimitive (int fieldIndex)
{
    Region R = updateableRegion;
    R.stride[3] = 1;

    if (fieldIndex != -1)
    {
        R.lower[3] = fieldIndex;
        R.upper[3] = fieldIndex + 1;
    }
    return P[R];
}

Cow::Array::Reference MeshData::getPrimitiveVector (int fieldIndex)
{
    Region R = updateableRegion;
    R.stride[3] = 1;
    R.lower[3] = fieldIndex;
    R.upper[3] = fieldIndex + 3;
    return P[R];
}

Cow::Array::Reference MeshData::getZoneHealth()
{
    return Z[Region()];
}
