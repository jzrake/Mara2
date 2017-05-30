#include "BlockDecomposition.hpp"
#include "CartesianMeshGeometry.hpp"

using namespace Cow;




// ============================================================================
BlockDecomposition::BlockDecomposition (const std::shared_ptr<MeshGeometry> globalGeometry) :
globalGeometry (globalGeometry)
{
    // The geometry instance needs to be cartesian.
    assert (dynamic_cast<const CartesianMeshGeometry*> (globalGeometry.get()));

    auto world = MpiCommunicator::world();
    communicator = world.createCartesian (3, globalGeometry->fleshedOutAxes());
}

Cow::Shape BlockDecomposition::getGlobalShape() const
{
    return globalGeometry->cellsShape();
}

std::shared_ptr<MeshGeometry> BlockDecomposition::decompose() const
{
    auto globalGeom = dynamic_cast<const CartesianMeshGeometry*> (globalGeometry.get());
    auto localGeom = std::make_shared<CartesianMeshGeometry>();
    auto localRegion = getPatchRegion();

    localGeom->assignPatchIndex (getPatchIndex());
    localGeom->shape = localRegion.shape();
    localGeom->lower = globalGeom->coordinateAtIndex (
        localRegion.lower[0] - 0.5,
        localRegion.lower[1] - 0.5,
        localRegion.lower[2] - 0.5);
    localGeom->upper = globalGeom->coordinateAtIndex (
        localRegion.upper[0] - 0.5,
        localRegion.upper[1] - 0.5,
        localRegion.upper[2] - 0.5);
    return localGeom;
}

Cow::Region BlockDecomposition::getPatchRegion (MeshGeometry::PatchIndex index) const
{
    auto globalGeom = dynamic_cast<const CartesianMeshGeometry*> (globalGeometry.get());
    auto cartCoords = communicator.getCoordinates (index[0]);
    auto globalDims = communicator.getDimensions();
    auto globalShape = globalGeom->cellsShape();
    auto R = Region();

    for (int n = 0; n < 3; ++n)
    {
        R.lower[n] = startIndex (globalShape[n], globalDims[n], cartCoords[n]);
        R.upper[n] = R.lower[n] + partition (globalShape[n], globalDims[n], cartCoords[n]);
    }
    return R.absolute (globalShape);
}

Cow::Region BlockDecomposition::getPatchRegion() const
{
    return getPatchRegion (getPatchIndex());
}

MeshGeometry::PatchIndex BlockDecomposition::getPatchIndex() const
{
    return {{communicator.rank(), 0, 0, 0, 0}};
}

int BlockDecomposition::partition (int numElements, int numPartitions, int whichPartition) const
{
    const int elementsPerPartition = numElements / numPartitions;
    const int elementsLeftOver = numElements % numPartitions;
    return elementsPerPartition + (whichPartition < elementsLeftOver ? 1 : 0);
}

int BlockDecomposition::startIndex (int numElements, int numPartitions, int whichPartition) const
{
    const int i = whichPartition;
    const int n = numElements / numPartitions;
    const int r = numElements % numPartitions;
    const int s = i - r; // number of partitions of size n
    const int t = i - s; // number of partitions of size n + 1
    return n * s + (n + 1) * t;
}
