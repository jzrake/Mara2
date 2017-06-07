#include <iostream>
#include "BlockDecomposition.hpp"
#include "CartesianMeshGeometry.hpp"

using namespace Cow;




// ============================================================================
class BlockDecomposedBC : public BoundaryCondition
{
public:
    BlockDecomposedBC (const BlockDecomposition& block, std::shared_ptr<BoundaryCondition> physicalBC)
    : block (block), physicalBC (physicalBC)
    {

    }

    void apply (
        Cow::Array& A,
        MeshLocation location,
        MeshBoundary boundary,
        int axis,
        int numGuard) const
    {
        auto send = Region(); // region to send (valid region)
        auto recv = Region(); // region to receive (guard region)

        switch (boundary)
        {
            // ----------------------------------------------------------------
            // Send to the right, receive from the left. Note the 'R' in the
            // shiftExchange call refers to the send direction, while the
            // boundary function argument (left in this case) refers to which
            // boundary of our mesh needs to be filled.
            // ----------------------------------------------------------------
            case MeshBoundary::left:
            {
                send.lower[axis] = -2 * numGuard;
                send.upper[axis] = -1 * numGuard;
                recv.lower[axis] = +0;
                recv.upper[axis] = +1 * numGuard;
                block.communicator.shiftExchange (A, axis, 'R', send, recv);
                break;
            }
            // ----------------------------------------------------------------
            // Send to the left, receive from the right. Note the 'L' in the
            // shiftExchange call refers to the send direction, while the
            // boundary function argument (right in this case) refers to which
            // boundary of our mesh needs to be filled.
            // ----------------------------------------------------------------
            case MeshBoundary::right:
            {
                send.lower[axis] = +1 * numGuard;
                send.upper[axis] = +2 * numGuard;
                recv.lower[axis] = -1 * numGuard;
                recv.upper[axis] = -0;
                block.communicator.shiftExchange (A, axis, 'L', send, recv);
                break;
            }
        }

        // --------------------------------------------------------------------
        // We apply physical boundary conditions if we are at the edge of the
        // block-decomposed domain.
        // --------------------------------------------------------------------
        auto C = block.communicator;
        bool isWallL = C.getCoordinates()[axis] == 0;
        bool isWallR = C.getCoordinates()[axis] == C.getDimensions()[axis] - 1;

        if (   (boundary == MeshBoundary::left  && isWallL)
            || (boundary == MeshBoundary::right && isWallR))
        {
            physicalBC->apply (A, location, boundary, axis, numGuard);
        }
    }

private:
    const BlockDecomposition& block;
    std::shared_ptr<BoundaryCondition> physicalBC;
};




// ============================================================================
BlockDecomposition::BlockDecomposition (const std::shared_ptr<MeshGeometry> globalGeometry) :
globalGeometry (globalGeometry)
{
    // The geometry instance needs to be cartesian.
    assert (dynamic_cast<const CartesianMeshGeometry*> (globalGeometry.get()));

    auto world = MpiCommunicator::world();
    communicator = world.createCartesian (3, globalGeometry->fleshedOutAxes());
    communicator.onMasterOnly ( [&] ()
    {
        std::cout
        << "[BlockDecomposition] running on "
        << communicator.size() << " "
        << "MPI processes ["
        << communicator.getDimensions()[0] << " "
        << communicator.getDimensions()[1] << " "
        << communicator.getDimensions()[2] << "]\n";
    });
}

std::shared_ptr<BoundaryCondition> BlockDecomposition::createBoundaryCondition (
    std::shared_ptr<BoundaryCondition> physicalBC) const
{
    return std::make_shared<BlockDecomposedBC> (*this, physicalBC);
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
    localGeom->setCellsShape (localRegion.shape());
    localGeom->setLowerUpper (
        globalGeom->coordinateAtIndex (
        localRegion.lower[0] - 0.5,
        localRegion.lower[1] - 0.5,
        localRegion.lower[2] - 0.5),
        globalGeom->coordinateAtIndex (
        localRegion.upper[0] - 0.5,
        localRegion.upper[1] - 0.5,
        localRegion.upper[2] - 0.5));
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
    const int s = i - r > 0 ? i - r : 0; // number of partitions of size n
    const int t = i - s;                 // number of partitions of size n + 1
    return n * s + (n + 1) * t;
}

std::vector<double> BlockDecomposition::volumeAverageOverPatches (const std::vector<double>& diagnostics) const
{
    auto result = communicator.sum (diagnostics);

    for (auto& x : result)
    {
        x /= globalGeometry->meshVolume();
    }
    return result;
}
