#include <iostream>
#include "BlockDecomposition.hpp"
#include "CartesianMeshGeometry.hpp"

using namespace Cow;




// ============================================================================
class BlockDecomposedBC : public BoundaryCondition
{
public:
    BlockDecomposedBC (const BlockDecomposition& block) : block (block) {}

    void apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const override
    {
        // if (P.size(0) > 1) applyToAxis (P, numGuard, 0);
        // if (P.size(1) > 1) applyToAxis (P, numGuard, 1);
        // if (P.size(2) > 1) applyToAxis (P, numGuard, 2);
    }

    void applyToCellCenteredB (Cow::Array& B, int numGuard) const override
    {
        // if (B.size(0) > 1) applyToAxis (B, numGuard, 0);
        // if (B.size(1) > 1) applyToAxis (B, numGuard, 1);
        // if (B.size(2) > 1) applyToAxis (B, numGuard, 2);
    }

    void applyToGodunovFluxes (Cow::Array& F, int numGuard, int axis) const override
    {
        // if (F.size (axis) > 1) applyToAxis (F, numGuard, axis);
    }

    void applyToAxis (Cow::Array& P, int numGuard, int axis) const
    {
        auto sendL = Region(); // region to send to the left
        auto sendR = Region(); // region to send to the right
        auto recvL = Region(); // region to receive from the left
        auto recvR = Region(); // region to receive from the right

        sendR.lower[axis] = -2 * numGuard;
        sendR.upper[axis] = -1 * numGuard;
        sendL.lower[axis] = +1 * numGuard;
        sendL.upper[axis] = +2 * numGuard;

        recvR.lower[axis] = -1 * numGuard;
        recvR.upper[axis] = -0;
        recvL.lower[axis] = +0;
        recvL.upper[axis] = +1 * numGuard;

        block.communicator.shiftExchange (P, axis, 'L', sendL, recvR);
        block.communicator.shiftExchange (P, axis, 'R', sendR, recvL);
    }
private:
    const BlockDecomposition& block;
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
        << communicator.size()
        << " MPI processes ["
        << communicator.getDimensions()[0] << " "
        << communicator.getDimensions()[1] << " "
        << communicator.getDimensions()[2] << "]\n";
    });
}

std::shared_ptr<BoundaryCondition> BlockDecomposition::createBoundaryCondition (
    std::shared_ptr<BoundaryCondition> physicalBC) const
{
    return std::make_shared<BlockDecomposedBC> (*this);
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
