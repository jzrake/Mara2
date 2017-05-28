#include "BlockDecomposition.hpp"

using namespace Cow;




BlockDecomposition::BlockDecomposition (std::shared_ptr<MeshGeometry> globalGeometry) :
globalGeometry (globalGeometry)
{
    auto world = MpiCommunicator::world();
    topology = world.createCartesian (3, globalGeometry->fleshedOutAxes());
}

std::shared_ptr<MeshGeometry> BlockDecomposition::decompose() const
{
    return globalGeometry;
}
