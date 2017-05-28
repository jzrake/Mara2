#include "CartesianMeshGeometry.hpp"




// ============================================================================
std::vector<bool> MeshGeometry::fleshedOutAxes() const
{
    auto shape = domainShape();
    return { shape[0] > 1, shape[1] > 1, shape[2] > 1};
}
