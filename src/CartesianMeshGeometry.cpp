#include "CartesianMeshGeometry.hpp"




// ============================================================================
MeshGeometry::MeshGeometry()
{
    patchIndex[0] = 0;
    patchIndex[1] = 0;
    patchIndex[2] = 0;
    patchIndex[3] = 0;
    patchIndex[4] = 0;
}

std::vector<bool> MeshGeometry::fleshedOutAxes() const
{
    auto shape = cellsShape();
    return { shape[0] > 1, shape[1] > 1, shape[2] > 1};
}

void MeshGeometry::assignPatchIndex (PatchIndex newPatchIndex)
{
    patchIndex = newPatchIndex;
}

MeshGeometry::PatchIndex MeshGeometry::getPatchIndex() const
{
    return patchIndex;
}
