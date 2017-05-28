#ifndef BlockDecomposition_hpp
#define BlockDecomposition_hpp

#include "Mara.hpp"
#include "MPI.hpp"




class BlockDecomposition : public MeshDecomposition
{
public:
    BlockDecomposition (std::shared_ptr<MeshGeometry> globalGeometry);
    std::shared_ptr<MeshGeometry> decompose () const;
private:
    std::shared_ptr<MeshGeometry> globalGeometry;
    Cow::MpiCartComm topology;
};

#endif
