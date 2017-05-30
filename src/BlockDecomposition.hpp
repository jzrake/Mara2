#ifndef BlockDecomposition_hpp
#define BlockDecomposition_hpp

#include "Mara.hpp"
#include "MPI.hpp"




class BlockDecomposition : public MeshDecomposition
{
public:

    BlockDecomposition() { }

    /**
    Constructor, initializes an MPI cartesian communicator.
    */
    BlockDecomposition (const std::shared_ptr<MeshGeometry> globalGeometry);

    /**
    Return a const reference to the communicator in use. Be sure the
    BlockDecomposition object lives longer than the returned reference.
    */
    const Cow::MpiCommunicator& getCommunicator() const { return communicator; }

    /**
    Return the shape of the global domain.
    */
    Cow::Shape getGlobalShape() const;

    /**
    Return a MeshGeometry::MeshGeometry object which corresponds to a patch in
    the global mesh. The patch index is assumed to be the MPI rank of this
    processes in the cartesian communicator.
    */
    std::shared_ptr<MeshGeometry> decompose() const;

    /**
    Return a MeshGeometry::PatchIndex object associated with this process's
    MPI rank in the cartesian communicator.
    */
    MeshGeometry::PatchIndex getPatchIndex() const;

    /**
    Return a Region object for the given patch index, with respect to the
    global mesh.
    */
    Cow::Region getPatchRegion (MeshGeometry::PatchIndex index) const;

    /**
    Convenience function for getPatchRegion (getPatchIndex())
    */
    Cow::Region getPatchRegion() const;

private:
    int partition (int numElements, int numPartitions, int whichPartition) const;
    int startIndex (int numElements, int numPartitions, int whichPartition) const;
    const std::shared_ptr<MeshGeometry> globalGeometry;
    Cow::MpiCartComm communicator;
};

#endif
