#ifndef Mara_hpp
#define Mara_hpp

#include <memory>
#include <vector>
#include <string>
#include <functional>
#include "Array.hpp"
#include "Matrix.hpp"
#include "Reconstruction.hpp"




class BoundaryCondition;
class ConservationLaw;
class ConstrainedTransport;
class IntercellFluxScheme;
class MaraSession;
class MeshGeometry;
class MeshDecomposition;
class RiemannSolver;
class SimulationSetup;
class SimulationStatus;

using InitialDataFunction = std::function<std::vector<double> (double, double, double)>;
using AreaElement = std::array<double, 3>;



class MaraSession
{
public:
    int launch (SimulationSetup& setup);
};




class SimulationSetup
{
public:
    SimulationSetup();

    // Run description
    double finalTime;
    double checkpointInterval;
    double vtkOutputInterval;
    double cflParameter;
    int rungeKuttaOrder;
    bool vtkUseBinary;
    std::string outputDirectory;
    std::string runName;
    std::string luaScript;
    std::string restartFile;
    
    // Solver options
    bool disableCT;

    // Algorithms
    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<IntercellFluxScheme> intercellFluxScheme;
    std::shared_ptr<ConstrainedTransport> constrainedTransport;
    std::shared_ptr<BoundaryCondition> boundaryCondition;
    std::shared_ptr<RiemannSolver> riemannSolver;
    InitialDataFunction initialDataFunction;
    InitialDataFunction vectorPotentialFunction;
};




class SimulationStatus
{
public:
    SimulationStatus();
    double simulationTime;
    int simulationIter;
    int checkpointsWrittenSoFar;
    int vtkOutputsWrittenSoFar;
};




/**
This class represents a decomposition of a global mesh into patches. There may
be one patch per MPI process, or the decomposition may be more general.
Currently there is only one implementation, which is a block decomposition.
*/
class MeshDecomposition
{
public:
    std::shared_ptr<MeshGeometry> decompose() const;
};




class MeshGeometry
{
public:

    /**
    A type to represent a three-component coordinate vector.
    */
    using Coordinate = std::array<double, 3>;

    /**
    A type used to identify a patch in a global composite mesh.
    */
    using PatchIndex = std::array<int, 5>;

    /** Default constructor, initializes private data members. */
    MeshGeometry();

    /**
    A convenience method that returns a vector of 3 bool's, indicating which
    of the axes have dimension greater than 1.
    */
    std::vector<bool> fleshedOutAxes() const;

    /**
    Assign a patch identifier to this mesh. A block decomposition object may
    use the patch index to identify this mesh as a subset of a larger
    composite mesh. By default, the patch index is all zeros.
    */
    void assignPatchIndex (PatchIndex newPatchIndex);

    /**
    Return the index of this patch in a global composite mesh.
    */
    PatchIndex getPatchIndex() const;

    /**
    Return an object that describes the number of cells (volumes) contained in
    the mesh. For cartesian topology, this is just the number of cells in each
    direction. This is not meant to include guard zone regions as may required
    by various solvers.
    */
    virtual Cow::Shape cellsShape() const = 0;

    /**
    Return the total number of cells in the mesh.
    */
    virtual unsigned long totalCellsInMesh() const = 0;

    /**
    Return the coordinates assocaited to a given index location in the mesh.
    Integer values are generally used for cell centroids.
    */
    virtual Coordinate coordinateAtIndex (double i, double j, double k) const = 0;

    /**
    Return the linear dimension of the cell with given index along a given axis.
    */
    virtual double cellLength (int i, int j, int k, int axis) const = 0;

    /**
    Return the surface area of the face at the given index and axis. Indexing
    convention is that a face with index i sits to the left of the cell with
    index i (so that face and cell 0 are both left-most along a given axis).
    */
    virtual double faceArea (int i, int j, int k, int axis) const = 0;

    /**
    Return the volume of the cell at the given index.
    */
    virtual double cellVolume (int i, int j, int k) const = 0;

private:
    PatchIndex patchIndex;
};




class BoundaryCondition
{
public:
    enum class MeshLocation { vert, edge, face, cell };
    enum class MeshBoundary { left, right };

    /**
    This is the only method that needs to be implemented by BoundaryCondition
    derived classes.

    @param A                The full array of data on which to apply the
                            boundary condition.

    @param location         MeshLocation flag, whether the supplied data
                            exists on mesh verts, edges, faces, or cells.
                            Implementations should throw std::runtime_error if
                            they cannot apply boundary conditions at the
                            requested mesh location.

    @param boundary         Which boundary of the mesh needs to be filled,
                            either left or right.

    @param axis             Which axis of the array to apply boundary
                            conditions to. For a rectilinear mesh this is 0,
                            1, or 2.

    @param numGuard         Width of the guard zone region the caller expects
                            to be replaced with valid data. The caller
                            promises that A contains a region of valid data,
                            of at least the same width, adjacent to the guard
                            zone region.
    */
    virtual void apply (
        Cow::Array& A,
        MeshLocation location,
        MeshBoundary boundary,
        int axis,
        int numGuard) const = 0;

    /**
    Assign a mesh geometry instance. Derived classes may re-implement
    this method and keep a shared pointer to the geometry if they need it.
    */
    virtual void setMeshGeometry (std::shared_ptr<MeshGeometry> geometry) {}

    /**
    Assign a conservation law instance. Derived classes may re-implement
    this method and keep a shared pointer to the geometry if they need it.
    */
    virtual void setConservationLaw (std::shared_ptr<ConservationLaw> law) {}
};




class ConstrainedTransport
{
public:
    enum class MeshLocation { vert, edge, face, cell };
    virtual void assignCellCenteredB (Cow::Array) = 0;
    virtual Cow::Array computeMonopole (MeshLocation) const = 0;
    virtual void setMeshGeometry (std::shared_ptr<MeshGeometry>) {}
    virtual void setBoundaryCondition (std::shared_ptr<BoundaryCondition>) {}
};




class ConservationLaw
{
public:
    enum class VariableType
    {
        density,
        velocity,
        pressure,
        magnetic,
    };

    struct State
    {
        std::array<double, 8> P; // Primitive quantities
        std::array<double, 8> U; // Conserved densities
        std::array<double, 8> F; // Fluxes in given direction
        std::array<double, 8> A; // Eigenvalues
        Cow::Matrix L; // Left eigenvector matrix
        Cow::Matrix R; // Right eigenvector matrix
        int healthFlag;
    };

    struct Request
    {
        Request();
        bool getPrimitive;
        bool getConserved;
        bool getFluxes;
        bool getEigenvalues;
        AreaElement areaElement;
    };

    class StateFailure : public std::exception
    {
    public:
        StateFailure (const State& failedState) : failedState (failedState) {}
        const char* what() const noexcept override;
        Cow::Index zoneIndex;
        State failedState;
    };

    using StateVector = std::vector<State>;

    /**
    Derived classes may override this method to implement a pressure floor.
    */
    virtual void setPressureFloor (double) {}

    /**
    Generate a state from the given information request, and a double pointer
    of conserved quantities. U must point to valid conserved quantity data
    with numConserved consecutive doubles.
    */
    virtual State fromConserved (const Request& request, const double* U) const = 0;

    /**
    Generate a state from the given information request, and a double pointer
    of primitive quantities. P must point to valid primitive quantity data
    with numConserved consecutive doubles.
    */
    virtual State fromPrimitive (const Request& request, const double* P) const = 0;

    /**
    Get the number of primitive and conserved quantities for this conservation
    law. For example, the Euler equation has 5 (density, total energy, and
    three components of momentum), and MHD has 8 (three magnetic field fluxes
    in addition). If there are passive scalars then those will increase that
    number.
    */
    virtual int getNumConserved() const = 0;

    /**
    Return the index at which the given primitive variable data resides. If
    this conservation law does not have that variable type, the return value
    will be -1.
    */
    virtual int getIndexFor (VariableType type) const = 0;

    /**
    Return a name for the primitive variable at the given index.
    */
    virtual std::string getPrimitiveName (int fieldIndex) const = 0;

    /**
    Return the numeric average the primitive quantities from the two given
    states, and return a state computed with from that average and the given
    request.
    */
    State averageStates (const Request& request, const State& L, const State& R) const;

    /**
    A helper function that returns a vector of states by calling fromPrimitive
    numStates times, where P.shape() = [numStates, numConserved].
    */
    StateVector fromPrimitive (const Request& request, const Cow::Array& P) const;

    /**
    A helper function that returns the maximum absolute value of all
    eigenvalues for the given state.
    */
    double maxEigenvalueMagnitude (const State& state) const;

    /**
    A helper function that the maximum eigenvalue magnitude, over a vector of
    states. This is useful for Lax-Friedrichs flux splitting.
    */
    double maxEigenvalueMagnitude (const StateVector& states) const;
};




class RiemannSolver
{
public:
    using State = ConservationLaw::State;
    virtual State solve (const State& L, const State& R, AreaElement dA) const = 0;
};




class IntercellFluxScheme
{
public:
    struct FaceData
    {
    public:
        Cow::Array stencilData;
        AreaElement areaElement;
        std::shared_ptr<ConservationLaw> conservationLaw;
    };

    virtual ConservationLaw::State intercellFlux (const FaceData&) const = 0;
    virtual int getStencilSize() const = 0;
};

#endif
