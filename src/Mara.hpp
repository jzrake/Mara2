#ifndef Mara_hpp
#define Mara_hpp


#include <memory>
#include <vector>
#include <string>
#include <functional>
#include "Array.hpp"
#include "Logger.hpp"
#include "Matrix.hpp"
#include "UnitVector.hpp"
#include "Variant.hpp"

#define MARA_NUM_FIELDS 10




/**
These are algorithm classes that may have inter-dependencies.
*/
class BoundaryCondition;
class ConservationLaw;
class ConstrainedTransport;
class FieldOperator;
class IntercellFluxScheme;
class ParticleData;
class MeshData;
class MeshGeometry;
class MeshOperator;
class RiemannSolver;
class SolutionScheme;
class SubProgram;
class TaskScheduler;
class TimeSeriesManager;




/**
These classes are here to support dependency injection.
*/
class MayUseBoundaryCondition    { public: virtual void setBoundaryCondition    (std::shared_ptr<BoundaryCondition>)    {} };
class MayUseConservationLaw      { public: virtual void setConservationLaw      (std::shared_ptr<ConservationLaw>)      {} };
class MayUseConstrainedTransport { public: virtual void setConstrainedTransport (std::shared_ptr<ConstrainedTransport>) {} };
class MayUseFieldOperator        { public: virtual void setFieldOperator        (std::shared_ptr<FieldOperator>)        {} };
class MayUseMeshOperator         { public: virtual void setMeshOperator         (std::shared_ptr<MeshOperator>)         {} };
class MayUseIntercellFluxScheme  { public: virtual void setIntercellFluxScheme  (std::shared_ptr<IntercellFluxScheme>)  {} };
class MayUseMeshGeometry         { public: virtual void setMeshGeometry         (std::shared_ptr<MeshGeometry>)         {} };
class MayUseRiemannSolver        { public: virtual void setRiemannSolver        (std::shared_ptr<RiemannSolver>)        {} };
class MayUseLogger
{
public:
    MayUseLogger() : logger (new Logger) {}
    void setLogger (std::shared_ptr<Logger> loggerToUse) { logger = loggerToUse; }
    std::shared_ptr<Logger> getLogger() { return logger; }
protected:
    std::shared_ptr<Logger> logger;
};




using StateArray = std::array<double, MARA_NUM_FIELDS>;
using InitialDataFunction = std::function<std::vector<double> (double x, double y, double z)>;
using SourceTermsFunction = std::function<StateArray (double x, double y, double z, StateArray primitive)>;
// using SourceTermsGeometry = std::function<StateArray (std::array<std::array<double, 2>, 3>, StateArray primitive)>;

using AreaElement = std::array<double, 3>;
using Coordinate = std::array<double, 3>;
enum class MeshLocation { vert, edge, face, cell };
enum class MeshBoundary { left, right };




class SubProgram
{
public:
    virtual int run (int argc, const char* argv[]) = 0;
};




class SolutionScheme :
public MayUseBoundaryCondition,
public MayUseFieldOperator,
public MayUseMeshOperator,
public MayUseIntercellFluxScheme
{
public:
    virtual int getStencilSize() const = 0;
    virtual void advance (MeshData& solution, double dt) const = 0;
    virtual void advance (MeshData& solution, ParticleData& particles, double dt) const {};
};




class SimulationStatus
{
public:
    SimulationStatus();
    void print (std::ostream& stream);
    void update (const Variant::NamedValues& values);
    Variant::NamedValues pack() const;
    double simulationTime = 0.0;
    int simulationIter = 0;
    long totalCellsInMesh = 0;
    double wallMinutes = 0.0;
};




/**
This class represents a decomposition of a global mesh into patches. There may
be one patch per MPI process, or the decomposition may be more general.
Currently there is only one implementation, which is called BlockDecomposition.
*/
class MeshDecomposition
{
public:
    virtual std::shared_ptr<MeshGeometry> decompose() const = 0;
};




class MeshGeometry
{
public:
    using Index = Cow::Index;

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
    Assign a patch identifier to this mesh. A mesh decomposition object may
    use the patch index to identify this mesh as a subset of a larger
    composite mesh. By default, the patch index is all zeros.
    */
    void assignPatchIndex (PatchIndex newPatchIndex);

    /**
    Set the starting index of this patch in a global geometry. This method
    should only be called by a MeshDecomposition object.
    */
    void assignStartIndex (Index index);

    /**
    Return the starting index of this patch in a global geometry.
    */
    Index getStartIndex() const;

    /**
    Return the index of this patch in a global composite mesh.
    */
    PatchIndex getPatchIndex() const;

    /**
    Convenience function which calls the derived class method with int
    arguments.
    */
    Coordinate coordinateAtIndex (Cow::Index index) const;

    /**
    Derived classes override this to set mesh resolution.
    */
    virtual void setCellsShape (Cow::Shape) = 0;

    /**
    Derived classes override this to set domain limits.
    */
    virtual void setLowerUpper (Coordinate, Coordinate) {}

    /**
    Return an object that describes the number of cells (volumes) contained in
    the mesh. For cartesian topology, this is just the number of cells in each
    direction. This is not meant to include guard zone regions as may required
    by various solvers.
    */
    virtual Cow::Shape cellsShape() const = 0;

    /**
    Return the index of the grid cell containing the given coordinate.
    */
    virtual Cow::Index indexAtCoordinate (Coordinate x) const = 0;

    /**
    Return the coordinates associated to a given index location in the mesh.
    Integer values are generally used for cell centroids.
    */
    virtual Coordinate coordinateAtIndex (double i, double j, double k) const = 0;

    /**
    Return the total number of cells in the mesh.
    */
    virtual unsigned long totalCellsInMesh() const = 0;

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
    Return the unit vector for the given face normal.
    */
    virtual UnitVector faceNormal (int i, int j, int k, int axis) const = 0;

    /**
    Return the length of the edge at the given index and axis.
    */
    virtual double edgeLength (int i, int j, int k, int axis) const = 0;

    /**
    Return the unit vector parallel to the edge at the given index and axis.
    */
    virtual UnitVector edgeVector (int i, int j, int k, int axis) const = 0;

    /**
    Return the volume of the cell at the given index.
    */
    virtual double cellVolume (int i, int j, int k) const = 0;

    /**
    Return the total volume of this grid patch.
    */
    virtual double meshVolume() const = 0;

    /**
    Return a 1D array of coordiate values along the given axis. This method
    only makes sense for logically cartesian, grid-aligned meshes.
    */
    virtual Cow::Array getPointCoordinates (int axis) const = 0;

    /**
    Return a copy of this mesh geometry as a new object.
    */
    virtual std::shared_ptr<MeshGeometry> duplicate() const = 0;

    /**
    Return a string identifier indicating the type of mesh represented.
    */
    virtual std::string getType() const = 0;

private:
    Index startIndex;
    PatchIndex patchIndex;
};




class BoundaryCondition :
public MayUseMeshGeometry,
public MayUseConservationLaw
{
public:
    /**
    This is the only method that needs to be implemented by BoundaryCondition
    derived classes.

    @param A                The full array of data on which to apply the
                            boundary condition.

    @param location         MeshLocation flag, whether the supplied data
                            exists on mesh verts, edges, faces, or cells.
                            Implementations should throw std::logic_error if
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
    A convenience function that calls the derived class's apply on both sides
    of each axis with size > 1 and the given boundaryShape used for the
    numGuard argument. Location is cells.
    */
    void applySimple (Cow::Array& A, Cow::Shape boundaryShape) const;

    /**
    To utilize a callback function, derived classes may override this function.
    */
    virtual void setBoundaryValueFunction (InitialDataFunction) {}

    /**
    Derived classes should override this function to indicate which of their
    axes are periodic. The result may be used, for example, by the patch
    (internal) BC object generated by a MeshDecomposition in determining
    whether it would be correct to apply physics BC's on a given patch.
    */
    virtual bool isAxisPeriodic (int axis) { return false; }

    virtual void setSimulationTime (double newSimulationTime) {}
};




class ConstrainedTransport :
public MayUseMeshGeometry,
public MayUseBoundaryCondition
{
public:
    virtual void assignCellCenteredB (Cow::Array) = 0;
    virtual Cow::Array computeMonopole (MeshLocation) const = 0;
    virtual Cow::Array computeCurrent (MeshLocation) const = 0;
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
        std::array<double, MARA_NUM_FIELDS> P; // Primitive quantities
        std::array<double, MARA_NUM_FIELDS> U; // Conserved densities
        std::array<double, MARA_NUM_FIELDS> F; // Fluxes in given direction
        std::array<double, MARA_NUM_FIELDS> A; // Eigenvalues
        Cow::Matrix L; // Left eigenvector matrix
        Cow::Matrix R; // Right eigenvector matrix
        int numFields = 0;
        int healthFlag = 0;
    };

    struct Request
    {
        Request();
        Request oriented (const UnitVector& nhat) const;
        bool getPrimitive;
        bool getConserved;
        bool getFluxes;
        bool getEigenvalues;
        AreaElement areaElement;
    };

    class StateFailure : public std::exception
    {
    public:
        StateFailure (const State& failedState);
        const char* what() const noexcept override;
        void setZoneIndex (Cow::Index I);
    private:
        void updateWhatMessage();
        std::string whatMessage;
        Cow::Index zoneIndex;
        State failedState;
    };

    using StateVector = std::vector<State>;

    virtual void setPressureFloor (double) {}
    virtual void setGammaLawIndex (double) {}
    virtual void setCoolingRate (double) {}
    virtual void setAdvectionSpeed (double, double, double) {}

    /**
    Derived classes may override this method to return a rich set of derived
    state quantities, such as entropy, Mach number, etc.
    */
    virtual std::vector<double> makeDiagnostics (const State& state) const
    {
        return std::vector<double>();
    }

    /**
    Return the names associated with each of the fields returned by the
    makeDiagnostics function.
    */
    virtual std::vector<std::string> getDiagnosticNames() const
    {
        return std::vector<std::string>();
    }

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
    This method may be called by the solution scheme, adjusting the time
    derivative L to account for geometrical or other source terms. If
    overridden, this method must *add* to (or subtract from) L, not replace
    its values.
    */
    virtual void addSourceTerms (const Cow::Array& P, Cow::Array& L) const {}

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
    Convenience function to return the field variable names. The vector is
    ordered the same way as the field indexes.
    */
    std::vector<std::string> getPrimitiveNames() const;

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




class IntercellFluxScheme : public MayUseRiemannSolver
{
public:
    struct FaceData
    {
    public:
        Cow::Array stencilData;
        AreaElement areaElement;
        std::shared_ptr<ConservationLaw> conservationLaw;
    };

    virtual void setPlmTheta (double plmTheta) {}
    virtual ConservationLaw::State intercellFlux (const FaceData&) const = 0;
    virtual int getStencilSize() const = 0;
};




int maraMainLoop (
    SimulationStatus& status,
    std::function<double ()> timestep,
    std::function<bool ()> condition,
    std::function<void (double)> advance,
    TaskScheduler& scheduler,
    Logger& logger);

#endif
