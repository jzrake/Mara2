#ifndef MeshData_hpp
#define MeshData_hpp

#include "Array.hpp"
#include "Mara.hpp"




/**
This class contains the actual solution data for a grid patch. It provides get
and assign methods which operate only on the non-boundary elements, but it
also exposes the raw data arrays publicly to simplify the interface. Region
objects are private, so that users are not aware of how this class retrieves
and assigns to its interior data. This way, if non-cartesian meshes are used
in the future, the get / assign methods may be overridden in a derived class.
*/
class MeshData
{
public:
	using Array = Cow::Array;
    using Shape = Cow::Shape;
    using Shape3D = Cow::Shape3D;
    using Region = Cow::Region;

    enum { includeGuard = 1 };

	MeshData (Shape baseShape, Shape boundaryShape, int numComponents);

    /**
    Set the index of velocity field in primitive data (-1 by default, indicating
    no velocity field data).
    */
    void setVelocityIndex (int velocityInd) { velocityIndex = velocityInd; }

    /**
    Set the index of magnetic field in primitive data (-1 by default, indicating
    no magnetic field data).
    */
    void setMagneticIndex (int magneticInd) { magneticIndex = magneticInd; }

    /**
    Assign a new set of primitive variables to the interior region.
    */
    void assignPrimitive (Array newP);

    /**
    Assign a new set of magnetic field variables to the interior region of
    either mesh or faces.
    */
    void assignMagneticField (Array newB, MeshLocation location, int flags=0);

    /**
    Allocate storage for diagnostic fields with given names.
    */
    void allocateDiagnostics (std::vector<std::string> diagnosticFields);

    /**
    Assign the given array to the diagnostic field at the given index.
    */
    void assignDiagnostic (Array newD, int index, int flags=0);

    /**
    Return a reference to the updateable region in the primitive variable
    array. All fields are returned, but if fieldIndex is non-negative then
    only that component of the primitive data is returned.
    */
    Array::Reference getPrimitive (int fieldIndex=-1, int flags=0);

    /**
    Return three contiguous components of the primitive data, starting at the
    given index.
    */
    Array::Reference getPrimitiveVector (int fieldIndex, int flags=0);

    /**
    Get the magnetic field at the given location (either cell or face).
    */
    Array::Reference getMagneticField (MeshLocation location, int flags=0);

    /**
    Get the velocity field at the given location (must be cell).
    */
    Array::Reference getVelocityField (MeshLocation location, int flags=0);

    /**
    Get the electric field values at the given location (must be cell). The
    velocity and magnetic field indexes must both be already set by the user.
    */
    Array getElectricField (MeshLocation location, int flags=0);

    /**
    Return an array indicating the health of zones. 0 means no errors occured
    there.
    */
    Array::Reference getZoneHealth(int flags=0);

    /**
    Return the diagnostic field with the given index.
    */
    Array::Reference getDiagnostic (int index, int flags=0);

    /**
    Return the number of diagnostic fields.
    */
    int getNumDiagnostics() const;

    /**
    Return the name of one of the diagnostic fields.
    */
    std::string getDiagnosticName (int index) const;

    /**
    Return the shape of the boundary zones. This is the number of cells, along
    each axis, by which the boundary region extends inwards from the array
    extent. It has the same value as the object with the same name that was
    passed to the constructor.
    */
    Shape3D getBoundaryShape() const;

    void applyBoundaryCondition (BoundaryCondition& bc);

    Array P; /**< Primitive variable array */
    Array B; /**< Face-centered magnetic field */
    Array Z; /**< Zone health array */
    Array D; /**< Diagnostics array (not allocated by default) */

private:
    Region getRegion (int flags) const;
    int velocityIndex;
    int magneticIndex;
    Region interior;
    Shape boundaryShape;
    std::vector<std::string> diagnosticFieldNames;
};

#endif
