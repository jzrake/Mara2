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
    using Region = Cow::Region;

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
    void assignMagneticField (Array newB, MeshLocation location);

    /**
    Return a reference to the updateable region in the primitive variable
    array. All fields are returned, but if fieldIndex is non-negative then
    only that component of the primitive data is returned.
    */
    Array::Reference getPrimitive (int fieldIndex=-1);

    /**
    Return three contiguous components of the primitive data, starting at the
    given index.
    */
    Array::Reference getPrimitiveVector (int fieldIndex);

    /**
    Get the magnetic field at the given location (either cell or face).
    */
    Array::Reference getMagneticField (MeshLocation location);

    /**
    Return an array indicating the health of zones. 0 means no errors occured
    there.
    */
    Array::Reference getZoneHealth();

    Array P; /**< Primitive variable array */
    Array B; /**< Face-centered magnetic field */
    Array Z; /**< Zone health array */

private:
    int velocityIndex;
    int magneticIndex;
    Region interior;
};

#endif
