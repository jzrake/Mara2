#ifndef MeshData_hpp
#define MeshData_hpp

#include "Array.hpp"




class MeshData
{
public:
	using Array = Cow::Array;
	using Shape = Cow::Shape;
    using Region = Cow::Region;

	MeshData (Shape baseShape, Shape guardZone, int numComponents);

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
    Return an array indicating the health of corresponding zones. 0 means no
    errors occured there.
    */
    Array::Reference getZoneHealth();

private:
	Array P; // Primitive variable array
	Array B; // Cell-centered magnetic field
	Array Z; // Zone health array

    Region updateableRegion;
    Region magneticIndices;
};

#endif
