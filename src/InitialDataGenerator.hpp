#ifndef InitialDataGenerator_hpp
#define InitialDataGenerator_hpp

#include "Mara.hpp"




class InitialDataGenerator
{
public:
    Cow::Array generatePrimitive (InitialDataFunction F, const MeshGeometry& geometry) const;
    Cow::Array generateMagnetic (InitialDataFunction A, const MeshGeometry& geometry) const;
};

#endif
