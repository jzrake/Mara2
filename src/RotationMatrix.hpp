#ifndef RotationMatrix_hpp
#define RotationMatrix_hpp

#include "UnitVector.hpp"




class RotationMatrix
{
public:
    static RotationMatrix aboutX (double angle);
    static RotationMatrix aboutY (double angle);
    static RotationMatrix aboutZ (double angle);
    UnitVector operator* (const UnitVector& nhat);
private:
    RotationMatrix() {}
    double M[3][3];    
};


#endif
