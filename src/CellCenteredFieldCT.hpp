#ifndef CellCenteredFieldCT_hpp
#define CellCenteredFieldCT_hpp

#include "Array.hpp"
#include "Mara.hpp"




class CellCenteredFieldCT
{
public:
    using Array = Cow::Array;
    void correctGodunovFluxes (Array& F, int magneticIndex) const;
    Array generateGodunovFluxes (const Array& F, int magneticIndex) const;
    Array monopole (Array& B, MeshLocation location) const;
    Array current (Array& B, MeshLocation location) const;
};

#endif