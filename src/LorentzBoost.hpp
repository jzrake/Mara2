#ifndef LorentzBoost_hpp
#define LorentzBoost_hpp


#include "FourVector.hpp"


// ============================================================================
class LorentzBoost
{
public:
    LorentzBoost();
    LorentzBoost (const FourVector& boostVector);

    LorentzBoost inverted();

    /**
    Return a four vector u' = L * u transformed by this matrix.
    */
    FourVector operator* (const FourVector& u) const;

private:
    double elements[4][4];
    FourVector boostVector;
};

#endif
