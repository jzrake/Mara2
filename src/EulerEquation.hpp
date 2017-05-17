#ifndef EulerEquation_hpp
#define EulerEquation_hpp

#include "Mara.hpp"




class EulerEquation : public ConservationLaw
{
public:
    EulerEquation();
    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    int getNumConserved() const override;
    std::string getPrimitiveName (int fieldIndex) const override;
private:
    double gammaLawIndex;
};

#endif
