#ifndef EulerEquation_hpp
#define EulerEquation_hpp

#include "Mara.hpp"




class ScalarAdvection : public ConservationLaw
{
public:
    ScalarAdvection (double waveSpeed=1.0);
    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    int getNumConserved() const override;
    int getIndexFor (VariableType type) const override;
    std::string getPrimitiveName (int fieldIndex) const override;
private:
    double waveSpeed;
};



class NewtonianHydro : public ConservationLaw
{
public:
    NewtonianHydro();
    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    int getNumConserved() const override;
    int getIndexFor (VariableType type) const override;
    std::string getPrimitiveName (int fieldIndex) const override;
private:
    double gammaLawIndex;
};



class NewtonianMHD : public ConservationLaw
{
public:
    NewtonianMHD();
    
    /**
    When the pressure floor value is set to something positive, then
    negative pressures are dealt with by setting p = pressureFloor * density.
    By default no floor is used.
    */
    void setPressureFloor (double floorValueToUse) override;

    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    int getNumConserved() const override;
    int getIndexFor (VariableType type) const override;
    std::string getPrimitiveName (int fieldIndex) const override;
private:
    double gammaLawIndex;
    double pressureFloor;
};

#endif
