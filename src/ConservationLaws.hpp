#ifndef EulerEquation_hpp
#define EulerEquation_hpp

#include "Mara.hpp"




class ScalarAdvection : public ConservationLaw
{
public:
    ScalarAdvection();
    void setAdvectionSpeed (double u, double v, double w) override { waveSpeed = u; }
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
    void setGammaLawIndex (double gm) override { gammaLawIndex = gm; }
    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    int getNumConserved() const override;
    int getIndexFor (VariableType type) const override;
    std::string getPrimitiveName (int fieldIndex) const override;
    std::vector<double> makeDiagnostics (const State& state) const override;
    std::vector<std::string> getDiagnosticNames() const override;
private:
    double gammaLawIndex;
};



class NewtonianMHD : public ConservationLaw
{
public:
    NewtonianMHD();

    void setCoolingRate (double) override;
    void setGammaLawIndex (double) override;

    /**
    When the pressure floor value is set to something positive, then
    negative pressures are dealt with by setting p = pressureFloor * density.
    By default no floor is used.
    */
    void setPressureFloor (double) override;

    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    void addSourceTerms (const Cow::Array& P, Cow::Array& L) const override;
    int getNumConserved() const override;
    int getIndexFor (VariableType type) const override;
    std::string getPrimitiveName (int fieldIndex) const override;
    std::vector<double> makeDiagnostics (const State& state) const override;
    std::vector<std::string> getDiagnosticNames() const override;
private:
    int numAuxiliary;
    double gammaLawIndex;
    double pressureFloor;
    double coolingRate;
};



class RelativisticMHD : public ConservationLaw
{
public:
    RelativisticMHD();
    void setGammaLawIndex (double) override;

    /**
    When the pressure floor value is set to something positive, then
    negative pressures are dealt with by setting p = pressureFloor * density.
    By default no floor is used.
    */
    void setPressureFloor (double) override;

    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    void addSourceTerms (const Cow::Array& P, Cow::Array& L) const override;
    int getNumConserved() const override;
    int getIndexFor (VariableType type) const override;
    std::string getPrimitiveName (int fieldIndex) const override;
    std::vector<double> makeDiagnostics (const State& state) const override;
    std::vector<std::string> getDiagnosticNames() const override;
private:
    double gammaLawIndex;
    double pressureFloor;
};
#endif
