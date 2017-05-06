#ifndef FluxConservativeSystem_hpp
#define FluxConservativeSystem_hpp

#include <functional>
#include "Array.hpp"




class MeshGeometry
{
public:
    using Coordinate = std::array<double, 3>;
    virtual Cow::Shape domainShape() const = 0;
    virtual Coordinate coordinateAtIndex (double i, double j, double k) const = 0;
    virtual double faceArea (int i, int j, int k, int axis) const = 0;
    virtual double cellVolume (int i, int j, int k) const = 0;
};




class ConservationLaw
{
public:
    struct State
    {
        std::vector<double> U; // Conserved densities
        std::vector<double> P; // Primitive quantities
        std::vector<double> F; // Fluxes in given direction
        std::vector<double> A; // Eigenvalues
    };

    struct Request
    {
        Request();
        bool getPrimitive;
        bool getConserved;
        bool getFluxes;
        bool getEigenvalues;
        double areaElement[3];
    };

    virtual State fromConserved (const Request& request, const double* U) const = 0;
    virtual State fromPrimitive (const Request& request, const double* P) const = 0;
    virtual int getNumConserved() const = 0;
};




class IntercellFluxEstimator
{
public:
    using StateVector = std::vector<ConservationLaw::State>;
    virtual ConservationLaw::State intercellFlux (StateVector stateVector) const = 0;
    virtual int getSchemeOrder() const = 0;
};




class FluxConservativeSystem
{
public:
    using InitialDataFunction = std::function<ConservationLaw::State (double, double, double)>;

    FluxConservativeSystem (
        MeshGeometry* meshGeometry,
        ConservationLaw* conservationLaw,
        IntercellFluxEstimator* intercellFluxEstimator);

    Cow::Array::Reference getPrimitive();

    void setInitialData (InitialDataFunction F);
    void computeIntercellFluxes();
    void computeTimeDerivative();
    void applyBoundaryConditions();
    void updateConserved (double dt, double rungeKuttaParameter=0);
    void recoverPrimitive();

private:
    int numDimensions;
    int numConserved;
    int schemeOrder;

    Cow::Shape domainShape;
    Cow::Shape domainShapeAndState;
    Cow::Region updateableRegion;

    Cow::Array U; // Conserved quantities
    Cow::Array P; // Primitive quantities
    Cow::Array L; // Time derivative of U
    Cow::Array F1; // Flux along axis 1
    Cow::Array F2; // Flux along axis 2
    Cow::Array F3; // Flux along axis 3

    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<IntercellFluxEstimator> intercellFluxEstimator;
};

#endif
