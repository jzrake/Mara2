#include <iostream>
#include "FluxConservativeSystem.hpp"




class CartesianMeshGeometry : public MeshGeometry
{
public:
    CartesianMeshGeometry()
    {
        shape = {{32, 1, 1, 1, 1}};
    }

    Cow::Shape domainShape() const override
    {
        return shape;
    }

    Coordinate coordinateAtIndex (double i, double j, double k) const override
    {
        return Coordinate ({{
            (i + 0.5) / shape[0],
            (j + 0.5) / shape[1],
            (k + 0.5) / shape[2]}});
    }

    double faceArea (int i, int j, int k, int axis) const override
    {
        return 1.0;
    }

    double cellVolume (int i, int j, int k) const override
    {
        return 1.0 / shape[0];
    }

    Cow::Shape shape;
};




class ScalarAdvection : public ConservationLaw
{
public:

    State fromConserved (const Request& request, const double* U) const override
    {
        State S;
        return S;
    }

    State fromPrimitive (const Request& request, const double* P) const override
    {
        State S;
        return S;
    }

    int getNumConserved() const override
    {
        return 1;
    }
};




class ScalarUpwind : public IntercellFluxEstimator
{
public:
    ConservationLaw::State intercellFlux (StateVector stateVector) const override
    {
        auto S = ConservationLaw::State();
        S.F.push_back (0);
        return S;
    }

    int getSchemeOrder() const override
    {
        return 1;
    }
};




int main()
{
    auto system = FluxConservativeSystem (new CartesianMeshGeometry, new ScalarAdvection, new ScalarUpwind);

    system.setInitialData ([] (double x, double y, double z)
        {
            std::cout << x << " " << y << " " << z << std::endl;
            auto S = ConservationLaw::State();
            S.P.push_back (0);
            return S;
        });

    system.computeIntercellFluxes();
    system.computeTimeDerivative();
    system.updateConserved (0.01);

	return 0;
}
