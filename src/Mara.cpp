#include <iostream>
#include <cmath>
#include "FluxConservativeSystem.hpp"
#include "HDF5.hpp"




class CartesianMeshGeometry : public MeshGeometry
{
public:
    CartesianMeshGeometry()
    {
        shape = {{256, 1, 1, 1, 1}};
    }

    Cow::Shape domainShape() const override
    {
        return shape;
    }

    Coordinate coordinateAtIndex (double i, double j, double k) const override
    {
        return Coordinate ({{
            -0.5 + (i + 0.5) / shape[0],
            -0.5 + (j + 0.5) / shape[1],
            -0.5 + (k + 0.5) / shape[2]}});
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
    ScalarAdvection()
    {
        waveSpeed = 1.0;
    }

    State fromConserved (const Request& request, const double* U) const override
    {
        double u = U[0];
        State S;
        S.P.push_back (u);
        S.U.push_back (u);
        S.A.push_back (waveSpeed);
        S.F.push_back (waveSpeed * u);
        return S;
    }

    State fromPrimitive (const Request& request, const double* P) const override
    {
        double u = P[0];
        State S;
        S.P.push_back (u);
        S.U.push_back (u);
        S.A.push_back (waveSpeed);
        S.F.push_back (waveSpeed * u);
        return S;
    }

    int getNumConserved() const override
    {
        return 1;
    }

private:
    double waveSpeed;
};




class ScalarUpwind : public IntercellFluxEstimator
{
public:
    ConservationLaw::State intercellFlux (StateVector stateVector) const override
    {
        auto S = ConservationLaw::State();
        const auto& L = stateVector[0];
        const auto& R = stateVector[1];

        if (L.A[0] > 0 && R.A[0] > 0)
        {
            S.F.push_back (L.F[0]);
        }
        else if (L.A[0] < 0 && R.A[0] < 0)
        {
            S.F.push_back (R.F[0]);
        }

        return S;
    }

    int getSchemeOrder() const override
    {
        return 1;
    }
};




int main()
{
    using namespace Cow;

    auto system = FluxConservativeSystem (new CartesianMeshGeometry, new ScalarAdvection, new ScalarUpwind);

    system.setInitialData ([] (double x, double y, double z)
        {
            auto S = ConservationLaw::State();
            S.P.push_back (std::exp (-x * x / 0.01));
            return S;
        });

    

    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0000.h5", "w");
        file.write ("primitive", P);
    }

    double t = 0.0;
    double dt = 0.0005;

    while (t < 0.05)
    {
        system.computeIntercellFluxes();
        system.computeTimeDerivative();
        system.updateConserved (dt);
        system.recoverPrimitive();
        t += dt;
    }

    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0001.h5", "w");
        file.write ("primitive", P);
    }

	return 0;
}
