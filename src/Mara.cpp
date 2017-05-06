#include <iostream>
#include <cmath>
#include "FluxConservativeSystem.hpp"
#include "HDF5.hpp"

#include "sol.hpp"


class CartesianMeshGeometry : public MeshGeometry
{
public:
    CartesianMeshGeometry()
    {
        shape = {{128, 1, 1, 1, 1}};
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
        waveSpeed = -1.0;
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


void testLua (std::string filename)
{
    sol::state lua;
    lua.open_libraries (sol::lib::base);

    lua.script_file (filename);


    for (auto entry : lua.globals())
    {
        std::cout << entry.first.as<std::string>() <<": " << entry.second.as<std::string>() << std::endl;
    }



    // std::cout << lua.get<std::string>("run_name") << std::endl;

    // for (auto entry : lua.get<sol::table>("resolution"))
    // {
    //     auto F = entry.second;
    //     std::cout << F.is<int>() << std::endl;
    // }

    // auto res = lua.get<sol::table>("resolution");
    // std::cout << res.size() << std::endl;



    //std::cout << initial_data() << std::endl;

    // auto initialDataFunction = [] (double x, double y, double z)
    // {
    //     auto S = ConservationLaw::State();
    //     S.P.push_back (std::exp (-x * x / 0.01));
    //     return S;
    // };
}


int main(int argc, const char* argv[])
{
    using namespace Cow;

    if (argc == 1)
    {
        std::cout << "usage: mara config.lua\n";
        return 0;
    }

    sol::state lua;
    lua.open_libraries (sol::lib::base);
    lua.script_file (argv[1]);

    auto makeInitialDataFunction = [] (sol::function& func)
    {
        auto F = [=] (double x, double y, double z)
        {
            auto S = ConservationLaw::State();
            sol::table result = func();

            for (int n = 0; n < result.size(); ++n)
            {
                S.P.push_back (result[n + 1]);
            }
            return S;
        };
        return F;
    };

    sol::function initial_data = lua["initial_data"];
    auto initialDataFunction = makeInitialDataFunction (initial_data);


    return 0;


    auto geometry = new CartesianMeshGeometry;
    geometry->shape[0] = lua["resolution"][1].get_or (128);
    geometry->shape[1] = lua["resolution"][2].get_or (1);
    geometry->shape[2] = lua["resolution"][3].get_or (1);


    auto system = FluxConservativeSystem (geometry, new ScalarAdvection, new ScalarUpwind);
    system.setInitialData (initialDataFunction);


    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0000.h5", "w");
        file.write ("primitive", P);
    }

    double t = 0.0;
    double dt = 0.0025;

    while (t < 0.5)
    {
        std::cout << "t=" << t << std::endl;
        system.computeIntercellFluxes();
        system.computeTimeDerivative();
        system.updateConserved (dt);
        system.recoverPrimitive();
        system.applyBoundaryConditions();
        t += dt;
    }

    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0001.h5", "w");
        file.write ("primitive", P);
    }

	return 0;
}
