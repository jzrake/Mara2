#include <iostream>
#include <iomanip>
#include <cmath>
#include "Mara.hpp"
#include "Configuration.hpp"
#include "FluxConservativeSystem.hpp"
#include "HDF5.hpp"




// ============================================================================
SimulationSetup::SimulationSetup()
{
    finalTime = 1.0;
    checkpointInterval = 1.0;
    cflParameter = 0.25;
}




// ============================================================================
SimulationStatus::SimulationStatus()
{
    simulationTime = 0.0;
    simulationIter = 0;
    outputsWrittenSoFar = 0;
}




// ============================================================================
int main(int argc, const char* argv[])
{
    using namespace Cow;

    if (argc == 1)
    {
        std::cout << "usage: mara config.lua\n";
        return 0;
    }
    auto configuration = Configuration();
    auto setup = configuration.fromLuaFile (argv[1]);


    // Setup lines specific to conservation law problems
    auto system = FluxConservativeSystem (setup);
    system.setInitialData (setup.initialDataFunction);


    auto status = SimulationStatus();
    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0000.h5", "w");
        file.write ("primitive", P);
    }


    while (status.simulationTime < setup.finalTime)
    {
        double dt = setup.cflParameter * system.getCourantTimestep();

        std::cout << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        std::cout << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        std::cout << "dt=" << std::setprecision (2) << std::scientific << dt << "\n";

        system.computeIntercellFluxes();
        system.computeTimeDerivative();
        system.updateConserved (dt);
        system.recoverPrimitive();
        system.applyBoundaryConditions();

        status.simulationTime += dt;
        status.simulationIter += 1;
    }

    {
        auto P = system.getPrimitive();
        auto file = H5::File ("chkpt.0001.h5", "w");
        file.write ("primitive", P);
    }

	return 0;
}
