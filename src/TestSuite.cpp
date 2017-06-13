#define CATCH_CONFIG_RUNNER
#include "Catch.hpp"
#include "Array.hpp"

#include "TestSuite.hpp"

using namespace Cow;




// ============================================================================
#include "Stencil.hpp"

SCENARIO ("Stencil operations should behave correctly")
{
    GIVEN ("A stencil with shape [2, 3, 4]")
    {
        auto stencilShape = Shape();
        auto stencil = Stencil();
        auto operation = [&] (const Array& sd, Array& od) { stencilShape = sd.shape(); };
        stencil.setFootprintLower (-2,-3,-4);
        stencil.setFootprintUpper ( 2, 3, 4);
        stencil.setCodomainRank (4, 5);

        WHEN ("The input array is [12, 12, 12] and the codomain rank is [4, 5]")
        {
            auto source = Array (12, 12, 12);
            auto result = stencil.evaluate (operation, source);

            THEN ("The output array shape is [8, 6, 4, 4, 5]")
            {
                CHECK (result.size(0) == 8);
                CHECK (result.size(1) == 6);
                CHECK (result.size(2) == 4);
                CHECK (result.size(3) == 4);
                CHECK (result.size(4) == 5);
            }

            THEN ("The array of stencil data has shape [5, 7, 9, 1, 1]")
            {
                CHECK (stencilShape[0] == 5);
                CHECK (stencilShape[1] == 7);
                CHECK (stencilShape[2] == 9);
                CHECK (stencilShape[3] == 1);
                CHECK (stencilShape[4] == 1);
            }
        }
    }

    GIVEN ("A stencil with shape [2, 2, 0]")
    {
        auto stencilShape = Shape();
        auto stencil = Stencil();
        auto operation = [&] (const Array& sd, Array& od) { stencilShape = sd.shape(); };
        auto source = Array (12, 12, 1, 8);
        stencil.setFootprintLower (-2,-2, 0);
        stencil.setFootprintUpper ( 2, 2, 0);
        stencil.setCodomainRank (3, 1);

        WHEN ("The input array is [12, 12, 1, 8] and the codomain rank is [3, 1]")
        {
            auto result = stencil.evaluate (operation, source);

            THEN ("The output array shape is [8, 8, 1, 3, 1]")
            {
                CHECK (result.size(0) == 8);
                CHECK (result.size(1) == 8);
                CHECK (result.size(2) == 1);
                CHECK (result.size(3) == 3);
                CHECK (result.size(4) == 1);
            }

            THEN ("The array of stencil data has shape [5, 5, 1, 8, 1]")
            {
                CHECK (stencilShape[0] == 5);
                CHECK (stencilShape[1] == 5);
                CHECK (stencilShape[2] == 1);
                CHECK (stencilShape[3] == 8);
                CHECK (stencilShape[4] == 1);
            }
        }

        WHEN ("The stencil footprint is changed to [2, 2, 1]")
        {
            stencil.setFootprintLower (-2,-2,-1);
            stencil.setFootprintUpper ( 2, 2, 1);

            THEN ("An exception is raised because the source data shape is too small")
            {
                CHECK_THROWS_AS (stencil.evaluate (operation, source), std::logic_error);
            }
        }
    }

    GIVEN ("A second order finite difference stencil")
    {
        auto stencil = Stencil();
        auto x = Array (256);
        auto y = Array (256);

        for (int i = 0; i < 256; ++i)
        {
            x[i] = (i + 0.5) / 256;
            y[i] = std::sin (x[i]);
        }
        stencil.setCodomainRank (2, 1);

        WHEN ("The stencil is centered, of the type [-1, 0, +1]")
        {
            stencil.setFootprintLower (-1, 0, 0);
            stencil.setFootprintUpper ( 1, 0, 0);

            auto operation = [] (const Array& x0, const Array& y0, Array& deriv)
            {
                deriv[0] = (y0[2] - y0[0]) / (x0[2] - x0[0]);
                deriv[1] = std::cos (x0[1]);
            };
            auto deriv = stencil.evaluate (operation, x, y);

            THEN ("The resulting array is smaller by two data points")
            {
                CHECK (deriv.size(0) == 254);
            }

            THEN ("The derivative approximation is correct to one part in 10^5")
            {
                CHECK (deriv (  0, 0, 0, 0) == Approx (deriv (  0, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (128, 0, 0, 0) == Approx (deriv (128, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (253, 0, 0, 0) == Approx (deriv (253, 0, 0, 1)).epsilon (1e-5));
            }
        }

        WHEN ("The stencil is centered, of the type [0, +1]")
        {
            stencil.setFootprintLower (0, 0, 0);
            stencil.setFootprintUpper (1, 0, 0);

            auto operation = [] (const Array& x0, const Array& y0, Array& deriv)
            {
                deriv[0] = (y0[1] - y0[0]) / (x0[1] - x0[0]);
                deriv[1] = std::cos (0.5 * (x0[0] + x0[1]));
            };
            auto deriv = stencil.evaluate (operation, x, y);

            THEN ("The resulting array is smaller by one data point")
            {
                CHECK (deriv.size(0) == 255);
            }

            THEN ("The derivative approximation is correct to one part in 10^5")
            {
                CHECK (deriv (  0, 0, 0, 0) == Approx (deriv (  0, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (128, 0, 0, 0) == Approx (deriv (128, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (254, 0, 0, 0) == Approx (deriv (254, 0, 0, 1)).epsilon (1e-5));
            }
        }
    }
}




// ============================================================================
#include "Mara.hpp"
#include "CartesianMeshGeometry.hpp"
#include "BoundaryConditions.hpp"
#include "ConservationLaws.hpp"
#include "ConstrainedTransport.hpp"
#include "IntercellFluxSchemes.hpp"
#include "RiemannSolvers.hpp"

SCENARIO ("Mara session should launch if given minimal setup")
{
    auto session = MaraSession();
    auto setup = SimulationSetup();

    setup.finalTime               = 0.001;
    setup.checkpointInterval      = 0.0;
    setup.vtkOutputInterval       = 0.0;
    setup.boundaryCondition       = std::make_shared<PeriodicBoundaryCondition>();
    setup.conservationLaw         = std::make_shared<NewtonianHydro>();
    setup.riemannSolver           = std::make_shared<HlleRiemannSolver>();
    setup.meshGeometry            = std::make_shared<CartesianMeshGeometry>();
    setup.intercellFluxScheme     = std::make_shared<MethodOfLinesPlm>();
    setup.constrainedTransport    = std::make_shared<UniformCartesianCT>();
    setup.initialDataFunction     = [] (double x, double, double)
    {
        return std::vector<double> { 1.0 + 0.1 * std::sin (2 * M_PI * x), 0., 0., 0., 1. };
    };

    session.getLogger()->setLogToNull();
    CHECK (session.launch (setup).simulationIter == 1);
}




// ============================================================================
#include "HDF5.hpp"
#include "TimeSeriesManager.hpp"

SCENARIO ("Time series manager should behave reasonably")
{
    GIVEN ("A time series manager populated with data series, with unequal lengths")
    {
        TimeSeriesManager timeSeriesManager;
        timeSeriesManager.getLogger()->setLogToNull();
        timeSeriesManager.append ("mean_energy", 1.1);
        timeSeriesManager.append ("num_unhealthy_zones", 12);
        timeSeriesManager.append ("mean_pressure", 2.2);
        timeSeriesManager.append ("num_unhealthy_zones", 13);

        WHEN ("Columns are queried")
        {
            THEN ("They have the expected size")
            {
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 2);
                CHECK (timeSeriesManager.getSeriesNamesInt().size() == 1);
                CHECK (timeSeriesManager.getSeriesDouble ("mean_energy").size() == 1);
                CHECK (timeSeriesManager.getSeriesDouble ("mean_pressure").size() == 1);
                CHECK (timeSeriesManager.getSeriesInt ("num_unhealthy_zones").size() == 2);
                CHECK_THROWS_AS (timeSeriesManager.getSeriesInt ("not_there"), std::out_of_range);
            }
        }

        WHEN ("Time series data is written to HDF5")
        {
            auto hdf5File = Cow::H5::File ("time_series_test.h5", "w");
            timeSeriesManager.write (hdf5File);

            THEN ("There are the expected number of columns in the file")
            {
                int numDatasets = 0;
                hdf5File.iterate ([&] (std::string) { numDatasets++; });
                CHECK (numDatasets == 3);
            }

            THEN ("Time series data can be cleared and restored from the file")
            {
                timeSeriesManager.clear();
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 0);

                timeSeriesManager.load (hdf5File);
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 2);
                CHECK (timeSeriesManager.getSeriesInt ("num_unhealthy_zones").size() == 2);
            }
        }
    }
}




// ============================================================================
// #include "InitialDataGenerator.hpp"
// #include "ConstrainedTransport.hpp"

// SCENARIO ("Constrained transport should work as expected")
// {
//     GIVEN ("An instance of UniformCartesianCT")
//     {
//         auto id = InitialDataGenerator();
//         auto ct = UniformCartesianCT();
//         auto bc = PeriodicBoundaryCondition();

//     }
// }




// ============================================================================
#include "InitialDataGenerator.hpp"
#include "BoundaryConditions.hpp"

SCENARIO ("Initial data generator and periodic boundaries should work as expected")
{
    GIVEN ("An initial data generator and cartesian geometry")
    {
        auto fn = [] (double x, double, double)
        {
            return std::vector<double> { std::sin (2 * M_PI * x), std::cos (2 * M_PI * x) };
        };
        auto id = InitialDataGenerator();
        auto mg = CartesianMeshGeometry();
        auto bc = PeriodicBoundaryCondition();

        WHEN ("Initial data is generated")
        {
            auto P = id.generatePrimitive (fn, mg);
    
            THEN ("The returned array has the expected size")
            {
                CHECK (P.size(0) == mg.cellsShape()[0]);
                CHECK (P.size(3) == 2);
            }

            THEN ("The returned array has the expected data")
            {
                CHECK (P (0, 0, 0, 0) == fn (mg.coordinateAtIndex(0, 0, 0)[0], 0., 0.)[0]);
                CHECK (P (0, 0, 0, 1) == fn (mg.coordinateAtIndex(0, 0, 0)[0], 0., 0.)[1]);
            }
        }

        WHEN ("Two guard zones are not set by the initial data function")
        {
            int numGuard = 2;
            auto fullShape = mg.cellsShape();
            fullShape[0] += 2 * numGuard;
            fullShape[3] = 2;

            auto interiorRegion = Region().withLower (0, numGuard).withUpper (0, -numGuard);
            auto P = Array (fullShape);

            P[interiorRegion] = id.generatePrimitive (fn, mg);
            bc.apply (P, BoundaryCondition::MeshLocation::cell, BoundaryCondition::MeshBoundary::left,  0, numGuard);
            bc.apply (P, BoundaryCondition::MeshLocation::cell, BoundaryCondition::MeshBoundary::right, 0, numGuard);

            THEN ("Periodic boundary conditions set the expected values to the guard zones")
            {
                CHECK (P (0, 0, 0, 0) == P (mg.cellsShape()[0], 0, 0, 0));
                CHECK (P (0, 0, 0, 1) == P (mg.cellsShape()[0], 0, 0, 1));
                CHECK (P (mg.cellsShape()[0] + numGuard, 0, 0, 0) == P (numGuard, 0, 0, 0));
                CHECK (P (mg.cellsShape()[0] + numGuard, 0, 0, 1) == P (numGuard, 0, 0, 1));
            }
        }
    }
}




// ============================================================================
int TestSuite::runAllTests (int argc, const char* argv[])
{
    int result = Catch::Session().run (argc, argv);
    return (result < 0xff ? result : 0xff);
}
