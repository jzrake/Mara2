#define CATCH_CONFIG_RUNNER
#include "Catch.hpp"

#include "Array.hpp"

#include "TestSuite.hpp"
#include "Stencil.hpp"
// #include "Mara.hpp"
// #include "TimeSeriesManager.hpp"
// #include "HDF5.hpp"


using namespace Cow;



// ============================================================================
SCENARIO ("Stencil operation returns array with correct sizes")
{
    GIVEN ("A stencil with shape [2, 3, 4]")
    {
        auto stencil = Stencil();
        stencil.setFootprint (2, 3, 4);
        stencil.setOperation ([] (const Array&, Array&) { });
        stencil.setCodomainRank (4, 5);

        auto source = Array (12, 12, 12);
        auto result = stencil.evaluate (source);

        REQUIRE(result.size(0) == source.size(0) - 2 * stencil.getFootprint(0));
        REQUIRE(result.size(1) == source.size(1) - 2 * stencil.getFootprint(1));
        REQUIRE(result.size(2) == source.size(2) - 2 * stencil.getFootprint(2));
    }
}


int TestSuite::runAllTests (int argc, const char* argv[])
{
    int result = Catch::Session().run (argc, argv);
    return (result < 0xff ? result : 0xff);


    // TimeSeriesManager timeSeriesManager;

    // auto measurement1 = Variant::NamedValues();
    // auto measurement2 = Variant::NamedValues();
    // auto status = SimulationStatus();
    // auto hdf5File = Cow::H5::File ("time_series.h5", "w");

    // measurement1["mean_energy"] = 1.1;
    // measurement1["num_unhealthy_zones"] = 12;
    // timeSeriesManager.append (status, measurement1);

    // measurement2["mean_pressure"] = 2.2;
    // measurement2["num_unhealthy_zones"] = 13;
    // timeSeriesManager.append (status, measurement2);
    // timeSeriesManager.write (hdf5File);

    // hdf5File.createGroup ("group1");
    // hdf5File.createGroup ("group2");

    // auto group1 = hdf5File.getGroup ("group1");
    // timeSeriesManager.write (group1);
    // timeSeriesManager.load (group1);
    // timeSeriesManager.load (hdf5File);

    return 0;
}
