#include "TestSuite.hpp"
#include "Mara.hpp"
#include "TimeSeriesManager.hpp"
#include "HDF5.hpp"




int TestSuite::runAllTests()
{
    TimeSeriesManager timeSeriesManager;

    auto measurement1 = Variant::NamedValues();
    auto measurement2 = Variant::NamedValues();
    auto status = SimulationStatus();
    auto hdf5File = Cow::H5::File ("time_series.h5", "w");

    measurement1["mean_energy"] = 1.1;
    measurement1["num_unhealthy_zones"] = 12;
    timeSeriesManager.append (status, measurement1);

    measurement2["mean_pressure"] = 2.2;
    measurement2["num_unhealthy_zones"] = 13;
    timeSeriesManager.append (status, measurement2);
    timeSeriesManager.write (hdf5File);

    hdf5File.createGroup ("group1");
    hdf5File.createGroup ("group2");

    auto group1 = hdf5File.getGroup ("group1");
    timeSeriesManager.write (group1);
    timeSeriesManager.load (group1);
    timeSeriesManager.load (hdf5File);

    return 0;
}
