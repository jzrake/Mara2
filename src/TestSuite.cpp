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

    measurement1["mean_energy"] = 1.0;
    measurement1["num_unhealthy_zones"] = 12;
    timeSeriesManager.append (status, measurement1);

    measurement2["mean_pressure"] = 2.0;
    measurement2["num_unhealthy_zones"] = 13;
    timeSeriesManager.append (status, measurement2);

    timeSeriesManager.write (hdf5File);

    return 0;
}
