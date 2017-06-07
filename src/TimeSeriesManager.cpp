#include "TimeSeriesManager.hpp"
#include "HDF5.hpp"

using namespace Cow;




void TimeSeriesManager::load (H5::DataSetCreator& location)
{
    location.iterate ([&] (std::string name)
    {
        if (location.hasDataSet (name))
        {
            auto dset = location.getDataSet (name);
            auto type = dset.getType();

            if (type == H5::DataType::nativeDouble())
            {
                seriesDoubles[name] = location.readVectorDouble (name);
                logger->log ("TimeSeriesManager")
                << "loading '"
                << name
                << "' [double, length " << seriesDoubles[name].size() << "]"
                << std::endl;
                return;
            }
            if (type == H5::DataType::nativeInt())
            {
                seriesInts[name] = location.readVectorInt (name);
                logger->log ("TimeSeriesManager")
                << "loading '"
                << name
                << "' [int, length " << seriesInts[name].size() << "]"
                << std::endl;
                return;
            }
        }
        logger->log ("TimeSeriesManager") << "skipping " << name << std::endl;
    });
}

void TimeSeriesManager::write (H5::DataSetCreator& location) const
{
    for (auto column : seriesDoubles)
    {
        location.writeVectorDouble (column.first, column.second);
    }
    for (auto column : seriesInts)
    {
        location.writeVectorInt (column.first, column.second);
    }
}
void TimeSeriesManager::append (std::string name, Variant value)
{
    switch (value.getType())
    {
        case 'd': seriesDoubles[name].push_back (value); break;
        case 'i': seriesInts[name].push_back (value); break;
        default: throw std::logic_error ("[TimeSeriesManager] only accepts "
            "entries of type int and double");
    }
}

void TimeSeriesManager::append (SimulationStatus status, Variant::NamedValues columns)
{
    append ("iteration", status.simulationIter);
    append ("time", status.simulationTime);

    for (auto entry : columns)
    {
        append (entry.first, entry.second);
    }
}

void TimeSeriesManager::clear()
{
    seriesDoubles.clear();
    seriesInts.clear();
}