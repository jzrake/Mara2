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
                << "' [type: double, size: " << seriesDoubles[name].size() << "]"
                << std::endl;
                return;
            }
            if (type == H5::DataType::nativeInt())
            {
                seriesInts[name] = location.readVectorInt (name);
                logger->log ("TimeSeriesManager")
                << "loading '"
                << name
                << "' [type: int, size: " << seriesInts[name].size() << "]"
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

void TimeSeriesManager::append (SimulationStatus status, Variant::NamedValues columns)
{
    for (auto entry : columns)
    {
        switch (entry.second.getType())
        {
            case 'd': seriesDoubles[entry.first].push_back (entry.second); break;
            case 'i': seriesInts[entry.first].push_back (entry.second); break;
            default: throw std::logic_error ("[TimeSeriesManager] can only accept "
                "entries of type int and double");
        }
    }
}

void TimeSeriesManager::clear()
{
    seriesDoubles.clear();
    seriesInts.clear();
}
