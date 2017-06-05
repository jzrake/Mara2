#include "TimeSeriesManager.hpp"
#include "HDF5.hpp"




void TimeSeriesManager::load (Cow::H5::File& file) const
{

}

void TimeSeriesManager::write (Cow::H5::File& file) const
{
    for (auto column : seriesDoubles)
    {
        std::cout << column.first << " (double) " << column.second.size() << std::endl;
        file.writeVectorDouble (column.first, column.second);
    }
    for (auto column : seriesInts)
    {
        std::cout << column.first << " (int) " << column.second.size() << std::endl;
        file.writeVectorInt (column.first, column.second);
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
            default: throw std::logic_error ("[TimeSeriesManager] can only accept int and double");
        }
    }
}

void TimeSeriesManager::clear()
{
    seriesDoubles.clear();
    seriesInts.clear();
}
