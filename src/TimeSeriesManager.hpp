#ifndef TimeSeriesManager_hpp
#define TimeSeriesManager_hpp

#include <map>
#include <string>
#include <vector>
#include "Mara.hpp"
#include "Variant.hpp"

namespace Cow { namespace H5 { class File; } }




/**
This class maintains a running 
*/
class TimeSeriesManager
{
public:
    /**
    Load time series data from an HDF5 file.
    */
    void load (Cow::H5::File& file) const;

    /**
    Write time series data into the given HDF5 file.
    */
    void write (Cow::H5::File& file) const;

    /**
    Add a new entry to the time series data. For each key in the columns
    parameter, a column of data with that name is ensured to exist, with the
    type of the variant's value (must be either int or double). Iteration and
    simulation time (and maybe others) are recorded from the status struct.
    */
    void append (SimulationStatus status, Variant::NamedValues columns);

    /**
    Clear existing time series data.
    */
    void clear();

private:
    std::map<std::string, std::vector<double>> seriesDoubles;
    std::map<std::string, std::vector<int>> seriesInts;
};


#endif