#ifndef TimeSeriesManager_hpp
#define TimeSeriesManager_hpp

#include <map>
#include <string>
#include <vector>
#include "Mara.hpp"
#include "Variant.hpp"

namespace Cow { namespace H5 { class DataSetCreator; } }




/**
This class handles logic associated with reading and writing time series
simulation data. It may be used as a table, with named columns, each
having homogeneous data type either int or double. However, there is no
requirement that the columns have uniform length, and you can have two
columns with the same name, one of each type, if you really want to. Calling
append with set of named Variant objects creates or appends to a column for
each entry.
*/
class TimeSeriesManager
{
public:
    /**
    Load time series data from an HDF5 location. This operation interprets
    all 1D data sets at the given HDF5 location, whose type is either int
    or double, as time series data and loads them into memory. If the name
    of the data set already exists as a column in memory then that column
    is overwritten with the contents of the data set.
    */
    void load (Cow::H5::DataSetCreator& location);

    /**
    Write time series data into the given HDF5 location.
    */
    void write (Cow::H5::DataSetCreator& location) const;

    /**
    Add a new entry to the time series data. For each key in the columns
    parameter, a column of data with that name is ensured to exist, with the
    type of the variant's value (must be either int or double). Iteration and
    simulation time (and maybe others) are recorded from the status struct.
    */
    void append (SimulationStatus status, Variant::NamedValues columns);

    /**
    Clear all existing time series data.
    */
    void clear();

private:
    std::map<std::string, std::vector<double>> seriesDoubles;
    std::map<std::string, std::vector<int>> seriesInts;
};


#endif