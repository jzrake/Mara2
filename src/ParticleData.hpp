#pragma once
#include <vector>

namespace Cow { namespace H5 { class Location; }}




class ParticleData
{
public:
    struct Particle
    {
        double position = 0.;
        double velocity = 0.;
        double opticalDepth = 0.;
    };

    ParticleData();

    /**
    Load particle data from the given HDF5 location.
    */
    void load (Cow::H5::Location& location);

    /**
    Write particle data into the given HDF5 location.
    */
    void write (Cow::H5::Location& location) const;

    std::vector<Particle> particles;
};
