#pragma once

#include <list>
#include "FourVector.hpp"




class GridLocator
{

};



class ParticleData
{
public:
    struct Particle
    {
        double position;
        double momentum;
        double velocity;
    };

    ParticleData();
    std::vector<Particle> particles;
};





// class ParticleData
// {
// public:
//     class Particle
//     {
//     public:
//         FourVector position;
//         FourVector momentum;
//         double charge, mass;
//     };

//     /**
//     Load particle data from the given HDF5 location.
//     */
//     void load (Cow::H5::Location& location);

//     /**
//     Write particle data into the given HDF5 location.
//     */
//     void write (Cow::H5::Location& location) const;

//     /**
//     Exchange particle data with other grids.
//     */
//     void exchangeParticles (GridLocator& locator);
// private:
//     std::list<Particle> data;
// };
