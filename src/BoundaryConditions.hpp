#ifndef BoudnaryConditions_hpp
#define BoudnaryConditions_hpp

#include "Mara.hpp"




/**
Simple periodic boundary condition to be used in single-core jobs. Works in
1D, 2D, and 3D.
*/
class PeriodicBoundaryCondition : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, int numGuard) const override;
};




/**
2D planar pipe flow with gas moving to the right: no-slip boundary condition
along top and bottom walls, inflow on the left and outflow on the right.
*/
class PlanarPipeFlow : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, int numGuard) const override;
};




/**
Implements driving of magnetic field foot points from the domain boundary.
Intended for use with MHD in 3D.
*/
class DrivenMHDBoundary : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, int numGuard) const override;
};


#endif
