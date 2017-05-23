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
    void apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const override;
    void applyToCellCenteredB (Cow::Array& B, int numGuard) const override;
    void applyToGodunovFluxes (Cow::Array& F, int numGuard, int axis) const override;
    void applyToAxis (Cow::Array& P, int numGuard, int axis) const;
};




/**
A simple outflow boundary condition, zero-gradient is imposed on all variables.
*/
class OutflowBoundaryCondition : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const override;
    void applyToAxis (Cow::Array& P, int numGuard, int axis) const;
};




/**
A simple reflecting boundary condition, solution is mirrored around the
boundary surface.
*/
class ReflectingBoundaryCondition : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const override;
    void applyToAxis (Cow::Array& P, const ConservationLaw& law, int numGuard, int axis) const;
};




/**
2D planar pipe flow with gas moving to the right: no-slip boundary condition
along top and bottom walls, inflow on the left and outflow on the right.
*/
class PlanarPipeFlow : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const override;
};




/**
Implements driving of magnetic field foot points from the domain boundary.
Intended for use with MHD in 3D.
*/
class DrivenMHDBoundary : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const override;
};


#endif
