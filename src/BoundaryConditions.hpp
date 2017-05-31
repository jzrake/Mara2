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
    void apply (
        Cow::Array& A,
        MeshLocation location,
        MeshBoundary boundary,
        int axis,
        int numGuard,
        const MeshGeometry& geometry,
        const ConservationLaw& law) const override;
};




/**
A simple outflow boundary condition, zero-gradient is imposed on all variables.
*/
class OutflowBoundaryCondition : public BoundaryCondition
{
public:
    void apply (
        Cow::Array& A,
        MeshLocation location,
        MeshBoundary boundary,
        int axis,
        int numGuard,
        const MeshGeometry& geometry,
        const ConservationLaw& law) const override;
};




/**
A simple reflecting boundary condition, solution is mirrored around the
boundary surface.
*/
class ReflectingBoundaryCondition : public BoundaryCondition
{
public:
    void apply (
        Cow::Array& A,
        MeshLocation location,
        MeshBoundary boundary,
        int axis,
        int numGuard,
        const MeshGeometry& geometry,
        const ConservationLaw& law) const override;
private:
    Cow::Array reflect (
        const Cow::Array::Reference& validData,
        const ConservationLaw& law,
        int axis) const;
};




// /**
// Implements driving of magnetic field foot points from the domain boundary.
// Intended for use with MHD in 3D.
// */
// class DrivenMHDBoundary : public BoundaryCondition
// {
// public:
//     DrivenMHDBoundary();
//     void setVelocityFunction (InitialDataFunction newVelocityFunction);
//     void apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const override;
//     void applyToCellCenteredB (Cow::Array& B, int numGuard) const override;
//     void applyToGodunovFluxes (Cow::Array& F, int numGuard, int axis) const override;
// private:
//     InitialDataFunction velocityFunction;
// };


#endif
