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
        int numGuard) const override;

    bool isAxisPeriodic (int axis) override { return true; }
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
        int numGuard) const override;

    bool isAxisPeriodic (int axis) override { return false; }
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
        int numGuard) const override;

    void setConservationLaw (std::shared_ptr<ConservationLaw> law) override;

    bool isAxisPeriodic (int axis) override { return false; }

private:
    Cow::Array reflect (const Cow::Array::Reference& validData, int axis) const;
    std::shared_ptr<ConservationLaw> conservationLaw;
};




/**
Implements driving of magnetic field foot points from the domain boundary.
Intended for use with MHD in 3D.
*/
class DrivenMHDBoundary : public BoundaryCondition
{
public:
    DrivenMHDBoundary();
    void setBoundaryValueFunction (InitialDataFunction) override;
    void setConservationLaw (std::shared_ptr<ConservationLaw>) override;
    void setMeshGeometry (std::shared_ptr<MeshGeometry>) override;
    void apply (
        Cow::Array& A,
        MeshLocation location,
        MeshBoundary boundary,
        int axis,
        int numGuard) const override;
    bool isAxisPeriodic (int axis) override { return axis == 0 || axis == 1; }

private:
    InitialDataFunction velocityFunction;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<MeshGeometry> meshGeometry;
};


#endif
