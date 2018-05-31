#include <cmath>
#include "BoundaryConditions.hpp"




// ============================================================================
void BoundaryCondition::applySimple (Cow::Array& A, Cow::Shape boundaryShape) const
{
    for (int axis = 0; axis < 3; ++axis)
    {
        if (boundaryShape[axis] > 0)
        {
            apply (A, MeshLocation::cell, MeshBoundary::left , axis, boundaryShape[axis]);
            apply (A, MeshLocation::cell, MeshBoundary::right, axis, boundaryShape[axis]);
        }
    }
}




// ============================================================================
void PeriodicBoundaryCondition::apply (
    Cow::Array& A,
    MeshLocation location,
    MeshBoundary boundary,
    int axis,
    int numGuard) const
{
    if (A.size (axis) == 1)
    {
        throw std::runtime_error ("Attempt to apply boundary "
            "condition on flattened axis " + std::to_string (axis));
    }

    auto guardZone = Cow::Region();
    auto validZone = Cow::Region();

    switch (boundary)
    {
        case MeshBoundary::left:
        {
            guardZone.lower[axis] =  0;
            guardZone.upper[axis] =  numGuard;
            validZone.lower[axis] = -numGuard * 2;
            validZone.upper[axis] = -numGuard;
            break;
        }
        case MeshBoundary::right:
        {
            guardZone.lower[axis] = -numGuard;
            guardZone.upper[axis] =  0;
            validZone.lower[axis] =  numGuard;
            validZone.upper[axis] =  numGuard * 2;
            break;
        }
    }
    A[guardZone] = A[validZone];
};




// ============================================================================
void OutflowBoundaryCondition::apply (
    Cow::Array& A,
    MeshLocation location,
    MeshBoundary boundary,
    int axis,
    int numGuard) const
{
    if (A.size (axis) == 1)
    {
        throw std::runtime_error ("Attempt to apply boundary "
            "condition on flattened axis " + std::to_string (axis));
    }

    for (int n = 0; n < numGuard; ++n)
    {
        auto guardZone = Cow::Region();
        auto validZone = Cow::Region();

        switch (boundary)
        {
            case MeshBoundary::left:
            {
                guardZone.lower[axis] = n;
                guardZone.upper[axis] = n + 1;
                validZone.lower[axis] = numGuard;
                validZone.upper[axis] = numGuard + 1;
                break;
            }
            case MeshBoundary::right:
            {
                guardZone.lower[axis] = -n - 1;
                guardZone.upper[axis] = -n;
                validZone.lower[axis] = -numGuard - 1;
                validZone.upper[axis] = -numGuard;
                break;
            }
        }
        A[guardZone] = A[validZone];
    }
};




// ============================================================================
void ReflectingBoundaryCondition::apply (
    Cow::Array& A,
    MeshLocation location,
    MeshBoundary boundary,
    int axis,
    int numGuard) const
{
    if (A.size (axis) == 1)
    {
        throw std::logic_error ("Attempt to apply boundary "
            "condition on flattened axis " + std::to_string (axis));
    }

    for (int n = 0; n < numGuard; ++n)
    {
        auto guardZone = Cow::Region();
        auto validZone = Cow::Region();

        switch (boundary)
        {
            case MeshBoundary::left:
            {
                guardZone.lower[axis] = n;
                guardZone.upper[axis] = n + 1;
                validZone.lower[axis] = 2 * numGuard - n - 1;
                validZone.upper[axis] = 2 * numGuard - n;
                break;
            }
            case MeshBoundary::right:
            {
                guardZone.lower[axis] = -n - 1;
                guardZone.upper[axis] = -n;
                validZone.lower[axis] = -2 * numGuard + n;
                validZone.upper[axis] = -2 * numGuard + n + 1;
                break;
            }
        }
        A[guardZone] = reflect (A[validZone], axis);
    }
};

void ReflectingBoundaryCondition::setConservationLaw (std::shared_ptr<ConservationLaw> law)
{
    conservationLaw = law;
}

Cow::Array ReflectingBoundaryCondition::reflect (const Cow::Array::Reference& validData, int axis) const
{
    if (conservationLaw == nullptr)
        throw std::logic_error ("ReflectingBoundaryCondition needs an instance of ConservationLaw");

    const int ivel = conservationLaw->getIndexFor (ConservationLaw::VariableType::velocity);

    auto reflectedData = Cow::Array (validData);
    auto fieldAxis = Cow::Region();
    fieldAxis.stride[3] = reflectedData.size(3);

    for (auto it = reflectedData[fieldAxis].begin();
        it != reflectedData[fieldAxis].end();
        ++it)
    {
        it[ivel + axis] *= -1;
    }

    return reflectedData;
}




// ============================================================================
#define SIGN(x) ((x > 0.) - (x < 0.))

DrivenMHDBoundary::DrivenMHDBoundary()
{
    setBoundaryValueFunction (nullptr);
}

void DrivenMHDBoundary::setConservationLaw (std::shared_ptr<ConservationLaw> law)
{
    conservationLaw = law;
}

void DrivenMHDBoundary::setMeshGeometry (std::shared_ptr<MeshGeometry> geometry)
{
    meshGeometry = geometry;
}

void DrivenMHDBoundary::setBoundaryValueFunction (InitialDataFunction newVelocityFunction)
{
    if (newVelocityFunction == nullptr)
    {
        velocityFunction = [] (double x, double y, double z)
        {
            const double vx = 0.2 * std::sin (4 * M_PI * y) * SIGN(z);
            const double vy = 0.2 * std::cos (4 * M_PI * x) * SIGN(z);
            return std::vector<double> {{vx, vy}};
        };
        return;
    }
    else if (newVelocityFunction (0, 0, 0).size() != 2)
    {
        throw std::runtime_error ("DrivenMHDBoundary: "
            "velocityFunction returned a vector of length != 2 (should be [vx, vy])");
    }

    velocityFunction = newVelocityFunction;
}

void DrivenMHDBoundary::apply (
    Cow::Array& A,
    MeshLocation location,
    MeshBoundary boundary,
    int axis,
    int numGuard) const
{
    if (! meshGeometry)    throw std::logic_error ("DrivenMHDBoundary has no mesh geometry");
    if (! conservationLaw) throw std::logic_error ("DrivenMHDBoundary has no conservation law");

    if (axis == 0 || axis == 1)
    {
        auto periodic = PeriodicBoundaryCondition();
        periodic.apply (A, location, boundary, axis, numGuard);
    }

    // If A has 3 components, then it's CT calling, and the data in A is
    // either E or B field. We use outflow BC's on them.
    if (A.size(3) == 3)
    {
        auto outflow = OutflowBoundaryCondition();
        outflow.setConservationLaw (conservationLaw);
        outflow.setMeshGeometry (meshGeometry);
        outflow.apply (A, location, boundary, axis, numGuard);
        return;
    }
    Cow::Array& P = A; // It's primitive data, so rename it.

    // We set reflecting BC's on all variables first.
    auto reflect = ReflectingBoundaryCondition();
    reflect.setConservationLaw (conservationLaw);
    reflect.setMeshGeometry (meshGeometry);
    reflect.apply (P, location, boundary, axis, numGuard);

    int ivel = conservationLaw->getIndexFor (ConservationLaw::VariableType::velocity);
    int ng = numGuard;
    int nk = P.size(2) - 2 * ng;
    int k0 = 0;
    int k1 = 0;

    switch (boundary)
    {
        case MeshBoundary::left: k0 = 0; k1 = ng; break;
        case MeshBoundary::right: k0 = nk + ng; k1 = nk + 2 * ng; break;
    }

    for (int i = 0; i < P.size(0); ++i)
    {
        for (int j = 0; j < P.size(1); ++j)
        {
            for (int k = k0; k < k1; ++k)
            {
                auto X = meshGeometry->coordinateAtIndex (i - ng, j - ng, k - ng);
                auto v = velocityFunction (X[0], X[1], X[2]);
                P (i, j, k, ivel + 0) = v[0];
                P (i, j, k, ivel + 1) = v[1];
            }
        }
    }
}
