#include <cmath>
#include "BoundaryConditions.hpp"





// ============================================================================
void PeriodicBoundaryCondition::apply (
    Cow::Array& A,
    MeshLocation location,
    MeshBoundary boundary,
    int axis,
    int numGuard,
    const MeshGeometry& geometry,
    const ConservationLaw& law) const
{

};




// ============================================================================
void OutflowBoundaryCondition::apply (
    Cow::Array& A,
    MeshLocation location,
    MeshBoundary boundary,
    int axis,
    int numGuard,
    const MeshGeometry& geometry,
    const ConservationLaw& law) const
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
    int numGuard,
    const MeshGeometry& geometry,
    const ConservationLaw& law) const
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

        A[guardZone] = reflect (A[validZone], law, axis);
    }
};

Cow::Array ReflectingBoundaryCondition::reflect (
    const Cow::Array::Reference& validData,
    const ConservationLaw& law,
    int axis) const
{
    const int ivel = law.getIndexFor (ConservationLaw::VariableType::velocity);

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




// // ============================================================================
// #define SIGN(x) ((x > 0.) - (x < 0.))

// DrivenMHDBoundary::DrivenMHDBoundary()
// {
//     setVelocityFunction (nullptr);
// }

// void DrivenMHDBoundary::setVelocityFunction (InitialDataFunction newVelocityFunction)
// {
//     if (newVelocityFunction == nullptr)
//     {
//         velocityFunction = [] (double x, double y, double z)
//         {
//             const double vx = 0.2 * std::sin (4 * M_PI * y) * SIGN(z);
//             const double vy = 0.2 * std::cos (4 * M_PI * x) * SIGN(z);
//             return std::vector<double> {{vx, vy}};
//         };
//         return;
//     }
//     else if (newVelocityFunction (0, 0, 0).size() != 2)
//     {
//         throw std::runtime_error ("DrivenMHDBoundary: "
//             "velocityFunction returned a vector of length != 2 (should be [vx, vy])");
//     }

//     velocityFunction = newVelocityFunction;
// }

// void DrivenMHDBoundary::apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const
// {
//     const auto periodic = PeriodicBoundaryCondition();
//     const auto vertical = ReflectingBoundaryCondition();

//     if (P.size(0) > 1) periodic.applyToAxis (P, numGuard, 0);
//     if (P.size(1) > 1) periodic.applyToAxis (P, numGuard, 1);
//     if (P.size(2) > 1) vertical.applyToAxis (P, law, numGuard, 2);

//     const int ivel = law.getIndexFor (ConservationLaw::VariableType::velocity);
//     const int ng = numGuard;
//     const int ni = P.size(0) - 2 * ng;
//     const int nj = P.size(1) - 2 * ng;
//     const int nk = P.size(2) - 2 * ng;

//     for (int i = 0; i < P.size(0); ++i)
//     {
//         for (int j = 0; j < P.size(1); ++j)
//         {
//             const double x = (i - ng + 0.5) / ni - 0.5;
//             const double y = (j - ng + 0.5) / nj - 0.5;
//             const auto vtop = velocityFunction (x, y, +1);
//             const auto vbot = velocityFunction (x, y, -1);

//             for (int k = 0; k < ng; ++k)
//             {
//                 P (i, j, k,           ivel + 0) = vbot[0];
//                 P (i, j, k + nk + ng, ivel + 0) = vtop[0];
//                 P (i, j, k,           ivel + 1) = vbot[1];
//                 P (i, j, k + nk + ng, ivel + 1) = vtop[1];
//             }
//         }
//     }
// }

// void DrivenMHDBoundary::applyToCellCenteredB (Cow::Array& B, int numGuard) const
// {
//     const auto periodic = PeriodicBoundaryCondition();
//     const auto outflow = OutflowBoundaryCondition();

//     periodic.applyToAxis (B, numGuard, 0);
//     periodic.applyToAxis (B, numGuard, 1);
//     outflow.applyToAxis (B, numGuard, 2);
// }

// void DrivenMHDBoundary::applyToGodunovFluxes (Cow::Array& F, int numGuard, int axis) const
// {
//     const auto periodic = PeriodicBoundaryCondition();
//     const auto outflow = OutflowBoundaryCondition();

//     switch (axis)
//     {
//         case 0: periodic.applyToAxis (F, numGuard, 0); break;
//         case 1: periodic.applyToAxis (F, numGuard, 1); break;
//         case 2: outflow.applyToAxis (F, numGuard, 2); break;
//     }
// }
