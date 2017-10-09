#include <cmath>
#include "LorentzBoost.hpp"

#define DOT(Ai, x) (Ai)[0] * (x)[0] + (Ai)[1] * (x)[1] + (Ai)[2] * (x)[2] + (Ai)[3] * (x)[3]


LorentzBoost::LorentzBoost() : LorentzBoost (FourVector (1, 0, 0, 0))
{

}

LorentzBoost::LorentzBoost (const FourVector& boostVector) : boostVector (boostVector)
{
    double gm = boostVector.components[0];
    double nx = boostVector.components[1];
    double ny = boostVector.components[2];
    double nz = boostVector.components[3];
    double nt = std::sqrt (nx * nx + ny * ny + nz * nz);

    if (std::fabs (nt) > 1e-12)
    {
        nx /= nt;
        ny /= nt;
        nz /= nt;
    }

    elements[0][0] = gm;
    elements[1][1] = 1 + (gm - 1) * nx * nx;
    elements[2][2] = 1 + (gm - 1) * ny * ny;
    elements[3][3] = 1 + (gm - 1) * nz * nz;

    elements[0][1] = elements[1][0] = -boostVector.components[1];
    elements[0][2] = elements[2][0] = -boostVector.components[2];
    elements[0][3] = elements[3][0] = -boostVector.components[3];

    elements[1][2] = elements[2][1] = (gm - 1) * nx * ny;
    elements[2][3] = elements[3][2] = (gm - 1) * ny * nz;
    elements[3][1] = elements[1][3] = (gm - 1) * nz * nx;
}

LorentzBoost LorentzBoost::inverted()
{
    return LorentzBoost (-boostVector);
}

FourVector LorentzBoost::operator* (const FourVector& other) const
{
    double v[4] = {
        DOT (elements[0], other.components),
        DOT (elements[1], other.components),
        DOT (elements[2], other.components),
        DOT (elements[3], other.components) };

    return FourVector (v);
}
