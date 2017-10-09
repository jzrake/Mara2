#include <cmath>
#include "RotationMatrix.hpp"

#define DOT(Mi, x) (Mi)[0] * (x)[0] + (Mi)[1] * (x)[1] + (Mi)[2] * (x)[2]

RotationMatrix RotationMatrix::aboutX (double angle)
{
    RotationMatrix R;
    double c = std::cos (angle);
    double s = std::sin (angle);
    R.M[0][0] =  1; R.M[0][1] =  0; R.M[0][2] =  0;
    R.M[1][0] =  0; R.M[1][1] =  c; R.M[1][2] = -s;
    R.M[2][0] =  0; R.M[2][1] =  s; R.M[2][2] =  c;
    return R;
}

RotationMatrix RotationMatrix::aboutY (double angle)
{
    RotationMatrix R;
    double c = std::cos (angle);
    double s = std::sin (angle);
    R.M[0][0] =  c; R.M[0][1] =  0; R.M[0][2] =  s;
    R.M[1][0] =  0; R.M[1][1] =  1; R.M[1][2] =  0;
    R.M[2][0] = -s; R.M[2][1] =  0; R.M[2][2] =  c;
    return R;
}

RotationMatrix RotationMatrix::aboutZ (double angle)
{
    RotationMatrix R;
    double c = std::cos (angle);
    double s = std::sin (angle);
    R.M[0][0] =  c; R.M[0][1] = -s; R.M[0][2] =  0;
    R.M[1][0] =  s; R.M[1][1] =  c; R.M[1][2] =  0;
    R.M[2][0] =  0; R.M[2][1] =  0; R.M[2][2] =  1;
    return R;
}

UnitVector RotationMatrix::operator* (const UnitVector& nhat)
{
    double n[3]; nhat.getCartesianComponents (n[0], n[1], n[2]);
    double v[3] = { DOT (M[0], n), DOT (M[1], n), DOT (M[2], n) };
    return UnitVector::normalizeFrom (v[0], v[1], v[2], true);
}
