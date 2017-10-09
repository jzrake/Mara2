#include <cmath>
#include <iostream>
#include "FourVector.hpp"
#include "RotationMatrix.hpp"
#include "LorentzBoost.hpp"




// ============================================================================
FourVector::FourVector()
{
    for (int n = 0; n < 4; ++n)
    {
        components[n] = 0.0;
    }
}

FourVector::FourVector (const FourVector& other)
{
    for (int n = 0; n < 4; ++n)
    {
        components[n] = other.components[n];
    }
}

FourVector::FourVector (double u[4])
{
    for (int n = 0; n < 4; ++n)
    {
        components[n] = u[n];
    }
}

FourVector::FourVector (double E, double px, double py, double pz)
{
    components[0] = E;
    components[1] = px;
    components[2] = py;
    components[3] = pz;
}

FourVector FourVector::fromThreeVelocity (double vx, double vy, double vz)
{
    double gm = 1.0 / std::sqrt (1 - (vx * vx + vy * vy + vz * vz));
    return FourVector (gm, gm * vx, gm * vy, gm * vz);
}

FourVector FourVector::fromFourVelocity (double ux, double uy, double uz)
{
    double gm = std::sqrt (1 + (ux * ux + uy * uy + uz * uz));
    return FourVector (gm, ux, uy, uz);
}

FourVector FourVector::nullWithUnitVector (UnitVector nhat)
{
    double u[4] = {1, 0, 0, 0};
    nhat.getCartesianComponents (u[1], u[2], u[3]);
    return FourVector (u);
}

FourVector FourVector::fromGammaBetaAndUnitVector (double gammaBeta, UnitVector nhat)
{
    if (gammaBeta == 0)
        return FourVector (1, 0, 0, 0);

    double g = std::sqrt (1 + gammaBeta * gammaBeta);
    double nx, ny, nz;
    nhat.getCartesianComponents (nx, ny, nz);
    double u[4] = {g, gammaBeta * nx, gammaBeta * ny, gammaBeta * nz};
    return FourVector (u);
}

FourVector FourVector::fromBetaAndUnitVector (double beta, UnitVector nhat)
{
    if (beta == 0)
        return FourVector (1, 0, 0, 0);
    
    return fromGammaBetaAndUnitVector (beta / std::sqrt (1 - beta * beta), nhat);
}

double FourVector::getTimeComponent() const
{
    return components[0];
}

const double& FourVector::operator[] (int index) const
{
    return components[index];
}

double& FourVector::operator[] (int index)
{
    return components[index];
}

double FourVector::getThreeVelocityMagnitude() const
{
    const double *u = components;
    return std::sqrt (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0];
}

double FourVector::getThreeVelocityAlong (const UnitVector& nhat) const
{
    const double *u = components;
    double nx, ny, nz;
    nhat.getCartesianComponents (nx, ny, nz);
    return (u[1] * nx + u[2] * ny + u[3] * nz) / u[0];
}

double FourVector::getSpatialComponentMagnitude() const
{
    const double *u = components;
    return std::sqrt (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]);
}

UnitVector FourVector::getUnitThreeVector() const
{
    const double *u = components;
    return UnitVector::normalizeFrom (u[1], u[2], u[3]);
}

FourVector FourVector::operator+ (const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return FourVector (u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}

FourVector FourVector::operator- (const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return FourVector (u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}

FourVector FourVector::operator-() const
{
    const double *u = components;
    return FourVector (u[0], -u[1], -u[2], -u[3]);
}

FourVector& FourVector::operator+= (const FourVector& other)
{
    for (int n = 0; n < 4; ++n)
        components[n] += other.components[n];

    return *this;
}

FourVector& FourVector::operator-= (const FourVector& other)
{
    for (int n = 0; n < 4; ++n)
        components[n] -= other.components[n];

    return *this;
}

FourVector FourVector::operator* (double scalar) const
{
    const double *u = components;
    const double s = scalar;
    return FourVector (u[0] * s, u[1] * s, u[2] * s, u[3] * s);
}

FourVector FourVector::operator/ (double scalar) const
{
    const double *u = components;
    const double s = scalar;
    return FourVector (u[0] / s, u[1] / s, u[2] / s, u[3] / s);
}

double FourVector::operator* (const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return -u[0] * v[0] + u[1] * v[1] + u[2] * v[2] + u[3] * v[3];
}

FourVector FourVector::transformedBy (const LorentzBoost& L) const
{
    return L * (*this);
}

bool FourVector::isNull (double tol) const
{
    return std::fabs ((*this) * (*this)) < tol;
}

bool FourVector::isFourVelocity (double tol) const
{
    return std::fabs ((*this) * (*this) + 1) < tol;
}

bool FourVector::isSpacelike() const
{
    return (*this) * (*this) < 0;
}

bool FourVector::isTimelike() const
{
    return (*this) * (*this) > 0;
}

void FourVector::printToStream (std::ostream& stream) const
{
    stream << "<four vector> ("
    << components[0] << ", "
    << components[1] << ", "
    << components[2] << ", "
    << components[3] << ")";
}

double FourVector::betaFromGammaBeta (double gammaBeta)
{
    return gammaBeta / std::sqrt (1 + gammaBeta * gammaBeta);
}

std::ostream& operator<< (std::ostream& os, const FourVector& u)
{
    u.printToStream (os);
    return os;
}
