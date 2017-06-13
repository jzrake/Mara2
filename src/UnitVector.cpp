#include <cmath>
#include "UnitVector.hpp"
#include "RotationMatrix.hpp"




// ============================================================================
UnitVector UnitVector::xhat (0, 0);
UnitVector UnitVector::yhat (0, M_PI / 2);
UnitVector UnitVector::zhat (1, 0);

UnitVector UnitVector::fromCartesian (double vx, double vy, double vz, bool normalized)
{
    double cosTheta = vz / (normalized ? 1 : std::sqrt (vx * vx + vy * vy + vz * vz));
    double phi = std::atan2 (vy, vx);
    return UnitVector (cosTheta, phi);
}

UnitVector::UnitVector (double pitchAngleMu, double azimuthalAnglePhi) :
pitchAngleMu (pitchAngleMu),
azimuthalAnglePhi (azimuthalAnglePhi)
{

}

void UnitVector::getCartesianComponents (double& nx, double& ny, double& nz) const
{
    const double cosTheta = pitchAngleMu;
    const double sinTheta = std::sqrt (1 - cosTheta * cosTheta);
    nx = sinTheta * std::cos (azimuthalAnglePhi);
    ny = sinTheta * std::sin (azimuthalAnglePhi);
    nz = cosTheta;
}

double UnitVector::getX() const
{
    double nx, ny, nz;
    getCartesianComponents (nx, ny, nz);
    return nx;
}

double UnitVector::getY() const
{
    double nx, ny, nz;
    getCartesianComponents (nx, ny, nz);
    return ny;
}

double UnitVector::getZ() const
{
    double nx, ny, nz;
    getCartesianComponents (nx, ny, nz);
    return nz;
}

double UnitVector::pitchAngleWith (const UnitVector& other) const
{
    double u[3];
    double v[3];
    this->getCartesianComponents (u[0], u[1], u[2]);
    other.getCartesianComponents (v[0], v[1], v[2]);
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

double UnitVector::project (double vx, double vy, double vz) const
{
    double nx, ny, nz;
    getCartesianComponents (nx, ny, nz);
    return vx * nx + vy * ny + vz * nz;
}

UnitVector UnitVector::withPolarAxis (const UnitVector& newPolarAxis)
{
    double theta = std::acos (newPolarAxis.pitchAngleMu);
    double phi = newPolarAxis.azimuthalAnglePhi;
    return RotationMatrix::aboutZ (phi) * (RotationMatrix::aboutY (theta) * (*this));
}
