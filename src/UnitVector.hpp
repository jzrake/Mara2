#ifndef UnitVector_hpp
#define UnitVector_hpp



/**
A class to encapsulate a 3D unit vector.
*/
class UnitVector
{
public:
    static UnitVector xhat;
    static UnitVector yhat;
    static UnitVector zhat;

    /**
    Construct a unit vector that points toward the given point in cartestian
    coordinates. If normalized in true, then this function assumes that (vx,
    vy, vz) is on the unit sphere.
    */
    static UnitVector fromCartesian (double vx, double vy, double vz, bool normalized=false);

    /**
    Construct a unit vector from its pitch angle, mu = cos (theta), and azimuth phi.
    */
    UnitVector (double pitchAngleMu, double azimuthalAnglePhi);

    /**
    Return the cartesian components of this unit vector.
    */
    void getCartesianComponents (double& nx, double& ny, double& nz) const;

    /** Return one of the cartesian components */
    double getX() const;

    /** Return one of the cartesian components */
    double getY() const;

    /** Return one of the cartesian components */
    double getZ() const;

    /**
    Return the cosine of the angle between two unit vectors.
    */
    double pitchAngleWith (const UnitVector& other) const;

    /**
    Return the projection of a cartesian vector onto this unit vector.
    */
    double project (double vx, double vy, double vz) const;

    /**
    Return the unit vector which results if this vector's pitch and azimuthal
    angles were with respect to newPolarAxis, rather than the z-axis. For
    example: zhat.withPolarAxis (xhat) = xhat.
    */
    UnitVector withPolarAxis (const UnitVector& newPolarAxis);

    double pitchAngleMu;
    double azimuthalAnglePhi;
};

#endif
