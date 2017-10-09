
#ifndef FourVector_hpp
#define FourVector_hpp

#include <ostream>
#include <random>
#include <functional>
//#include "RandomVariable.hpp"
#include "UnitVector.hpp"

class LorentzBoost;




/**
A class to encapsulate four-vector operations. The (-,+,+,+) metric is assumed.
*/
class FourVector
{
public:
    typedef std::function<FourVector (FourVector)> Field;

    FourVector();
    FourVector (const FourVector& other);
    FourVector (double u[4]);
    FourVector (double E, double px, double py, double pz);
    static FourVector fromThreeVelocity (double vx, double vy, double vz);
    static FourVector fromFourVelocity (double ux, double uy, double uz);
    static FourVector nullWithUnitVector (UnitVector nhat);
    static FourVector fromGammaBetaAndUnitVector (double gammaBeta, UnitVector nhat);
    static FourVector fromBetaAndUnitVector (double beta, UnitVector nhat);

    /**
    Return the magnitude of the vector's associated three-velocity, sqrt (1 -
    1 / u0^2), where u0 is the Lorentz factor. This assumes isFourVelocity()
    would evaluate to true.
    */
    double getThreeVelocityMagnitude() const;

    /**
    Compute beta.nhat where beta is the three-velocity vector associated with
    this four-vector.
    */
    double getThreeVelocityAlong (const UnitVector& nhat) const;

    /**
    Return the norm of the vector's spatial components. This does not assume
    anything about the four-vector.
    */
    double getSpatialComponentMagnitude() const;

    /**
    Return a unit vector from this four vector's spatial components.
    */
    UnitVector getUnitThreeVector() const;

    /**
    Return the time component u[0] of the four vector.
    */
    double getTimeComponent() const;
    double getLorentzFactor() const { return getTimeComponent(); }

    /**
    Return one of the components. This function does not do a range check on
    the index.
    */
    const double& operator[] (int index) const;
    double& operator[] (int index);

    /**
    Return four vector addition.
    */
    FourVector operator+ (const FourVector& other) const;

    /**
    Return four vector subtraction.
    */
    FourVector operator- (const FourVector& other) const;

    /**
    Modify this four vector by adding another one.
    */
    FourVector& operator+= (const FourVector& other);

    /**
    Modify this four vector by subtracting another one.
    */
    FourVector& operator-= (const FourVector& other);

    /**
    Return four vector whose spatial components are negated.
    */
    FourVector operator-() const;

    /**
    Return four vector with all components multiplied by the given scalar.
    */
    FourVector operator* (double scalar) const;

    /**
    Return four vector with all components divided by the given scalar.
    */
    FourVector operator/ (double scalar) const;

    /**
    Return the four vector u' = L * u.
    */
    FourVector transformedBy (const LorentzBoost& L) const;

    /**
    Return the contraction of this with another four-vector u.
    */
    double operator* (const FourVector& other) const;

    /**
    Return true if this four vector has zero contraction with itself.
    */
    bool isNull (double tol=1e-10) const;

    /**
    Return true if this four vector is a four-velocity, that is u.u == -1.
    */
    bool isFourVelocity (double tol=1e-10) const;

    /**
    Return true if this four vector has positive contraction with itself.
    */
    bool isSpacelike() const;

    /**
    Return true if this four vector has negative contraction with itself.
    */
    bool isTimelike() const;

    /**
    Print the four vector in a format like <four vector>: (ut, ux, uy, uz).
    */
    void printToStream (std::ostream& stream) const;

    /**
    Return v = u / sqrt (1 + u^2) where u = gammaBeta is the four velocity and
    v is the three velocity.
    */
    static double betaFromGammaBeta (double gammaBeta);

private:
    friend class LorentzBoost;
    double components[4];
};

std::ostream& operator<< (std::ostream& os, const FourVector& u);


#endif
