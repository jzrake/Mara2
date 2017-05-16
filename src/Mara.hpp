#ifndef Mara_hpp
#define Mara_hpp

#include <memory>
#include <vector>
#include <string>
#include <functional>
#include "Array.hpp"
#include "Matrix.hpp"
#include "Reconstruction.hpp"




class SimulationSetup;
class SimulationStatus;
class MeshGeometry;
class BoundaryCondition;
class ConservationLaw;
class IntercellFluxScheme;

using InitialDataFunction = std::function<std::vector<double> (double, double, double)>;
using AreaElement = std::array<double, 3>;



class SimulationSetup
{
public:
    SimulationSetup();

    // Run description
    double finalTime;
    double checkpointInterval;
    double cflParameter;
    int rungeKuttaOrder;
    std::string outputDirectory;
    std::string runName;

    // Algorithms
    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<IntercellFluxScheme> intercellFluxScheme;
    std::shared_ptr<BoundaryCondition> boundaryCondition;
    InitialDataFunction initialDataFunction;
};




class SimulationStatus
{
public:
    SimulationStatus();

    double simulationTime;
    int simulationIter;
    int outputsWrittenSoFar;
};




class MeshGeometry
{
public:
    using Coordinate = std::array<double, 3>;
    virtual Cow::Shape domainShape() const = 0;
    virtual Coordinate coordinateAtIndex (double i, double j, double k) const = 0;
    virtual double cellLength (int i, int j, int k, int axis) const = 0;
    virtual double faceArea (int i, int j, int k, int axis) const = 0;
    virtual double cellVolume (int i, int j, int k) const = 0;
};




class BoundaryCondition
{
public:
    virtual void apply (Cow::Array& P, int numGuard) const = 0;
};




class ConservationLaw
{
public:
    struct State
    {
        std::vector<double> U; // Conserved densities
        std::vector<double> P; // Primitive quantities
        std::vector<double> F; // Fluxes in given direction
        std::vector<double> A; // Eigenvalues
        Cow::Matrix L; // Left eigenvector matrix
        Cow::Matrix R; // Right eigenvector matrix
    };

    struct Request
    {
        Request();
        bool getPrimitive;
        bool getConserved;
        bool getFluxes;
        bool getEigenvalues;
        AreaElement areaElement;
    };

    using StateVector = std::vector<State>;

    /**
    Generate a state from the given information request, and a double pointer
    of conserved quantities. U must point to valid conserved quantity data
    with numConserved consecutive doubles.
    */
    virtual State fromConserved (const Request& request, const double* U) const = 0;

    /**
    Generate a state from the given information request, and a double pointer
    of primitive quantities. P must point to valid primitive quantity data
    with numConserved consecutive doubles.
    */
    virtual State fromPrimitive (const Request& request, const double* P) const = 0;

    /**
    Get the number of primitive and conserved quantities for this conservation
    law. For example, the Euler equation has 5 (density, total energy, and
    three components of momentum), and MHD has 8 (three magnetic field fluxes
    in addition). If there are passive scalars then those will increase that
    number.
    */
    virtual int getNumConserved() const = 0;

    /**
    Return the numeric average the primitive quantities from the two given
    states, and return a state computed with from that average and the given
    request.
    */
    State averageStates (const Request& request, const State& L, const State& R) const;

    /**
    A helper function that returns a vector of states by calling fromPrimitive
    numStates times, where P.shape() = [numStates, numConserved].
    */
    StateVector fromPrimitive (const Request& request, const Cow::Array& P) const;

    /**
    A helper function that returns the maximum absolute value of all
    eigenvalues for the given state.
    */
    double maxEigenvalueMagnitude (const State& state) const;

    /**
    A helper function that the maximum eigenvalue magnitude, over a vector of
    states. This is useful for Lax-Friedrichs flux splitting.
    */
    double maxEigenvalueMagnitude (const StateVector& states) const;
};




class RiemannSolver
{
public:
    using State = ConservationLaw::State;
    virtual State solve (const State& L, const State& R, AreaElement dA) = 0;
};




class IntercellFluxScheme
{
public:
    struct FaceData
    {
    public:
        Cow::Array stencilData;
        AreaElement areaElement;
        std::shared_ptr<ConservationLaw> conservationLaw;
    };

    virtual ConservationLaw::State intercellFlux (const FaceData&) const = 0;
    virtual int getStencilSize() const = 0;
};




// Classes below will be moved to implementation files soon
// ============================================================================
#include <iostream>



class UpwindRiemannSolver : public RiemannSolver
{
public:
    State solve (const State& L, const State& R, AreaElement dA) override
    {
        auto S = ConservationLaw::State();

        if (L.A[0] > 0 && R.A[0] > 0)
        {
            S.F = L.F;
        }
        else if (L.A[0] < 0 && R.A[0] < 0)
        {
            S.F = R.F;
        }
        else
        {
            S.F = std::vector<double> (L.F.size(), 0.0);
        }
        return S;
    }
};




class ScalarAdvection : public ConservationLaw
{
public:
    ScalarAdvection (double waveSpeed);
    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    int getNumConserved() const override;
private:
    double waveSpeed;
};




// ============================================================================
class ScalarUpwind : public IntercellFluxScheme
{
public:
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override;
    int getStencilSize() const override;
};




// ============================================================================
class MethodOfLines : public IntercellFluxScheme
{
public:
    MethodOfLines (double plmTheta);
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override;
    int getStencilSize() const override;
private:
    Reconstruction plm;
};




// ============================================================================
class MethodOfLinesWeno : public IntercellFluxScheme
{
public:
    MethodOfLinesWeno (double shenZhaA=50);
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override;
    int getStencilSize() const override;
private:
    Reconstruction weno;
};




// ============================================================================
class PeriodicBoundaryCondition : public BoundaryCondition
{
public:
    void apply (Cow::Array& P, int numGuard) const override
    {
        int ng = numGuard;
        int n1 = P.size(0) - 2 * ng;
        int nq = P.size(3);

        for (int i = 0; i < ng; ++i)
        {
            for (int q = 0; q < nq; ++q)
            {
                P (i, 0, 0, q) = P (i + n1, 0, 0, q);
                P (i + ng + n1, 0, 0, q) = P (i + ng, 0, 0, q);
            }
        }
    }
};


#endif
