#ifndef Mara_hpp
#define Mara_hpp

#include <memory>
#include <vector>
#include <functional>
#include "Array.hpp"




class SimulationSetup;
class SimulationStatus;
class MeshGeometry;
class BoundaryCondition;
class ConservationLaw;
class IntercellFluxScheme;

using InitialDataFunction = std::function<std::vector<double> (double, double, double)>;




class SimulationSetup
{
public:
    SimulationSetup();

    double finalTime;
    double checkpointInterval;
    double cflParameter;

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
    };

    struct Request
    {
        Request();
        bool getPrimitive;
        bool getConserved;
        bool getFluxes;
        bool getEigenvalues;
        double areaElement[3];
    };

    virtual State fromConserved (const Request& request, const double* U) const = 0;
    virtual State fromPrimitive (const Request& request, const double* P) const = 0;
    virtual int getNumConserved() const = 0;
};




class IntercellFluxScheme
{
public:
    using StateVector = std::vector<ConservationLaw::State>;
    virtual ConservationLaw::State intercellFlux (StateVector stateVector) const = 0;
    virtual int getSchemeOrder() const = 0;
};




// Classes below are only here temporarily:
#include <iostream>


// ============================================================================
class ScalarAdvection : public ConservationLaw
{
public:
    ScalarAdvection()
    {
        waveSpeed = -1.0;
    }

    State fromConserved (const Request& request, const double* U) const override
    {
        double u = U[0];
        State S;
        S.P.push_back (u);
        S.U.push_back (u);
        S.A.push_back (waveSpeed);
        S.F.push_back (waveSpeed * u);
        return S;
    }

    State fromPrimitive (const Request& request, const double* P) const override
    {
        double u = P[0];
        State S;
        S.P.push_back (u);
        S.U.push_back (u);
        S.A.push_back (waveSpeed);
        S.F.push_back (waveSpeed * u);
        return S;
    }

    int getNumConserved() const override
    {
        return 1;
    }

private:
    double waveSpeed;
};




// ============================================================================
class ScalarUpwind : public IntercellFluxScheme
{
public:
    ConservationLaw::State intercellFlux (StateVector stateVector) const override
    {
        auto S = ConservationLaw::State();
        const auto& L = stateVector[0];
        const auto& R = stateVector[1];

        if (L.A[0] > 0 && R.A[0] > 0)
        {
            S.F.push_back (L.F[0]);
        }
        else if (L.A[0] < 0 && R.A[0] < 0)
        {
            S.F.push_back (R.F[0]);
        }

        return S;
    }

    int getSchemeOrder() const override
    {
        return 1;
    }
};




// ============================================================================
class PeriodicBoundaryCondition : public BoundaryCondition
{
public:

    void apply (Cow::Array& P, int numGuard) const override
    {
        int ng = numGuard;
        int n1 = P.shape()[0] - 2 * ng;
        int nq = P.shape()[3];

        for (int i = 0; i < nq; ++i)
        {
            for (int q = 0; q < nq; ++q)
            {
                P (i, 0, 0, q) = P (n1 - 2 * ng + i, 0, 0, q);
                P (n1 - ng + i, 0, 0, q) = P (ng + i, 0, 0, q);
            }
        }
    }
};

#endif
