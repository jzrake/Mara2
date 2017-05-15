#ifndef Mara_hpp
#define Mara_hpp

#include <memory>
#include <vector>
#include <string>
#include <functional>
#include "Array.hpp"




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

    virtual State fromConserved (const Request& request, const double* U) const = 0;
    virtual State fromPrimitive (const Request& request, const double* P) const = 0;
    virtual int getNumConserved() const = 0;
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
    using StateVector = std::vector<ConservationLaw::State>; // DEP

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
    ScalarAdvection(double waveSpeed) : waveSpeed (waveSpeed)
    {

    }

    State fromConserved (const Request& request, const double* U) const override
    {
        double u = U[0];
        State S;
        S.P = {u, 0};
        S.U = {u, 0};
        S.A = {waveSpeed, 0};
        S.F = {waveSpeed * u, 0};
        return S;
    }

    State fromPrimitive (const Request& request, const double* P) const override
    {
        double u = P[0];
        State S;
        S.P = {u, 0};
        S.U = {u, 0};
        S.A = {waveSpeed, 0};
        S.F = {waveSpeed * u, 0};
        return S;
    }

    int getNumConserved() const override
    {
        return 2;
    }

private:
    double waveSpeed;
};




// ============================================================================
class ScalarUpwind : public IntercellFluxScheme
{
public:
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override
    {
        ConservationLaw::Request request;
        request.areaElement = faceData.areaElement;

        const double* PL = &faceData.stencilData (0);
        const double* PR = &faceData.stencilData (1);

        auto L = faceData.conservationLaw->fromPrimitive (request, PL);
        auto R = faceData.conservationLaw->fromPrimitive (request, PR);

        UpwindRiemannSolver riemannSolver;
        return riemannSolver.solve (L, R, faceData.areaElement);
    }

    int getStencilSize() const override
    {
        return 1;
    }
};




// ============================================================================
#include "Reconstruction.hpp"

class MethodOfLines : public IntercellFluxScheme
{
public:
    MethodOfLines (double plmTheta)
    {
        plm.setPlmTheta (plmTheta);
    }

    ConservationLaw::State intercellFlux (const FaceData& faceData) const override
    {
        ConservationLaw::Request request;
        request.areaElement = faceData.areaElement;

        const double* P0 = &faceData.stencilData (0);
        const double* P1 = &faceData.stencilData (1);
        const double* P2 = &faceData.stencilData (2);
        const double* P3 = &faceData.stencilData (3);
        const int nq = faceData.conservationLaw->getNumConserved();

        std::vector<double> PL (nq);
        std::vector<double> PR (nq);

        for (int q = 0; q < nq; ++q)
        {
            const double ps[4] = {P0[q], P1[q], P2[q], P3[q]};
            PL[q] = plm.reconstruct (&ps[1], Reconstruction::PLM_C2R);
            PR[q] = plm.reconstruct (&ps[2], Reconstruction::PLM_C2L);
        }
        auto L = faceData.conservationLaw->fromPrimitive (request, &PL[0]);
        auto R = faceData.conservationLaw->fromPrimitive (request, &PR[0]);

        UpwindRiemannSolver riemannSolver;
        return riemannSolver.solve (L, R, faceData.areaElement);
    }

    int getStencilSize() const override
    {
        return 2;
    }

private:
    Reconstruction plm;
};




// ============================================================================
class MethodOfLinesWeno : public IntercellFluxScheme
{
public:
    MethodOfLinesWeno (double shenZhaA=50)
    {
        weno.setSmoothnessIndicator (Reconstruction::ImprovedShenZha10);
    }

    ConservationLaw::State intercellFlux (const FaceData& faceData) const override
    {
        ConservationLaw::Request request;
        request.areaElement = faceData.areaElement;

        const double* P0 = &faceData.stencilData (0);
        const double* P1 = &faceData.stencilData (1);
        const double* P2 = &faceData.stencilData (2);
        const double* P3 = &faceData.stencilData (3);
        const double* P4 = &faceData.stencilData (4);
        const double* P5 = &faceData.stencilData (5);

        const auto S0 = faceData.conservationLaw->fromPrimitive (request, P0);
        const auto S1 = faceData.conservationLaw->fromPrimitive (request, P1);
        const auto S2 = faceData.conservationLaw->fromPrimitive (request, P2);
        const auto S3 = faceData.conservationLaw->fromPrimitive (request, P3);
        const auto S4 = faceData.conservationLaw->fromPrimitive (request, P4);
        const auto S5 = faceData.conservationLaw->fromPrimitive (request, P5);

        const auto A = 1; // !!! hard-coded to 1 for now

        const double fp[6] = {
            S0.F[0] + A * S0.U[0],
            S1.F[0] + A * S1.U[0],
            S2.F[0] + A * S2.U[0],
            S3.F[0] + A * S3.U[0],
            S4.F[0] + A * S4.U[0],
            S5.F[0] + A * S5.U[0]};

        const double fm[6] = {
            S0.F[0] - A * S0.U[0],
            S1.F[0] - A * S1.U[0],
            S2.F[0] - A * S2.U[0],
            S3.F[0] - A * S3.U[0],
            S4.F[0] - A * S4.U[0],
            S5.F[0] - A * S5.U[0]};

        double fhatp = weno.reconstruct (fp + 2, Reconstruction::WENO5_FD_C2R);
        double fhatm = weno.reconstruct (fm + 3, Reconstruction::WENO5_FD_C2L);
        auto S = ConservationLaw::State();
        S.F = {0.5 * (fhatp + fhatm)};
        return S;
    }

    int getStencilSize() const override
    {
        return 3;
    }

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
        int n1 = P.shape()[0] - 2 * ng;
        int nq = P.shape()[3];

        for (int i = 0; i < ng; ++i)
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
