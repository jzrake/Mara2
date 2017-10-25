#include <algorithm>
#include "RiemannSolvers.hpp"
#define MIN2(a, b) ((a) < (b) ? a : b)
#define MIN3(a, b, c) ((a) < (b) ? MIN2(a, c) : MIN2(b, c))
#define MAX2(a, b) ((a) > (b) ? a : b)
#define MAX3(a, b, c) ((a) > (b) ? MAX2(a, c) : MAX2(b, c))

using State = ConservationLaw::State;




// ============================================================================
State UpwindRiemannSolver::solve (const State& L, const State& R, AreaElement dA) const
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
        S.F[0] = 0.0;
        S.F[1] = 0.0;
    }
    return S;
}




// ============================================================================
State HlleRiemannSolver::solve (const State& L, const State& R, AreaElement dA) const
{
    const int nq = L.U.size();

    const double epl = *std::max_element (L.A.begin(), L.A.end());
    const double eml = *std::min_element (L.A.begin(), L.A.end());
    const double epr = *std::max_element (R.A.begin(), R.A.end());
    const double emr = *std::min_element (R.A.begin(), R.A.end());
    const double ap = MAX3(0.0, epl, epr);
    const double am = MIN3(0.0, eml, emr);

    auto S = ConservationLaw::State();

    for (int q = 0; q < nq; ++q)
    {
        S.U[q] = (ap * R.U[q] - am * L.U[q] + (L.F[q] - R.F[q])) / (ap - am);
        S.F[q] = (ap * L.F[q] - am * R.F[q] + ap * am * (R.U[q] - L.U[q])) / (ap - am);
    }
    return S;
}




// ============================================================================
State HllcNewtonianHydroRiemannSolver::solve (const State& L, const State& R, AreaElement dA) const
{
    enum { ddd, px, py, pz, nrg }; // Conserved
    enum { rho, vx, vy, vz, pre }; // Primitive

    const double epl = *std::max_element (L.A.begin(), L.A.end());
    const double eml = *std::min_element (L.A.begin(), L.A.end());
    const double epr = *std::max_element (R.A.begin(), R.A.end());
    const double emr = *std::min_element (R.A.begin(), R.A.end());
    const double ap = MAX3 (0.0, epl, epr);
    const double am = MIN3 (0.0, eml, emr);

    const auto& Pl = L.P;
    const auto& Pr = R.P;
    const auto& Ul = L.U;
    const auto& Ur = R.U;

    const double vnr = R.P[vx] * dA[0] + R.P[vy] * dA[1] + R.P[vz] * dA[2];
    const double vnl = L.P[vx] * dA[0] + L.P[vy] * dA[1] + L.P[vz] * dA[2];

    double Ul_[5], Ur_[5]; // The star states
    double lc = (
        + (Pr[pre] - Pr[rho] * vnr * (ap - vnr))
        - (Pl[pre] - Pl[rho] * vnl * (am - vnl))) / (Pl[rho] * (am - vnl) - Pr[rho] * (ap - vnr)); // eqn 10.58

    double fl = Pl[rho] * (am - vnl) / (am - lc); // eqn 10.33
    double fr = Pr[rho] * (ap - vnr) / (ap - lc); // eqn 10.33

    Ul_[ddd] = fl;
    Ul_[nrg] = fl * (Ul[nrg] / Pl[rho] + (lc - vnl) * (lc + Pl[pre] / (Pl[rho] * (am - vnl))));
    Ul_[px ] = fl * ((lc - vnl) * dA[0] + L.P[vx]);
    Ul_[py ] = fl * ((lc - vnl) * dA[1] + L.P[vy]);
    Ul_[pz ] = fl * ((lc - vnl) * dA[2] + L.P[vz]);

    Ur_[ddd] = fr;
    Ur_[nrg] = fr * (Ur[nrg] / Pr[rho] + (lc - vnr) * (lc + Pr[pre] / (Pr[rho] * (ap - vnr))));
    Ur_[px ] = fr * ((lc - vnr) * dA[0] + R.P[vx]);
    Ur_[py ] = fr * ((lc - vnr) * dA[1] + R.P[vy]);
    Ur_[pz ] = fr * ((lc - vnr) * dA[2] + R.P[vz]);


    const double s = 0.0; // s := x / t, set to zero for the interface
    auto S = ConservationLaw::State();

    if      (         s<=am ) for (int i=0; i<5; ++i) S.U[i] = Ul [i];
    else if ( am<s && s<=lc ) for (int i=0; i<5; ++i) S.U[i] = Ul_[i];
    else if ( lc<s && s<=ap ) for (int i=0; i<5; ++i) S.U[i] = Ur_[i];
    else if ( ap<s          ) for (int i=0; i<5; ++i) S.U[i] = Ur [i];

    if      (         s<=am ) for (int i=0; i<5; ++i) S.F[i] = L.F[i];
    else if ( am<s && s<=lc ) for (int i=0; i<5; ++i) S.F[i] = L.F[i] + am * (Ul_[i] - Ul[i]);
    else if ( lc<s && s<=ap ) for (int i=0; i<5; ++i) S.F[i] = R.F[i] + ap * (Ur_[i] - Ur[i]);
    else if ( ap<s          ) for (int i=0; i<5; ++i) S.F[i] = R.F[i];

    return S;
}
