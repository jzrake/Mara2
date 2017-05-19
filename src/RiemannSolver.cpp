#include "RiemannSolver.hpp"
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
        //S.F = std::vector<double> (L.F.size(), 0.0);
        S.F[0] = 0.0;
        S.F[1] = 0.0;
    }
    return S;
}




// ============================================================================
State HlleRiemannSolver::solve (const State& L, const State& R, AreaElement dA) const
{
    int nq = L.U.size();

    const double epl = *std::max_element (L.A.begin(), L.A.end());
    const double eml = *std::min_element (L.A.begin(), L.A.end());
    const double epr = *std::max_element (R.A.begin(), R.A.end());
    const double emr = *std::min_element (R.A.begin(), R.A.end());
    const double ap = MAX3(0.0, epl, epr);
    const double am = MIN3(0.0, eml, emr);

    auto S = ConservationLaw::State();
    // S.U.resize (nq);
    // S.F.resize (nq);

    for (int q = 0; q < nq; ++q)
    {
        S.U[q] = (ap * R.U[q] - am * L.U[q] + (L.F[q] - R.F[q])) / (ap - am);
        S.F[q] = (ap * L.F[q] - am * R.F[q] + ap * am * (R.U[q] - L.U[q])) / (ap - am);
    }
    return S;
}
