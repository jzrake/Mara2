#include <cmath>
#include <sstream>
#include <cassert>
#include "ConservationLaws.hpp"

// Indexes to primitive quanitites P
#define RHO 0
#define V11 1
#define V22 2
#define V33 3
#define PRE 4
#define B11 5
#define B22 6
#define B33 7

// Indexes to conserved quanitites U
#define DDD 0
#define S11 1
#define S22 2
#define S33 3
#define NRG 4
#define H11 5
#define H22 6
#define H33 7




// ============================================================================
ConservationLaw::StateFailure::StateFailure (const State& failedState) : failedState (failedState)
{
    zoneIndex[0] = 0;
    zoneIndex[1] = 0;
    zoneIndex[2] = 0;
    zoneIndex[3] = 0;
    zoneIndex[4] = 0;
    updateWhatMessage();
}

const char* ConservationLaw::StateFailure::what() const noexcept
{
    return whatMessage.c_str();
}

void ConservationLaw::StateFailure::setZoneIndex (Cow::Index I)
{
    zoneIndex = I;
    updateWhatMessage();
}

void ConservationLaw::StateFailure::updateWhatMessage()
{
    auto stream = std::ostringstream();
    stream << "at zone index [" << zoneIndex[0] << " " << zoneIndex[1] << " " << zoneIndex[2] << "]\n";
    stream << "P = {";

    for (int q = 0; q < failedState.numFields; ++q)
    {
        stream << failedState.P[q] << " ";
    }
    stream << "}\n";
    stream << "U = {";

    for (int q = 0; q < failedState.numFields; ++q)
    {
        stream << failedState.U[q] << " ";
    }
    stream << "}";

    whatMessage = stream.str();
}




// ============================================================================
std::vector<std::string> ConservationLaw::getPrimitiveNames() const
{
    auto names = std::vector<std::string>();

    for (int q = 0; q < getNumConserved(); ++q)
    {
        names.push_back (getPrimitiveName (q));
    }
    return names;
}

ConservationLaw::State ConservationLaw::averageStates (const Request& request,
    const State& L, const State& R) const
{
    int nq = getNumConserved();
    auto Paverage = std::vector<double> (nq);

    for (int q = 0; q < nq; ++q)
    {
        Paverage[q] = 0.5 * (L.P[q] + R.P[q]);
    }
    return fromPrimitive (request, &Paverage[0]);
}

std::vector<ConservationLaw::State> ConservationLaw::fromPrimitive
(const Request& request, const Cow::Array& P) const
{
    assert (P.size (1) == getNumConserved());
    auto states = std::vector<ConservationLaw::State>(P.size (0));

    for (int n = 0; n < P.size (0); ++n)
    {
        states[n] = fromPrimitive (request, &P(n));
    }
    return states;
}

double ConservationLaw::maxEigenvalueMagnitude (const State& state) const
{
    double maxLambda = 0.0;
    int nq = getNumConserved();

    for (int n = 0; n < nq; ++n)
    {
        if (maxLambda < std::fabs (state.A[n]))
        {
            maxLambda = std::fabs (state.A[n]);
        }
    }
    return maxLambda;
}

double ConservationLaw::maxEigenvalueMagnitude (const StateVector& states) const
{
    double maxLambda = 0.0;

    for (unsigned int n = 0; n < states.size(); ++n)
    {
        double A = maxEigenvalueMagnitude (states[n]);

        if (maxLambda < A) maxLambda = A;
    }
    return maxLambda;
}




// ============================================================================
ConservationLaw::Request::Request()
{
    getPrimitive = false;
    getConserved = false;
    getFluxes = false;
    getEigenvalues = false;
    areaElement[0] = 1.0;
    areaElement[1] = 0.0;
    areaElement[2] = 0.0;
}

ConservationLaw::Request ConservationLaw::Request::oriented (const UnitVector& nhat) const
{
    auto R = *this;
    R.areaElement = nhat.cartesian();
    return R;
}




// ============================================================================
ScalarAdvection::ScalarAdvection() : waveSpeed (1.0)
{

}

ConservationLaw::State ScalarAdvection::fromConserved (const Request& request, const double* U) const
{
    double u = U[0];
    double v = U[1];
    State S;
    S.P = {{u, v}};
    S.U = {{u, v}};
    S.A = {{waveSpeed, waveSpeed}};
    S.F = {{waveSpeed * u, waveSpeed * v}};
    S.L = Cow::Matrix (2, 2); // identity
    S.R = Cow::Matrix (2, 2);
    S.numFields = 2;
    return S;
}

ConservationLaw::State ScalarAdvection::fromPrimitive (const Request& request, const double* P) const
{
    double u = P[0];
    double v = P[1];
    State S;
    S.P = {{u, v}};
    S.U = {{u, v}};
    S.A = {{waveSpeed, waveSpeed}};
    S.F = {{waveSpeed * u, waveSpeed * v}};
    S.L = Cow::Matrix (2, 2); // identity
    S.R = Cow::Matrix (2, 2);
    S.numFields = 2;
    return S;
}

int ScalarAdvection::getNumConserved() const
{
    return 2;
}

int ScalarAdvection::getIndexFor (VariableType type) const
{
    switch (type)
    {
        case VariableType::density: return 0;
        default: return -1;
    }
}

std::string ScalarAdvection::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case 0: return "u";
        case 1: return "v";
        default: return "";
    }
}




// ============================================================================
NewtonianHydro::NewtonianHydro() : gammaLawIndex (5./3)
{

}

ConservationLaw::State NewtonianHydro::fromConserved (const Request& request, const double* U) const
{
    const double gm1 = gammaLawIndex - 1.0;
    const double pp = U[S11] * U[S11] + U[S22] * U[S22] + U[S33] * U[S33];
    double P[5];

    P[RHO] =  U[DDD];
    P[PRE] = (U[NRG] - 0.5 * pp / U[DDD]) * gm1;
    P[V11] =  U[S11] / U[DDD];
    P[V22] =  U[S22] / U[DDD];
    P[V33] =  U[S33] / U[DDD];

    return fromPrimitive (request, P);
}

ConservationLaw::State NewtonianHydro::fromPrimitive (const Request& request, const double* P) const
{
    const auto dAA = request.areaElement;
    const double gm0 = gammaLawIndex;
    const double gm1 = gammaLawIndex - 1.0;
    const double cs = std::sqrt (gm0 * P[PRE] / P[RHO]);
    const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];

    auto S = State();
    S.numFields = 5;

    S.P[RHO] = P[RHO];
    S.P[V11] = P[V11];
    S.P[V22] = P[V22];
    S.P[V33] = P[V33];
    S.P[PRE] = P[PRE];

    S.U[DDD] = P[RHO];
    S.U[S11] = P[RHO] * P[V11];
    S.U[S22] = P[RHO] * P[V22];
    S.U[S33] = P[RHO] * P[V33];
    S.U[NRG] = P[RHO] * 0.5 * vv + P[PRE] / gm1;

    S.F[DDD] = vn * S.U[DDD];
    S.F[S11] = vn * S.U[S11] + P[PRE] * dAA[0];
    S.F[S22] = vn * S.U[S22] + P[PRE] * dAA[1];
    S.F[S33] = vn * S.U[S33] + P[PRE] * dAA[2];
    S.F[NRG] = vn * S.U[NRG] + P[PRE] * vn;

    S.A[0] = vn - cs;
    S.A[1] = vn;
    S.A[2] = vn;
    S.A[3] = vn;
    S.A[4] = vn + cs;

    return S;
}

int NewtonianHydro::getNumConserved() const
{
    return 5;
}

int NewtonianHydro::getIndexFor (VariableType type) const
{
    switch (type)
    {
        case VariableType::density: return 0;
        case VariableType::velocity: return 1;
        case VariableType::pressure: return 4;
        default: return -1;
    }
}

std::string NewtonianHydro::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case RHO: return "density";
        case V11: return "velocity1";
        case V22: return "velocity2";
        case V33: return "velocity3";
        case PRE: return "pressure";
        default: return "";
    }
}

std::vector<double> NewtonianHydro::makeDiagnostics (const State& state) const
{
    auto D = std::vector<double> (5);
    D[0 ] = state.U[RHO];
    D[1 ] = state.U[S11];
    D[2 ] = state.U[RHO] * std::log(state.P[PRE] / std::pow(state.U[RHO], 5.0 / 3.0)); //use this to get entropic density; originally state.U[S22]
    D[3 ] = state.U[S33];
    D[4 ] = state.U[NRG];
    return D;
}

std::vector<std::string> NewtonianHydro::getDiagnosticNames() const
{
    auto N = std::vector<std::string>(5);
    N[0 ] = "mass";
    N[1 ] = "momentum1";
    N[2 ] = "total_entropy"; //originally "momentum2"
    N[3 ] = "momentum3";
    N[4 ] = "total_energy";
    return N;
}




// ============================================================================
NewtonianMHD::NewtonianMHD()
{
    setCoolingRate (-1.0);
    setGammaLawIndex (5. / 3);
    setPressureFloor (-1.0);
}

void NewtonianMHD::setCoolingRate (double cr)
{
    coolingRate = cr;
    numAuxiliary = cr > 0.0 ? 1 : 0;
}

void NewtonianMHD::setGammaLawIndex (double gm)
{
    gammaLawIndex = gm;
}

void NewtonianMHD::setPressureFloor (double pf)
{
    pressureFloor = pf;
}

ConservationLaw::State NewtonianMHD::fromConserved (const Request& request, const double* U) const
{
    const double gm1 = gammaLawIndex - 1.0;
    const double pp = U[S11] * U[S11] + U[S22] * U[S22] + U[S33] * U[S33];
    const double BB = U[H11] * U[H11] + U[H22] * U[H22] + U[H33] * U[H33];
    double P[MARA_NUM_FIELDS];

    P[RHO] =  U[DDD];
    P[PRE] = (U[NRG] - 0.5 * pp / U[DDD] - 0.5 * BB) * gm1;
    P[V11] =  U[S11] / U[DDD];
    P[V22] =  U[S22] / U[DDD];
    P[V33] =  U[S33] / U[DDD];
    P[B11] =  U[H11];
    P[B22] =  U[H22];
    P[B33] =  U[H33];

    for (int n = 0; n < numAuxiliary; ++n)
    {
        P[8 + n] = U[8 + n];
    }

    // ------------------------------------------------------------------------
    // Bad state detection
    // ------------------------------------------------------------------------
    if (P[PRE] < 0.0 || P[RHO] < 0.0 || U[DDD] < 0.0 || U[NRG] < 0.0)
    {
        if (P[PRE] < 0.0 && pressureFloor > 0.0)
        {
            P[PRE] = pressureFloor * P[RHO];
            auto S = fromPrimitive (request, P);
            S.healthFlag = 1;
            return S;
        }
        auto S = State();
        S.numFields = getNumConserved();

        for (int q = 0; q < getNumConserved(); ++q)
        {
            S.U[q] = U[q];
            S.P[q] = P[q];
        }
        throw ConservationLaw::StateFailure (S);
    }
    return fromPrimitive (request, P);
}

ConservationLaw::State NewtonianMHD::fromPrimitive (const Request& request, const double* P) const
{
    const auto dAA = request.areaElement;
    const double gm0 = gammaLawIndex;
    const double gm1 = gammaLawIndex - 1.0;
    const double cs2 = gm0 * P[PRE] / P[RHO];
    const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
    const double BB = P[B11] * P[B11] + P[B22] * P[B22] + P[B33] * P[B33];
    const double Bv = P[B11] * P[V11] + P[B22] * P[V22] + P[B33] * P[V33];
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];
    const double Bn = P[B11] * dAA[0] + P[B22] * dAA[1] + P[B33] * dAA[2];
    const double ps = P[PRE] + 0.5 * BB; // total pressure

    auto S = State();
    S.numFields = getNumConserved();

    S.P[RHO] = P[RHO];
    S.P[V11] = P[V11];
    S.P[V22] = P[V22];
    S.P[V33] = P[V33];
    S.P[PRE] = P[PRE];
    S.P[B11] = P[B11];
    S.P[B22] = P[B22];
    S.P[B33] = P[B33];

    S.U[DDD] = P[RHO];
    S.U[S11] = P[RHO] * P[V11];
    S.U[S22] = P[RHO] * P[V22];
    S.U[S33] = P[RHO] * P[V33];
    S.U[NRG] = P[RHO] * 0.5 * vv + 0.5 * BB + P[PRE] / gm1;
    S.U[H11] = P[B11];
    S.U[H22] = P[B22];
    S.U[H33] = P[B33];

    S.F[DDD] = vn * S.U[DDD];
    S.F[S11] = vn * S.U[S11] - Bn * P[B11] + ps * dAA[0];
    S.F[S22] = vn * S.U[S22] - Bn * P[B22] + ps * dAA[1];
    S.F[S33] = vn * S.U[S33] - Bn * P[B33] + ps * dAA[2];
    S.F[NRG] = vn * S.U[NRG] - Bn * Bv     + ps * vn;
    S.F[H11] = vn * S.U[H11] - Bn * P[V11];
    S.F[H22] = vn * S.U[H22] - Bn * P[V22];
    S.F[H33] = vn * S.U[H33] - Bn * P[V33];

    // ------------------------------------------------------------------------
    // See Antony Jameson's notes at
    // http://aero-comlab.stanford.edu/Papers/jameson.mhd.pdf
    // ------------------------------------------------------------------------
    const double Bn2 = Bn * Bn;
    const double ca2 = BB / P[RHO]; /* Alfven */
    const double cn2 = Bn2 / P[RHO]; /* Alfven (in field direction) */
    const double cw4 = (cs2 + ca2) * (cs2 + ca2);
    const double cF2 = 0.5 * (cs2 + ca2 + std::sqrt (cw4 - 4 * cs2 * cn2)); /* fast */
    const double cS2 = 0.5 * (cs2 + ca2 - std::sqrt (cw4 - 4 * cs2 * cn2)); /* slow */

    S.A[0] = vn - std::sqrt (cF2);
    S.A[1] = vn - std::sqrt (cn2);
    S.A[2] = vn - std::sqrt (cS2);
    S.A[3] = vn;
    S.A[4] = vn;
    S.A[5] = vn + std::sqrt (cS2);
    S.A[6] = vn + std::sqrt (cn2);
    S.A[7] = vn + std::sqrt (cF2);

    for (int n = 0; n < numAuxiliary; ++n)
    {
        S.P[8 + n] = P[8 + n];
        S.U[8 + n] = P[8 + n];
        S.F[8 + n] = 0.0;
        S.A[8 + n] = 0.0;
    }
    return S;
}

void NewtonianMHD::addSourceTerms (const Cow::Array& P, Cow::Array& L) const
{
    if (coolingRate < 0.0) return;

    const double gm1 = gammaLawIndex - 1.0;

    P.shape3D().deploy ([&] (int i, int j, int k)
    {
        const double p = P (i, j, k, PRE);
        L (i, j, k, NRG) -= p / gm1 * coolingRate;
        L (i, j, k, 8  ) += p / gm1 * coolingRate;
    });
}

int NewtonianMHD::getNumConserved() const
{
    return 8 + numAuxiliary;
}

int NewtonianMHD::getIndexFor (VariableType type) const
{
    switch (type)
    {
        case VariableType::density: return 0;
        case VariableType::velocity: return 1;
        case VariableType::pressure: return 4;
        case VariableType::magnetic: return 5;
        default: return -1;
    }
}

std::string NewtonianMHD::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case RHO: return "density";
        case V11: return "velocity1";
        case V22: return "velocity2";
        case V33: return "velocity3";
        case PRE: return "pressure";
        case B11: return "magnetic1";
        case B22: return "magnetic2";
        case B33: return "magnetic3";
        default: return "auxiliary" + std::to_string (fieldIndex - 7);
    }
}

std::vector<double> NewtonianMHD::makeDiagnostics (const State& state) const
{
    const double gm = gammaLawIndex;
    const double* P = state.P.begin();
    const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
    const double BB = P[B11] * P[B11] + P[B22] * P[B22] + P[B33] * P[B33];
    const double Bv = P[B11] * P[V11] + P[B22] * P[V22] + P[B33] * P[V33];

    const double cs = std::sqrt (gm * state.P[PRE] / state.P[RHO]);
    const double ca = std::sqrt (BB / state.P[RHO]);
    const double s0 = state.P[PRE] / std::pow (state.P[RHO], gm);
    const double Ms = std::sqrt (vv) / cs;
    const double Ma = std::sqrt (vv) / ca;

    auto D = std::vector<double> (17 + numAuxiliary);
    D[0 ] = state.U[RHO];
    D[1 ] = state.U[S11];
    D[2 ] = state.U[S22];
    D[3 ] = state.U[S33];
    D[4 ] = state.U[NRG];
    D[5 ] = state.U[B11];
    D[6 ] = state.U[B22];
    D[7 ] = state.U[B33];
    D[8 ] = 0.5 * vv * state.P[RHO];   // kinetic energy
    D[9 ] = state.P[PRE] / (gm - 1.0); // internal energy
    D[10] = 0.5 * BB;                  // magnetic energy
    D[11] = cs;                        // sound speed
    D[12] = ca;                        // Alfven speed
    D[13] = Ms;                        // sonic Mach number
    D[14] = Ma;                        // Alfven Mach number
    D[15] = s0;                        // specific entropy
    D[16] = Bv;                        // cross helicity

    for (int n = 0; n < numAuxiliary; ++n)
    {
        D[17 + n] = state.U[8 + n];
    }
    return D;
}

std::vector<std::string> NewtonianMHD::getDiagnosticNames() const
{
    auto N = std::vector<std::string>(17 + numAuxiliary);
    N[0 ] = "mass";
    N[1 ] = "momentum1";
    N[2 ] = "momentum2";
    N[3 ] = "momentum3";
    N[4 ] = "total_energy";
    N[5 ] = "magnetic_flux1";
    N[6 ] = "magnetic_flux2";
    N[7 ] = "magnetic_flux3";
    N[8 ] = "kinetic_energy";
    N[9 ] = "internal_energy";
    N[10] = "magnetic_energy";
    N[11] = "sound_speed";
    N[12] = "alfven_speed";
    N[13] = "sonic_mach_number";
    N[14] = "alfven_mach_number";
    N[15] = "specific_entropy";
    N[16] = "cross_helicity";

    for (int n = 0; n < numAuxiliary; ++n)
    {
        N[17 + n] = "auxiliary" + std::to_string (n + 1);
    }
    return N;
}



// ============================================================================
#include "QuarticPolynomial.hpp"

RelativisticMHD::RelativisticMHD()
{
    setGammaLawIndex (5. / 3);
    setPressureFloor (-1.0);
}

void RelativisticMHD::setGammaLawIndex (double gm)
{
    gammaLawIndex = gm;
    solver.set_gamma (gm);
}

void RelativisticMHD::setPressureFloor (double pf)
{
    pressureFloor = pf;
    solver.set_pressure_floor (pf);
}

ConservationLaw::State RelativisticMHD::fromConserved (const Request& request, const double* U) const
{
    double P[MARA_NUM_FIELDS];
    solver.new_state(U);
    solver.estimate_from_cons();

    if (int error = solver.solve_anton2dzw(P))
    {
        std::cout << "U = ["
        << U[0] << " "
        << U[1] << " "
        << U[2] << " "
        << U[3] << " "
        << U[4] << " "
        << U[5] << " "
        << U[6] << " "
        << U[7] << "]\n";

        throw std::runtime_error (solver.get_error (error));
    }
    return fromPrimitive (request, P);
}

ConservationLaw::State RelativisticMHD::fromPrimitive (const Request& request, const double* P) const
{
    const auto dAA = request.areaElement;
    const double gm0 = gammaLawIndex;
    const double gm1 = gammaLawIndex - 1.0;
    const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
    const double BB = P[B11] * P[B11] + P[B22] * P[B22] + P[B33] * P[B33];
    const double Bv = P[B11] * P[V11] + P[B22] * P[V22] + P[B33] * P[V33];
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];
    const double Bn = P[B11] * dAA[0] + P[B22] * dAA[1] + P[B33] * dAA[2];
    const double W2 = 1.0 / (1.0 - vv);
    const double W0 = std::sqrt (W2);
    const double b0 = W0 * Bv;
    const double bb = (BB + b0 * b0) / W2;
    const double b1 = (P[B11] + b0 * W0 * P[V11]) / W0;
    const double b2 = (P[B22] + b0 * W0 * P[V22]) / W0;
    const double b3 = (P[B33] + b0 * W0 * P[V33]) / W0;
    const double bn = b1 * dAA[0] + b2 * dAA[1] + b3 * dAA[2];
    const double e0 = (P[PRE] / P[RHO]) / gm1;
    const double es = e0     + 0.5 * bb / P[RHO];
    const double p0 = P[PRE];
    const double ps = p0 + 0.5 * bb;
    const double h0 = 1.0 + e0 + p0 / P[RHO];
    const double hs = 1.0 + es + ps / P[RHO];

    auto S = State();
    S.numFields = getNumConserved();


    /**
    Primitive
    ---------------------------------------------------------------------------
    */
    S.P[RHO] = P[RHO];
    S.P[V11] = P[V11];
    S.P[V22] = P[V22];
    S.P[V33] = P[V33];
    S.P[PRE] = P[PRE];
    S.P[B11] = P[B11];
    S.P[B22] = P[B22];
    S.P[B33] = P[B33];


    /**
    Conserved
    ---------------------------------------------------------------------------
    */
    S.U[DDD] = P[RHO] * W0;
    S.U[NRG] = P[RHO] * hs * W2 - ps     - b0 * b0 - S.U[DDD];
    S.U[S11] = P[RHO] * hs * W2 * P[V11] - b0 * b1;
    S.U[S22] = P[RHO] * hs * W2 * P[V22] - b0 * b2;
    S.U[S33] = P[RHO] * hs * W2 * P[V33] - b0 * b3;
    S.U[H11] = P[B11];
    S.U[H22] = P[B22];
    S.U[H33] = P[B33];


    /**
    Fluxes
    ---------------------------------------------------------------------------
    */
    S.F[DDD] = vn * S.U[DDD];
    S.F[S11] = vn * S.U[S11] - Bn * b1 / W0 + ps * dAA[0];
    S.F[S22] = vn * S.U[S22] - Bn * b2 / W0 + ps * dAA[1];
    S.F[S33] = vn * S.U[S33] - Bn * b3 / W0 + ps * dAA[2];
    S.F[NRG] = vn * S.U[NRG] - Bn * b0 / W0 + ps * vn;
    S.F[H11] = vn * S.U[H11] - Bn * P[V11];
    S.F[H22] = vn * S.U[H22] - Bn * P[V22];
    S.F[H33] = vn * S.U[H33] - Bn * P[V33];


    /**
    Eigenvalues
    ---------------------------------------------------------------------------
    */
    const double cs2 = gm0 * p0 / (P[RHO] * h0);
    const double W4 = W2 * W2;
    const double v2 = vn * vn;
    const double v3 = vn * v2;
    const double v4 = vn * v3;

    const double K  =  W4 * (P[RHO] * h0 * (1. / cs2 - 1.));
    const double L  = -W2 * (P[RHO] * h0 + bb / cs2);
    const double A4 =      K      - L             -     b0 * b0;
    const double A3 = -4 * K * vn + L * vn * 2    + 2 * b0 * bn;
    const double A2 =  6 * K * v2 + L * (1. - v2) +     b0 * b0 - bn * bn;
    const double A1 = -4 * K * v3 - L * vn * 2    - 2 * b0 * bn;
    const double A0 =      K * v4 + L * v2        +     bn * bn;

    QuarticPolynomial quartic (A4, A3, A2, A1, A0);
    double roots[4];
    const int nr = quartic.solve (roots);
    const double C = std::sqrt (h0 + bb); /* constant in Alfven wave expression */

    if (nr == 4)
    {
        S.A[0] = roots[0];
        S.A[1] = (bn - C * vn * W0) / (b0 - C * W0);
        S.A[2] = roots[1];
        S.A[3] = vn;
        S.A[4] = vn;
        S.A[5] = roots[2];
        S.A[6] = (bn + C * vn * W0) / (b0 + C * W0);
        S.A[7] = roots[3];
    }
    else if (nr == 2) // This generally happens when B=0
    {
        S.A[0] = roots[0];
        S.A[1] = (bn - C * vn * W0) / (b0 - C * W0);
        S.A[2] = vn;
        S.A[3] = vn;
        S.A[4] = vn;
        S.A[5] = vn;
        S.A[6] = (bn + C * vn * W0) / (b0 + C * W0);
        S.A[7] = roots[1];
    }
    else
    {
        std::cout << "nr = " << nr << std::endl;

        std::cout << "R = ["
        << roots[0] << " "
        << roots[1] << " "
        << roots[2] << " "
        << roots[3] << "]\n";

        std::cout << "C = ["
        << A0 << " "
        << A1 << " "
        << A2 << " "
        << A3 << " "
        << A4 << "]\n";

        std::cout << "P = ["
        << S.P[0] << " "
        << S.P[1] << " "
        << S.P[2] << " "
        << S.P[3] << " "
        << S.P[4] << " "
        << S.P[5] << " "
        << S.P[6] << " "
        << S.P[7] << "]\n";

        std::cout << "A = ["
        << S.A[0] << " "
        << S.A[1] << " "
        << S.A[2] << " "
        << S.A[3] << " "
        << S.A[4] << " "
        << S.A[5] << " "
        << S.A[6] << " "
        << S.A[7] << "]\n";

        throw std::runtime_error ("got bad eigenvalues");
    }
    return S;
}

void RelativisticMHD::addSourceTerms (const Cow::Array& P,  Cow::Array& L) const
{
}

int RelativisticMHD::getNumConserved() const
{
    return 8;
}

int RelativisticMHD::getIndexFor (VariableType type) const
{
    switch (type)
    {
        case VariableType::density: return 0;
        case VariableType::velocity: return 1;
        case VariableType::pressure: return 4;
        case VariableType::magnetic: return 5;
        default: return -1;
    }
}

std::string RelativisticMHD::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case RHO: return "density";
        case V11: return "velocity1";
        case V22: return "velocity2";
        case V33: return "velocity3";
        case PRE: return "pressure";
        case B11: return "magnetic1";
        case B22: return "magnetic2";
        case B33: return "magnetic3";
        default: throw std::logic_error ("Invalid field index");
    }
}

std::vector<double> RelativisticMHD::makeDiagnostics (const State& state) const
{
    auto D = std::vector<double> ();
    return D;
}

std::vector<std::string> RelativisticMHD::getDiagnosticNames() const
{
    auto N = std::vector<std::string>();
    return N;
}
