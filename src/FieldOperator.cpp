#include <limits>
#include <algorithm>
#include "FieldOperator.hpp"
#define MIN2(a, b) ((a) < (b) ? a : b)

using namespace Cow;




// ============================================================================
const char* FieldOperator::InversionFailure::what() const noexcept
{
    return whatMessage.c_str();
}

void FieldOperator::InversionFailure::updateWhatMessage()
{
    std::stringstream stream;

    for (auto& e : failedStates)
    {
        stream << e.what() << std::endl;
    }
    whatMessage = stream.str();
}




// ============================================================================
FieldOperator::FieldOperator (std::shared_ptr<ConservationLaw> law) : law (law)
{
}

void FieldOperator::setConservationLaw (std::shared_ptr<ConservationLaw> cl)
{
    law = cl;
}

Array FieldOperator::recoverPrimitive (Array::Reference U, Array::Reference P) const
{
    int numConserved = law->getNumConserved();
    auto request = ConservationLaw::Request();
    // auto shape = P.shape();
    auto zoneHealth = Array (U.getArray().shape3D());
    auto inversionFailure = InversionFailure();

    auto Treg = P.getRegion().withStride (3, law->getNumConserved());
    auto Preg = P.getArray()[Treg];
    auto Ureg = U.getArray()[Treg];
    auto pit = Preg.begin();
    auto uit = Ureg.begin();

    for ( ; uit != U.end(); ++pit, ++uit)
    {
        try
        {
            auto S = law->fromConserved (request, uit);

            for (int q = 0; q < numConserved; ++q)
            {
                pit[q] = S.P[q];
            }
            const auto index = uit.index();
            zoneHealth (index[0], index[1], index[2]) = S.healthFlag;
        }
        catch (ConservationLaw::StateFailure& stateFailure)
        {
            stateFailure.setZoneIndex (uit.relativeIndex());
            inversionFailure.failedStates.push_back (stateFailure);
        }
    }
    if (inversionFailure.failedStates.size() > 0)
    {
        inversionFailure.updateWhatMessage();
        throw inversionFailure;
    }
    return zoneHealth;
}

Array FieldOperator::recoverPrimitive (Array::Reference U) const
{
    auto P = Array (U.shape());
    recoverPrimitive (U, P);
    return P;
}

void FieldOperator::generateConserved (Array::Reference P, Array::Reference U) const
{
    int numConserved = law->getNumConserved();
    auto request = ConservationLaw::Request();

    // This ensures that we are stepping over the axis (3) with the field
    // components.
    // ------------------------------------------------------------------------
    auto Treg = P.getRegion().withStride (3, law->getNumConserved());
    auto Preg = P.getArray()[Treg];
    auto Ureg = U.getArray()[Treg];
    auto pit = Preg.begin();
    auto uit = Ureg.begin();

    for ( ; uit != Ureg.end(); ++pit, ++uit)
    {
        auto S = law->fromPrimitive (request, pit);

        for (int q = 0; q < numConserved; ++q)
        {
            uit[q] = S.U[q];
        }
    }
}

Array FieldOperator::generateConserved (Array::Reference P) const
{
    auto U = Array (P.shape());
    generateConserved (P, U);
    return U;
}

double FieldOperator::getCourantTimestep (
    Array::Reference P,
    Array::Reference L) const
{
    if (P.size(0) != L.size(0)
     || P.size(1) != L.size(1)
     || P.size(2) != L.size(2))
    {
        throw std::logic_error ("Primitive and linear cell dimension arrays have different sizes");
    }

    auto request0 = ConservationLaw::Request().oriented (UnitVector::xhat);
    auto request1 = ConservationLaw::Request().oriented (UnitVector::yhat);
    auto request2 = ConservationLaw::Request().oriented (UnitVector::zhat);
    double courantTimestep = std::numeric_limits<double>::max();

    // This ensures that we are stepping over the axis (3) with the field
    // components.
    // ------------------------------------------------------------------------
    auto Treg = P.getRegion().withStride (3, law->getNumConserved());
    auto Preg = P.getArray()[Treg];
    auto pit = Preg.begin();
    auto lit = L.begin();

    for ( ; pit != Preg.end(); ++pit, ++lit)
    {
        auto S0 = law->fromPrimitive (request0, pit);
        auto S1 = law->fromPrimitive (request1, pit);
        auto S2 = law->fromPrimitive (request2, pit);

        auto lambdas = std::array<double, 3>
        {{
            law->maxEigenvalueMagnitude (S0),
            law->maxEigenvalueMagnitude (S1),
            law->maxEigenvalueMagnitude (S2)
        }};
        const double dt = *lit / *std::max_element (lambdas.begin(), lambdas.end());

        if (dt < courantTimestep)
        {
            courantTimestep = dt;
        }
    }
    return courantTimestep;
}

std::vector<double> FieldOperator::volumeIntegratedDiagnostics (
    Array::Reference P,
    Array::Reference V) const
{
    auto request = ConservationLaw::Request();
    auto diagnostics = std::vector<double> (law->getDiagnosticNames().size());

    if (P.size(0) != V.size(0)
     || P.size(1) != V.size(1)
     || P.size(2) != V.size(2))
    {
        throw std::logic_error ("Primitive and cell volume arrays have different sizes");
    }

    // This ensures that we are stepping over the axis (3) with the field
    // components.
    // ------------------------------------------------------------------------
    auto Treg = P.getRegion().withStride (3, law->getNumConserved());
    auto Preg = P.getArray()[Treg];
    auto pit = Preg.begin();
    auto vit = V.begin();

    for ( ; pit != Preg.end(); ++pit, ++vit)
    {
        auto state = law->fromPrimitive (request, pit);
        auto celld = law->makeDiagnostics (state);
        const double Vol = *vit;

        for (unsigned int n = 0; n < celld.size(); ++n)
        {
            diagnostics[n] += celld[n] * Vol;
        }
    }
    return diagnostics;
}
