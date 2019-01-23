#include "Mara.hpp"
#include "Array.hpp"




/**
Class responsible for deploying point-wise (non-stencil-based) operations over
field variables. Services include converting arrays of conserved quanitites to
primitives and vice-versa, and determining the Courant condition for a patch
of data. Methods in this class must be agnostic to the mesh geometry, and thus
to the shape of the input data. However (at least for now) the input arrays
are assumed to have a layout where axes 0, 1, and 2 are spatial, while axis 3
has components of the primitive variable data.
*/
class FieldOperator
{
public:
    using Array = Cow::Array;

    /**
    An exception class that holds a list of zones that have failed primitive
    variable inversion.
    */
    class InversionFailure : public std::exception
    {
    public:
        const char* what() const noexcept override;
        void updateWhatMessage();
        std::string whatMessage;
        std::vector<ConservationLaw::StateFailure> failedStates;
    };

    FieldOperator (std::shared_ptr<ConservationLaw> law=nullptr);

    /**
    Set the conservation law to be used.
    */
    void setConservationLaw (std::shared_ptr<ConservationLaw>);

    std::shared_ptr<ConservationLaw> getConservationLaw() { return law; }

    /**
    Deploy the ConservationLaw::fromConserved function over the array of
    conserved variables U. Note that this function may throw an exception if
    the primitive variable recovery fails for any state. Even if the recovery
    formally succeeded, the ConservationLaw might have set a health flag on a
    given state, indicating that some tricks (e.g. pressure floor) were
    invoked. The array of health flags is returned in the PrimitiveRecovery
    struct.
    */
    Array recoverPrimitive (Array::Reference U, Array::Reference P, std::array<double, 8> t={}) const;

    /**
    Convenience function; Returns P rather than zone health array.
    */
    Array recoverPrimitive (Array::Reference U, std::array<double, 8> t={}) const;

    /**
    Convert an array of primitive variables to conserved variables.
    */
    void generateConserved (Array::Reference P, Array::Reference U, std::array<double, 8> t={}) const;

    /**
    Convenience function.
    */
    Array generateConserved (Array::Reference P, std::array<double, 8> t={}) const;

    /**
    Determine the shortest wave-propagation time for the given patch of
    primitive data. The array of linear cell dimensions L are prepared in
    advance so that this class remains agnostic to the mesh geometry.

    @param P                An array of primitive data, shape (ni, nj, nk, nq).

    @param L                An array of linear cell dimensions, shape is same
                            as P on axis 0, 1, and 2, but has only one
                            component on axis 3.
    */
    double getCourantTimestep (Array::Reference P, Array::Reference L) const;

    /**
    Return the volume integral (not volume average) of named diagnostics
    defined by the ConservationLaw.

    @param P                An array of primitive data, shape (ni, nj, nk, nq).

    @param V                An array of linear cell volumes, shape is same
                            as P on axis 0, 1, and 2, but has only one
                            component on axis 3.
    */
    std::vector<double> volumeIntegratedDiagnostics (Array::Reference P, Array::Reference V) const;

private:
    std::shared_ptr<ConservationLaw> law;
};
