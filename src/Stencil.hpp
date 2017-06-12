#include "Array.hpp"
#include <functional>



/**
A class to encapulate colocated, symmetric stencil operations on logically
cartesian data (the output is located at the same mesh locations as the input
data). For these operations the footprint size is the number of mesh locations
to either side of the target location. So, including the target location the
stencil size is odd on each axis.
*/
class Stencil
{
public:
    using StencilOperation = std::function<void (const Cow::Array& stencilData, Cow::Array& output)>;

    /**
    At construction the stencil object is in a deterministic bun unusable
    state. Each of the three set functions below need to be called.
    */
    Stencil();

    /**
    Set the shape of the stencil footprint. For example, if evaluating a
    fourth order 3D gradient, the footprint might be (2, 2, 2). This ensures
    that 2 mesh locations exist to either side of the target location, i.e.
    there the stencil data will have a shape (5, 5, 5).
    */
    void setFootprint (int f1, int f2, int f3);

    /**
    Set the shape of the output array data, that is the vector space rank of
    the operation's target domain, R^N by R^M.
    */
    void setCodomainRank (int M, int N);

    /**
    Assign an operation to be called on the stencil data. The operation's
    first argument (the stencil data array) will have size:

    (2 * F[0] + 1, 2 * F[1] + 1, 2 * F[2] + 1, R, S)

    where F is the footprint and R and S depend are size(3) and size(4) of the
    source array supplied to the evaluate function.
    */
    void setOperation (StencilOperation operationToUse);

    /**
    Return the footprint size on a given axis.
    */
    int getFootprint (int axis) const { return footprint[axis]; }

    /**
    Execute the stencil operation over the given source array. The output
    array will have size

    (S[0] - 2 * F[0], S[1] - 2 * F[1], S[2] - 2 * F[2], M, N)

    where S is the shape of the source array, and M and N are the codomain
    ranks.
    */
    Cow::Array evaluate (const Cow::Array& source) const;

private:
    int tensorM;
    int tensorN;
    Cow::Shape footprint;
    StencilOperation operation;
};
