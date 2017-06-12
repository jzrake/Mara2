#include "Array.hpp"
#include <functional>



/**
A class to facilitate general stencil operations on logically cartesian data.
Such operations are a mapping from a rectangular "footprint" within one or
more arrays of (identically sized) source data to a single point in an array
of output data. The dimensions of source data are determined by the size of
the last two axes in the input array, while the dimensions of output data are
set using the setCodomainRank() method. The returned array is smaller than the
input arrays by the size of the footprint.
*/
class Stencil
{
public:
    using Array = Cow::Array;
    using Shape = Cow::Shape;
    using Region = Cow::Region;
    using StencilOperation1 = std::function<void (const Array&, Array&)>;
    using StencilOperation2 = std::function<void (const Array&, const Array&, Array&)>;
    using StencilOperation3 = std::function<void (const Array&, const Array&, Array&, Array&)>;

    /**
    At construction the stencil object has zero footprint and codomain rank
    (1, 1).
    */
    Stencil();

    /**
    Set the lower bounds of the stencil's footprint. For example, if i0 is set
    to -1, then the stencil extends one point to the left of the target point.
    */
    void setFootprintLower (int i0, int j0, int k0);

    /**
    Set the upper bounds of the stencil's footprint. For example, if i1 is set
    to +1, then the stencil extends one point to the right of the target
    point. Note: the upper bound is inclusive; the stencil width on the first
    axis is i1 - i0 + 1.
    */
    void setFootprintUpper (int i1, int j1, int k1);

    /**
    Set the shape of the output array data, that is the vector space rank of
    the operation's target domain, R^N by R^M.
    */
    void setCodomainRank (int M, int N);

    /**
    Execute the stencil operation over the given source array. The operation's
    first argument (the array of local stencil data) will have size:

    (2 * F[0] + 1, 2 * F[1] + 1, 2 * F[2] + 1, R, S)

    where F is the footprint and R and S depend are size(3) and size(4) of the
    source array x supplied to this function. The output array will have size

    (S[0] - 2 * F[0], S[1] - 2 * F[1], S[2] - 2 * F[2], M, N)

    where S is the shape of the source array, and M and N are the codomain
    ranks.
    */
    Array evaluate (StencilOperation1 op, const Array& x) const;

    /**
    Execute the stencil over two input arrays.
    */
    Array evaluate (StencilOperation2 op, const Array& x, const Array& y) const;

    /**
    Execute the stencil over three input arrays.
    */
    Array evaluate (StencilOperation3 op, const Array& x, const Array& y, const Array& z) const;

private:
    /* @internal */
    Array runStencil (Shape, std::function<void (const Region&, Array&)>) const;
    std::vector<Array> makeStencilDataArrays (std::vector<const Array*>, Shape&) const;
    std::array<int, 3> footprintLower;
    std::array<int, 3> footprintUpper;
    int tensorM;
    int tensorN;
};
