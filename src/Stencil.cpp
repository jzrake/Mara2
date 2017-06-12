#include "Stencil.hpp"
#include "Array.hpp"
#include <functional>

using namespace Cow;




// ============================================================================
Stencil::Stencil()
{
    tensorM = 1;
    tensorN = 1;
    footprint[0] = 1;
    footprint[1] = 1;
    footprint[2] = 1;
}

void Stencil::setFootprint (int f1, int f2, int f3)
{
    footprint[0] = f1;
    footprint[1] = f2;
    footprint[2] = f3;
}

void Stencil::setCodomainRank (int M, int N)
{
    tensorM = M;
    tensorN = N;
}

void Stencil::setOperation (StencilOperation operationToUse)
{
    operation = operationToUse;
}

Array Stencil::evaluate (const Array& source) const
{
    auto resultShape = Shape
    {{
        source.size(0) - 2 * footprint[0],
        source.size(1) - 2 * footprint[1],
        source.size(2) - 2 * footprint[2],
        tensorM,
        tensorN
    }};

    auto stdataShape = Shape {{
        2 * footprint[0] + 1,
        2 * footprint[1] + 1,
        2 * footprint[2] + 1,
        tensorM,
        tensorN
    }};

    auto result = Array (resultShape);
    auto stdata = Array (stdataShape);
    auto locres = Array (tensorM, tensorN);
    auto footprintRegion = Region();

    for (int i = footprint[0]; i < source.size(0) - 2 * footprint[0]; ++i)
    {
        footprintRegion.lower[0] = i - footprint[0];
        footprintRegion.upper[0] = i + footprint[0];

        for (int j = footprint[1]; j < source.size(1) - 2 * footprint[1]; ++j)
        {
            footprintRegion.lower[1] = j - footprint[1];
            footprintRegion.upper[1] = j + footprint[1];

            for (int k = footprint[2]; k < source.size(2) - 2 * footprint[2]; ++k)
            {
                footprintRegion.lower[2] = k - footprint[2];
                footprintRegion.upper[2] = k + footprint[2];

                stdata = source.extract (footprintRegion);
                operation (stdata, locres);

                for (int m = 0; m < tensorM; ++m)
                {
                    for (int n = 0; n < tensorN; ++n)
                    {
                        result (i, j, k, m, n) = locres (m, n);
                    }
                }
            }
        }
    }
    return result;
}
