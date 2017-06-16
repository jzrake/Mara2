#include <iostream> // DEBUG
#include "Stencil.hpp"
#include "Array.hpp"

using namespace Cow;




// ============================================================================
Stencil::Stencil()
{
    tensorM = 1;
    tensorN = 1;
    footprintLower[0] = 0;
    footprintLower[1] = 0;
    footprintLower[2] = 0;
    footprintUpper[0] = 0;
    footprintUpper[1] = 0;
    footprintUpper[2] = 0;
}

void Stencil::setFootprintLower (int i0, int j0, int k0)
{
    footprintLower[0] = i0;
    footprintLower[1] = j0;
    footprintLower[2] = k0;
}

void Stencil::setFootprintUpper (int i1, int j1, int k1)
{
    footprintUpper[0] = i1;
    footprintUpper[1] = j1;
    footprintUpper[2] = k1;
}

void Stencil::setCodomainRank (int M, int N)
{
    tensorM = M;
    tensorN = N;
}

Array Stencil::evaluate (StencilOperation1 operation, const Array& x) const
{
    auto inputDataArrays = std::vector<const Array*> { &x };
    auto resultShape = Shape();
    auto stencilDataArrays = makeStencilDataArrays (inputDataArrays, resultShape);

    auto F = [&] (const Region& footprintRegion, Array& localResult)
    {
        stencilDataArrays[0].copyFrom (x, Region(), footprintRegion);
        operation (stencilDataArrays[0], localResult);
    };
    return runStencil (resultShape, F);
}

Array Stencil::evaluate (StencilOperation2 operation, const Array& x, const Array& y) const
{
    auto inputDataArrays = std::vector<const Array*> { &x, &y };
    auto resultShape = Shape();
    auto stencilDataArrays = makeStencilDataArrays (inputDataArrays, resultShape);

    auto F = [&] (const Region& footprintRegion, Array& localResult)
    {
        stencilDataArrays[0].copyFrom (x, Region(), footprintRegion);
        stencilDataArrays[1].copyFrom (y, Region(), footprintRegion);
        operation (stencilDataArrays[0], stencilDataArrays[1], localResult);
    };
    return runStencil (resultShape, F);
}

Array Stencil::evaluate (StencilOperation3 operation, const Array& x, const Array& y, const Array& z) const
{
    auto inputDataArrays = std::vector<const Array*> { &x, &y, &z };
    auto resultShape = Shape();
    auto stencilDataArrays = makeStencilDataArrays (inputDataArrays, resultShape);

    auto F = [&] (const Region& footprintRegion, Array& localResult)
    {
        stencilDataArrays[0].copyFrom (x, Region(), footprintRegion);
        stencilDataArrays[1].copyFrom (y, Region(), footprintRegion);
        stencilDataArrays[2].copyFrom (z, Region(), footprintRegion);
        operation (stencilDataArrays[0], stencilDataArrays[1], stencilDataArrays[2], localResult);
    };
    return runStencil (resultShape, F);
}

std::vector<Array> Stencil::makeStencilDataArrays (std::vector<const Array*> inputDataArrays, Shape &resultShape) const
{
    auto stencilDataArrays = std::vector<Array>();

    for (const auto& inputDataArray : inputDataArrays)
    {
        if (   inputDataArray->size(0) != inputDataArrays[0]->size(0)
            || inputDataArray->size(1) != inputDataArrays[0]->size(1)
            || inputDataArray->size(2) != inputDataArrays[0]->size(2))
        {
            throw std::logic_error ("input arrays do not have the same shape");
        }
        stencilDataArrays.push_back (Array (Shape {{
            footprintUpper[0] - footprintLower[0] + 1,
            footprintUpper[1] - footprintLower[1] + 1,
            footprintUpper[2] - footprintLower[2] + 1,
            inputDataArray->size(3),
            inputDataArray->size(4)
        }}));
    }

    resultShape = Shape
    {{
        inputDataArrays[0]->size(0) - (footprintUpper[0] - footprintLower[0]),
        inputDataArrays[0]->size(1) - (footprintUpper[1] - footprintLower[1]),
        inputDataArrays[0]->size(2) - (footprintUpper[2] - footprintLower[2]),
        tensorM,
        tensorN
    }};

    if (   resultShape[0] <= 0
        || resultShape[1] <= 0
        || resultShape[2] <= 0)
    {
        throw std::logic_error ("stencil footprint too large for input array data");
    }
    return stencilDataArrays;
}

Array Stencil::runStencil (
    Shape resultShape,
    std::function<void (
        const Region& footprintRegion,
        Array& localResult)> loadDataAndCall) const
{
    auto result = Array (resultShape);
    auto locres = Array (tensorM, tensorN);
    auto footprintRegion = Region();

    for (int i = 0; i < resultShape[0]; ++i)
    {
        footprintRegion.lower[0] = i;
        footprintRegion.upper[0] = i + (footprintUpper[0] - footprintLower[0]) + 1;

        for (int j = 0; j < resultShape[1]; ++j)
        {
            footprintRegion.lower[1] = j;
            footprintRegion.upper[1] = j + (footprintUpper[1] - footprintLower[1]) + 1;

            for (int k = 0; k < resultShape[2]; ++k)
            {
                footprintRegion.lower[2] = k;
                footprintRegion.upper[2] = k + (footprintUpper[2] - footprintLower[2]) + 1;

                loadDataAndCall (footprintRegion, locres);

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
