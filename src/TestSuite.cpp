#include <iostream>
#define CATCH_CONFIG_FAST_COMPILE
#include "Catch.hpp"
#include "Array.hpp"

using namespace Cow;




// ============================================================================
#include "Array.hpp"

SCENARIO ("Heap allocations can be move-constructed and move-assigned", "[HeapAllocation]")
{
    GIVEN ("A HeapAllocation H(10) and another I(0)")
    {
        auto H = HeapAllocation(10);
        auto I = HeapAllocation();
        auto J = HeapAllocation(100);

        WHEN ("H is copy-assigned to I and J")
        {
            I = H;
            J = H;
            THEN ("H, I, and J all have size 10")
            {
                CHECK (H.size() == 10);
                CHECK (I.size() == 10);
                CHECK (J.size() == 10);
            }
        }

        WHEN ("H is std::move'd to I")
        {
            I = std::move(H);

            THEN ("H has zero size and I has size 10")
            {
                CHECK (H.size() == 0);
                CHECK (I.size() == 10);
            }
        }

        WHEN ("H is std::move'd to J")
        {
            J = std::move(H);

            THEN ("H has zero size and J has size 10")
            {
                CHECK (H.size() == 0);
                CHECK (J.size() == 10);
            }
        }

        WHEN ("H is std::move constructed to a new allocation J")
        {
            auto J = std::move(H);

            THEN ("H has zero size and J has size 10")
            {
                CHECK (H.size() == 0);
                CHECK (J.size() == 10);
            }
        }
    }
}

SCENARIO ("Arrays should work with slicing and indexing", "[Array]")
{
    GIVEN ("An array with 100 elements")
    {
        auto A = Array (100);

        THEN ("The array is zero-initialized")
        {
            for (auto x : A)
            {
                CHECK (x == 0.0);
            }
        }

        WHEN ("The array is indexed out-of-bounds")
        {
            THEN ("We get an exception")
            {
                if (! Array::isBoundsCheckDisabled())
                {
                    CHECK_THROWS (A (-1));
                    CHECK_THROWS (A (101));
                }
            }
        }
    }

    GIVEN ("An array A with 100 elements, B with 10, a contiguous region R, and a strided one S")
    {
        auto A = Array (100);
        auto B = Array (10);
        auto R = Region (A.shape()).withRange (0, 0, 10, 1);
        auto S = Region (A.shape()).withRange (0, 0, 100, 10);
        auto T = Region (A.shape()).withRange (0, 0, 100, 20);

        for (int i = 0; i < B.size(); ++i)
        {
            B(i) = i;
        }

        THEN ("R and S have the same shape")
        {
            CHECK (R.shape() == S.shape());
        }

        THEN ("S and T do not have the same shape")
        {
            CHECK (S.shape() != T.shape());
        }

        THEN ("B can be inserted into A[R]")
        {
            CHECK_NOTHROW (A[R] = B);

            for (int i = 0; i < 10; ++i)
            {
                CHECK (A(i) == B(i));
            }
        }

        THEN ("B can be inserted into A[S]")
        {
            CHECK_NOTHROW (A[S] = B);

            for (int i = 0; i < 100; i += 10)
            {
                CHECK (A(i) == B(i / 10));
            }
        }

        THEN ("B cannot be inserted into A[T]")
        {
            CHECK_THROWS (A[T] = B);
        }
    }

    GIVEN ("An array A with 100 elements, B with 100, and two non-overlapping regions R and S")
    {
        auto A = Array (100);
        auto B = Array (100);
        auto R = Region().withRange (0, 0, 10);
        auto S = Region().withRange (0, 90, 100);
        auto T = Region().withRange (0, 90, 101);

        for (int i = 0; i < B.size(); ++i)
        {
            B(i) = i;
        }

        THEN ("The last 10 elements of B can be loaded into the first of A")
        {
            CHECK_NOTHROW (A[R] = B[S]);

            for (int i = 0; i < 10; ++i)
            {
                CHECK (A(i) == B(90 + i));
            }
        }

        THEN ("An out-of-bounds region cannot be loaded")
        {
            CHECK_THROWS (A[R] = B[T]);
        }
    }

    GIVEN ("An array A with [8,8,8,1,3] elements, B with [8,8,8,1,1]")
    {
        auto A = Array (8, 8, 8, 1, 3);
        auto B = Array (8, 8, 8, 1, 1);

        for (auto& x : B)
        {
            x = 1.;
        }

        auto R0 = Region(A.shape()).withRange (4, 0, 1);
        auto R1 = Region(A.shape()).withRange (4, 1, 2);
        auto R2 = Region(A.shape()).withRange (4, 2, 3);
        auto R3 = Region(A.shape()).withRange (4, 3, 4); // bad region
        auto R4 = Region(A.shape()).withRange (3, 0, -2); // bad region

        THEN ("B can be inserted into axis 0 of A")
        {
            CHECK_NOTHROW (A[R0] = B);

            for (auto it = A[R0].begin(); it != A.end(); ++it)
            {
                if (it.index()[4] == 0)
                {
                    CHECK (*it == 1.0);
                }
                else
                {
                    CHECK (*it == 0.0);
                }
            }
        }

        THEN ("B can be inserted into axis 1 of A")
        {
            CHECK_NOTHROW (A[R1] = B);

            for (auto it = A[R1].begin(); it != A.end(); ++it)
            {
                if (it.index()[4] == 1)
                {
                    CHECK (*it == 1.0);
                }
                else
                {
                    CHECK (*it == 0.0);
                }
            }
        }

        THEN ("B can be inserted into axis 2 of A")
        {
            CHECK_NOTHROW (A[R2] = B);

            for (auto it = A[R2].begin(); it != A.end(); ++it)
            {
                if (it.index()[4] == 2)
                {
                    CHECK (*it == 1.0);
                }
                else
                {
                    CHECK (*it == 0.0);
                }
            }
        }

        THEN ("B cannot be inserted into axis 3 of A (doesn't exist)")
        {
            if (! Array::isBoundsCheckDisabled())
            {
                CHECK_THROWS (A[R3] = B);
            }
        }

        THEN ("A cannot be referenced with a bad region")
        {
            if (! Array::isBoundsCheckDisabled())
            {
                CHECK_THROWS (A[R4] = B);
            }
        }
    }
}

SCENARIO ("Arrays can be move-constructed and move-assigned", "[Array]")
{
    GIVEN ("An array A(10) and another B(0)")
    {
        auto A = Array(10);
        auto B = Array();

        WHEN ("A is copy-constructed to B")
        {
            B = A;

            THEN ("A has size 10 and B has size 10")
            {
                CHECK (A.size() == 10);
                CHECK (B.size() == 10);
            }
        }

        WHEN ("A is copy-assigned to C")
        {
            auto C = A;

            THEN ("A has size 10 and C has size 10")
            {
                CHECK (A.size() == 10);
                CHECK (C.size() == 10);
            }
        }

        WHEN ("A is std::move'd to B")
        {
            B = std::move(A);
            THEN ("A has zero size and B has size 10")
            {
                CHECK (A.size() == 0);
                CHECK (B.size() == 10);
            }
        }

        WHEN ("A is std::move constructed to a new array C")
        {
            auto C = std::move(A);

            THEN ("A has zero size and C has size 10")
            {
                CHECK (A.size() == 0);
                CHECK (C.size() == 10);
            }
        }
    }
}

SCENARIO ("Shape3D objects should work correctly", "[Shape3D]")
{
    GIVEN ("Three shapes A, B, C where A contains B but not C")
    {
        auto A = Shape3D (4, 2, 8);
        auto B = Shape3D (4, 0, 8);
        auto C = Shape3D (1, 1, 9);

        THEN ("The contains member function works as expected")
        {
            CHECK (A.contains(A));
            CHECK (A.contains(B));
            CHECK_FALSE (A.contains(C));
            CHECK_FALSE (B.contains(C));
            CHECK_FALSE (C.contains(A));
        }

        WHEN ("C is reduced by one on the last axis")
        {
            C = C.reduced (2, 1);

            THEN ("A now contains C")
            {
                CHECK (A.contains(C));
            }
        }
    }
}




// ============================================================================
#include "Variant.hpp"

SCENARIO ("Variant works as expected", "[Variant]")
{
    GIVEN ("A Variant::NamedValues instance and some mock command line arguments")
    {
        auto user = Variant::NamedValues();
        user["thing1"] = 1;
        user["thing2"] = "2";

        const char* argv[] = {"thing1=2", "thing2=3"};
        Variant::updateFromCommandLine (user, 2, argv);

        THEN ("User parameters have the correct values")
        {
            REQUIRE (user["thing1"].getType() == 'i');
            REQUIRE (int (user["thing1"]) == 2);
            REQUIRE (user["thing2"].getType() == 's');
            REQUIRE (std::string (user["thing2"]) == "3");
        }
    }
}




// ============================================================================
#include "UnitVector.hpp"
#include "RotationMatrix.hpp"

SCENARIO ("Rotation matrices should behave appropriately with unit vectors", "[UnitVector]")
{
    GIVEN ("The unit vector zhat")
    {
        UnitVector xhat = UnitVector::normalizeFrom (1, 0, 0);
        UnitVector yhat = UnitVector::normalizeFrom (0, 1, 0);
        UnitVector zhat = UnitVector::normalizeFrom (0, 0, 1);

        THEN ("Y (pi / 2) zhat = xhat")
        {
            REQUIRE ((RotationMatrix::aboutY (M_PI / 2) * zhat).pitchAngleMu == Approx (xhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutY (M_PI / 2) * zhat).azimuthalAnglePhi == Approx (xhat.azimuthalAnglePhi));
        }

        THEN ("X (-pi / 2) zhat = yhat")
        {
            REQUIRE ((RotationMatrix::aboutX (-M_PI / 2) * zhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutX (-M_PI / 2) * zhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("Z (pi / 2) xhat = yhat")
        {
            REQUIRE ((RotationMatrix::aboutZ (M_PI / 2) * xhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutZ (M_PI / 2) * xhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (xhat) = xhat")
        {
            REQUIRE (zhat.withPolarAxis (xhat).pitchAngleMu == Approx (xhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (xhat).azimuthalAnglePhi == Approx (xhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (yhat) = yhat")
        {
            REQUIRE (zhat.withPolarAxis (yhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (yhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (zhat) = zhat")
        {
            REQUIRE (zhat.withPolarAxis (zhat).pitchAngleMu == Approx (zhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (zhat).azimuthalAnglePhi == Approx (zhat.azimuthalAnglePhi));
        }
    }

    REQUIRE (UnitVector::xhat != UnitVector::yhat);
    REQUIRE (UnitVector::xhat == UnitVector::xhat);
    REQUIRE (UnitVector::yhat == UnitVector::yhat);
    REQUIRE (UnitVector::zhat == UnitVector::zhat);
}




// ============================================================================
#include "QuarticPolynomial.hpp"

SCENARIO ("Quartic solver should work correctly", "[QuarticPolynomial]")
{
    GIVEN ("A polynomial with 4 easy roots -2, -1, 1, 2")
    {
        QuarticPolynomial q (1, 0, -5, 0, 4);

        THEN ("The roots are found correctly to one part in 10^14")
        {
            double roots[4];
            const int nr = q.solve (roots);
            CHECK (nr == 4);
            CHECK (roots[0] == Approx(-2).epsilon (1e-14));
            CHECK (roots[1] == Approx(-1).epsilon (1e-14));
            CHECK (roots[2] == Approx(+1).epsilon (1e-14));
            CHECK (roots[3] == Approx(+2).epsilon (1e-14));
        }
    }
    GIVEN ("A polynomial with 4 challenging roots -1, -1e-5, 1e-5, 1")
    {
        QuarticPolynomial q (1., 0., -1., 0., 1e-10);

        THEN ("The roots are found correctly to one part in 10^10")
        {
            double roots[4];
            const int nr = q.solve (roots);
            CHECK (nr == 4);
            CHECK (roots[0] == Approx(-1).epsilon (1e-10));
            CHECK (roots[1] == Approx(-1e-5).epsilon (1e-10));
            CHECK (roots[2] == Approx(+1e-5).epsilon (1e-10));
            CHECK (roots[3] == Approx(+1).epsilon (1e-10));
        }
    }
    GIVEN ("A polynomial with 2 easy roots -1, 1")
    {
        QuarticPolynomial q (-1, 0, 0, 0, 1);

        THEN ("The roots are found correctly to one part in 10^10")
        {
            double roots[4];
            const int nr = q.solve (roots);
            CHECK (nr == 2);
            CHECK (roots[0] == Approx(-1).epsilon (1e-14));
            CHECK (roots[1] == Approx(+1).epsilon (1e-14));
        }
    }
}




// ============================================================================
#include "Stencil.hpp"

SCENARIO ("Stencil operations should behave correctly", "[Stencil]")
{
    GIVEN ("A stencil with shape [2, 3, 4]")
    {
        auto stencilShape = Shape();
        auto stencil = Stencil();
        auto operation = [&] (const Array& sd, Array& od) { stencilShape = sd.shape(); };
        stencil.setFootprintLower (-2,-3,-4);
        stencil.setFootprintUpper ( 2, 3, 4);
        stencil.setCodomainRank (4, 5);

        WHEN ("The input array is [12, 12, 12] and the codomain rank is [4, 5]")
        {
            auto source = Array (12, 12, 12);
            auto result = stencil.evaluate (operation, source);

            THEN ("The output array shape is [8, 6, 4, 4, 5]")
            {
                CHECK (result.size(0) == 8);
                CHECK (result.size(1) == 6);
                CHECK (result.size(2) == 4);
                CHECK (result.size(3) == 4);
                CHECK (result.size(4) == 5);
            }

            THEN ("The array of stencil data has shape [5, 7, 9, 1, 1]")
            {
                CHECK (stencilShape[0] == 5);
                CHECK (stencilShape[1] == 7);
                CHECK (stencilShape[2] == 9);
                CHECK (stencilShape[3] == 1);
                CHECK (stencilShape[4] == 1);
            }
        }
    }

    GIVEN ("A stencil with shape [2, 2, 0]")
    {
        auto stencilShape = Shape();
        auto stencil = Stencil();
        auto operation = [&] (const Array& sd, Array& od) { stencilShape = sd.shape(); };
        auto source = Array (12, 12, 1, 8);
        stencil.setFootprintLower (-2,-2, 0);
        stencil.setFootprintUpper ( 2, 2, 0);
        stencil.setCodomainRank (3, 1);

        WHEN ("The input array is [12, 12, 1, 8] and the codomain rank is [3, 1]")
        {
            auto result = stencil.evaluate (operation, source);

            THEN ("The output array shape is [8, 8, 1, 3, 1]")
            {
                CHECK (result.size(0) == 8);
                CHECK (result.size(1) == 8);
                CHECK (result.size(2) == 1);
                CHECK (result.size(3) == 3);
                CHECK (result.size(4) == 1);
            }

            THEN ("The array of stencil data has shape [5, 5, 1, 8, 1]")
            {
                CHECK (stencilShape[0] == 5);
                CHECK (stencilShape[1] == 5);
                CHECK (stencilShape[2] == 1);
                CHECK (stencilShape[3] == 8);
                CHECK (stencilShape[4] == 1);
            }
        }

        WHEN ("The stencil footprint is changed to [2, 2, 1]")
        {
            stencil.setFootprintLower (-2,-2,-1);
            stencil.setFootprintUpper ( 2, 2, 1);

            THEN ("An exception is raised because the source data shape is too small")
            {
                CHECK_THROWS_AS (stencil.evaluate (operation, source), std::logic_error);
            }
        }
    }

    GIVEN ("A second order finite difference stencil")
    {
        auto stencil = Stencil();
        auto x = Array (256);
        auto y = Array (256);

        for (int i = 0; i < 256; ++i)
        {
            x[i] = (i + 0.5) / 256;
            y[i] = std::sin (x[i]);
        }
        stencil.setCodomainRank (2, 1);

        WHEN ("The stencil is centered, of the type [-1, 0, +1]")
        {
            stencil.setFootprintLower (-1, 0, 0);
            stencil.setFootprintUpper ( 1, 0, 0);

            auto operation = [] (const Array& x0, const Array& y0, Array& deriv)
            {
                deriv[0] = (y0[2] - y0[0]) / (x0[2] - x0[0]);
                deriv[1] = std::cos (x0[1]);
            };
            auto deriv = stencil.evaluate (operation, x, y);

            THEN ("The resulting array is smaller by two data points")
            {
                CHECK (deriv.size(0) == 254);
            }

            THEN ("The derivative approximation is correct to one part in 10^5")
            {
                CHECK (deriv (  0, 0, 0, 0) == Approx (deriv (  0, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (128, 0, 0, 0) == Approx (deriv (128, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (253, 0, 0, 0) == Approx (deriv (253, 0, 0, 1)).epsilon (1e-5));
            }
        }

        WHEN ("The stencil is centered, of the type [0, +1]")
        {
            stencil.setFootprintLower (0, 0, 0);
            stencil.setFootprintUpper (1, 0, 0);

            auto operation = [] (const Array& x0, const Array& y0, Array& deriv)
            {
                deriv[0] = (y0[1] - y0[0]) / (x0[1] - x0[0]);
                deriv[1] = std::cos (0.5 * (x0[0] + x0[1]));
            };
            auto deriv = stencil.evaluate (operation, x, y);

            THEN ("The resulting array is smaller by one data point")
            {
                CHECK (deriv.size(0) == 255);
            }

            THEN ("The derivative approximation is correct to one part in 10^5")
            {
                CHECK (deriv (  0, 0, 0, 0) == Approx (deriv (  0, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (128, 0, 0, 0) == Approx (deriv (128, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (254, 0, 0, 0) == Approx (deriv (254, 0, 0, 1)).epsilon (1e-5));
            }
        }
    }
}




// ============================================================================
#include "ConservationLaws.hpp"

SCENARIO ("Relativistic MHD class produces reasonable results", "[ConservationLaws]")
{
    auto nrmhd = NewtonianMHD();
    auto srmhd = RelativisticMHD();
    ConservationLaw::Request R;

    GIVEN ("A state with cs << c and B=0")
    {
        double P[8] = {1., 0., 0., 0., 1e-6, 0, 0, 0};

        WHEN ("The eigenvalues from Newtonian and relativistic MHD are compared")
        {
            THEN ("They are approximately equal")
            {
                for (int d = 0; d < 3; ++d)
                {
                    R.areaElement[0] = 0;
                    R.areaElement[1] = 0;
                    R.areaElement[2] = 0;
                    R.areaElement[d] = 1;

                    auto S1 = nrmhd.fromPrimitive (R, P);
                    auto S2 = srmhd.fromPrimitive (R, P);

                    for (int i = 0; i < 8; ++i)
                    {
                        CHECK (S1.A[i] == Approx (S2.A[i]));
                    }
                }
            }
        }
        WHEN ("We try to convert from prim to cons and back again")
        {
            auto S1 = srmhd.fromPrimitive (R, P);
            auto S2 = srmhd.fromConserved (R, S1.U.begin());

            for (int i = 0; i < 8; ++i)
            {
               CHECK (P[i] == Approx (S2.P[i]));
            }
        }
    }
    GIVEN ("A state with cs << c, v << c, and B << rho c^2")
    {
        double P[8] = {1., 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};

        WHEN ("The eigenvalues from Newtonian and relativistic MHD are compared")
        {
            auto S1 = nrmhd.fromPrimitive (R, P);
            auto S2 = srmhd.fromPrimitive (R, P);

            THEN ("They are approximately equal")
            {
                for (int d = 0; d < 3; ++d)
                {
                    R.areaElement[0] = 0;
                    R.areaElement[1] = 0;
                    R.areaElement[2] = 0;
                    R.areaElement[d] = 1;

                    auto S1 = nrmhd.fromPrimitive (R, P);
                    auto S2 = srmhd.fromPrimitive (R, P);

                    for (int i = 0; i < 8; ++i)
                    {
                        CHECK (S1.A[i] == Approx (S2.A[i]));
                    }
                }
            }
        }
        WHEN ("We try to convert from prim to cons and back again")
        {
            auto S1 = srmhd.fromPrimitive (R, P);
            auto S2 = srmhd.fromConserved (R, S1.U.begin());

            for (int i = 0; i < 8; ++i)
            {
               CHECK (P[i] == Approx (S2.P[i]).epsilon (1e-12));
            }
        }
    }
    GIVEN ("An state with cs ~ c, v ~ c, and B ~ rho c^2")
    {
        double P[8] = {1., 0.5, 0.4, 0.3, 1, 1, 1, 1};

        WHEN ("We try to convert from prim to cons and back again")
        {
            auto S1 = srmhd.fromPrimitive (R, P);
            auto S2 = srmhd.fromConserved (R, S1.U.begin());

            for (int i = 0; i < 8; ++i)
            {
               CHECK (P[i] == Approx (S2.P[i]).epsilon (1e-12));
            }
        }
    }
}




// ============================================================================
#include "HDF5.hpp"
#include "TimeSeriesManager.hpp"

SCENARIO ("Time series manager should behave reasonably", "[TimeSeriesManager]")
{
    GIVEN ("A time series manager populated with data series, with unequal lengths")
    {
        TimeSeriesManager timeSeriesManager;
        timeSeriesManager.getLogger()->setLogToNull();
        timeSeriesManager.append ("mean_energy", 1.1);
        timeSeriesManager.append ("num_unhealthy_zones", 12);
        timeSeriesManager.append ("mean_pressure", 2.2);
        timeSeriesManager.append ("num_unhealthy_zones", 13);

        WHEN ("Columns are queried")
        {
            THEN ("They have the expected size")
            {
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 2);
                CHECK (timeSeriesManager.getSeriesNamesInt().size() == 1);
                CHECK (timeSeriesManager.getSeriesDouble ("mean_energy").size() == 1);
                CHECK (timeSeriesManager.getSeriesDouble ("mean_pressure").size() == 1);
                CHECK (timeSeriesManager.getSeriesInt ("num_unhealthy_zones").size() == 2);
                CHECK_THROWS_AS (timeSeriesManager.getSeriesInt ("not_there"), std::out_of_range);
            }
        }

        WHEN ("Time series data is written to HDF5")
        {
            auto hdf5File = Cow::H5::File ("time_series_test.h5", "w");
            timeSeriesManager.write (hdf5File);

            THEN ("There are the expected number of columns in the file")
            {
                int numDatasets = 0;
                hdf5File.iterate ([&] (std::string) { numDatasets++; });
                CHECK (numDatasets == 3);
            }

            THEN ("Time series data can be cleared and restored from the file")
            {
                timeSeriesManager.clear();
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 0);

                timeSeriesManager.load (hdf5File);
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 2);
                CHECK (timeSeriesManager.getSeriesInt ("num_unhealthy_zones").size() == 2);
            }
        }
    }
}




// ============================================================================
#include "SphericalMeshGeometry.hpp"

SCENARIO ("Spherical polar mesh geometry should do its geometry right", "[SphericalMeshGeometry]")
{
    GIVEN ("An instance of SphericalMeshGeometry with shape [100, 1, 1]")
    {
        SphericalMeshGeometry geometry (100, 1, 1);
        geometry.setLowerUpper ({{1.0, 0.0, 0.0}}, {{10.0, M_PI, 2 * M_PI}});

        WHEN ("The geometry is queried for basic values")
        {
            THEN ("The results are correct")
            {
                CHECK (geometry.cellsShape()[0] == 100);
                CHECK (geometry.cellsShape()[1] == 1);
                CHECK (geometry.cellsShape()[2] == 1);
                CHECK (geometry.coordinateAtIndex (-0.5, 0.0, 0.0)[0] == Approx ( 1.0).epsilon (1e-14));
                CHECK (geometry.coordinateAtIndex (99.5, 0.0, 0.0)[0] == Approx (10.0).epsilon (1e-14));
                CHECK (geometry.meshVolume() == 4. / 3 * M_PI * 999);
            }
        }
    }
}




// ============================================================================
#include "MeshOperator.hpp"
#include "BoundaryConditions.hpp"
#include "CartesianMeshGeometry.hpp"

SCENARIO ("Mesh operator should generate data, take div and curl", "[MeshOperator]")
{
    GIVEN ("A 1D mesh operator and cartesian geometry")
    {
        auto fn = [] (double x, double, double)
        {
            return std::vector<double> { std::sin (2 * M_PI * x), std::cos (2 * M_PI * x) };
        };
        auto mo = std::make_shared<MeshOperator>();
        auto mg = std::make_shared<CartesianMeshGeometry>();
        auto bc = std::make_shared<PeriodicBoundaryCondition>();

        mo->setMeshGeometry (mg);

        WHEN ("Data is generated in volumes")
        {
            auto P = mo->generate (fn, MeshLocation::cell);
    
            THEN ("The returned array has the expected size")
            {
                CHECK (P.size(0) == mg->cellsShape()[0]);
                CHECK (P.size(3) == 2);
            }

            THEN ("The returned array has the expected data")
            {
                CHECK (P (0, 0, 0, 0) == fn (mg->coordinateAtIndex(0, 0, 0)[0], 0., 0.)[0]);
                CHECK (P (0, 0, 0, 1) == fn (mg->coordinateAtIndex(0, 0, 0)[0], 0., 0.)[1]);
            }
        }

        WHEN ("Two guard zones are left untouched by the initial data function")
        {
            int numGuard = 2;
            auto fullShape = mg->cellsShape();
            fullShape[0] += 2 * numGuard;
            fullShape[3] = 2;

            auto interiorRegion = Region().withLower (0, numGuard).withUpper (0, -numGuard);
            auto P = Array (fullShape);

            P[interiorRegion] = mo->generate (fn, MeshLocation::cell);
            bc->apply (P, MeshLocation::cell, MeshBoundary::left,  0, numGuard);
            bc->apply (P, MeshLocation::cell, MeshBoundary::right, 0, numGuard);

            THEN ("Periodic boundary conditions assign the expected values to the guard zones")
            {
                CHECK (P (0, 0, 0, 0) == P (mg->cellsShape()[0], 0, 0, 0));
                CHECK (P (0, 0, 0, 1) == P (mg->cellsShape()[0], 0, 0, 1));
                CHECK (P (mg->cellsShape()[0] + numGuard, 0, 0, 0) == P (numGuard, 0, 0, 0));
                CHECK (P (mg->cellsShape()[0] + numGuard, 0, 0, 1) == P (numGuard, 0, 0, 1));
            }
        }

        WHEN ("Data is generated on faces")
        {
            auto B = mo->generate (fn, MeshLocation::face);
            auto M = mo->divergence (B);

            THEN ("The returned array has the expected shape")
            {
                CHECK (B.size(0) == mg->cellsShape()[0] + 1);
                CHECK (B.size(3) == 2);
                CHECK (B.size(4) == 3);
            }

            THEN ("The array of corresponding divergences has the correct shape")
            {
                CHECK (M.size(0) == mg->cellsShape()[0]);
            }
        }
    }

    GIVEN ("A 3D mesh operator and a [16, 32, 8] cartesian geometry")
    {
        auto abcField = [] (double x, double y, double z)
        {
            const double k = 2 * M_PI;
            const double A = 1.0;
            const double B = 1.0;
            const double C = 1.0;
            const double b1 = C * std::cos (k * z) - B * std::sin (k * y);
            const double b2 = A * std::cos (k * x) - C * std::sin (k * z);
            const double b3 = B * std::cos (k * y) - A * std::sin (k * x);
            return std::vector<double> { b1, b2, b3 };
        };
        auto divLikeField = [] (double x, double y, double z)
        {
            const double k = 2 * M_PI;
            const double b1 = std::cos (k * x);
            const double b2 = std::cos (k * y);
            const double b3 = std::cos (k * z);
            return std::vector<double> { b1, b2, b3 };
        };
        auto mo = std::make_shared<MeshOperator>();
        auto mg = std::make_shared<CartesianMeshGeometry>();

        mg->setCellsShape ({{16, 32, 8}});
        mo->setMeshGeometry (mg);

        WHEN ("An ABC field is generated on faces in flux mode")
        {
            auto B = mo->generate (abcField, MeshLocation::face, MeshOperator::VectorMode::fluxish);
            auto M = mo->divergence (B);

            THEN ("The B field has the expected shape [17, 33, 9] with one rank (1, 3)")
            {
                CHECK (B.size(0) == 17);
                CHECK (B.size(1) == 33);
                CHECK (B.size(2) == 9);
                CHECK (B.size(3) == 1);
                CHECK (B.size(4) == 3);
            }

            THEN ("The monopole field has the expected shape [16, 32, 8] with rank (1, 1)")
            {
                CHECK (M.size(0) == 16);
                CHECK (M.size(1) == 32);
                CHECK (M.size(2) == 8);
                CHECK (M.size(3) == 1);
                CHECK (M.size(4) == 1);
            }

            THEN ("The monopole field is zero at some randomly checked locations")
            {
                CHECK (M (2, 6, 6) == Approx (0.0));
                CHECK (M (6, 0, 6) == Approx (0.0));
                CHECK (M (6, 6, 2) == Approx (0.0));
                CHECK (M (0, 1, 0) == Approx (0.0));
            }
        }

        WHEN ("A non-solenoidal field is generated on faces in flux mode")
        {
            auto V = mo->generate (divLikeField, MeshLocation::face, MeshOperator::VectorMode::fluxish);
            auto D = mo->divergence (V);

            THEN ("The divergence field is non-zero at some randomly checked locations")
            {
                CHECK (D (2, 6, 6) != Approx (0.0));
                CHECK (D (6, 0, 6) != Approx (0.0));
                CHECK (D (6, 6, 2) != Approx (0.0));
                CHECK (D (0, 1, 0) != Approx (0.0));
            }
        }

        WHEN ("An ABC field is generated on edges in scalars mode")
        {
            auto A = mo->generate (abcField, MeshLocation::edge, MeshOperator::VectorMode::scalars);

            THEN ("The A field has the expected shape [17, 33, 9] with rank (3, 3)")
            {
                CHECK (A.size(0) == 17);
                CHECK (A.size(1) == 33);
                CHECK (A.size(2) == 9);
                CHECK (A.size(3) == 3);
                CHECK (A.size(4) == 3);
            }
        }

        WHEN ("An ABC field is generated on edges in EMF mode")
        {
            auto A = mo->generate (abcField, MeshLocation::edge, MeshOperator::VectorMode::emflike);
            auto B = mo->curl (A);
            auto M = mo->divergence (B);

            THEN ("The A field has the expected shape [17, 33, 9] with rank (1, 3)")
            {
                CHECK (A.size(0) == 17);
                CHECK (A.size(1) == 33);
                CHECK (A.size(2) == 9);
                CHECK (A.size(3) == 1);
                CHECK (A.size(4) == 3);
            }

            THEN ("The curl of A has the expected shape [17, 33, 9] with rank (1, 3)")
            {
                CHECK (B.size(0) == 17);
                CHECK (B.size(1) == 33);
                CHECK (B.size(2) == 9);
                CHECK (B.size(3) == 1);
                CHECK (B.size(4) == 3);
            }

            THEN ("The div of B has the expected shape [16, 32, 8] with rank (1, 1)")
            {
                CHECK (M.size(0) == 16);
                CHECK (M.size(1) == 32);
                CHECK (M.size(2) == 8);
                CHECK (M.size(3) == 1);
                CHECK (M.size(4) == 1);
            }

            THEN ("The div of B is zero at some randomly checked locations")
            {
                CHECK (M (2, 6, 6) == Approx (0.0));
                CHECK (M (6, 0, 6) == Approx (0.0));
                CHECK (M (6, 6, 2) == Approx (0.0));
                CHECK (M (0, 1, 0) == Approx (0.0));
            }
        }
    }
}




// ============================================================================
#include "MeshOperator.hpp"

SCENARIO ("Mesh operator should run flux sweeps", "[MeshOperator]")
{
    GIVEN ("A 1D mesh operator, cartesian geometry, and a trivial Godunov flux function")
    {
        auto mg = std::make_shared<CartesianMeshGeometry>(8, 1, 1);
        auto mo = std::make_shared<MeshOperator>(mg);

        auto cd = [] (double x, double, double) { return std::vector<double> { std::sin (M_PI * x) }; };
        auto fd = [] (double x, double, double) { return std::vector<double> { 0., 1., 2. }; };
        auto gd = [] (GodunovStencil& stencil)
        {
            if (stencil.faceNormal == UnitVector::xhat)
            {
                CHECK (stencil.cellData.size(0) == 2);
                CHECK (stencil.cellData.size(1) == 1);
                CHECK (stencil.faceData.size(0) == 3); // is a 1D array of data on a single face
                CHECK (stencil.faceData[0] == 0.);
                CHECK (stencil.faceData[1] == 1.);
                CHECK (stencil.faceData[2] == 2.);
            }
            else
            {
                CHECK (stencil.cellData.size(0) == 0);
                CHECK (stencil.cellData.size(1) == 1);
                CHECK (stencil.faceData.size(0) == 3);
            }
        };

        auto cellData = mo->generate (cd, MeshLocation::cell);
        auto faceData = mo->generate (fd, MeshLocation::face);

        THEN ("Initial data (face and cell) has the expected shape")
        {
            CHECK (cellData.size(0) == 8);
            CHECK (cellData.size(3) == 1);
            CHECK (faceData.size(0) == 9);
            CHECK (faceData.size(3) == 3);
        }

        auto footprint = Shape {{ 2, 0, 0 }};
        auto faceFlux = mo->godunov (gd, cellData, faceData, footprint);

        THEN ("Godunov fluxes have the expected shape")
        {
            CHECK (faceFlux.size(0) == 9);
            CHECK (faceFlux.size(1) == 2);
            CHECK (faceFlux.size(2) == 2);
            CHECK (faceFlux.size(3) == 1);
        }
    }

    GIVEN ("A 3D mesh operator, cartesian geometry, and a trivial Godunov flux function")
    {
        auto mg = std::make_shared<CartesianMeshGeometry>(8, 8, 8);
        auto mo = std::make_shared<MeshOperator>(mg);

        auto id = [] (double x, double y, double z) { return std::vector<double> { 0., 1., 2. }; };
        auto cd = [] (double x, double y, double z) { return std::vector<double> { }; };
        auto gd = [] (GodunovStencil& stencil)
        {
            CHECK (stencil.cellData.size(0) == 2);
            CHECK (stencil.cellData.size(1) == 3);
            CHECK (stencil.cellCoords.size() == 2);
            CHECK (stencil.godunovFlux.size(0) == 3);
            stencil.godunovFlux[0] = stencil.cellCoords[0][0]; // Fhat = x of cell to left of face
            stencil.godunovFlux[1] = stencil.cellCoords[0][1]; // Fhat = y of cell to left of face
            stencil.godunovFlux[2] = stencil.cellCoords[0][2]; // Fhat = z of cell to left of face
        };

        auto cellData = mo->generate (id, MeshLocation::cell);
        auto faceData = mo->generate (cd, MeshLocation::face);

        THEN ("Initial data (face and cell) has the expected shape")
        {
            CHECK (cellData.size(0) == 8);
            CHECK (cellData.size(1) == 8);
            CHECK (cellData.size(2) == 8);
            CHECK (cellData.size(3) == 3);
            CHECK (faceData.size(0) == 9);
            CHECK (faceData.size(1) == 9);
            CHECK (faceData.size(2) == 9);
            CHECK (faceData.size(3) == 0);
        }

        auto footprint = Shape {{ 2, 2, 2 }};
        auto start = Shape {{ -1, 0, +1 }}; // choose an arbitrary start location
        auto godunovFlux = mo->godunov (gd, cellData, faceData, footprint, start, nullptr);

        THEN ("Godunov fluxes have the expected shape")
        {
            CHECK (godunovFlux.size(0) == 9);
            CHECK (godunovFlux.size(1) == 9);
            CHECK (godunovFlux.size(2) == 9);
            CHECK (godunovFlux.size(3) == 3);
        }

        THEN ("Godunov fluxes have values consistent with start index")
        {
            CHECK (godunovFlux (1, 0, 0, 0, 0) == mg->coordinateAtIndex (start[0], 0, 0)[0]);
            CHECK (godunovFlux (0, 1, 0, 1, 1) == mg->coordinateAtIndex (0, start[1], 0)[1]);
            CHECK (godunovFlux (0, 0, 1, 2, 2) == mg->coordinateAtIndex (0, 0, start[2])[2]);
        }
    }
}




// ============================================================================
#include "FieldOperator.hpp"

SCENARIO ("Field operator should convert prim -> cons", "[FieldOperator]")
{
    GIVEN ("A field operator and a 1D mesh")
    {
        auto mg = std::make_shared<CartesianMeshGeometry>();
        auto mo = std::make_shared<MeshOperator>(mg);
        auto cl = std::make_shared<NewtonianHydro>();
        auto fo = std::make_shared<FieldOperator>(cl);
        auto id = [] (double x, double, double)
        {
            return std::vector<double> { 1. + 0.5 * std::cos (M_PI * x), 1., 2., 3., 1. };
        };

        WHEN ("The mesh operator is used to generate primitives and measures")
        {
            auto P = mo->generate (id, MeshLocation::cell);
            auto V = mo->measure (MeshLocation::cell);
            auto L = mo->linearCellDimension();

            THEN ("The cell linear dimension is correct")
            {
                CHECK (L.shape() == mg->cellsShape());
                CHECK (L (0, 0, 0) == mg->cellLength (0, 0, 0, 0));
            }

            THEN ("The cell volume is correct")
            {
                CHECK (V.shape() == mg->cellsShape());
                CHECK (V (0, 0, 0) == mg->cellVolume (0, 0, 0));
            }
        }

        WHEN ("The field operator is used to generate conserved variables")
        {
            auto P = mo->generate (id, MeshLocation::cell);
            auto X = mo->cellCentroidCoordinates();
            auto U = fo->generateConserved (P);
            auto Q = fo->recoverPrimitive (U);

            THEN ("The conserved variables are as expected")
            {
                for (int i = 0; i < U.size(0); ++i)
                {
                    auto R = ConservationLaw::Request();
                    auto Pi = id (X (i), 0, 0);
                    auto Si = cl->fromPrimitive (R, &Pi[0]);

                    CHECK (U (i, 0, 0, 0) == Si.U[0]);
                    CHECK (U (i, 0, 0, 1) == Si.U[1]);
                    CHECK (U (i, 0, 0, 2) == Si.U[2]);
                    CHECK (U (i, 0, 0, 3) == Si.U[3]);
                    CHECK (U (i, 0, 0, 4) == Si.U[4]);
                }
            }

            THEN ("The re-generated primatives are the same as the originals")
            {
                CHECK (P.size() == Q.size());
                
                for (int n = 0; n < P.size(); ++n)
                {
                    CHECK (P[n] == Approx (Q[n]));
                }
            }
        }

        WHEN ("The and mesh field operators are used to generate the Courant time and diagnostics")
        {
            auto id = [] (double, double, double) { return std::vector<double> { 1., 1., 0., 0., 0. }; };
            auto P = mo->generate (id, MeshLocation::cell);
            auto V = mo->measure (MeshLocation::cell);
            auto L = mo->linearCellDimension();
            auto D = fo->volumeIntegratedDiagnostics (P, V);
            auto N = cl->getDiagnosticNames();

            THEN ("The Courant time is as expected")
            {
                CHECK (fo->getCourantTimestep (P, L) == mg->cellLength (0, 0, 0, 0));
            }

            THEN ("The diagnostics are reasonable")
            {
                CHECK (N[0] == "mass");
                CHECK (D[0] == mg->meshVolume());
            }
        }
    }
}




// ============================================================================
#include "MeshData.hpp"

SCENARIO ("Mesh data class works OK", "[MeshData]")
{
    GIVEN ("A mesh data class with [8, 12, 16] interior and 2 guard zones on all boundaries")
    {
        auto cs = Shape {{ 8, 12, 16 }}; // cells shape
        auto bs = Shape {{ 2, 2, 2 }};   // boundary shape
        auto md = std::make_shared<MeshData>(cs, bs, 3);
        auto mg = std::make_shared<CartesianMeshGeometry>(cs);
        auto mo = std::make_shared<MeshOperator>(mg);

        WHEN ("Magnetic field index is not set")
        {
            THEN ("The mesh primitive data has the correct shape")
            {
                CHECK (md->P.size(0) == cs[0] + 2 * bs[0]);
                CHECK (md->P.size(1) == cs[1] + 2 * bs[1]);
                CHECK (md->P.size(2) == cs[2] + 2 * bs[2]);
                CHECK (md->P.size(3) == 3);
                CHECK (md->P.size(4) == 1);
            }

            THEN ("The mesh magnetic data has the correct shape")
            {
                CHECK (md->B.size(0) == cs[0] + 2 * bs[0] + 1);
                CHECK (md->B.size(1) == cs[1] + 2 * bs[1] + 1);
                CHECK (md->B.size(2) == cs[2] + 2 * bs[2] + 1);
                CHECK (md->B.size(3) == 1);
                CHECK (md->B.size(4) == 3);
            }

            THEN ("We get an exception if retrieving non-existent magnetic field")
            {
                CHECK_THROWS (md->getMagneticField (MeshLocation::cell));
                CHECK_THROWS (md->getMagneticField (MeshLocation::face));
            }
        }

        WHEN ("Magnetic field index is set")
        {
            md->setMagneticIndex(0);

            auto P = md->getPrimitive();
            auto D = md->getPrimitive(2);
            auto H = md->getMagneticField (MeshLocation::cell);
            auto B = md->getMagneticField (MeshLocation::face);

            THEN ("No exception is raised on query magnetic field")
            {
                CHECK_NOTHROW (md->getMagneticField (MeshLocation::cell));
                CHECK_NOTHROW (md->getMagneticField (MeshLocation::face));
            }

            THEN ("The data returned by get methods has the correct shape")
            {
                CHECK (P.shape() == Shape ({{ cs[0], cs[1], cs[2], 3, 1 }}));
                CHECK (D.shape() == Shape ({{ cs[0], cs[1], cs[2], 1, 1 }}));
                CHECK (H.shape() == Shape ({{ cs[0], cs[1], cs[2], 3, 1 }}));
                CHECK (B.shape() == Shape ({{ cs[0] + 1, cs[1] + 1, cs[2] + 1, 1, 3 }}));
            }

            THEN ("Retrieved data can be sent back to assign methods")
            {
                CHECK_NOTHROW (md->assignPrimitive(P));
                CHECK_NOTHROW (md->assignMagneticField (H, MeshLocation::cell));
                CHECK_NOTHROW (md->assignMagneticField (B, MeshLocation::face));
            }

            THEN ("An exception is raised if we mixed up the order of field and cell data")
            {
                CHECK_THROWS (md->assignMagneticField (B, MeshLocation::cell)); // should be H
                CHECK_THROWS (md->assignMagneticField (H, MeshLocation::face)); // should be B
            }
        }

        WHEN ("Data is prepared with the mesh operator")
        {
            md->setMagneticIndex(0);

            auto id = [] (double, double, double) { return std::vector<double> { 0., 1., 2. }; };
            auto P = mo->generate (id, MeshLocation::cell);
            auto A = mo->generate (id, MeshLocation::edge, MeshOperator::VectorMode::emflike);
            auto B = mo->curl (A);

            THEN ("It is assigned to the mesh data interior as expected")
            {
                CHECK_NOTHROW (md->assignPrimitive (P));
                CHECK_NOTHROW (md->assignMagneticField (B, MeshLocation::face));
            }
        }
    }
}




// ============================================================================
#include "CellCenteredFieldCT.hpp"

SCENARIO ("Cell-centered field CT should behave reasonalby", "[CellCenteredFieldCT]")
{
    GIVEN ("An instance of CellCenteredFieldCT and random 2D Godunov flux (EMF)")
    {
        auto cs = Shape {{ 16, 16, 1 }};
        auto mg = std::make_shared<CartesianMeshGeometry>(cs);
        auto mo = std::make_shared<MeshOperator>(mg);
        auto ct = std::make_shared<CellCenteredFieldCT>();
        auto F = Array (17, 17, 2, 3, 3);

        for (int i = 0; i < 17; ++i)
        for (int j = 0; j < 17; ++j)
        {
            F(i, j, 0, 0, 0) = F(i, j, 1, 0, 0) = rand();
            F(i, j, 0, 0, 1) = F(i, j, 1, 0, 1) = rand();
            F(i, j, 0, 0, 2) = F(i, j, 1, 0, 2) = rand();
            F(i, j, 0, 1, 0) = F(i, j, 1, 1, 0) = rand();
            F(i, j, 0, 1, 1) = F(i, j, 1, 1, 1) = rand();
            F(i, j, 0, 1, 2) = F(i, j, 1, 1, 2) = rand();
            F(i, j, 0, 2, 0) = F(i, j, 1, 2, 0) = rand();
            F(i, j, 0, 2, 1) = F(i, j, 1, 2, 1) = rand();
            F(i, j, 0, 2, 2) = F(i, j, 1, 2, 2) = rand();
        }

        auto G = ct->generateGodunovFluxes (F, 0);
        auto Bf = mo->divergence (F);
        auto Bg = mo->divergence (G);
        auto Mf = ct->monopole (Bf, MeshLocation::cell);
        auto Nf = ct->monopole (Bf, MeshLocation::vert);
        auto Mg = ct->monopole (Bg, MeshLocation::cell);
        auto Ng = ct->monopole (Bg, MeshLocation::vert);

        CHECK (Mf.size(0) == 16);
        CHECK (Mf.size(1) == 16);
        CHECK (Mf.size(2) == 1);
        CHECK (Mf.size(3) == 1);
        CHECK (Mf.size(4) == 1);

        CHECK (Nf.size(0) == 17);
        CHECK (Nf.size(1) == 17);
        CHECK (Nf.size(2) == 2);
        CHECK (Nf.size(3) == 1);
        CHECK (Nf.size(4) == 1);

        THEN ("Cell-centered div-B's are zero")
        {
            CHECK (Mf (8, 8, 0) != Approx (0.0));
            CHECK (Mg (8, 8, 0) == Approx (0.0));
        }
    }

    GIVEN ("An instance of CellCenteredFieldCT and random 3D Godunov flux (EMF)")
    {
        auto cs = Shape {{ 16, 16, 16 }};
        auto mg = std::make_shared<CartesianMeshGeometry>(cs);
        auto mo = std::make_shared<MeshOperator>(mg);
        auto ct = std::make_shared<CellCenteredFieldCT>();
        auto F = Array (17, 17, 17, 3, 3);

        for (int i = 0; i < 17; ++i)
        for (int j = 0; j < 17; ++j)
        for (int k = 0; k < 17; ++k)
        {
            F(i, j, k, 0, 0) = rand();
            F(i, j, k, 0, 1) = rand();
            F(i, j, k, 0, 2) = rand();
            F(i, j, k, 1, 0) = rand();
            F(i, j, k, 1, 1) = rand();
            F(i, j, k, 1, 2) = rand();
            F(i, j, k, 2, 0) = rand();
            F(i, j, k, 2, 1) = rand();
            F(i, j, k, 2, 2) = rand();
        }

        auto G = ct->generateGodunovFluxes (F, 0);
        auto interior = Region().withRange (0, 2, -2).withRange (1, 2, -2).withRange (2, 2, -2);
        auto Bf = mo->divergence (F);
        auto Bg = mo->divergence (G);
        auto Mf = Array (ct->monopole (Bf, MeshLocation::cell)[interior]);
        auto Mg = Array (ct->monopole (Bg, MeshLocation::cell)[interior]);

        THEN ("Cell-centered div-B's are zero excluding a boudary with shape (2, 2, 2)")
        {
            Mg.shape3D().deploy ([&] (int i, int j, int k)
            {
                INFO (i << " " << j << " " << k);
                CHECK (Mf (i, j, k) != Approx (0.0));
                CHECK (Mg (i, j, k) == Approx (0.0));
            });
        }
    }
}




// ============================================================================
#include "SolutionSchemes.hpp"

SCENARIO ("Method of lines TVD scheme behaves reasonably", "[SolutionScheme]")
{
    GIVEN ("An instance of MethodOflinesTVD")
    {
        auto cs = Shape {{ 8, 1, 1 }}; // cells shape
        auto bs = Shape {{ 1, 0, 0 }}; // boundary shape

        auto mg = std::make_shared<CartesianMeshGeometry>(cs);
        auto mo = std::make_shared<MeshOperator>(mg);
        auto cl = std::make_shared<NewtonianHydro>();
        auto fo = std::make_shared<FieldOperator>(cl);
        auto ss = std::make_shared<MethodOfLinesTVD>();
        auto bc = std::make_shared<PeriodicBoundaryCondition>();
        auto md = std::make_shared<MeshData>(cs, bs, 5);
        auto id = [] (double x, double, double) { return std::vector<double> {1., 0., 0., 0., 1.}; };

        WHEN ("No field operator has been set")
        {
            THEN ("Calling advance() raises an exception as it should")
            {
                CHECK_THROWS (ss->advance (*md, 0.0, 0.0));
            }
        }

        WHEN ("All necessary operators have been set")
        {
            ss->setFieldOperator (fo);
            ss->setMeshOperator (mo);
            ss->setBoundaryCondition (bc);
            md->assignPrimitive (mo->generate (id, MeshLocation::cell));
            md->applyBoundaryCondition (*bc);

            THEN ("Calling advance() returns successfully")
            {
                CHECK_NOTHROW (ss->advance (*md, 0.0, 0.0));
            }

            THEN ("The solution is unchanged after update (it was uniform)")
            {
                auto P0 = md->P;

                ss->advance (*md, 0.0, 1.0);

                CHECK (md->P.shape() == P0.shape());
                
                for (int n = 0; n < md->P.size(); ++n)
                {
                    CHECK (md->P[n] == P0[n]);
                }
            }
        }

        WHEN ("The solution data guard zone region is too small for the scheme")
        {
            ss->setFieldOperator (fo);
            ss->setMeshOperator (mo);
            ss->setBoundaryCondition (bc);

            auto bs2 = Shape {{ 0, 0, 0 }};
            auto md2 = std::make_shared<MeshData>(cs, bs2, 5);
            md2->assignPrimitive (mo->generate (id, MeshLocation::cell));
            md2->applyBoundaryCondition (*bc);

            THEN ("The scheme works with the old mesh data")
            {
                CHECK_NOTHROW (ss->advance (*md, 0.0, 1.0));
            }

            THEN ("But it raises an exception with the new mesh data")
            {
                CHECK_THROWS (ss->advance (*md2, 0.0, 1.0));
            }
        }
    }
}




// ============================================================================
#include "TaskScheduler.hpp"

class Task1 : public TaskScheduler::Task
{
public:
    Task1() : numberOfRuns (0) {}
    void run (SimulationStatus status, int rep) override
    {
        CHECK (status.simulationIter % 12 == 0);
        CHECK (numberOfRuns == rep);
        ++numberOfRuns;
    }
    int numberOfRuns;
};

class Task2 : public TaskScheduler::Task
{
public:
    Task2() : numberOfRuns (0) {}
    void run (SimulationStatus status, int rep) override
    {
        CHECK (status.simulationTime == Approx (numberOfRuns * 0.1));
        CHECK (numberOfRuns == rep);
        ++numberOfRuns;
    }
    int numberOfRuns;
};

SCENARIO ("Task scheduler works reasonably", "[TaskScheduler]")
{
    GIVEN ("A scheduler with two tasks")
    {
        auto ts = TaskScheduler();
        auto t1 = std::make_shared<Task1>();
        auto t2 = std::make_shared<Task2>();
        auto status = SimulationStatus();

        ts.schedule (t1, TaskScheduler::Recurrence (0.0, 0.0, 12));
        ts.schedule (t2, TaskScheduler::Recurrence (0.1, 0.0));

        for (int i = 0; i < 144; ++i)
        {
            ts.dispatch (status);
            status.simulationIter += 1;
            status.simulationTime += 0.01;
        }
        ts.dispatch (status);

        THEN ("The tasks were run the expected number of times")
        {
            CHECK (t1->numberOfRuns == 13);
            CHECK (t2->numberOfRuns == 15);
        }
    }
}
