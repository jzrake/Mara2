#include <cmath>
#include <algorithm>
#include "QuarticPolynomial.hpp"




// ============================================================================
QuarticPolynomial::QuarticPolynomial (double d4, double d3, double d2, double d1, double d0) :
d4 (d4),
d3 (d3),
d2 (d2),
d1 (d1),
d0 (d0)
{

}

int QuarticPolynomial::solve (double roots[4]) const
{
    const double a3 = d3 / d4;
    const double a2 = d2 / d4;
    const double a1 = d1 / d4;
    const double a0 = d0 / d4;
    const double au2 = -a2;
    const double au1 = (a1 * a3 - 4 * a0) ;
    const double au0 = 4.0 * a0 * a2 - a1 * a1 - a0 * a3 * a3;

    int nr, nr12, nr34;
    double r1, r2, r3, r4;
    double x1, x2, x3, u1;
    solve_cubic_equation (1.0, au2, au1, au0, x1, x2, x3, nr);

    if (nr == 1)
    {
        u1 = x1;
    }
    else
    {
        u1 = (x1 > x3) ? x1 : x3;
    }

    double R2 = 0.25 * a3 * a3 + u1 - a2;
    double R = (R2 > 0.0) ? std::sqrt (R2) : 0.0;
    double D2, E2;

    if (R != 0.0)
    {
        const double foo1 = 0.75 * a3 * a3 - R2 - 2.0 * a2;
        const double foo2 = 0.25 * (4.0 * a3 * a2 - 8.0 * a1 - a3 * a3 * a3) / R;
        D2 = foo1 + foo2;
        E2 = foo1 - foo2;
    }
    else
    {
        const double foo1 = 0.75 * a3 * a3 - 2.0 * a2;
        const double foo2 = 2.00 * std::sqrt (u1 * u1 - 4.0 * a0);
        D2 = foo1 + foo2;
        E2 = foo1 - foo2;
    }

    if (D2 >= 0.0)
    {
        const double D = std::sqrt (D2);
        r1 = -0.25 * a3 + 0.5 * R - 0.5 * D;
        r2 = -0.25 * a3 + 0.5 * R + 0.5 * D;
        nr12 = 2;
    }
    else
    {
        r1 = r2 = -0.25 * a3 + 0.5 * R;
        nr12 = 0;
    }

    if (E2 >= 0.0)
    {
        const double E = std::sqrt (E2);
        r3 = -0.25 * a3 - 0.5 * R - 0.5 * E;
        r4 = -0.25 * a3 - 0.5 * R + 0.5 * E;
        nr34 = 2;
    }
    else
    {
        r3 = r4 = -0.25 * a3 - 0.5 * R;
        nr34 = 0;
    }

    int i = 0;

    if (nr12 != 0)
    {
        roots[i++] = r1;
        roots[i++] = r2;
    }
    if (nr34 != 0)
    {
        roots[i++] = r3;
        roots[i++] = r4;
    }

    std::sort (roots, roots + nr12 + nr34);
    return nr12 + nr34;
}

void QuarticPolynomial::solve_cubic_equation (
    double c3, double c2, double c1, double c0,
    double &x1, double &x2, double &x3, int &nr) const
{
    const double a2 = c2 / c3;
    const double a1 = c1 / c3;
    const double a0 = c0 / c3;
    const double q = a1 / 3.0 - a2 * a2 / 9.0;
    const double r = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 * a2 / 27.0;
    const double delta = q * q * q + r * r;

    if (delta > 0.0)
    {
        double s1 = r + std::sqrt (delta);
        double s2 = r - std::sqrt (delta);
        s1 = (s1 >= 0.0) ? std::pow (s1, 1./3) : -std::pow (-s1, 1./3);
        s2 = (s2 >= 0.0) ? std::pow (s2, 1./3) : -std::pow (-s2, 1./3);
        x1 = (s1 + s2) - a2 / 3.0;
        x2 = x3 = -0.5 * (s1 + s2) - a2 / 3.0;
        nr = 1;
    }
    else if (delta < 0.0)
    {
        const double theta = std::acos (r / std::sqrt (-q * q * q)) / 3.0;
        const double costh = std::cos (theta);
        const double sinth = std::sin (theta);
        const double sq = std::sqrt (-q);
        x1 = 2.0 * sq*costh - a2 / 3.0;
        x2 = -sq*costh - a2 / 3.0 - std::sqrt(3) * sq * sinth;
        x3 = -sq*costh - a2 / 3.0 + std::sqrt(3) * sq * sinth;
        nr = 3;
    }
    else
    {
        const double s = (r >= 0.0) ? std::pow (r, 1./3) : -std::pow (-r, 1./3);
        x1 = 2.0 * s - a2 / 3.0;
        x2 = x3 = -s - a2 / 3.0;
        nr = 3;
    }
}
