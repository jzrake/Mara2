#pragma once




class QuarticPolynomial
{
public:
    QuarticPolynomial (double d4, double d3, double d2, double d1, double d0);
    int solve (double roots[4]) const;
private:
    void solve_cubic_equation (
        double c3, double c2, double c1, double c0,
        double &x1, double &x2, double &x3, int &nr) const;
    const double d4, d3, d2, d1, d0;
};
