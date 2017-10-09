#pragma once




class ConsToPrimRMHD
{
public:
    enum
    {
        SRMHD_C2P_SUCCESS,
        SRMHD_C2P_CONS_CONTAINS_NAN,
        SRMHD_C2P_CONS_NEGATIVE_DENSITY,
        SRMHD_C2P_CONS_NEGATIVE_ENERGY,
        SRMHD_C2P_PRIM_CONTAINS_NAN,
        SRMHD_C2P_PRIM_NEGATIVE_PRESSURE,
        SRMHD_C2P_PRIM_NEGATIVE_RESTMASS,
        SRMHD_C2P_PRIM_SUPERLUMINAL,
        SRMHD_C2P_MAXITER
    } ;

    ConsToPrimRMHD();
    void set_gamma (double adiabatic_gamma);
    void new_state (const double *U);
    void estimate_from_cons();
    void set_starting_prim (const double *P);
    void get_starting_prim (double *P);
    int set_pressure_floor (double pf);
    int did_use_pressure_floor();
    int reconstruct_prim (double Z, double W, double *P);
    int solve_anton2dzw (double *P);
    int solve_noble1dw (double *P);
    int get_iterations();
    int check_cons (const double *U);
    int check_prim (const double *P);
    const char *get_error (int error);        

private:
    int MaxIterations;
    double Tolerance;
    double bigZ;
    double bigW;
    double smlZ;
    double smlW;

    int AppliedPressureFloor;
    int Iterations;
    double AdiabaticGamma;
    double gamf;
    double D,Tau;
    double S2,B2,BS,BS2;
    double Cons[8];
    double Z_start, W_start;
    double PressureFloor;
};
