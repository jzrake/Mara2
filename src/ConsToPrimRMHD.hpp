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

    /**
    Set the adiabatic index (defaults to 4/3)
    */
    void set_gamma (double adiabatic_gamma);

    /**
    Provide a new conserved state in memory. Before the solver is executed, the
    user must provide a guess for the initial primitive state, by calling either
    estimate_from_cons() or set_starting_prim(P).
    */
    void new_state (const double *U);

    /**
    This estimate becomes exact for no magnetic field in the NR limit.
    */
    void estimate_from_cons();

    /**
    Explicitly provide a guess to be used at the initial iteration.
    */
    void set_starting_prim (const double *P);

    /**
    Get the starting values being used for root finding.
    */
    void get_starting_prim (double *P);

    /**
    Set a new value to use as the pressure floor. If the value is less than or
    equal to zero, then no pressure floor is used.
    */
    int set_pressure_floor (double pf);

    /**
    Return 1 if the previous call to reconstruct_prim applied a pressure
    floor.
    */
    int did_use_pressure_floor();

    /**
    Using Z=rho*h*W^2, and W, get the primitive variables.
    */
    int reconstruct_prim (double Z, double W, double *P);

    /**
    Solution based on Anton & Zanotti (2006), equations 84 and 85.
    */
    int solve_anton2dzw (double *P);

    /**
    Solution based on Noble et. al. (2006), using Z = rho h W^2 as the single
    unkown. Unfortunately, Noble uses 'W' for what I call Z. This function
    should really be called '1dz', but I use this name to reflect the name of
    the section in which it appears.
    */
    int solve_noble1dw (double *P);

    /**
    Get number of iterations used on the last execution
    */
    int get_iterations();

    /**
    Verify the health of a conserved state. Return one of the enum values
    named above.
    */
    static int check_cons (const double *U);

    /**
    Verify the health of a primitive state. Return one of the enum values
    named above.
    */
    static int check_prim (const double *P);

    /**
    Return an error message for one of the enum values listed above.
    */
    static const char *get_error (int error);        

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
