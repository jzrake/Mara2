#ifndef Reconstruction_hpp
#define Reconstruction_hpp





class Reconstruction
{
public:
    enum Operation
    {
        PLM_C2L, PLM_C2R,
        WENO5_FD_C2R, WENO5_FD_C2L, // Finite-difference, center-to-right/left
        WENO5_FV_C2R, WENO5_FV_C2L, // Finite-volume, center-to-right/left
        WENO5_FV_C2A, WENO5_FV_A2C, // Finite-volume, average-to-center/center-to-average
    };

    enum SmoothnessIndicator
    {
        OriginalJiangShu96,
        ImprovedBorges08,
        ImprovedShenZha10,
    };

    Reconstruction();
    double reconstruct (const double* v, enum Operation type) const;
    void setSmoothnessIndicator (enum SmoothnessIndicator IS);
    void setPlmTheta (double theta);
    void setShenZha10A (double A);

private:
    double plmTheta;    // [1 -> 2 (most aggressive)]
    double shenZha10A;  // [0 -> ~100 (most aggressive)]
    enum SmoothnessIndicator modeIS;

    double minmod (double ul, double u0, double ur) const;
    double plm (const double *v, double sgn) const;
    double weno5 (const double *v, const double c[3][3], const double d[3]) const;


    static const double CeesA2C_FV[3][3];
    static const double CeesC2A_FV[3][3];
    static const double CeesC2L_FV[3][3];
    static const double CeesC2R_FV[3][3];
    static const double CeesC2L_FD[3][3];
    static const double CeesC2R_FD[3][3];
    static const double DeesA2C_FV[3];
    static const double DeesC2A_FV[3];
    static const double DeesC2L_FV[3];
    static const double DeesC2R_FV[3];
    static const double DeesC2L_FD[3];
    static const double DeesC2R_FD[3];
};

#endif
