#ifndef IntercellFluxSchemes_hpp
#define IntercellFluxSchemes_hpp

#include "Mara.hpp"




class ScalarUpwind : public IntercellFluxScheme
{
public:
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override;
    int getStencilSize() const override;
};




class MethodOfLines : public IntercellFluxScheme
{
public:
    MethodOfLines();
    void setRiemannSolver (std::shared_ptr<RiemannSolver> solverToUse);
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override;
    int getStencilSize() const override;
private:
    std::shared_ptr<RiemannSolver> riemannSolver;
};




class MethodOfLinesPlm : public IntercellFluxScheme
{
public:
    MethodOfLinesPlm (double plmTheta=1.5);
    void setRiemannSolver (std::shared_ptr<RiemannSolver> solverToUse);
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override;
    int getStencilSize() const override;
private:
    std::shared_ptr<RiemannSolver> riemannSolver;
    Reconstruction plm;
};




class MethodOfLinesWeno : public IntercellFluxScheme
{
public:
    MethodOfLinesWeno (double shenZhaA=50);
    ConservationLaw::State intercellFlux (const FaceData& faceData) const override;
    int getStencilSize() const override;
private:
    Reconstruction weno;
};

#endif
