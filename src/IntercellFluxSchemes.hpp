#ifndef IntercellFluxSchemes_hpp
#define IntercellFluxSchemes_hpp

#include "Mara.hpp"
#include "Reconstruction.hpp"




class ScalarUpwind : public IntercellFluxScheme
{
public:
    virtual ~ScalarUpwind() {}
    ConservationLaw::State intercellFlux (const FaceData& faceData, double t) const override;
    int getStencilSize() const override;
};




class MethodOfLines : public IntercellFluxScheme
{
public:
    MethodOfLines();
    virtual ~MethodOfLines() {}
    void setRiemannSolver (std::shared_ptr<RiemannSolver> solverToUse) override;
    ConservationLaw::State intercellFlux (const FaceData& faceData, double t) const override;
    int getStencilSize() const override;
private:
    std::shared_ptr<RiemannSolver> riemannSolver;
};




class MethodOfLinesPlm : public IntercellFluxScheme
{
public:
    MethodOfLinesPlm();
    virtual ~MethodOfLinesPlm() {}
    void setPlmTheta (double plmTheta) override;
    void setRiemannSolver (std::shared_ptr<RiemannSolver> solverToUse) override;
    ConservationLaw::State intercellFlux (const FaceData& faceData, double t) const override;
    int getStencilSize() const override;
private:
    std::shared_ptr<RiemannSolver> riemannSolver;
    Reconstruction plm;
};




class MethodOfLinesWeno : public IntercellFluxScheme
{
public:
    MethodOfLinesWeno();
    virtual ~MethodOfLinesWeno() {}
    ConservationLaw::State intercellFlux (const FaceData& faceData, double t) const override;
    int getStencilSize() const override;
private:
    Reconstruction weno;
};

#endif
