#ifndef SolutionSchemes_hpp
#define SolutionSchemes_hpp

#include "Mara.hpp"




class SchemeHelpers
{
public:
    static void makeFootprint (
        int stencilSize,
        Cow::Shape3D arrayShape,
        Cow::Shape3D boundaryShape,
        Cow::Shape3D& footprint,
        Cow::Index& startIndex,
        Cow::Region& interior);
};




class GenericSolutionScheme : public SolutionScheme
{
public:
    using Shape3D = Cow::Shape3D;
    using Index = Cow::Index;
    using Region = Cow::Region;

    int getStencilSize() const override;
    void setBoundaryCondition   (std::shared_ptr<BoundaryCondition> bc)   override { boundaryCondition = bc; }
    void setMeshOperator        (std::shared_ptr<MeshOperator> mo)        override { meshOperator = mo; }
    void setFieldOperator       (std::shared_ptr<FieldOperator> fo)       override { fieldOperator = fo; }
    void setIntercellFluxScheme (std::shared_ptr<IntercellFluxScheme> fs) override;
    void setSourceTermsFunction (SourceTermsFunction sourceTermsFunctionToUse) { sourceTermsFunction = sourceTermsFunctionToUse; }

protected:
    std::shared_ptr<BoundaryCondition>    boundaryCondition;
    std::shared_ptr<FieldOperator>        fieldOperator;
    std::shared_ptr<MeshOperator>         meshOperator;
    std::shared_ptr<IntercellFluxScheme>  fluxScheme;
    SourceTermsFunction                   sourceTermsFunction;
};




class MethodOfLinesTVD : public GenericSolutionScheme
{
public:
    MethodOfLinesTVD();
    virtual ~MethodOfLinesTVD() {}
    void setRungeKuttaOrder (int rungeKuttaOrderToUse);
    void setDisableFieldCT (bool shouldDisableFieldCT);
    void setViscousFluxFunction (std::function<void(const Cow::Array&, Cow::Array&, double)> viscousFluxToUse);
    void setStarParticleDerivatives (std::function<std::vector<double>(const Cow::Array&, const std::vector<double>&)> starParticleDerivativesToUse);
    void setSourceTermsWithParticles (SourceTermsWithParticles sourceTermsWithParticlesToUse);
    void advance (MeshData& solution, double t0, double dt) const override;

    Cow::Array computeAdvectiveFluxes (MeshData& solution, double t0) const;
    Cow::Array computeViscousFluxes (MeshData& solution, double t0) const;

private:
    void check_valid() const;
    int rungeKuttaOrder;
    bool disableFieldCT;
    std::function<void(const Cow::Array&, Cow::Array&, double)> viscousFlux = nullptr;
    std::function<std::vector<double>(const Cow::Array&, const std::vector<double>&)> starParticleDerivatives = nullptr;
    SourceTermsWithParticles sourceTermsWithParticles = nullptr;
};



#endif
