#ifndef SolutionSchemes_hpp
#define SolutionSchemes_hpp

#include "Mara.hpp"




class MethodOfLinesTVD : public SolutionScheme
{
public:
    using Shape3D = Cow::Shape3D;
    using Index = Cow::Index;

    MethodOfLinesTVD();
    int getStencilSize() const override;
    void advance (MeshData& solution, double dt) const override;
    void setBoundaryCondition (std::shared_ptr<BoundaryCondition> bc) override { boundaryCondition = bc; }
    void setMeshOperator      (std::shared_ptr<MeshOperator> mo)      override { meshOperator = mo; }
    void setFieldOperator     (std::shared_ptr<FieldOperator> fo)     override { fieldOperator = fo; }

private:
    std::shared_ptr<BoundaryCondition>    boundaryCondition;
    std::shared_ptr<FieldOperator>        fieldOperator;
    std::shared_ptr<MeshOperator>         meshOperator;
    std::shared_ptr<IntercellFluxScheme>  fluxScheme;
    Shape3D footprint;
    Index startIndex;
};

#endif
