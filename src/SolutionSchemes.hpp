#ifndef SolutionSchemes_hpp
#define SolutionSchemes_hpp

#include "Mara.hpp"




class MethodOfLinesTVD : public SolutionScheme
{
public:
    MethodOfLinesTVD();
    void advance (double dt, MeshData& solution) const override;
    void applyBoundaryCondition (MeshData& solution) const override;
    void setBoundaryCondition (std::shared_ptr<BoundaryCondition> bc) override { boundaryCondition = bc; }
    void setMeshOperator      (std::shared_ptr<MeshOperator> mo)      override { meshOperator = mo; }
    void setFieldOperator     (std::shared_ptr<FieldOperator> fo)     override { fieldOperator = fo; }

private:
    std::shared_ptr<BoundaryCondition>    boundaryCondition;
    std::shared_ptr<FieldOperator>        fieldOperator;
    std::shared_ptr<MeshOperator>         meshOperator;
};

#endif
