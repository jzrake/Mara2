#ifndef SolutionSchemes_hpp
#define SolutionSchemes_hpp

#include "Mara.hpp"




class MeshHelpers
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
protected:
    std::shared_ptr<BoundaryCondition>    boundaryCondition;
    std::shared_ptr<FieldOperator>        fieldOperator;
    std::shared_ptr<MeshOperator>         meshOperator;
    std::shared_ptr<IntercellFluxScheme>  fluxScheme;
};





class MethodOfLinesTVD : public GenericSolutionScheme
{
public:
    MethodOfLinesTVD();
    void setRungeKuttaOrder (int rungeKuttaOrderToUse);
    void setDisableFieldCT (bool shouldDisableFieldCT);
    void advance (MeshData& solution, double dt) const override;
private:
    int rungeKuttaOrder;
    bool disableFieldCT;
};




/**
A place-holder class for when we're ready to begin implementing radiation
hydro.
*/
class BasicRadHydroScheme : public GenericSolutionScheme
{
public:
    BasicRadHydroScheme();
    void advance (MeshData& solution, double dt) const override;
private:
};

#endif
