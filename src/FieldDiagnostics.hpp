#include "Mara.hpp"

namespace Cow { namespace H5 { class Location; } }




class FieldDiagnostics
{
public:
	FieldDiagnostics();
	void setMeshDecomposition (std::shared_ptr<MeshDecomposition>);
	void setBoundaryCondition (std::shared_ptr<BoundaryCondition>);
	void setConservationLaw (std::shared_ptr<ConservationLaw>);
	void setMeshGeometry (std::shared_ptr<MeshGeometry>);

	void computeAndCache (Cow::Array::Reference primitiveData);
	void write (Cow::H5::Location& location) const;
	void load (Cow::H5::Location& location);
	void clear();
private:
	std::shared_ptr<BoundaryCondition> boundaryCondition;
	std::shared_ptr<ConservationLaw> conservationLaw;
	std::shared_ptr<MeshGeometry> meshGeometry;
};
