#ifndef DiagnosticsMHD_hpp
#define DiagnosticsMHD_hpp

#include "Mara.hpp"

class CartesianMeshGeometry;
class BlockDecomposition;




class DiagnosticsMHD
{
public:
	std::vector<double> getDiagnostics (
		ConservationLaw& law,
		CartesianMeshGeometry& geometry,
		BlockDecomposition& block) const;
	std::vector<std::string> getDiagnosticNames() const;
};

#endif
