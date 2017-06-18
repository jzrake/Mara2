#ifndef Problems_hpp
#define Problems_hpp

#include "Mara.hpp"




class SimpleTestProgram
{
public:
	struct Problem;
	int run (int argc, const char* argv[]);

private:
	void setup (const Problem&);
	std::shared_ptr<MeshGeometry>        mg;
	std::shared_ptr<MeshOperator>        mo;
	std::shared_ptr<ConservationLaw>     cl;
	std::shared_ptr<FieldOperator>       fo;
	std::shared_ptr<SolutionScheme>      ss;
	std::shared_ptr<BoundaryCondition>   bc;
	std::shared_ptr<MeshData>            md;
};

#endif
