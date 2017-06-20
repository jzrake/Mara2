#ifndef Problems_hpp
#define Problems_hpp

#include "Mara.hpp"




class SimpleTestProgram
{
public:
	struct Problem;
	struct Scheme;
	int run (int argc, const char* argv[]);
	void run (const Problem& problem, const Scheme& scheme);
};

#endif
