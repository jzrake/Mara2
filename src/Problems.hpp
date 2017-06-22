#ifndef Problems_hpp
#define Problems_hpp

#include "Mara.hpp"




class Hydro1DTestProgram
{
public:
	struct Problem;
	struct Scheme;
	int run (int argc, const char* argv[]);
	void run (const Problem& problem, const Scheme& scheme);
};




class Hydro2DTestProgram
{
public:
    struct Problem;
    struct Scheme;
    int run (int argc, const char* argv[]);
    void run (const Problem& problem, const Scheme& scheme);
};



class NewtonianMHD2DTestProgram
{
public:
    struct Problem;
    struct Scheme;
    int run (int argc, const char* argv[]);
    void run (const Problem& problem, const Scheme& scheme);
};

#endif
