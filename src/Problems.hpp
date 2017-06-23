#ifndef Problems_hpp
#define Problems_hpp

#include "Mara.hpp"




class Hydro1DTestProgram : public SubProgram
{
public:
	struct Problem;
	struct Scheme;
	int run (int argc, const char* argv[]) override;
private:
	void runProblem (const Problem& problem, const Scheme& scheme);
};




class Hydro2DTestProgram : public SubProgram
{
public:
    struct Problem;
    struct Scheme;
    int run (int argc, const char* argv[]) override;
private:
    void runProblem (const Problem& problem, const Scheme& scheme);
};



class NewtonianMHD2DTestProgram : public SubProgram
{
public:
    struct Problem;
    struct Scheme;
    int run (int argc, const char* argv[]) override;
private:
    void runProblem (const Problem& problem, const Scheme& scheme);
};



class MagneticBraidingProgram : public SubProgram
{
public:
    int run (int argc, const char* argv[]) override;
};


#endif
