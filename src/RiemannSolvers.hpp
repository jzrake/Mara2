#ifndef RiemannSolver_hpp
#define RiemannSolver_hpp

#include "Mara.hpp"
#include "ConservationLaws.hpp"



/**
An extremely simple Riemann solver appropriate for scalar conservation laws
only. This is here to illustrate a prototypical implementation of the base
class.
*/
class UpwindRiemannSolver : public RiemannSolver
{
public:
	virtual ~UpwindRiemannSolver() {}
    State solve (const State& L, const State& R, AreaElement dA) const override;
};




/**
The classic HLL (or HLLE) Riemann solver.
*/
class HlleRiemannSolver : public RiemannSolver
{
public:
	virtual ~HlleRiemannSolver() {}

    State solve (const State& L, const State& R, AreaElement dA) const override;
};




/**
The HLLC Riemann solver of Toro for the non-relativistic Euler equations.
*/
class HllcNewtonianHydroRiemannSolver : public RiemannSolver
{
public:
	virtual ~HllcNewtonianHydroRiemannSolver() {}

    State solve (const State& L, const State& R, AreaElement dA) const override;
};


#endif
