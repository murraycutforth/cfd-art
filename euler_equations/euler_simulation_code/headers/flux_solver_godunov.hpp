/*
 *	DESCRIPTION:	An implementation of Godunov's first order
 * 			flux solver applied to different governing equations.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		25/07/2017
 */


#ifndef GODUNOV_ALLAIRE_H
#define GODUNOV_ALLAIRE_H

#include "flux_solver_base.hpp"
#include <cassert>
#include <cmath>


class GodunovEulerEquations : public FluxSolverBaseEulerEquations {
	/* Godunov's first order flux solver for Euler equations
	*/
	
public:
	
	GodunovEulerEquations (std::shared_ptr<riemann_solver_base_euler_equations> RS_ptr, sim_info params, double gamma, double pinf)
	:
		FluxSolverBaseEulerEquations(RS_ptr, params, gamma, pinf)
	{}
	
	void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double dt, double dx, double& u_star, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr)
	{
		assert(stencil.size() == 2);
			
		if (p_star_ptr) flux = RS_ptr->solve_RP(stencil[0], stencil[1], &u_star, p_star_ptr, v_star_ptr);
		else flux = RS_ptr->solve_RP(stencil[0], stencil[1], &u_star);
	}
	
	std::shared_ptr<FluxSolverBaseEulerEquations> clone ()
	{
		return std::make_shared<GodunovEulerEquations>(RS_ptr->clone(), params, eosparams.gamma, eosparams.pinf);
	}
};

#endif
