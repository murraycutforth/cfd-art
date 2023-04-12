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

class flux_solver_godunov : public flux_solver_base {
	
public:
	
	flux_solver_godunov (std::shared_ptr<riemann_solver_base> RS_ptr, sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		flux_solver_base(RS_ptr, params, gamma1, gamma2, pinf1, pinf2)
	{}
	
	void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double dt, double dx, double& u_star, double& z_star, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr)
	{
		assert(stencil.size() == 2);
			
		if (p_star_ptr) flux = RS_ptr->solve_RP(stencil[0], stencil[1], &u_star, p_star_ptr, v_star_ptr);
		else flux = RS_ptr->solve_RP(stencil[0], stencil[1], &u_star);
				
		z_star = ((1.0 + std::copysign(1.0, u_star)) / 2.0 ) * stencil[0](5)
			 + ((1.0 - std::copysign(1.0, u_star)) / 2.0 ) * stencil[1](5);
	}
	
	std::shared_ptr<flux_solver_base> clone ()
	{
		return std::make_shared<flux_solver_godunov>(RS_ptr->clone(), params, eosparams.gamma1, eosparams.gamma2, eosparams.pinf1, eosparams.pinf2);
	}
};


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
