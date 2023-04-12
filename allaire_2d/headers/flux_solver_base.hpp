/*
 *	DESCRIPTION:	An abstract base class for solvers used to find
 * 			the flux used in the conservative update
 *			in the hyperbolic five equation system of Allaire.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef FS_ALLAIRE_H
#define FS_ALLAIRE_H

#include "typedefs.hpp"
#include "riemann_solver_base.hpp"
#include "stiffened_gas_eos.hpp"
#include <vector>
#include <memory>

class flux_solver_base {
	
protected:

	std::shared_ptr<riemann_solver_base> RS_ptr;	// Riemann solver
	
	const sim_info params;	// Constants used in computation
	binarySGparams eosparams;
	
	
public:
	
	flux_solver_base (std::shared_ptr<riemann_solver_base> RS_ptr, sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		RS_ptr (RS_ptr),
		params (params),
		eosparams (gamma1, gamma2, pinf1, pinf2)
	{}
	
	virtual void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double dt, double dx, double& u_star, double& z_star, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr) =0;
	
	virtual std::shared_ptr<flux_solver_base> clone () =0;
};


class FluxSolverBaseEulerEquations {
	/* Abstract base class for flux solvers for the 2D Euler equations
	*/
	
protected:

	std::shared_ptr<riemann_solver_base_euler_equations> RS_ptr;	// Riemann solver
	const sim_info params;	// Constants used in computation
	SGparams eosparams;
	
public:
	
	FluxSolverBaseEulerEquations (std::shared_ptr<riemann_solver_base_euler_equations> RS_ptr, sim_info params, double gamma, double pinf)
	:
		RS_ptr (RS_ptr),
		params (params),
		eosparams (gamma, pinf)
	{}
	
	virtual void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double dt, double dx, double& u_star, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr) =0;
	
	virtual std::shared_ptr<FluxSolverBaseEulerEquations> clone () =0;
};

#endif
