/*
 *	DESCRIPTION:	An abstract base class for Riemann solvers used
 *			in the hyperbolic five equation system of Allaire.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef RS_ALLAIRE_2D
#define RS_ALLAIRE_2D

#include "sim_info.hpp"
#include "typedefs.hpp"
#include "stiffened_gas_eos.hpp"
#include <memory>

class riemann_solver_base {
	
protected:
	
	const sim_info params;	// Constants used in computation
	binarySGparams eosparams;


public:
	
	riemann_solver_base (sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		params (params),
		eosparams (gamma1, gamma2, pinf1, pinf2)
	{}
	
	virtual vectype solve_RP (const vectype& UL, const vectype& UR, double* u_star_ptr = nullptr, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr) =0;
	
	virtual std::shared_ptr<riemann_solver_base> clone () =0;
	
};	


class riemann_solver_base_euler_equations {
	
protected:
	
	const sim_info params;	// Constants used in computation
	SGparams eosparams;


public:
	
	riemann_solver_base_euler_equations (sim_info params, double gamma, double pinf)
	:
		params (params),
		eosparams (gamma, pinf)
	{}
	
	virtual vectype solve_RP (const vectype& UL, const vectype& UR, double* u_star_ptr = nullptr, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr) =0;
	
	virtual std::shared_ptr<riemann_solver_base_euler_equations> clone () =0;
	
};	

#endif
