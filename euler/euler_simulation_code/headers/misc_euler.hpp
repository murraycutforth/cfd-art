/*
 *	DESCRIPTION:	Miscellaneous useful functions for the Allaire system
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef MISC_EULER_H
#define MISC_EULER_H

#include "typedefs.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <cassert>



inline vectype flux_conserved_var_euler_equations (const SGparams& eosparams, const vectype& U)
{
	/*
	 * The flux of conserved variables in the x-direction for Euler equations
	 */
	 
	vectype flux (4);
	double rho = U(0);
	double u = U(1) / rho;
	double v = U(2) / rho;
	double e = U(3) / rho - 0.5 * (u * u + v * v);
	double p = eos::pressure(eosparams.gamma, eosparams.pinf, e, rho);
	
	flux(0) = U(0) * u;
	flux(1) = U(1) * u + p;
	flux(2) = U(2) * u;
	flux(3) = (U(3) + p) * u;
	
	return flux;
}




inline vectype flux_primitive_var_euler_equations (const SGparams& eosparams, const vectype& W)
{
	vectype flux (4);
	double rho = W(0);
	double p = W(3);
	double u = W(1);
	double v = W(2);
	double e = eos::specific_ie(eosparams.gamma, eosparams.pinf, p, rho);
	double E = rho * e + 0.5 * rho * (W(2) * W(2) + W(1) * W(1));
	
	flux(0) = rho * u;
	flux(2) = rho * u * u + p;
	flux(3) = rho * u * v;
	flux(4) = (E + p) * u;
	
	return flux;
}


inline vectype conserved_to_primitives_euler_equations (const SGparams& eosparams, const vectype& U)
{
	vectype prims (4);
	double rho = U(0);
	double u = U(1) / rho;
	double v = U(2) / rho;
	double e = U(3) / rho - 0.5 * (u * u + v * v);
	double p = eos::pressure(eosparams.gamma, eosparams.pinf, e, rho);
	
	prims(0) = rho;
	prims(1) = u;
	prims(2) = v;
	prims(3) = p;
	
	return prims;
}



inline vectype primitives_to_conserved_euler_equations (const SGparams& eosparams, const vectype& W)
{
	vectype conserved (4);
	double rho = W(0);
	double u = W(1);
	double v = W(2);
	double p = W(3);
	double e = eos::specific_ie(eosparams.gamma, eosparams.pinf, p, rho);
	
	conserved(0) = rho;
	conserved(1) = rho * u;
	conserved(2) = rho * v;
	conserved(3) = rho * e + 0.5 * rho * (u * u + v * v);
	
	return conserved;
}




inline void A_primitive_vars_euler_equations (const SGparams& eosparams, const vectype& W, Matrix4d& A)
{
	/*
	 * Use vector of primitive variables (W) to set the value of the
	 * matrix A which is the Jacobian of the system in primitive
	 * variable form.
	 */
	 
	double rho = W(0);
	double u = W(1);
	double p = W(3);
	double c = eos::soundspeed(eosparams.gamma, eosparams.pinf, p, rho);
	
	A = u * Eigen::Matrix<double, 4, 4>::Identity();
	A(0,1) = rho;
	A(1,3) = 1.0 / rho;
	A(3,1) = rho * c * c;
}



inline bool is_physical_state_euler_equations (const SGparams& eosparams, const vectype& U)
{
	double rho = U(0);
	double u = U(1) / rho;
	double v = U(2) / rho;
	double e = U(3) / rho - 0.5 * (u * u + v * v);
		
	return rho >= 0.0 && e >= 0.0;
}



#endif
