/*
 * 	DESCRIPTION:	Functions related to the stiffened gas equation
 * 			of state, defined by the equation:
 * 				
 * 			p = (gamma - 1) * rho * e - gamma * pinf
 * 
 * 			where p, rho, e are pressure, density, and
 * 			specific internal energy respectively, while
 * 			gamma and pinf are material-dependent constants.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		21/07/2017
 */

#ifndef STIFFENEDGASEOS_H
#define STIFFENEDGASEOS_H

#include <cmath>

struct SGparams {
	/* Store EOS parameters for single stiffened gas material
	*/
	double gamma;
	double pinf;

	SGparams() : gamma (0.0), pinf (0.0) {}
	SGparams(double gamma, double pinf) : gamma (gamma), pinf (pinf) {}
	SGparams(const SGparams& other) : gamma (other.gamma), pinf (other.pinf) {}
};


struct binarySGparams {
	
	double gamma1;
	double gamma2;
	double pinf1;
	double pinf2;
	
	binarySGparams ()
	:
		gamma1	(0.0),
		gamma2	(0.0),
		pinf1	(0.0),
		pinf2	(0.0)
	{}
	
	binarySGparams (double gamma1, double gamma2, double pinf1, double pinf2)
	:
		gamma1	(gamma1),
		gamma2	(gamma2),
		pinf1	(pinf1),
		pinf2	(pinf2)
	{}
	
	binarySGparams (const binarySGparams& other)
	:
		gamma1	(other.gamma1),
		gamma2	(other.gamma2),
		pinf1	(other.pinf1),
		pinf2	(other.pinf2)
	{}
		
};

namespace eos
{

	inline double pressure (double gamma, double pinf, double e, double rho)
	{
		return (gamma - 1.0) * rho * e - gamma * pinf;
	}
	
	inline double specific_ie (double gamma, double pinf, double p, double rho)
	{
		return (p + gamma * pinf) / ((gamma - 1.0) * rho);
	}
	
	inline double soundspeed (double gamma, double pinf, double p, double rho)
	{
		return sqrt((gamma * (p + pinf)) / rho);
	}
	
	inline double bulkmodulus (double gamma, double pinf, double p, double rho)
	{
		double c = soundspeed(gamma, pinf, p, rho);
		return c * c * rho;
	}
	
	inline double isentropic_extrapolation (double gamma, double pinf, double rho_old, double p_old, double p_new)
	{
		/* 
		 * Find the density rho_new such that the state (rho_new, p_new)
		 * has the same entropy as the state (rho_old, p_old) using the
		 * isentropic law for stiffened gases.
		 */
		
		return rho_old * pow((p_new + pinf) / (p_old + pinf), 1.0 / gamma);
	}
	
	inline double shockspeed_increase (double gamma, double pinf, double rho, double p, double p_star)
	{
		double A = 2.0 / ((gamma + 1.0) * rho);
		double B = ((gamma - 1.0) / (gamma + 1.0)) * (p + pinf);
		return sqrt((B + pinf + p_star) / A) / rho;
	}

}

#endif
