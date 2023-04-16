/*
 *	DESCRIPTION:	The primitive variable linearised Riemann solver
 * 			for ALlaire's five equation system.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		02/08/2017
 */

#ifndef PVRS_H
#define PVRS_H

#include "riemann_solver_base.hpp"
#include "typedefs.hpp"
#include "misc_euler.hpp"
#include "sim_info.hpp"
#include <vector>



class PVRS_riemann_solver_euler_equations : public riemann_solver_base_euler_equations {
	
private:

	vectype flux;
	vectype WR;
	vectype WL;
	
public:
	
	PVRS_riemann_solver_euler_equations (sim_info params, double gamma, double pinf)
	:
		riemann_solver_base_euler_equations (params, gamma, pinf),
		flux (4),
		WR (4),
		WL (4)
	{}
	
	PVRS_riemann_solver_euler_equations (const PVRS_riemann_solver_euler_equations& other)
	:
		riemann_solver_base_euler_equations (other.params, other.eosparams.gamma, other.eosparams.pinf),
		flux (4),
		WR (4),
		WL (4)
	{}
	
	vectype solve_RP (const vectype& UL, const vectype& UR, double* u_star_ptr = nullptr, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr)
	{
		assert(UL.rows() == UR.rows());
		assert(UL.rows() == 4);
		
		double rhoL = UL(0);
		double uL = UL(1) / rhoL;
		double vL = UL(2) / rhoL;
		double eL = UL(3) / rhoL - 0.5 * (uL * uL + vL * vL);
		double pL = eos::pressure(eosparams.gamma, eosparams.pinf, eL, rhoL);
		double cL = eos::soundspeed(eosparams.gamma, eosparams.pinf, pL, rhoL);
		
		double rhoR = UR(0);
		double uR = UR(1) / rhoR;
		double vR = UR(2) / rhoR;
		double eR = UR(3) / rhoR - 0.5 * (uR * uR + vR * vR);
		double pR = eos::pressure(eosparams.gamma, eosparams.pinf, eR, rhoR);
		double cR = eos::soundspeed(eosparams.gamma, eosparams.pinf, pR, rhoR);
		
		
		// Averaged quantities in linear problem
		
		double rhoavg = 0.5 * (rhoL + rhoR);
		double cavg = 0.5 * (cL + cR);
		double uavg = 0.5 * (uL + uR);
		
		
		// Star states
		
		double pstar = 0.5 * (pL + pR) + 0.5 * (uL - uR) * rhoavg * cavg;
		double ustar = uavg + 0.5 * (pL - pR) / (rhoavg * cavg);
		double rhostarL = rhoL + rhoavg * (uL - ustar) / cavg;
		double rhostarR = rhoR + rhoavg * (ustar - uR) / cavg;

		if (u_star_ptr) *u_star_ptr = ustar;
		if (p_star_ptr) *p_star_ptr = pstar;
		if (v_star_ptr)
		{
			if (ustar >= 0.0) *v_star_ptr = vL;
			else *v_star_ptr = vR;
		}
		
		
		if (uavg - cavg > 0.0)
		{
			// Flux is left state
			
			flux = flux_conserved_var_euler_equations(eosparams, UL);
		}
		else if (uavg > 0.0)
		{
			// Flux is left star state
			
			WR(0) = rhostarL;
			WR(1) = ustar;
			WR(2) = vL;
			WR(3) = pstar;
			
			flux = flux_primitive_var_euler_equations(eosparams, WR);
		}
		else if (uavg + cavg > 0.0)
		{
			// Flux is right star state
			
			WR(0) = rhostarR;
			WR(1) = ustar;
			WR(2) = vR;
			WR(3) = pstar;
			
			flux = flux_primitive_var_euler_equations(eosparams, WR);
		}
		else
		{
			// Flux is right state
			
			flux = flux_conserved_var_euler_equations(eosparams, UR);
		}
		
		return flux;
	}
	
	
	std::shared_ptr<riemann_solver_base_euler_equations> clone ()
	{	
		return std::make_shared<PVRS_riemann_solver_euler_equations>(*this);
	}
	
};

#endif
