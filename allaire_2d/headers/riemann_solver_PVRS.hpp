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
#include "misc.hpp"
#include "sim_info.hpp"
#include "mixturemodel.hpp"
#include <vector>


class PVRS_riemann_solver : public riemann_solver_base {
	
private:

	vectype flux;
	vectype WR;
	vectype WL;
	
public:
	
	PVRS_riemann_solver (sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		riemann_solver_base (params, gamma1, gamma2, pinf1, pinf2),
		flux (6),
		WR (6),
		WL (6)
	{}
	
	PVRS_riemann_solver (const PVRS_riemann_solver& other)
	:
		riemann_solver_base (other.params, other.eosparams.gamma1, other.eosparams.gamma2, other.eosparams.pinf1, other.eosparams.pinf2),
		flux (6),
		WR (6),
		WL (6)
	{}
	
	vectype solve_RP (const vectype& UL, const vectype& UR, double* u_star_ptr = nullptr, double* p_star_ptr = nullptr, double* v_star_ptr = nullptr)
	{
		assert(UL.rows() == UR.rows());
		assert(UL.rows() == 6);
		
		double zL = UL(5);
		double rhoL = UL(0) + UL(1);
		double rhoz1L = UL(0);
		double rhoz2L = UL(1);
		double uL = UL(2) / rhoL;
		double vL = UL(3) / rhoL;
		double eL = UL(4) / rhoL - 0.5 * (uL * uL + vL * vL);
		double pL = allairemodel::mixture_pressure(eosparams, rhoL, eL, zL);
		double cL = allairemodel::mixture_soundspeed(eosparams, rhoL, pL, zL);
		
		double zR = UR(5);
		double rhoR = UR(0) + UR(1);
		double rhoz1R = UR(0);
		double rhoz2R = UR(1);
		double uR = UR(2) / rhoR;
		double vR = UR(3) / rhoR;
		double eR = UR(4) / rhoR - 0.5 * (uR * uR + vR * vR);
		double pR = allairemodel::mixture_pressure(eosparams, rhoR, eR, zR);
		double cR = allairemodel::mixture_soundspeed(eosparams, rhoR, pR, zR);
		
		
		// Averaged quantities in linear problem
		
		double rhoz1avg = 0.5 * (rhoz1L + rhoz1R);
		double rhoz2avg = 0.5 * (rhoz2L + rhoz2R);
		double rhoavg = 0.5 * (rhoL + rhoR);
		double cavg = 0.5 * (cL + cR);
		double uavg = 0.5 * (uL + uR);
		
		
		
		
		// Star states
		
		double pstar = 0.5 * (pL + pR) + 0.5 * (uL - uR) * rhoavg * cavg;
		double ustar = uavg + 0.5 * (pL - pR) / (rhoavg * cavg);
		double rhoz1starL = rhoz1L + rhoz1avg * (uL - ustar) / cavg;
		double rhoz1starR = rhoz1R + rhoz1avg * (ustar - uR) / cavg;
		double rhoz2starL = rhoz2L + rhoz2avg * (uL - ustar) / cavg;
		double rhoz2starR = rhoz2R + rhoz2avg * (ustar - uR) / cavg;
		
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
			
			flux = flux_conserved_var(eosparams, UL);
		}
		else if (uavg > 0.0)
		{
			// Flux is left star state
			
			WR(0) = rhoz1starL;
			WR(1) = rhoz2starL;
			WR(2) = ustar;
			WR(3) = vL;
			WR(4) = pstar;
			WR(5) = zL;
			
			flux = flux_primitive_var(eosparams, WR);
		}
		else if (uavg + cavg > 0.0)
		{
			// Flux is right star state
			
			WR(0) = rhoz1starR;
			WR(1) = rhoz2starR;
			WR(2) = ustar;
			WR(3) = vR;
			WR(4) = pstar;
			WR(5) = zR;
			
			flux = flux_primitive_var(eosparams, WR);
		}
		else
		{
			// Flux is right state
			
			flux = flux_conserved_var(eosparams, UR);
		}
		
		return flux;
	}
	
	void compute_primitive_jumps (const vectype& UL, const vectype& UR, std::vector<vectype>& jumps)
	{
		assert(UL.rows() == UR.rows());
		assert(UL.rows() == 6);
		assert(jumps.size() == 3);
		
		double zL = UL(5);
		double rhoL = UL(0) + UL(1);
		double rhoz1L = UL(0);
		double rhoz2L = UL(1);
		double uL = UL(2) / rhoL;
		double vL = UL(3) / rhoL;
		double eL = UL(4) / rhoL - 0.5 * (uL * uL + vL * vL);
		double pL = allairemodel::mixture_pressure(eosparams, rhoL, eL, zL);
		double cL = allairemodel::mixture_soundspeed(eosparams, rhoL, pL, zL);
		
		double zR = UR(5);
		double rhoR = UR(0) + UR(1);
		double rhoz1R = UR(0);
		double rhoz2R = UR(1);
		double uR = UR(2) / rhoR;
		double vR = UR(3) / rhoR;
		double eR = UR(4) / rhoR - 0.5 * (uR * uR + vR * vR);
		double pR = allairemodel::mixture_pressure(eosparams, rhoR, eR, zR);
		double cR = allairemodel::mixture_soundspeed(eosparams, rhoR, pR, zR);
		
		
		// Averaged quantities in linear problem
		
		double rhoz1avg = 0.5 * (rhoz1L + rhoz1R);
		double rhoz2avg = 0.5 * (rhoz2L + rhoz2R);
		double rhoavg = 0.5 * (rhoL + rhoR);
		double cavg = 0.5 * (cL + cR);
		double uavg = 0.5 * (uL + uR);
		
		
		// Star states
		
		double pstar = 0.5 * (pL + pR) + 0.5 * (uL - uR) * rhoavg * cavg;
		double ustar = uavg + 0.5 * (pL - pR) / (rhoavg * cavg);
		double rhoz1starL = rhoz1L + rhoz1avg * (uL - ustar) / cavg;
		double rhoz1starR = rhoz1R + rhoz1avg * (ustar - uR) / cavg;
		double rhoz2starL = rhoz2L + rhoz2avg * (uL - ustar) / cavg;
		double rhoz2starR = rhoz2R + rhoz2avg * (ustar - uR) / cavg;
		
		
		// Set jump over wave 1
		
		WL = conserved_to_primitives(eosparams, UL);
		
		WR(0) = rhoz1starL;
		WR(1) = rhoz2starL;
		WR(2) = ustar;
		WR(3) = vL;
		WR(4) = pstar;
		WR(5) = zL;
		
		jumps[0] = WR - WL;
		
		
		// Set jump over wave 2
		
		WL = WR;
		
		WR(0) = rhoz1starR;
		WR(1) = rhoz2starR;
		WR(2) = ustar;
		WR(3) = vR;
		WR(4) = pstar;
		WR(5) = zR;
		
		jumps[1] = WR - WL;
		
		
		// Set jump over wave 3
		
		WL = WR;
		
		WR = conserved_to_primitives(eosparams, UR);
		
		jumps[2] = WR - WL;
	}
	
	std::shared_ptr<riemann_solver_base> clone ()
	{	
		return std::make_shared<PVRS_riemann_solver>(*this);
	}
	
};


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
