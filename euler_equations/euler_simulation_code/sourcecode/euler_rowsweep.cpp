#include "euler_problem.hpp"
#include "sim_info.hpp"
#include "misc_euler.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <string>


void Euler :: update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t)
{
	// Storage for flux and velocity across each internal cell edge
	
	rowtype fluxes (params.Nx + 2 * params.numGC, vectype(4));
	std::vector<double> u_stars (params.Nx + 2 * params.numGC, 0.0);
	// std::vector<double> z_stars (params.Nx + 2 * params.numGC, 0.0);
		
	
	// Compute cell edge fluxes across each edge of this row
	
	std::vector<vectype> stencil (2*params.stclsize);
		
	for (int j=params.numGC; j < params.Nx + params.numGC + 1; j++)
	{
		for (int l = j - params.stclsize; l <= j + params.stclsize - 1; l++)
		{
			stencil[l - j + params.stclsize] = grid[i][l];
		}

		fs_ptr->flux_computation(stencil, fluxes[j], dt, params.dx, u_stars[j]); //z_stars[j]);
	}
		
	
	// Conservative update formula
	
	for (int j=params.numGC; j < params.Nx + params.numGC; j++)
	{
		future_grid[i][j] = grid[i][j] + (dt / params.dx) * (fluxes[j] - fluxes[j+1]);
	}

	
	// Check for and correct any unphysical states
	// This will break conservation if it is triggered - but for the purposes of an art project we don't mind.
	// I have found this to be necessary for very strong shocks with very fine grids
	// This may point to a subtle bug in my code, I'm not sure if MUSCL is theoretically a positivity-preserving scheme?
	
	for (int j=params.numGC; j < params.Nx + params.numGC; j++)
	{
		if (! is_physical_state_euler_equations(eosparams, future_grid[i][j]))
		{
			std::cout << "Found unphysical state" << std::endl;

			// Fix by reverting to density or internal energy from previous step
			//
			double rho = future_grid[i][j](0);
			double u = future_grid[i][j](1) / rho;
			double v = future_grid[i][j](2) / rho;
			double e = future_grid[i][j](3) / rho - 0.5 * (u * u + v * v);

			double rho_old = grid[i][j](0);
			double u_old = grid[i][j](1) / rho_old;
			double v_old = grid[i][j](2) / rho_old;
			double e_old = grid[i][j](3) / rho_old - 0.5 * (u_old * u_old + v_old * v_old);

			if (rho <= 1e-12){
				rho = rho_old;
			}

			if (e <= 1e-12){
				e = e_old;
			}

			double p = eos::pressure(eosparams.gamma, eosparams.pinf, e, rho);

			vectype W (4);
			W(0) = rho;
			W(1) = u;
			W(2) = v;
			W(3) = p;
			future_grid[i][j] = primitives_to_conserved_euler_equations(eosparams, W);
		}
	}
	
	
	// Volume fraction update formula
	
	//for (int j=params.numGC; j < params.Nx + params.numGC; j++)
	//{
	//	future_grid[i][j](5) = zupdate_ptr->zupdate(params.dx, dt, grid[i][j-1](5), grid[i][j](5), grid[i][j+1](5), 
	//						    u_stars[j], u_stars[j+1], z_stars[j], z_stars[j+1]);
	//						    
	//	assert(is_physical_state(eosparams, future_grid[i][j]));
	//}
	
}
