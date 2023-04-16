#include "euler_problem.hpp"
#include "stiffened_gas_eos.hpp"
#include "flux_solver_godunov.hpp"
#include "flux_solver_MUSCLhancock.hpp"
#include "riemann_solver_HLLC.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <memory>
#include <fstream>
#include <string>

void Euler :: set_parameters (std::string test_case, sim_info& params, SGparams& eosparams)
{
	if (test_case == "circular_explosion")
	{
		eosparams.gamma = 1.4;
		eosparams.pinf = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 2.0/params.Nx;
		params.dy = 2.0/params.Ny;
		params.T = 1.0;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "reflective";
	}
	else if (test_case == "multiple_circles_1" || test_case == "multiple_circles_2" || test_case == "cross" || test_case == "internet" || test_case == "totem_pole_2")
	{
		eosparams.gamma = 1.4;
		eosparams.pinf = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 1.0;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (test_case == "totem_pole_1")
	{
		eosparams.gamma = 1.4;
		eosparams.pinf = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 1.0;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (test_case == "totem_pole_3")
	{
		eosparams.gamma = 1.4;
		eosparams.pinf = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 0.5/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 1.0;
		params.BC_L = "reflective";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	//else if (test_case == "underwater_shocked_bubble")
	//{
	//	// See "Practical techniques in ghost fluid method..", Xu, Communications in computational physics, 2016
	//	
	//	eosparams.gamma1 = 7.15;
	//	eosparams.gamma2 = 1.4;
	//	eosparams.pinf1 = 3309.0;
	//	eosparams.pinf2 = 0.0;
	//	
	//	params.x0 = 0.0;
	//	params.y0 = 0.0;
	//	params.dx = 12.0/params.Nx;
	//	params.dy = 12.0/params.Ny;
	//	params.T = 0.05;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "transmissive";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "transmissive";
	//}
	//else if (test_case == "underwater_explosion")
	//{
	//	eosparams.gamma1 = 1.4;
	//	eosparams.gamma2 = 7.15;
	//	eosparams.pinf1 = 0.0;
	//	eosparams.pinf2 = 3.309e8;
	//	
	//	params.x0 = -5.0;
	//	params.y0 = -5.0;
	//	params.dx = 10.0/params.Nx;
	//	params.dy = 10.0/params.Ny;
	//	params.T = 0.005;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "transmissive";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "transmissive";
	//}
	//else if (test_case == "underwater_explosion_modified")
	//{
	//	eosparams.gamma1 = 1.4;
	//	eosparams.gamma2 = 7.15;
	//	eosparams.pinf1 = 0.0;
	//	eosparams.pinf2 = 3.309e8;
	//	
	//	params.x0 = -5.0;
	//	params.y0 = -5.0;
	//	params.dx = 10.0/params.Nx;
	//	params.dy = 10.0/params.Ny;
	//	params.T = 0.010;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "reflective";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "reflective";
	//}
	//else if (test_case == "shocked_helium_bubble")
	//{
	//	eosparams.gamma1 = 1.4;
	//	eosparams.gamma2 = 1.667;
	//	eosparams.pinf1 = 0.0;
	//	eosparams.pinf2 = 0.0;
	//	
	//	params.x0 = 0.0;
	//	params.y0 = 0.0;
	//	params.dx = 325.0/params.Nx;
	//	params.dy = 89.0/params.Ny;
	//	params.T = 280.0;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "reflective";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "reflective";
	//}
	//else if (test_case == "shocked_R22_bubble")
	//{
	//	eosparams.gamma1 = 1.4;
	//	eosparams.gamma2 = 1.249;

	//	params.x0 = 0.0;
	//	params.y0 = 0.0;
	//	params.dx = 0.445/params.Nx;
	//	params.dy = 0.089/params.Ny;
	//	params.T = 0.001080;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "reflective";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "reflective";
	//}
	//else if (test_case == "shocked_SF6")
	//{
	//	eosparams.gamma1 = 1.4;
	//	eosparams.gamma2 = 1.076;
	//	eosparams.pinf1 = 0.0;
	//	eosparams.pinf2 = 0.0;
	//	
	//	params.x0 = 0.0;
	//	params.y0 = 0.0;
	//	params.dx = 0.45/params.Nx;
	//	params.dy = 0.2/params.Ny;
	//	params.T = 0.25;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "reflective";
	//	params.BC_R = "reflective";
	//	params.BC_B = "reflective";
	//}
	//else if (test_case == "RMI_SF6")
	//{
	//	// From "A volume of fluid method based ghost fluid method for compressible multi-fluid flows" - Computers & Fluids - 2014
	//	
	//	eosparams.gamma1 = 1.4;
	//	eosparams.gamma2 = 1.093;
	//	eosparams.pinf1 = 0.0;
	//	eosparams.pinf2 = 0.0;
	//	
	//	params.x0 = 0.0;
	//	params.y0 = 0.0;
	//	params.dx = 4.0/params.Nx;
	//	params.dy = 0.5/params.Ny;
	//	params.T = 10.0;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "reflective";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "reflective";
	//}
	//else if (test_case == "tin_air_implosion")
	//{
	//	eosparams.gamma1 = 3.27;
	//	eosparams.gamma2 = 1.4;
	//	eosparams.pinf1 = 149500.0;
	//	eosparams.pinf2 = 0.0;
	//	
	//	params.x0 = 0.0;
	//	params.y0 = -25.0;
	//	params.dx = 25.0/params.Nx;
	//	params.dy = 50.0/params.Ny;
	//	params.T = 0.08;
	//	params.BC_L = "reflective";
	//	params.BC_T = "transmissive";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "transmissive";
	//}
	//else if (test_case == "TSTM")
	//{
	//	eosparams.gamma1 = 1.5;
	//	eosparams.gamma2 = 1.4;
	//	eosparams.pinf1 = 0.0;
	//	eosparams.pinf2 = 0.0;
	//	
	//	params.x0 = 0.0;
	//	params.y0 = 0.0;
	//	params.dx = 7.0/params.Nx;
	//	params.dy = 3.0/params.Ny;
	//	params.T = 8.0;
	//	params.BC_L = "transmissive";
	//	params.BC_T = "transmissive";
	//	params.BC_R = "transmissive";
	//	params.BC_B = "transmissive";
	//}
	else
	{
		assert(!"[Euler] Invalid test_case in settings file.");
	}
}

//void allaire_diffuse :: set_halfspace_IC (const vectype& U_under, const double a, const double b, const double c, gridtype& grid, const sim_info& params)
//{
//	/*
//	 * Set every cell where a*x + b*y <= c to the state U_under.
//	 */
//	 
//	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
//	{
//		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
//		{
//			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
//			
//			if (cc(0)*a + cc(1)*b <= c)
//			{
//				grid[i][j] = U_under;
//			}
//		}
//	}
//}
//
//void allaire_diffuse :: set_planar_IC (const vectype& U_under, const vectype& U_over, const double a, const double b, const double c, gridtype& grid, const sim_info& params)
//{
//	/*
//	 * Set every cell where a*x + b*y <= c to the state U_under,
//	 * otherwise U_over. Useful for planar initial conditions.
//	 */
//	 
//	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
//	{
//		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
//		{
//			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
//			
//			if (cc(0)*a + cc(1)*b <= c)
//			{
//				grid[i][j] = U_under;
//			}
//			else
//			{
//				grid[i][j] = U_over;
//			}
//		}
//	}
//}


void Euler :: add_rectangle_to_IC (const vectype& W_in, const Eigen::Vector2d& bl, const double width, const double height, gridtype& grid, const sim_info& params, const int N)
{
	/*
	 * Set state using rectangular boundary conditions
	 */
	vectype U_in = primitives_to_conserved_euler_equations(eosparams, W_in);
	 
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			int totalnumsamples = N*N;
			int numinside = 0;
			double delx = params.dx/N;
			double dely = params.dy/N;
			
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			Eigen::Vector2d BL;
			BL(0) = cc(0) - 0.5 * params.dx;
			BL(1) = cc(1) - 0.5 * params.dy;
			Eigen::Vector2d samplepos;
			
			for (int a=0; a<N; a++)
			{
				for (int b=0; b<N; b++)
				{
					samplepos(0) = BL(0) + (a + 0.5) * delx;
					samplepos(1) = BL(1) + (b + 0.5) * dely;
					
					samplepos -= bl;
					
					if ((0.0 <= samplepos(0)) & 
					     (samplepos(0) <= width) &
					     (0.0 <= samplepos(1)) & 
					     (samplepos(1) <= height)
					   )
					{
						numinside++;
					}
				}
			}

			double frac = double(numinside) / totalnumsamples;
			
			vectype U_cell = frac * U_in;
			grid[i][j] += U_cell;
		}
	}
}
	



void Euler :: add_circle_to_IC (const vectype& W_in, const Eigen::Vector2d& centre, const double R, gridtype& grid, const sim_info& params, const int N)
{
	/*
	 * Set state according to weighting of volume fraction inside/outside
	 * the circle, according to Monte Carlo estimate using N^2 samples
	 */

	vectype U_in = primitives_to_conserved_euler_equations(eosparams, W_in);
	 
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			double frac;

			if ((cc - centre).norm() > R + 2.0 * params.dx)  // ASSUMPTION: dx and dy are similar
			{
				frac = 0.0;
			}
			else
			{

				int totalnumsamples = N*N;
				int numinside = 0;
				double delx = params.dx/N;
				double dely = params.dy/N;
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				Eigen::Vector2d samplepos;
				
				for (int a=0; a<N; a++)
				{
					for (int b=0; b<N; b++)
					{
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						samplepos -= centre;
						
						if (samplepos.norm() <= R) numinside++;
					}
				}

				frac = double(numinside) / totalnumsamples;
			}
			
			vectype U_cell = frac * U_in;
			grid[i][j] += U_cell;
		}
	}
}


void Euler :: set_circular_IC (const vectype& W_in, const vectype& W_out, const Eigen::Vector2d& centre, const double R, gridtype& grid, const sim_info& params, const int N)
{
	/*
	 * Set state according to weighting of volume fraction inside/outside
	 * the circle, according to Monte Carlo estimate using N^2 samples
	 */
	 
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			int totalnumsamples = N*N;
			int numinside = 0;
			double delx = params.dx/N;
			double dely = params.dy/N;
			
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			Eigen::Vector2d BL;
			BL(0) = cc(0) - 0.5 * params.dx;
			BL(1) = cc(1) - 0.5 * params.dy;
			Eigen::Vector2d samplepos;
			
			for (int a=0; a<N; a++)
			{
				for (int b=0; b<N; b++)
				{
					samplepos(0) = BL(0) + (a + 0.5) * delx;
					samplepos(1) = BL(1) + (b + 0.5) * dely;
					
					samplepos -= centre;
					
					if (samplepos.norm() <= R) numinside++;
				}
			}

			double frac = double(numinside)/totalnumsamples;
			
			vectype W = frac * W_in + (1.0 - frac) * W_out;
			grid[i][j] = primitives_to_conserved_euler_equations(eosparams, W);
		}
	}
}


std::shared_ptr<gridtype> Euler :: set_ICs (settings_file SF, sim_info& params)
{
	params.Nx = SF.Nx;
	params.Ny = SF.Ny;
	params.CFL = SF.CFL;
	params.outputname = SF.basename;
	params.output_freq = SF.output_freq;
	set_parameters(SF.test_case, params, eosparams);
	
	
	// Number of ghost cells needed
	
	if (SF.flux_solver == "Godunov")
	{
		params.numGC = 1;
		params.stclsize = 1;
	}
	else if (SF.flux_solver == "MUSCL")
	{
		params.numGC = 2;
		params.stclsize = 2;
	}
	else
	{
		assert(!"[Euler] Invalid flux_solver in settings file.");
	}
	
	
	// Riemann solver object
	
	std::shared_ptr<riemann_solver_base_euler_equations> RS_ptr = nullptr;	
	
	if (SF.riemann_solver == "HLLC")
	{
		RS_ptr = std::make_shared<HLLC_riemann_solver_euler_equations>(params, eosparams.gamma, eosparams.pinf);
	}
	else
	{
		assert(!"[Euler] Invalid riemann_solver in settings file.");
	}
	
	
	// Flux solver object
	
	if (SF.flux_solver == "Godunov")
	{
		fs_ptr = std::make_shared<GodunovEulerEquations>(RS_ptr, params, eosparams.gamma, eosparams.pinf);
	}
	else if (SF.flux_solver == "MUSCL")
	{
		fs_ptr = std::make_shared<MusclEulerEquations>(RS_ptr, params, eosparams.gamma, eosparams.pinf);
	}
	else
	{
		assert(!"[Euler] Invalid flux_solver in settings file.");
	}
			
	
	gridtype ICgrid (SF.Ny + 2 * params.numGC, rowtype(SF.Nx + 2 * params.numGC, vectype(4)));
	
	
	// Initial states
		
	if (SF.test_case == "circular_explosion")
	{
		vectype Win (4);
		Win << 1.0, 0.0, 0.0, 1.0;
		vectype Wout (4);
		Wout << 0.125, 0.0, 0.0, 0.1;
		
		Eigen::Vector2d centre;
		centre << 1.0, 1.0;
		
		double R = 0.4;
		
		set_circular_IC(Win, Wout, centre, R, ICgrid, params, 10);		
	}
	else if (SF.test_case == "multiple_circles_1")
	{
		vectype W_background (4);
		W_background << 0.125, 0.0, 0.0, 0.1;

		vectype U_background (4);
		U_background = primitives_to_conserved_euler_equations(eosparams, W_background);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				ICgrid[i][j] = U_background;
			}
		}

		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.5, 0.4), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.3, 0.7), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.2, 0.7), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.6, 0.5), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.1, 0.3), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.3, 0.5), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.4, 0.6), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.5, 0.4), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.5, 0.3), 0.1, ICgrid, params);		
	}
	else if (SF.test_case == "multiple_circles_2")
	{
		vectype W_background (4);
		W_background << 0.125, 0.0, 0.0, 0.1;

		vectype U_background (4);
		U_background = primitives_to_conserved_euler_equations(eosparams, W_background);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				ICgrid[i][j] = U_background;
			}
		}

		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(0.0, 1.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(1.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), Eigen::Vector2d(1.0, 1.0), 0.2, ICgrid, params);		

		add_circle_to_IC(Eigen::Vector4d(0.4, 0.0, 0.0, 0.0), Eigen::Vector2d(0.8, 0.5), 0.1, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(0.3, 0.0, 0.0, 0.0), Eigen::Vector2d(0.7, 0.5), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(0.2, 0.0, 0.0, 0.0), Eigen::Vector2d(0.6, 0.5), 0.3, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(0.1, 0.0, 0.0, 0.0), Eigen::Vector2d(0.5, 0.5), 0.4, ICgrid, params);		
	}
	else if (SF.test_case == "cross")
	{
		vectype W_background (4);
		W_background << 0.125, 0.0, 0.0, 0.1;

		vectype U_background (4);
		U_background = primitives_to_conserved_euler_equations(eosparams, W_background);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				ICgrid[i][j] = U_background;
			}
		}

		// Incoming shock from each corner
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.0, 1.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(1.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(1.0, 1.0), 0.2, ICgrid, params);		

		add_circle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.5, 0.1), 0.11, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.5, 0.9), 0.11, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.3, 0.7), 0.11, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.7, 0.7), 0.11, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.45, 0.1), 0.1, 0.8, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.3, 0.65), 0.4, 0.1, ICgrid, params);		
	}
	else if (SF.test_case == "totem_pole_1")
	{
		vectype W_background (4);
		W_background << 0.125, 0.0, 0.0, 0.1;

		vectype U_background (4);
		U_background = primitives_to_conserved_euler_equations(eosparams, W_background);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				ICgrid[i][j] = U_background;
			}
		}

		// Incoming shock from each corner
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.0, 1.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(1.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(1.0, 1.0), 0.2, ICgrid, params);		

		// Totem pole
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.45, 0.1), 0.1, 0.8, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.3, 0.65), 0.4, 0.1, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.2, 0.0, 0.0, 0.0), Eigen::Vector2d(0.4, 0.5), 0.2, 0.05, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.2, 0.0, 0.0, 0.0), Eigen::Vector2d(0.4, 0.375), 0.2, 0.075, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.2, 0.0, 0.0, 0.0), Eigen::Vector2d(0.4, 0.25), 0.2, 0.1, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(0.45, 0.1), 0.1, 0.1, ICgrid, params);		
	}
	else if (SF.test_case == "totem_pole_3")
	{
		vectype W_background (4);
		W_background << 0.125, 0.0, 0.0, 0.1;

		vectype U_background (4);
		U_background = primitives_to_conserved_euler_equations(eosparams, W_background);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				ICgrid[i][j] = U_background;
			}
		}

		// Incoming shock from each corner
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.5, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.5, 1.0), 0.2, ICgrid, params);		

		// Totem pole
		// Main vertical
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(-0.05, 0.0), 0.1, 0.9, ICgrid, params);		
		// Big cross bar at 0.65 to 0.75
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(-0.2, 0.65), 0.4, 0.1, ICgrid, params);		
		// Small cross bar at 0.5 to 0.55
		add_rectangle_to_IC(Eigen::Vector4d(0.2, 0.0, 0.0, 0.0), Eigen::Vector2d(-0.1, 0.5), 0.2, 0.1, ICgrid, params);		
		// Small cross bar at 0.35 to 0.4
		add_rectangle_to_IC(Eigen::Vector4d(0.4, 0.0, 0.0, 0.0), Eigen::Vector2d(-0.15, 0.375), 0.3, 0.05, ICgrid, params);		
		// Small cross bar at 0.2 to 0.25
		add_rectangle_to_IC(Eigen::Vector4d(0.2, 0.0, 0.0, 0.0), Eigen::Vector2d(-0.1, 0.2), 0.2, 0.1, ICgrid, params);		
		// Big cross bar at 0 to 0.1
		add_rectangle_to_IC(Eigen::Vector4d(0.5, 0.0, 0.0, 0.0), Eigen::Vector2d(-0.15, 0.0), 0.3, 0.15, ICgrid, params);		
	}
	else if (SF.test_case == "totem_pole_2")
	{
		vectype W_background (4);
		W_background << 0.125, 0.0, 0.0, 0.1;

		vectype U_background (4);
		U_background = primitives_to_conserved_euler_equations(eosparams, W_background);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				ICgrid[i][j] = U_background;
			}
		}

		// Incoming shock from each corner
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(0.0, 1.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(1.0, 0.0), 0.2, ICgrid, params);		
		add_circle_to_IC(Eigen::Vector4d(1.0, 0.0, 0.0, 5.0), Eigen::Vector2d(1.0, 1.0), 0.2, ICgrid, params);		

		// Totem pole
		add_rectangle_to_IC(Eigen::Vector4d(0.7, 0.0, 0.0, 0.0), Eigen::Vector2d(0.45, 0.1), 0.1, 0.8, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.7, 0.0, 0.0, 0.0), Eigen::Vector2d(0.3, 0.65), 0.4, 0.1, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.4, 0.0, 0.0, 0.0), Eigen::Vector2d(0.45, 0.575), 0.2, 0.05, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.4, 0.0, 0.0, 0.0), Eigen::Vector2d(0.35, 0.5), 0.2, 0.05, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.4, 0.0, 0.0, 0.0), Eigen::Vector2d(0.4, 0.375), 0.2, 0.075, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.4, 0.0, 0.0, 0.0), Eigen::Vector2d(0.35, 0.25), 0.3, 0.1, ICgrid, params);		
		add_rectangle_to_IC(Eigen::Vector4d(0.7, 0.0, 0.0, 0.0), Eigen::Vector2d(0.45, 0.1), 0.1, 0.1, ICgrid, params);		
	}
	else if (SF.test_case == "internet")
	{
		vectype W_background (4);
		W_background << 0.125, 0.0, 0.0, 0.1;

		vectype U_background (4);
		U_background = primitives_to_conserved_euler_equations(eosparams, W_background);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				ICgrid[i][j] = U_background;
			}
		}

		int num_people = 16;
		double dx = 1.0 / num_people;
		std::random_device rd;  // Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dis(0.0, 1.0);

		for (int i=0; i<num_people; i++)
		{
			for (int j=0; j<num_people; j++)
			{
				//Eigen::Vector2d centre (i * dx + dis(gen) * dx, j * dx + dis(gen) * dx);
				Eigen::Vector2d centre ((i + 0.5) * dx + dis(gen) * 0.0, (j + 0.5) * dx + dis(gen) * 0.0);
				double radius (0.02);
				double density (5.0 * dis(gen));
				double pressure  = density;
				Eigen::Vector4d W (density, 0.0, 0.0, pressure);

				add_circle_to_IC(W, centre, radius, ICgrid, params);
			}
		}

	}

	//else if (SF.test_case == "shocked_helium_bubble")
	//{
	//	vectype W_helium (6);
	//	W_helium << 0.0, 0.138, 0.0, 0.0, 1.0, 0.0;
	//	
	//	vectype W_preshock (6);
	//	W_preshock << 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
	//	
	//	vectype W_postshock (6);
	//	W_postshock << 1.3764, 0.0, -0.394, 0.0, 1.5698, 1.0;
	//	
	//	Eigen::Vector2d centre;
	//	centre << 175.0, 44.5;
	//	
	//	double R = 25.0;
	//	
	//	set_circular_IC(W_helium, W_preshock, centre, R, ICgrid, params, 10);
	//	
	//	set_halfspace_IC(primitives_to_conserved(eosparams, W_postshock), -1.0, 0.0, -225.0, ICgrid, params);	
	//}
	//else if (SF.test_case == "shocked_R22_bubble")
	//{
	//	vectype W_R22 (6);
	//	W_R22 << 0.0, 3.863, 0.0, 0.0, 1.01325e5, 0.0;
	//	
	//	vectype W_preshock (6);
	//	W_preshock << 1.225, 0.0, 0.0, 0.0, 1.01325e5, 1.0;
	//	
	//	vectype W_postshock (6);
	//	W_postshock << 1.686, 0.0, -113.5, 0.0, 1.59e5, 1.0;
	//	
	//	Eigen::Vector2d centre;
	//	centre << 0.225, 0.0445;
	//	
	//	double R = 0.025;
	//	
	//	set_circular_IC(W_R22, W_preshock, centre, R, ICgrid, params, 25);
	//	
	//	set_halfspace_IC(primitives_to_conserved(eosparams, W_postshock), -1.0, 0.0, -0.275, ICgrid, params);	
	//}
	//else if (SF.test_case == "underwater_shocked_bubble")
	//{
	//	vectype W_air (6);
	//	W_air << 0.0, 0.0012, 0.0, 0.0, 1.0, 0.0;
	//	
	//	vectype W_preshock (6);
	//	W_preshock << 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
	//	
	//	vectype W_postshock (6);
	//	W_postshock << 1.31, 0.0, 67.32, 0.0, 19000.0, 1.0;
	//	
	//	Eigen::Vector2d centre;
	//	centre << 6.0, 6.0;
	//	
	//	double R = 3.0;
	//	
	//	set_circular_IC(W_air, W_preshock, centre, R, ICgrid, params, 10);
	//	
	//	set_halfspace_IC(primitives_to_conserved(eosparams, W_postshock), 1.0, 0.0, 2.4, ICgrid, params);
	//}
	//else if (SF.test_case == "underwater_explosion" || SF.test_case == "underwater_explosion_modified")
	//{
	//	vectype W_water (6);
	//	W_water << 0.0, 1000.0, 0.0, 0.0, 1.0e5, 0.0;
	//	
	//	vectype W_air (6);
	//	W_air << 1.0, 0.0, 0.0, 0.0, 1.0e5, 1.0;
	//	
	//	vectype W_airbubble (6);
	//	W_airbubble << 1270.0, 0.0, 0.0, 0.0, 8.29e8, 1.0;
	//	
	//	Eigen::Vector2d centre;
	//	centre << 0.0, 0.0;
	//	
	//	double R = 1.0;
	//	
	//	set_circular_IC(W_airbubble, W_water, centre, R, ICgrid, params, 10);
	//	
	//	set_halfspace_IC(primitives_to_conserved(eosparams, W_air), 0.0, -1.0, -2.5, ICgrid, params);
	//}
	//else if (SF.test_case == "shocked_SF6")
	//{
	//	double rho_preshock = 1.153;
	//	double p_preshock = 9.6856;
	//	double u_preshock = 0.0;
	//	double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
	//	double rho_postshock = 1.6672;
	//	double p_postshock = 16.3256;
	//	double u_postshock = 1.33273;
	//	double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
	//	double rho2 = 5.805;
	//	double p2 = 9.6856;
	//	double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
	//	double v = 0.0;
	//	double z;
	//	
	//	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	//	{
	//		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
	//		{
	//			Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
	//			
	//			double rho1, u, e1;
	//			
	//			if (cc(0) > 0.05)
	//			{
	//				rho1 = rho_preshock;
	//				u = u_preshock;
	//				e1 = e_preshock;
	//			}
	//			else
	//			{
	//				rho1 = rho_postshock;
	//				u = u_postshock;
	//				e1 = e_postshock;
	//			}
	//			
	//			
	//			// Set z as fraction of area inside rectangular region [0.1, 0.25] x [0.0, 0.1]
	//			
	//			int numsamples = 10;
	//			int totalnumsamples = numsamples*numsamples;
	//			int numinside = 0;
	//			double delx = params.dx/numsamples;
	//			double dely = params.dy/numsamples;
	//			
	//			Eigen::Vector2d BL;
	//			BL(0) = cc(0) - 0.5 * params.dx;
	//			BL(1) = cc(1) - 0.5 * params.dy;
	//			
	//			for (int a=0; a<numsamples; a++)
	//			{
	//				for (int b=0; b<numsamples; b++)
	//				{
	//					Eigen::Vector2d samplepos;
	//					samplepos(0) = BL(0) + (a + 0.5) * delx;
	//					samplepos(1) = BL(1) + (b + 0.5) * dely;
	//					
	//					if (samplepos(0) >= 0.1 && samplepos(0) <= 0.25
	//						&& samplepos(1) >= 0.0 && samplepos(1) <= 0.1) 
	//					{
	//						numinside++;
	//					}
	//				}
	//			}

	//			z = 1.0 - double(numinside)/totalnumsamples;
	//			
	//			ICgrid[i][j](0) = z * rho1;
	//			ICgrid[i][j](1) = (1.0 - z) * rho2;
	//			ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
	//			ICgrid[i][j](5) = z;
	//		}
	//	}		
	//}
	//else if (SF.test_case == "TSTM")
	//{
	//	double rho_preshock = 0.125;
	//	double p_preshock = 0.1;
	//	double u_preshock = 0.0;
	//	double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
	//	double rho_postshock = 1.0;
	//	double p_postshock = 1.0;
	//	double u_postshock = 0.0;
	//	double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
	//	double rho2 = 1.0;
	//	double p2 = 0.1;
	//	double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
	//	double v = 0.0;
	//	double z;
	//	
	//	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	//	{
	//		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
	//		{
	//			Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
	//			
	//			double rho1, u, e1;
	//			
	//			if (cc(0) > 1.0)
	//			{
	//				rho1 = rho_preshock;
	//				u = u_preshock;
	//				e1 = e_preshock;
	//			}
	//			else
	//			{
	//				rho1 = rho_postshock;
	//				u = u_postshock;
	//				e1 = e_postshock;
	//			}
	//			
	//			
	//			// Set z as fraction of area inside rectangular region [1.0, 8.0] x [0.0, 1.5]
	//			
	//			int numsamples = 10;
	//			int totalnumsamples = numsamples*numsamples;
	//			int numinside = 0;
	//			double delx = params.dx/numsamples;
	//			double dely = params.dy/numsamples;
	//			
	//			Eigen::Vector2d BL;
	//			BL(0) = cc(0) - 0.5 * params.dx;
	//			BL(1) = cc(1) - 0.5 * params.dy;
	//			
	//			for (int a=0; a<numsamples; a++)
	//			{
	//				for (int b=0; b<numsamples; b++)
	//				{
	//					Eigen::Vector2d samplepos;
	//					samplepos(0) = BL(0) + (a + 0.5) * delx;
	//					samplepos(1) = BL(1) + (b + 0.5) * dely;
	//					
	//					if (samplepos(0) >= 1.0 && samplepos(0) <= 8.0
	//						&& samplepos(1) >= -1.0 && samplepos(1) <= 1.5) 
	//					{
	//						numinside++;
	//					}
	//				}
	//			}

	//			z = 1.0 - double(numinside)/totalnumsamples;
	//			
	//			ICgrid[i][j](0) = z * rho1;
	//			ICgrid[i][j](1) = (1.0 - z) * rho2;
	//			ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
	//			ICgrid[i][j](5) = z;
	//		}
	//	}		
	//}
	//else if (SF.test_case == "RMI_SF6")
	//{
	//	double rho_preshock = 1.0;
	//	double p_preshock = 1.0;
	//	double u_preshock = 0.0;
	//	double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
	//	double rho_postshock = 1.411;
	//	double p_postshock = 1.628;
	//	double u_postshock = -0.39;
	//	double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
	//	double rho2 = 5.04;
	//	double p2 = 1.0;
	//	double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
	//	double v = 0.0;
	//	double z;
	//	
	//	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	//	{
	//		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
	//		{
	//			Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
	//			
	//			double rho1, u, e1;
	//			
	//			if (cc(0) < 3.2)
	//			{
	//				rho1 = rho_preshock;
	//				u = u_preshock;
	//				e1 = e_preshock;
	//			}
	//			else
	//			{
	//				rho1 = rho_postshock;
	//				u = u_postshock;
	//				e1 = e_postshock;
	//			}
	//			
	//			
	//			// Set z as fraction of area inside sinusoidal interface
	//			
	//			double x1 = 2.9;
	//			double eps = 0.2;
	//			int numsamples = 10;
	//			int totalnumsamples = numsamples*numsamples;
	//			int numinside = 0;
	//			double delx = params.dx/numsamples;
	//			double pi = atan(1.0) * 4.0;
	//			double dely = params.dy/numsamples;
	//			
	//			Eigen::Vector2d BL;
	//			BL(0) = cc(0) - 0.5 * params.dx;
	//			BL(1) = cc(1) - 0.5 * params.dy;
	//			
	//			for (int a=0; a<numsamples; a++)
	//			{
	//				for (int b=0; b<numsamples; b++)
	//				{
	//					Eigen::Vector2d samplepos;
	//					samplepos(0) = BL(0) + (a + 0.5) * delx;
	//					samplepos(1) = BL(1) + (b + 0.5) * dely;
	//					
	//					double interfacex = x1 - eps * sin(2.0 * pi * (samplepos(1) + 0.25));
	//					
	//					
	//					if (interfacex > samplepos(0)) 
	//					{
	//						numinside++;
	//					}
	//				}
	//			}

	//			z = 1.0 - double(numinside)/totalnumsamples;
	//			
	//			ICgrid[i][j](0) = z * rho1;
	//			ICgrid[i][j](1) = (1.0 - z) * rho2;
	//			ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
	//			ICgrid[i][j](5) = z;
	//		}
	//	}		
	//}
	//else if (SF.test_case == "tin_air_implosion")
	//{
	//	double rho_preshock = 7.28;
	//	double p_preshock = 1.0;
	//	double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
	//	double rho_postshock = 11.84;
	//	double p_postshock = 1000000.0;
	//	double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
	//	double rho2 = 0.001;
	//	double p2 = 1.0;
	//	double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
	//	double v = 0.0;
	//	double z;
	//	
	//	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	//	{
	//		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
	//		{
	//			// Set tin state as fraction of area inside circle of radius 24
	//			
	//			int numsamples = 10;
	//			int totalnumsamples = numsamples*numsamples;
	//			int numinside = 0;
	//			double delx = params.dx/numsamples;
	//			double dely = params.dy/numsamples;
	//			
	//			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
	//			Eigen::Vector2d BL;
	//			BL(0) = cc(0) - 0.5 * params.dx;
	//			BL(1) = cc(1) - 0.5 * params.dy;
	//			
	//			for (int a=0; a<numsamples; a++)
	//			{
	//				for (int b=0; b<numsamples; b++)
	//				{
	//					Eigen::Vector2d samplepos;
	//					samplepos(0) = BL(0) + (a + 0.5) * delx;
	//					samplepos(1) = BL(1) + (b + 0.5) * dely;
	//					
	//					if (samplepos.norm() <= 24.0) numinside++;
	//				}
	//			}
	//			
	//			double insideratio = double(numinside) / totalnumsamples;
	//			
	//			double rho1, u = 0.0, e1;
	//			rho1 = insideratio * rho_preshock + (1.0 - insideratio) * rho_postshock;
	//			e1 = insideratio * e_preshock + (1.0 - insideratio) * e_postshock;
	//			
	//			
	//			// Set z as fraction of area inside circle of radius 20 at (0, 0)
	//			
	//			numsamples = 10;
	//			totalnumsamples = numsamples*numsamples;
	//			numinside = 0;
	//			delx = params.dx/numsamples;
	//			dely = params.dy/numsamples;
	//			
	//			cc = params.cellcentre_coord(i, j);
	//			BL(0) = cc(0) - 0.5 * params.dx;
	//			BL(1) = cc(1) - 0.5 * params.dy;
	//			
	//			for (int a=0; a<numsamples; a++)
	//			{
	//				for (int b=0; b<numsamples; b++)
	//				{
	//					Eigen::Vector2d samplepos;
	//					samplepos(0) = BL(0) + (a + 0.5) * delx;
	//					samplepos(1) = BL(1) + (b + 0.5) * dely;
	//					
	//					double theta = atan2(samplepos(1), samplepos(0));
	//					theta -= atan(1) * 2;
	//					double r_interface = 20.0 + 0.4 * cos(22 * theta) + 0.4 * cos(17 * theta) + 0.3 * cos(29 * theta);
	//					
	//					if (samplepos.norm() <= r_interface) numinside++;
	//				}
	//			}

	//			z = 1.0 - double(numinside)/totalnumsamples;
	//			
	//			ICgrid[i][j](0) = z * rho1;
	//			ICgrid[i][j](1) = (1.0 - z) * rho2;
	//			ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
	//			ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
	//			ICgrid[i][j](5) = z;
	//		}
	//	}		
	//}
	else
	{
		assert(!"[allaire_diffuse] Invalid test_case in settings file.");
	}
	
	return std::make_shared<gridtype>(ICgrid);
}


void Euler :: gnuplot_schlieren (const gridtype& grid, const sim_info& params, int n, double t)
{

	std::string filename2 = params.outputname + "-schlieren-" + std::to_string(t) + ".dat";
	std::ofstream outfile2;
	outfile2.open(filename2);

	std::vector<double> allgrads;
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho_L = grid[i][j-1](0);
			double rho_R = grid[i][j+1](0);
			double rho_T = grid[i+1][j](0);
			double rho_B = grid[i-1][j](0);

			double gradx = (rho_R - rho_L) / (2.0 * params.dx);
			double grady = (rho_T - rho_B) / (2.0 * params.dy);

			allgrads.push_back(sqrt(gradx*gradx + grady*grady));
		}
	}

	double maxgrad = *std::max_element(allgrads.begin(), allgrads.end());

	int counter = 0;

	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);

			outfile2 << CC(0) << " " << CC(1) << " " << allgrads[counter] / maxgrad << std::endl;
			counter++;
		}
		outfile2 << std::endl;
	}

	outfile2.close();
	std::cout << "[Euler] Schlieren output to gnuplot complete" << std::endl;
}


void Euler :: gnuplot_output (const gridtype& grid, const sim_info& params, int n, double t)
{
	std::string filename = params.outputname + "-state-" + std::to_string(t) + ".dat";
	std::ofstream outfile;
	outfile.open(filename);
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);
			
			double rho = grid[i][j](0);
			double u = grid[i][j](1) / rho;
			double v = grid[i][j](2) / rho;
			double E = grid[i][j](3);
			double e = E / rho - 0.5 * (u * u + v * v);
			double p = eos::pressure(eosparams.gamma, eosparams.pinf, rho, e);
			
			outfile << CC(0) << " " << CC(1) << " " << rho << " " << u << " " << v << " " << e << " " << p << std::endl;
		}
		outfile << std::endl;
	}
	
	outfile.close();
	std::cout << "[Euler] State output to gnuplot complete" << std::endl;
	
	
}


void Euler :: gnuplot_masschange (const sim_info& params)
{
	std::string filename = params.outputname + "-masschange.dat";
	std::ofstream outfile;
	outfile.open(filename);
		
	for (unsigned int k=0; k<time.size(); k++)
	{
		outfile << time[k] << " " << mass[k] - mass[0] << std::endl;
	}
	
	outfile.close();
	
	std::cout << "[Euler] Mass change output complete" << std::endl;
}
	
void Euler :: output (const gridtype& grid, const sim_info& params, int n, double t)
{
	const double outputinterval = params.T / params.output_freq;
	static double lastoutputtime = 0.0;


	if (n==0 || t==params.T)
	{
		gnuplot_output(grid, params, n, t);
		gnuplot_schlieren(grid, params, n, t);
		
		if (t == params.T)
		{
			gnuplot_masschange(params);
		}
			
	}
	else if (params.output_freq != 0.0)
	{
		if (t - lastoutputtime > outputinterval)
		{
			gnuplot_output(grid, params, n, t);
			gnuplot_schlieren(grid, params, n, t);
			lastoutputtime = t;
		}
	}
}
