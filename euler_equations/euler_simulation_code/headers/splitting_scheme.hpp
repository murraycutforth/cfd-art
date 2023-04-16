/*
 *	DESCRIPTION:	Class definition for the dimensional splitting
 * 			schemes used to solve initial-value 
 * 			hyperbolic conservation law problems in 2D.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		19/07/2017
 */

#ifndef SPLITTING_SCHEME_H
#define SPLITTING_SCHEME_H

#include "typedefs.hpp"
#include "settings_file.hpp"
#include "sim_info.hpp"
#include "problem_base.hpp"
#include <memory>
#include <cassert>
#include <string>
#include <omp.h>

class splitting_scheme_base {
	
	public:
	
	virtual void advance_timestep(std::shared_ptr<problem_base> problem, gridtype& grid, gridtype& future_grid, const sim_info& params, const double dt, const double t) =0;
};

class naive_splitting : public splitting_scheme_base {
	
	/*
	 * First order accurate dimensional splitting
	 */
	
	public:
	
	void advance_timestep(std::shared_ptr<problem_base> problem, gridtype& grid, gridtype& future_grid, const sim_info& params, const double dt, const double t)
	{
		problem->pre_sweep(grid, params);
				
		#pragma omp parallel
		{ 
			std::shared_ptr<problem_base> threadproblem = problem->clone();
			
			#pragma omp for schedule(dynamic) 
			for (int i=params.numGC; i<params.Ny + params.numGC; ++i)
			{
				threadproblem->update_row(grid, future_grid, params, i, dt, t);
			}
		}
		
		problem->post_sweep(grid, future_grid, params);
		problem->pre_sweep(grid, params);
		
		#pragma omp parallel
		{ 
			std::shared_ptr<problem_base> threadproblem = problem->clone();
			
			#pragma omp for schedule(dynamic)
			for (int j=params.numGC; j<params.Nx + params.numGC; ++j)
			{
				threadproblem->update_col(grid, future_grid, params, j, dt, t);
			}
		}
		
		problem->post_sweep(grid, future_grid, params);
	}
};

class strang_splitting : public splitting_scheme_base {
	
	/*
	 * Second order accurate strang splitting approach
	 */
	
	public:
	
	void advance_timestep(std::shared_ptr<problem_base> problem, gridtype& grid, gridtype& future_grid, const sim_info& params, const double dt, const double t)
	{
		problem->pre_sweep(grid, params);
				
		#pragma omp parallel
		{ 
			std::shared_ptr<problem_base> threadproblem = problem->clone();
			
			#pragma omp for schedule(dynamic) 
			for (int i=params.numGC; i<params.Ny + params.numGC; ++i)
			{
				threadproblem->update_row(grid, future_grid, params, i, 0.5 * dt, t);
			}
		}
		
		problem->post_sweep(grid, future_grid, params);
		problem->pre_sweep(grid, params);
		
		#pragma omp parallel
		{ 
			std::shared_ptr<problem_base> threadproblem = problem->clone();
			
			#pragma omp for schedule(dynamic)
			for (int j=params.numGC; j<params.Nx + params.numGC; ++j)
			{
				threadproblem->update_col(grid, future_grid, params, j, dt, t);
			}
		}
		
		problem->post_sweep(grid, future_grid, params);
		problem->pre_sweep(grid, params);
		
		#pragma omp parallel
		{ 
			std::shared_ptr<problem_base> threadproblem = problem->clone();
			
			#pragma omp for schedule(dynamic) 
			for (int i=params.numGC; i<params.Ny + params.numGC; ++i)
			{
				threadproblem->update_row(grid, future_grid, params, i, 0.5 * dt, t);
			}
		}
		
		problem->post_sweep(grid, future_grid, params);
	}
};

std::shared_ptr<splitting_scheme_base> set_splitting_scheme (settings_file SF)
{
	std::shared_ptr<splitting_scheme_base> ss = nullptr;
	
	if (SF.splitting_scheme == "naive_splitting")
	{
		ss = std::make_shared<naive_splitting>();
	}
	else if (SF.splitting_scheme == "strang_splitting")
	{
		ss = std::make_shared<strang_splitting>();
	}
	else
	{
		assert(!"Unknown splitting scheme specified in settingsfile");
	}
	
	return ss;
}

#endif
