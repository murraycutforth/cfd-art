/*
 *	DESCRIPTION:	Class definition for a simulation object. This
 *			is constructed and called externally in order to fully
 *			run a problem as specified in a given settings
 *			file.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		19/07/2017
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include "sim_info.hpp"
#include "settings_file.hpp"
#include "typedefs.hpp"
#include "splitting_scheme.hpp"

class simulation {
	
	public:
	
	settings_file SF;
	std::shared_ptr<problem_base> problem;
	
	simulation (settings_file SF, std::shared_ptr<problem_base> problem)
	:
		SF (SF),
		problem (problem)
	{}
	
	int run_simulation ()
	{
		sim_info params;
		std::shared_ptr<gridtype> grid (problem->set_ICs(SF, params));
		std::shared_ptr<gridtype> future_grid (problem->set_ICs(SF, params));
		std::shared_ptr<splitting_scheme_base> splitting_scheme (set_splitting_scheme(SF));
		
		int n = 0;
		double t = 0.0, dt;
		
		std::cout << std::endl << "[simulation] Beginning time steps.." << std::endl;
		
		while (t < params.T)
		{
			problem->output(*grid, params, n, t);
			
			dt = problem->compute_dt(*grid, params, n, t);
			
			splitting_scheme->advance_timestep(problem, *grid, *future_grid, params, dt, t);
						
			t += dt;
			n++;
			if (fabs(t - params.T) < 1e-12) t = params.T;
			
			std::cout << "[simulation] Time step " << n << " complete. Time t = " << t << "." << std::endl;
		}
		
		std::cout << "[simulation] Time stepping complete." << std::endl;
		
		problem->output(*grid, params, n, t);
		return 0;
	}
};

#endif
			
