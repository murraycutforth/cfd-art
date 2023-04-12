/*
 *	DESCRIPTION:	Class definition for the abstract base "problem"
 * 			class. This defines how the grid is updated and
 * 			how time steps are calculated.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		19/07/2017
 */

#ifndef PROBLEM_BASE_H
#define PROBLEM_BASE_H

#include "typedefs.hpp"
#include "settings_file.hpp"
#include "sim_info.hpp"
#include <memory>

class problem_base {
	
public:
	
	virtual std::shared_ptr<gridtype> set_ICs (settings_file SF, sim_info& params) =0;
	
	virtual void output (const gridtype& grid, const sim_info& params, int n, double t) =0;
	
	virtual double compute_dt (const gridtype& grid, const sim_info& params, int n, double t) =0;
	
	virtual void pre_sweep (gridtype& grid, const sim_info& params) =0;
	
	virtual void update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t) =0;
	
	virtual void update_col (const gridtype& grid, gridtype& future_grid, const sim_info& params, int j, double dt, double t) =0;
	
	virtual void post_sweep (gridtype& grid, gridtype& future_grid, const sim_info& params) =0;
	
	
	// Clone function which returns a shared_ptr<problem_base> object for deep copy of problem between threads
	
	virtual std::shared_ptr<problem_base> clone () =0;
};

#endif
