/*
 *	DESCRIPTION:	An implementation of the 2D Euler equations
 *
 *	In this model the conserved variables are U = [rho, rho * u, rho * v, E]^T
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		12/04/2023
 */


#ifndef EULERPROBLEM_H
#define EULERPROBLEM_H

// From PhD-2D-HCLsolver-framework
#include "problem_base.hpp"

// From PhD-Common
#include "sim_info.hpp"
#include "typedefs.hpp"
#include "settings_file.hpp"

// From this project
#include "flux_solver_base.hpp"
#include "stiffened_gas_eos.hpp"

// From STL
#include <vector>

class Euler : public problem_base {
	
protected:
	
	SGparams eosparams;			// EOS parameters of the stiffened gas fluid
	std::shared_ptr<FluxSolverBaseEulerEquations> fs_ptr;	// Algorithm for updating conservative variables
	std::vector<double> time;			// Storage for mass o fluid at each time step
	std::vector<double> mass;
	
	
	// Functions specific to this problem
	
	void set_parameters (std::string test_case, sim_info& params, SGparams& eosparams);

	void add_circle_to_IC (const vectype& W_in, const Eigen::Vector2d& centre, const double R, gridtype& grid, const sim_info& params, const int N = 10);

	void add_rectangle_to_IC (const vectype& W_in, const Eigen::Vector2d& bl, const double width, const double height, gridtype& grid, const sim_info& params, const int N = 10);
	
	void set_circular_IC (const vectype& W_in, const vectype& W_out, const Eigen::Vector2d& centre, const double R, gridtype& grid, const sim_info& params, const int N);
	
	void set_boundary_conditions (gridtype& grid, const sim_info& params);
	
	void gnuplot_output (const gridtype& grid, const sim_info& params, int n, double t);
	
	void gnuplot_schlieren (const gridtype& grid, const sim_info& params, int n, double t);
	
	void gnuplot_masschange (const sim_info& params);
	
	double get_rho (const vectype& U)
	{
		return U(0);
	}
	
	double get_u (const vectype& U)
	{
		return U(1) / get_rho(U);
	}
	
	double get_v (const vectype& U)
	{
		return U(2) / get_rho(U);
	}
	
	double get_e (const vectype& U)
	{
		return U(3) / get_rho(U) - 0.5 * (get_u(U) * get_u(U) + get_v(U) * get_v(U));
	}
	
	
public:
	
	Euler ()
	{}
	
	Euler (double gamma, double pinf, std::shared_ptr<FluxSolverBaseEulerEquations> FS_ptr)
	:
		eosparams (gamma, pinf),
		fs_ptr (FS_ptr)
	{}
	
	
	// Over-ride all pure virtual member functions of problem_base
	
	std::shared_ptr<gridtype> set_ICs (settings_file SF, sim_info& params);
	
	void output (const gridtype& grid, const sim_info& params, int n, double t);
	
	double compute_dt (const gridtype& grid, const sim_info& params, int n, double t);
	
	void pre_sweep (gridtype& grid, const sim_info& params);
	
	void update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t);
	
	void update_col (const gridtype& grid, gridtype& future_grid, const sim_info& params, int j, double dt, double t);
	
	void unsplit_update (const gridtype& grid, gridtype& future_grid, const sim_info& params, int j, double dt, double t){}
	
	void post_sweep (gridtype& grid, gridtype& future_grid, const sim_info& params);
	
	
	// Clone function which returns a shared_ptr<problem_base> object for deep copy of problem between threads
	
	std::shared_ptr<problem_base> clone ()
	{
		return std::make_shared<Euler>(eosparams.gamma, eosparams.pinf, fs_ptr->clone());
	}
};

#endif
