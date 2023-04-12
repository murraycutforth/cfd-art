#include "euler_problem.hpp"
#include "stiffened_gas_eos.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>

double Euler :: compute_dt (const gridtype& grid, const sim_info& params, int n, double t)
{
	double maxu = 0.0, maxv = 0.0, masssum = 0.0;
	
	
	#pragma omp parallel for schedule(dynamic) reduction(+:masssum)
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho = grid[i][j](0);
			double u = grid[i][j](1) / rho;
			double v = grid[i][j](2) / rho;
			double e = grid[i][j](3) / rho - 0.5 * (u * u + v * v);
			double p = eos::pressure(eosparams.gamma, eosparams.pinf, e, rho);
			double c = eos::soundspeed(eosparams.gamma, eosparams.pinf, p, rho);
			maxu = std::max(maxu, fabs(u) + c);
			maxv = std::max(maxv, fabs(v) + c);
			masssum += rho;
		}
	}
	
	
	// Record fluid masses at this time step
	
	time.push_back(t);
	mass.push_back(masssum * params.dx * params.dy);
	
	double CFL = params.CFL;
	
	if (n < 5) CFL = std::min(CFL, 0.2);
	
	double dt = CFL * std::min(params.dx / maxu, params.dy / maxv);
	
	if (t + dt > params.T) dt = params.T - t;
	
	return dt;	
}
