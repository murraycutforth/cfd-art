#include "sim_info.hpp"
#include "typedefs.hpp"
#include <cmath>
#include <algorithm>

inline double fastsweep_positivecellupdate (const sim_info& params, griddoubletype& phi, int i, int j)
{
	double phi_ij = phi[i][j];
	
	double a = std::min(phi[i][j-1], phi[i][j+1]);
	double b = std::min(phi[i-1][j], phi[i+1][j]);
	double phi_bar;
	
	if (a - b <= - params.dx)
	{
		phi_bar = params.dx + a;
	}
	else if (a - b >= params.dy)
	{
		phi_bar = params.dy + b;
	}
	else
	{
		double dxdy = params.dx * params.dy;
		double dxsq = params.dx * params.dx;
		double dysq = params.dy * params.dy;
		
		phi_bar = (a * dysq + b * dxsq + dxdy * sqrt(dxsq + dysq - (a - b) * (a - b))) / (dxsq + dysq);
	}
	
	return fabs(phi_ij) < fabs(phi_bar) ? phi_ij : phi_bar;
}


inline double fastsweep_negativecellupdate (const sim_info& params, griddoubletype& phi, int i, int j)
{
	double phi_ij = phi[i][j];
	
	double a = std::max(phi[i][j-1], phi[i][j+1]);
	double b = std::max(phi[i-1][j], phi[i+1][j]);
	double phi_bar;


	if (b - a <= - params.dx)
	{
		phi_bar = a - params.dx;
	}
	else if (b - a >= params.dy)
	{
		phi_bar = b - params.dy;
	}
	else
	{
		double dxdy = params.dx * params.dy;
		double dxsq = params.dx * params.dx;
		double dysq = params.dy * params.dy;
		
		phi_bar = (a * dysq + b * dxsq - dxdy * sqrt(dxsq + dysq - (a - b) * (a - b))) / (dxsq + dysq);
	}
	
	return fabs(phi_ij) < fabs(phi_bar) ? phi_ij : phi_bar;
}

inline void fastsweep_sdf (const sim_info& params, const gridbooltype& gridfrozencells, griddoubletype& phi)
{
	// Assume that the sign of phi in unfrozen cells is consistent with the side of the interface which they're on
		
	for (int j=params.numGC; j<params.Nx + params.numGC; j++)
	{
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			if (! gridfrozencells[i][j])
			{
				if (phi[i][j] > 0.0)
				{
					phi[i][j] = fastsweep_positivecellupdate(params, phi, i, j);
				}
				else
				{
					phi[i][j] = fastsweep_negativecellupdate(params, phi, i, j);
				}
			}
		}
	}
	
	for (int j=params.Nx + params.numGC - 1; j>=params.numGC; j--)
	{
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			if (! gridfrozencells[i][j])
			{
				if (phi[i][j] > 0.0)
				{
					phi[i][j] = fastsweep_positivecellupdate(params, phi, i, j);
				}
				else
				{
					phi[i][j] = fastsweep_negativecellupdate(params, phi, i, j);
				}
			}
		}
	}
	
	for (int j=params.Nx + params.numGC - 1; j>=params.numGC; j--)
	{
		for (int i=params.Ny + params.numGC - 1; i>=params.numGC; i--)
		{
			if (! gridfrozencells[i][j])
			{
				if (phi[i][j] > 0.0)
				{
					phi[i][j] = fastsweep_positivecellupdate(params, phi, i, j);
				}
				else
				{
					phi[i][j] = fastsweep_negativecellupdate(params, phi, i, j);
				}
			}
		}
	}
	
	for (int j=params.numGC; j<params.Nx + params.numGC; j++)
	{
		for (int i=params.Ny + params.numGC - 1; i>=params.numGC; i--)
		{
			if (! gridfrozencells[i][j])
			{
				if (phi[i][j] > 0.0)
				{
					phi[i][j] = fastsweep_positivecellupdate(params, phi, i, j);
				}
				else
				{
					phi[i][j] = fastsweep_negativecellupdate(params, phi, i, j);
				}
			}
		}
	}
}
