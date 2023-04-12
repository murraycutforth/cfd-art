/*
 *	DESCRIPTION:	Class definition for a sim_info object. This
 * 			contains all information needed to define a 2D
 * 			Cartesian grid, as well as other necessary data.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		19/07/2017
 */

#ifndef SIM_INFO_H
#define SIM_INFO_H

#include <Eigen/Dense>
#include <utility>

struct sim_info {
	
	int Nx;
	int Ny;
	double x0;
	double y0;
	double dx;
	double dy;
	int numGC;
	std::string BC_T;
	std::string BC_L;
	std::string BC_B;
	std::string BC_R;
	
	double T;
	double CFL;
	std::string outputname;
	int output_freq;
	int stclsize;

	
	sim_info ()
	{}
	
	Eigen::Vector2d cellcentre_coord (int i, int j) const
	{
		Eigen::Vector2d cc;
		cc(0) = x0 + double(j - numGC) * dx + 0.5 * dx;
		cc(1) = y0 + double(i - numGC) * dy + 0.5 * dy;
		return cc;
	}
	
	Eigen::Vector2d cellBL_coord (int i, int j) const
	{
		Eigen::Vector2d cc;
		cc(0) = x0 + double(j - numGC) * dx;
		cc(1) = y0 + double(i - numGC) * dy;
		return cc;
	}
	
	std::pair<int, int> cellindex (const Eigen::Vector2d& pos) const
	{
		int i, j;
		i = static_cast<int>(floor((pos(1) - y0) / dy) + numGC);
		j = static_cast<int>(floor((pos(0) - x0) / dx) + numGC);
		return std::make_pair(i, j);
	}
};

#endif
