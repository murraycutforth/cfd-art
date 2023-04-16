/*
 *	DESCRIPTION:	Definition of bilinear interpolation on 2D grid.
 * 			Templated over grid element type.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		07/12/2017
 */

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "sim_info.hpp"
#include <Eigen/Dense>
#include <vector>




template <typename T>
inline T linear_interpolate
(
	const T& fL, 
	const T& fR, 
	const double alpha
)
{
	return fL + alpha * (fR - fL);
}




template <typename T>
inline T bilinear_interpolate
(
	const T& fBL, 
	const T& fBR, 
	const T& fTL, 
	const T& fTR, 
	const double alpha, 
	const double beta
)
{
	T fB = linear_interpolate<T>(fBL, fBR, alpha);
	T fT = linear_interpolate<T>(fTL, fTR, alpha);
	return linear_interpolate<T>(fB, fT, beta);
}




template <typename T>
inline T grid_bilinear_interpolate
(
	const sim_info& params, 
	const std::vector<std::vector<T>>& grid, 
	const Eigen::Vector2d& pos
)
{
	Eigen::Vector2d displaced_pos;
	
	displaced_pos(0) = pos(0) + 0.5 * params.dx;
	displaced_pos(1) = pos(1) + 0.5 * params.dy;
	
	std::pair<int, int> ind = params.cellindex(displaced_pos);
	int i = ind.first;
	int j = ind.second;
	
	Eigen::Vector2d local_pos = pos - params.cellcentre_coord(i-1, j-1);
	double alpha = local_pos(0) / params.dx;
	double beta = local_pos(1) / params.dy;
	
	return bilinear_interpolate<T>(grid[i-1][j-1], grid[i-1][j], grid[i][j-1], grid[i][j], alpha, beta);
}

#endif
