/**
 * dist-inv-solver.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "dist-inv-solver.hpp"

namespace dragon8
{

DistInvSolver::DistInvSolver(const std::shared_ptr<Shape> shape, const uint32_t n)
	: Solver(shape, n)
{
}

PointsState DistInvSolver::solve() const
{
	// TODO implement
	return PointsState();
}

}