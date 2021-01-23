/**
 * dist-inv-solver.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#ifndef DIST_INV_SOLVER
#define DIST_INV_SOLVER

#include <cstdint>
#include <memory>

#include "points-state.hpp"
#include "solver.hpp"

namespace dragon8
{

class DistInvSolver : public Solver
{
public:
	DistInvSolver(const std::shared_ptr<Shape> shape, const uint32_t n);
	
	PointsState solve() const override;
};

}

#endif // DIST_INV_SOLVER