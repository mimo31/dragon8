/**
 * min-dist-solver.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/26
 */
#ifndef MIN_DIST_SOLVER
#define MIN_DIST_SOLVER

#include <memory>

#include "points-state.hpp"
#include "solver.hpp"

namespace dragon8
{

class MinDistSolver : public Solver
{
private:
	void jumpSolve(PointsState& ps, const uint32_t iters, RGen& rgen) const;
	void shiftSolve(PointsState& ps, const uint32_t iters, RGen& rgen) const;

public:
	MinDistSolver(const std::shared_ptr<Shape> shape, const uint32_t n);

	PointsState solve() const override;
};

}

#endif // MIN_DIST_SOLVER