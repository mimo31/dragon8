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
	void jump_solve(PointsState& ps, const uint32_t iters, RGen& rgen) const;
	void shift_solve(PointsState& ps, const uint32_t iters, RGen& rgen) const;

	bool fit_enough(const vec2d p0, const double d) const;
	PointsState get_fit(const vec2d p0, const double d) const;
	PointsState regular_fit(const vec2d p0, double& score) const;
	PointsState search_regular_fits(const uint32_t iters, RGen& rgen) const;

public:
	MinDistSolver(const std::shared_ptr<Shape> shape, const uint32_t n);

	PointsState solve() const override;
};

}

#endif // MIN_DIST_SOLVER