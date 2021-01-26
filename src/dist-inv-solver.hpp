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

struct CompState
{
	PointsState ps;
	double lm;
	double score;

	CompState(const ShapePtr shape, const uint32_t n, RGen& rgen);

	bool operator<(const CompState& other) const;
};

class DistInvSolver : public Solver
{
private:
	CompState singleSolve(const uint32_t iters) const;
	CompState selectSolve(const uint32_t pool_size, const uint32_t init_iters, const uint32_t final_iters) const;
	void iter(CompState& cs, vec2d *const fvec, vec2d *const tempstate) const;

public:
	DistInvSolver(const std::shared_ptr<Shape> shape, const uint32_t n);
	
	PointsState solve() const override;
};

}

#endif // DIST_INV_SOLVER