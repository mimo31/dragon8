/**
 * dist-inv-solver.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "dist-inv-solver.hpp"

#include <algorithm>
#include <iostream>
#include <random>

#include "evaluator.hpp"

namespace dragon8
{

CompState::CompState(const ShapePtr shape, const uint32_t n, RGen& rgen)
	: ps(shape->gen_state(n, rgen)),
	lm(1.0),
	score(DistInvEvaluator().evalv(ps))
{
}

bool CompState::operator<(const CompState& other) const
{
	return score < other.score;
}

DistInvSolver::DistInvSolver(const std::shared_ptr<Shape> shape, const uint32_t n)
	: Solver(shape, n)
{
}

void DistInvSolver::iter(CompState& cs, vec2d *const fvec, vec2d *const tempstate) const
{
	std::fill_n(fvec, n, vec2d(0, 0));
	for (uint32_t j = 0; j < n; j++)
	{
		for (uint32_t k = 0; k < n; k++)
		{
			if (j == k)
				continue;
			const vec2d rv = cs.ps[j] - cs.ps[k];
			fvec[j] += rv.get_unit() / rv.len2();
		}
	}
	std::copy(cs.ps.begin(), cs.ps.end(), tempstate);
	for (uint32_t j = 0; j < n; j++)
		tempstate[j] = shape->bound(cs.ps[j], cs.ps[j] + fvec[j] * cs.lm);
	const double val = DistInvEvaluator().eval(tempstate, n);
	if (val < cs.score)
	{
		std::copy_n(tempstate, n, cs.ps.begin());
		cs.score = val;
	}
	else
	{
		cs.lm /= 2;
	}
}

CompState DistInvSolver::singleSolve(const uint32_t iters) const
{
	std::random_device rd;
	std::mt19937 rgen(rd());

	vec2d *fvec = new vec2d[n];
	vec2d *tempstate = new vec2d[n];

	CompState cs(shape, n, rgen);

	for (uint32_t i = 0; i < iters; i++)
		iter(cs, fvec, tempstate);

	delete[] fvec;
	delete[] tempstate;

	return cs;
}

CompState DistInvSolver::selectSolve(const uint32_t pool_size, const uint32_t init_iters, const uint32_t final_iters) const
{
	std::random_device rd;
	std::mt19937 rgen(rd());

	vec2d *fvec = new vec2d[n];
	vec2d *tempstate = new vec2d[n];

	std::vector<CompState> states;
	for (uint32_t i = 0; i < pool_size; i++)
		states.push_back(CompState(shape, n, rgen));

	for (uint32_t i = 0; i < pool_size; i++)
		for (uint32_t j = 0; j < init_iters; j++)
			iter(states[i], fvec, tempstate);

	uint32_t best = 0;
	for (uint32_t i = 0; i < pool_size; i++)
	{
		if (states[i] < states[best])
			best = i;		
	}

	for (uint32_t i = 0; i < final_iters; i++)
		iter(states[best], fvec, tempstate);

	delete[] fvec;
	delete[] tempstate;

	return states[best];
}

PointsState DistInvSolver::solve() const
{
	if (n < 100)
		return selectSolve(300, 500, 5000).ps;
	else if (n < 1000)
		return selectSolve(300, 50, 5000).ps;
	else if (n < 2000)
		return singleSolve(5000).ps;
	else
		return singleSolve(0).ps;
}

}