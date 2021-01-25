/**
 * dist-inv-solver.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "dist-inv-solver.hpp"

#include <iostream>
#include <random>

#include "evaluator.hpp"

namespace dragon8
{

DistInvSolver::DistInvSolver(const std::shared_ptr<Shape> shape, const uint32_t n)
	: Solver(shape, n)
{
}

PointsState DistInvSolver::solve() const
{
	std::random_device rd;
	std::mt19937 rgen(rd());

	PointsState state = shape->gen_state(n, rgen);

	std::cout << "starting: E = " << DistInvEvaluator().eval(state) << '\n';

	vec2d *fvec = new vec2d[n];

	constexpr uint32_t iters = 3000;
	constexpr double lm = .01;

	for (uint32_t i = 0; i < iters; i++)
	{
		std::fill_n(fvec, n, vec2d(0, 0));
		for (uint32_t j = 0; j < n; j++)
		{
			for (uint32_t k = 0; k < n; k++)
			{
				if (j == k)
					continue;
				const vec2d rv = state[j] - state[k];
				fvec[j] += rv.get_unit() / rv.len2();
			}
		}
		for (uint32_t j = 0; j < n; j++)
			state[j] = shape->bound(state[j], state[j] + fvec[j] * lm);
		std::cout << "after " << i << ": E = " << DistInvEvaluator().eval(state) << '\n';
	}

	delete[] fvec;

	return state;
}

}