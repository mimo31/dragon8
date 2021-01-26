/**
 * min-dist-solver.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/26
 */
#include "min-dist-solver.hpp"

#include <iostream>
#include <limits>

#include "evaluator.hpp"
#include "graphics.hpp"

namespace dragon8
{

MinDistSolver::MinDistSolver(const std::shared_ptr<Shape> shape, const uint32_t n)
	: Solver(shape, n)
{
}

PointsState MinDistSolver::solve() const
{
	std::random_device rd;
	std::mt19937 rgen(rd());

	constexpr uint32_t iters = 500;

	PointsState ps = shape->gen_state(n, rgen);

	for (uint32_t it = 0; it < iters; it++)
	{
		double leastd2 = std::numeric_limits<double>::max();
		uint32_t li, lj;
		for (uint32_t i = 0; i < n; i++)
		{
			for (uint32_t j = i + 1; j < n; j++)
			{
				const double d2 = (ps[i] - ps[j]).len2();
				if (d2 < leastd2)
				{
					leastd2 = d2;
					li = i;
					lj = j;
				}
			}
		}
		std::uniform_int_distribution<> dist(0, 1);
		if (dist(rgen))
			std::swap(li, lj);
		double mxmove = std::numeric_limits<double>::max();
		int32_t limiting = -1;
		const vec2d s0 = ps[li], s1 = ps[lj], s = s1 - s0;
		for (uint32_t i = 0; i < n; i++)
		{
			const vec2d p = ps[i];
			if (s.dot(p - s1) <= 0)
				continue;
			const double moveavail = .5 * (p - s0).len2() / s.dot(p - s0) - 1;
			if (moveavail < mxmove)
			{
				mxmove = moveavail;
				limiting = i;
			}
		}
		if (limiting == -1)
			mxmove = .5;
		ps[lj] = shape->bound(s1, s1 + s * mxmove);
		if (limiting != -1)
		{
			mxmove = std::numeric_limits<double>::max();
			const vec2d sp0 = (ps[li] + ps[limiting]) / 2, sp1 = ps[lj], sp = sp1 - sp0;
			for (uint32_t i = 0; i < n; i++)
			{
				const vec2d p = ps[i];
				if (sp.dot(p - sp1) <= 0)
					continue;
				const double moveavail = ((p + s0) / 2 - sp1).dot(p - s0) / sp.dot(p - s0);
				mxmove = std::min(mxmove, moveavail);
			}
			if (mxmove == std::numeric_limits<double>::max())
				mxmove = .5;
			ps[lj] = shape->bound(sp1, sp1 + sp * mxmove);
		}
		std::cout << "score after it " << it << ": " << MinDistEvaluator().evalv(ps) << std::endl;
		//write_image(shape, ps, "it" + std::to_string(it) + ".png");
	}

	return ps;
}

}