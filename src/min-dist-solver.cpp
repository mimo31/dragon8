/**
 * min-dist-solver.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/26
 */
#include "min-dist-solver.hpp"

#include <iostream>
#include <limits>
#include <queue>
#include <set>

#include "evaluator.hpp"
#include "fortune.hpp"
#include "graphics.hpp"

template <typename T>
using vec = std::vector<T>;

using std::cout;
using std::endl;

namespace dragon8
{

MinDistSolver::MinDistSolver(const std::shared_ptr<Shape> shape, const uint32_t n)
	: Solver(shape, n)
{
}

void MinDistSolver::shift_solve(PointsState& ps, const uint32_t iters, RGen& rgen) const
{
	double bestscore = MinDistEvaluator().evalv(ps);
	PointsState beststate = ps;
	for (uint32_t it = 0; it < iters; it++)
	{
		/*double leastd2 = std::numeric_limits<double>::max();
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
			std::swap(li, lj);*/
		
		std::uniform_int_distribution<> dist(0, n - 1);
		const uint32_t lj = dist(rgen);
		double mdist2 = std::numeric_limits<double>::max();
		uint32_t li;
		for (uint32_t i = 0; i < n; i++)
		{
			if (i == lj)
				continue;
			const double dist2 = (ps[lj] - ps[i]).len2();
			if (dist2 < mdist2)
			{
				mdist2 = dist2;
				li = i;
			}
		}

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
		/*if (limiting != -1)
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
		const double curscore = MinDistEvaluator().evalv(ps);
		if (curscore > bestscore)
		{
			bestscore = curscore;
			beststate = ps;
		}*/
	}
	//ps = beststate;
}

void MinDistSolver::jump_solve(PointsState& ps, const uint32_t iters, RGen& rgen) const
{
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
		while (it < iters)
		{
			if (dist(rgen))
				std::swap(li, lj);
			const vec2d np = shape->gen_point(rgen);
			bool better = true;
			for (uint32_t i = 0; i < n; i++)
			{
				if (i == li)
					continue;
				if ((np - ps[i]).len2() <= leastd2)
				{
					better = false;
					break;
				}
			}
			if (better)
			{
				ps[li] = np;
				break;
			}
			it++;
		}
	}
}

bool MinDistSolver::fit_enough(const vec2d p0, const double d) const
{
	const rectangle2d box = shape->get_box();
	const double h = sqrt(3) / 2 * d;
	uint32_t cou = 0;
	int32_t yi0 = 0;
	while ((yi0 - 1) * h + p0.y >= box.c0.y)
		yi0--;
	int32_t xi0 = 0;
	while ((xi0 - 1) * d + p0.x >= box.c0.x)
		xi0--;
	xi0--;
	for (int32_t yi = yi0; yi * h + p0.y <= box.c1.y; yi++)
	{
		const double sh = yi % 2 ? 0 : d / 2;
		for (int32_t xi = xi0; xi * d + sh + p0.x <= box.c1.x; xi++)
		{
			const vec2d p(xi * d + sh + p0.x, yi * h + p0.y);
			if (shape->is_inside(p))
			{
				cou++;
				if (cou >= n)
				{
					return true;
				}
			}
		}
	}
	return false;
}

PointsState MinDistSolver::get_fit(const vec2d p0, const double d) const
{
	const rectangle2d box = shape->get_box();
	const double h = sqrt(3) / 2 * d;
	PointsState points;
	int32_t yi0 = 0;
	while ((yi0 - 1) * h + p0.y >= box.c0.y)
		yi0--;
	int32_t xi0 = 0;
	while ((xi0 - 1) * d + p0.x >= box.c0.x)
		xi0--;
	xi0--;
	for (int32_t yi = yi0; yi * h + p0.y <= box.c1.y; yi++)
	{
		const double sh = yi % 2 ? 0 : d / 2;
		for (int32_t xi = xi0; xi * d + sh + p0.x <= box.c1.x; xi++)
		{
			const vec2d p(xi * d + sh + p0.x, yi * h + p0.y);
			if (shape->is_inside(p))
			{
				points.push_back(p);
				if (points.size() >= n)
					return points;
			}
		}
	}
	return points;
}

PointsState MinDistSolver::regular_fit(const vec2d p0, double& score) const
{
	double d = shape->get_box().get_width();
	double dm, dp;
	if (fit_enough(p0, d))
	{
		while (fit_enough(p0, d * 2))
			d *= 2;
		dp = d * 2;
		dm = d;
	}
	else
	{
		while (!fit_enough(p0, d / 2))
			d /= 2;
		dp = d;
		dm = d / 2;
	}
	for (uint32_t i = 0; i < 50; i++)
	{
		const double mid = (dm + dp) / 2;
		if (fit_enough(p0, mid))
			dm = mid;
		else
			dp = mid;
	}
	score = dm;
	return get_fit(p0, dm);
}

PointsState MinDistSolver::search_regular_fits(const uint32_t iters, RGen& rgen) const
{
	PointsState best;
	double score = 0;
	for (uint32_t i = 0; i < iters; i++)
	{
		const vec2d p0 = shape->gen_point(rgen);
		double sc;
		const PointsState state = regular_fit(p0, sc);
		if (sc > score)
		{
			score = sc;
			best = state;
		}
	}
	return best;
}

PointsState MinDistSolver::solve() const
{
	std::random_device rd;
	std::mt19937 rgen(rd());

	PointsState ps = shape->gen_state(n, rgen);
	//PointsState ps{ vec2d(.49, .1), vec2d(.3, .4), vec2d(.8, .35), vec2d(.5, .7) };
	//PointsState ps{ vec2d(.49, .1), vec2d(.3, .4), vec2d(.8, .35) };
	//PointsState ps{ vec2d(.49, .8), vec2d(.3, .4), vec2d(.8, .35) };
	//PointsState ps{ vec2d(.6, .38), /*vec2d(.49, .1), */vec2d(.3, .4), vec2d(.8, .35), vec2d(.5, .7) };
	const vec<LineSegment> unitsquare{
		LineSegment(vec2d(0, 0), vec2d(1, 0)),
		LineSegment(vec2d(1, 0), vec2d(1, 1)),
		LineSegment(vec2d(1, 1), vec2d(0, 1)),
		LineSegment(vec2d(0, 1), vec2d(0, 0))
	};

	constexpr uint32_t iters = 10000;

	std::uniform_int_distribution<> dist(0, n - 1);
	for (uint32_t i = 0; i < iters; i++)
	{
		const uint32_t movei = dist(rgen);

		vec<vec2d> rempoints;
		rempoints.reserve(n - 1);
		for (uint32_t j = 0; j < n; j++)
		{
			if (j != movei)
				rempoints.push_back(ps[j]);
		}

		Fortune frt(shape, rempoints, unitsquare);

		vec2d bpoi = frt.furthest;
		double bd2 = frt.dist2;

		const vec<vec2d> corners{ vec2d(0, 0), vec2d(0, 1), vec2d(1, 0), vec2d(1, 1) };
		for (const vec2d c : corners)
		{
			double mndist2 = -1;
			for (const vec2d p : rempoints)
			{
				const double dst2 = p.dist2(c);
				if (mndist2 == -1 || dst2 < mndist2)
					mndist2 = dst2;
			}
			if (mndist2 > bd2)
			{
				bpoi = c;
				bd2 = mndist2;
			}
		}

		ps[movei] = shape->bound(ps[movei], bpoi);
		//cout << "moving to " << frt.furthest << endl;

		if (i && i % 512 == 0)
			write_image(shape, ps, "fortune" + std::to_string(i) + ".png", 0, frt.pots);
	}

	return ps;

	/*

	PointsState ps = shape->gen_state(n, rgen);
	if (n < 200)
	{
		jump_solve(ps, 500'000, rgen);
		shift_solve(ps, 50'000, rgen);
	}
	else if (n < 1000)
	{
		jump_solve(ps, 500'000, rgen);
		shift_solve(ps, 1'000, rgen);
	}
	else if (n < 2000)
	{
		jump_solve(ps, 10'000, rgen);
		shift_solve(ps, 100, rgen);
	}
	else
	{
		jump_solve(ps, 100, rgen);
		shift_solve(ps, 10, rgen);
	}
	uint32_t fit_iters;
	if (n < 50)
	{
		fit_iters = 50'000;
	}
	else if (n < 200)
	{
		fit_iters = 7'000;
	}
	else if (n < 1000)
	{
		fit_iters = 1'000;
	}
	else
	{
		fit_iters = 100;
	}

	PointsState rfit = search_regular_fits(fit_iters, rgen);
	const double fit_score = MinDistEvaluator().evalv(rfit);
	std::cout << "pure fit: " << fit_score << std::endl;
	shift_solve(rfit, 2000, rgen);
	const double shift_score = MinDistEvaluator().evalv(rfit);
	std::cout << "after shift: " << shift_score << std::endl;
	
	return fit_score > MinDistEvaluator().evalv(ps) ? rfit : ps;*/
}

}