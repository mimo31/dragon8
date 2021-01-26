/**
 * evaluator.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "evaluator.hpp"

#include <limits>

namespace dragon8
{

double Evaluator::evalv(const PointsState& state) const
{
	if (state.size() == 0)
		return eval(nullptr, 0);
	else
		return eval(&state[0], state.size());
}

double DistInvEvaluator::eval(const vec2d *const ar, const uint32_t n) const
{
	double sm = 0;
	for (uint32_t i = 0; i < n; i++)
		for (uint32_t j = i + 1; j < n; j++)
			sm += 1 / (ar[i] - ar[j]).len();
	return sm;
}

double MinDistEvaluator::eval(const vec2d *const ar, const uint32_t n) const
{
	double mn = std::numeric_limits<double>::max();
	for (uint32_t i = 0; i < n; i++)
		for (uint32_t j = i + 1; j < n; j++)
			mn = std::min(mn, (ar[i] - ar[j]).len2());
	return sqrt(mn);
}

}