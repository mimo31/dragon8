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

double DistInvEvaluator::eval(const PointsState& state) const
{
	double sm = 0;
	const uint32_t n = state.size();
	for (uint32_t i = 0; i < n; i++)
		for (uint32_t j = i + 1; j < n; j++)
			sm += 1 / (state[i] - state[j]).len();
	return sm;
}

double MinDistEvaluator::eval(const PointsState& state) const
{
	double mn = std::numeric_limits<double>::max();
	const uint32_t n = state.size();
	for (uint32_t i = 0; i < n; i++)
		for (uint32_t j = i + 1; j < n; j++)
			mn = std::min(mn, (state[i] - state[j]).len2());
	return sqrt(mn);
}

}