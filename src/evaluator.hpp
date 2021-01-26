/**
 * evaluator.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP

#include "points-state.hpp"

namespace dragon8
{

class Evaluator
{
public:
	virtual double eval(const vec2d *const ar, const uint32_t n) const = 0;
	double evalv(const PointsState&) const;
};

class DistInvEvaluator : public Evaluator
{
public:
	double eval(const vec2d *const ar, const uint32_t n) const override;
};

class MinDistEvaluator : public Evaluator
{
public:
	double eval(const vec2d *const ar, const uint32_t n) const override;
};

}

#endif // EVALUATOR_HPP