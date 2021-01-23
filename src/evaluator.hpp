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
	virtual double eval(const PointsState&) const = 0;
};

class DistInvEvaluator : public Evaluator
{
public:
	double eval(const PointsState&) const override;
};

class MinDistEvaluator : public Evaluator
{
public:
	double eval(const PointsState&) const override;
};

}

#endif // EVALUATOR_HPP