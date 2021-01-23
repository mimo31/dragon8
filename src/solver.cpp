/**
 * solver.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "solver.hpp"

namespace dragon8
{

Solver::Solver(const std::shared_ptr<Shape> shape, const uint32_t n)
	: shape(shape), n(n)
{
}

}