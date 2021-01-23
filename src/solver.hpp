/**
 * solver.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cstdint>
#include <memory>

#include "points-state.hpp"
#include "shape.hpp"

namespace dragon8
{

class Solver
{
private:
	uint32_t n;
	std::shared_ptr<Shape> shape;

protected:
	Solver(const std::shared_ptr<Shape> shape, const uint32_t n);

public:
	virtual PointsState solve() const = 0;
	
};

}


#endif // SOLVER_HPP