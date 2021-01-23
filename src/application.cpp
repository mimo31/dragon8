/**
 * application.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "application.hpp"

#include <iostream>
#include <memory>

#include "dist-inv-solver.hpp"
#include "points-state.hpp"
#include "shape.hpp"

namespace dragon8
{

void Application::run()
{
	const std::shared_ptr<Shape> shape = std::make_shared<ShapeCircle>(1);
	const uint32_t n = 100;

	const PointsState ps = DistInvSolver(shape, n).solve();

	for (const vec2d v : ps)
	{
		std::cout << "(" << v.x << ", " << v.y << ")" << '\n';
	}
}

}