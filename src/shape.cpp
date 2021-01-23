/**
 * shape.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "shape.hpp"

#include <cmath>

namespace dragon8
{

PointsState Shape::gen_state(const uint32_t n, RGen& rgen) const
{
	PointsState state;
	for (uint32_t i = 0; i < n; i++)
		state.push_back(gen_point(rgen));
	return state;
}

ShapeCircle::ShapeCircle(const double r)
	: r(r)
{
}

vec2d ShapeCircle::gen_point(RGen& rgen) const
{
	std::uniform_real_distribution<> dist(-r, r);
	double x, y;
	do
	{
		x = dist(rgen);
		y = dist(rgen);
	} while (x * x + y * y > r * r);
	return vec2d(x, y);
}

vec2d ShapeCircle::bound(const vec2d vfrom, const vec2d vto) const
{
	if (vto.len2() <= r * r)
		return vto;
	return vto.get_unit() * r;
	/*const vec2d vmove = vto - vfrom;
	const double vmlen2 = vmove.len2();
	const double b = vfrom.dot(vmove);
	const double t = sqrt(b * b + (r * r - vfrom.len2()) * vmlen2) / vmlen2;
	return vfrom + vmove * t;*/
}

}