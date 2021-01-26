/**
 * rectangle2d.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#ifndef RECTAGLE2D_HPP
#define RECTAGLE2D_HPP

#include "vec2d.hpp"

namespace dragon8
{

struct rectangle2d
{
	vec2d c0, c1;

	rectangle2d(const vec2d c0, const vec2d c1) : c0(c0), c1(c1)
	{
	}

	rectangle2d(const double x0, const double y0, const double x1, const double y1)
		: c0(vec2d(x0, y0)), c1(vec2d(x1, y1))
	{
	}

	double get_width() const
	{
		return c1.x - c0.x;
	}

	double get_height() const
	{
		return c1.y - c0.y;
	}

};

}

#endif // RECTAGLE2D_HPP