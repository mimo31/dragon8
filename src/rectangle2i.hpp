/**
 * rectangle2i.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#ifndef RECTANGLE2I
#define RECTANGLE2I

#include "point.hpp"

namespace dragon8
{

struct rectangle2i
{
	point c0, c1;

	rectangle2i(const point c0, const point c1) : c0(c0), c1(c1)
	{
	}

	rectangle2i(const int32_t x0, const int32_t y0, const int32_t x1, const int32_t y1)
		: c0(point(x0, y0)), c1(point(x1, y1))
	{
	}

	int32_t get_width() const
	{
		return c1.x - c0.x;
	}

	int32_t get_height() const
	{
		return c1.y - c0.y;
	}
};

}

#endif // RECTANGLE2I