/**
 * point.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#ifndef POINT_HPP
#define POINT_HPP

#include <cstdint>

namespace dragon8
{

struct point
{
	int32_t x, y;

	point(const int32_t x, const int32_t y) : x(x), y(y)
	{
	}
};

}

#endif // POINT_HPP