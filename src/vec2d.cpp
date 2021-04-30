/**
 * vec2d.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "vec2.hpp"

namespace dragon8
{

vec2d operator*(const double a, const vec2d& v)
{
	return v * a;
}

std::ostream& operator<<(std::ostream& os, const vec2d& v)
{
	os << "vec2d(" << v.x << ", " << v.y << ")";
	return os;
}

}