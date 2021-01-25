/**
 * rescaling.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#include "rescaling.hpp"

namespace dragon8
{

Rescaling::Rescaling(const vec2d origin, const double sx, const double sy)
	: origin(origin), sx(sx), sy(sy)
{
}

vec2d Rescaling::map(const vec2d v) const
{
	return vec2d(v.x * sx, v.y * sy) + origin;
}

}