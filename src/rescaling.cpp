/**
 * rescaling.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#include "rescaling.hpp"

namespace dragon8
{

Rescaling::Rescaling(const double sx, const double sy, const vec2d origin)
	: sx(sx), sy(sy), origin(origin)
{
}

Rescaling::Rescaling(const rectangle2d& mapof, const rectangle2d& mapto)
	: sx((mapto.c1.x - mapto.c0.x) / (mapof.c1.x - mapof.c0.x)),
	sy((mapto.c1.y - mapto.c0.y) / (mapof.c1.y - mapof.c0.y)),
	origin(mapto.c0.x - mapof.c0.x * sx, mapto.c0.y - mapof.c0.y * sy)
{
}

vec2d Rescaling::map(const vec2d v) const
{
	return vec2d(v.x * sx, v.y * sy) + origin;
}

}