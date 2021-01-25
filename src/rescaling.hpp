/**
 * rescaling.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#ifndef RESCALING_HPP
#define RESCALING_HPP

#include "vec2d.hpp"

namespace dragon8
{

class Rescaling
{
public:
	vec2d origin;
	double sx, sy;

	Rescaling(const vec2d origin, const double sx, const double sy);

	vec2d map(const vec2d v) const;
};

}

#endif // RESCALING_HPP