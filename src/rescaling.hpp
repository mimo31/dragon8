/**
 * rescaling.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#ifndef RESCALING_HPP
#define RESCALING_HPP

#include "rectangle2d.hpp"
#include "vec2.hpp"

namespace dragon8
{

class Rescaling
{
public:
	double sx, sy;
	vec2d origin;

	Rescaling() = default;

	Rescaling(const double sx, const double sy, const vec2d origin);
	Rescaling(const rectangle2d& mapof, const rectangle2d& mapto);

	vec2d map(const vec2d v) const;
};

}

#endif // RESCALING_HPP