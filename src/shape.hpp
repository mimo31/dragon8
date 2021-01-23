/**
 * shape.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#ifndef SHAPE_HPP
#define SHAPE_HPP

#include <random>

#include "vec2d.hpp"

namespace dragon8
{

typedef std::mt19937 RGen;

class Shape
{
public:
	virtual vec2d gen_point(RGen& rgen) const = 0;
	virtual vec2d bound(const vec2d vfrom, const vec2d vto) const = 0;
	// TODO: add some draw method here
};

class ShapeCircle : public Shape
{
private:
	double r;

public:
	ShapeCircle(const double r);

	vec2d gen_point(RGen& rgen) const override;
	vec2d bound(const vec2d vfrom, const vec2d vto) const override;
};

class ShapePolygon : public Shape
{
private:
	std::vector<vec2d> verts;

public:
	ShapePolygon(const std::vector<vec2d>& verts);

	vec2d gen_point(RGen& rgen) const override;
	vec2d bound(const vec2d vfrom, const vec2d vto) const override;
};

}

#endif // SHAPE_HPP