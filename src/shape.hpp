/**
 * shape.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#ifndef SHAPE_HPP
#define SHAPE_HPP

#include <memory>
#include <random>

#include "CImg.h"

#include "points-state.hpp"
#include "rectangle2d.hpp"
#include "rectangle2i.hpp"
#include "rescaling.hpp"
#include "vec2d.hpp"

#include "test/test.hpp"

namespace dragon8
{

bool intersects(const vec2d v0, const vec2d v1, const vec2d w0, const vec2d w1);

typedef std::mt19937 RGen;

class Shape
{
public:
	virtual bool is_inside(const vec2d p) const = 0;
	virtual vec2d gen_point(RGen& rgen) const = 0;
	PointsState gen_state(const uint32_t n, RGen& rgen) const;
	virtual vec2d bound(const vec2d vfrom, const vec2d vto) const = 0;
	virtual rectangle2d get_box() const = 0;
	virtual Rescaling draw(cimg_library::CImg<unsigned char>& img, const rectangle2i box) const = 0; // a transformation of internal coordinates to image coordinates
};

typedef std::shared_ptr<Shape> ShapePtr;

class ShapeCircle : public Shape
{
private:
	double r;

public:
	ShapeCircle(const double r);

	bool is_inside(const vec2d p) const override;
	vec2d gen_point(RGen& rgen) const override;
	vec2d bound(const vec2d vfrom, const vec2d vto) const override;
	rectangle2d get_box() const override;
	Rescaling draw(cimg_library::CImg<unsigned char>& img, const rectangle2i box) const override;
};

class ShapePolygon : public Shape
{
private:
	std::vector<vec2d> verts;

public:
	ShapePolygon(const std::vector<vec2d>& verts);

	bool is_inside(const vec2d p) const override;
	vec2d gen_point(RGen& rgen) const override;
	vec2d bound(const vec2d vfrom, const vec2d vto) const override;
	rectangle2d get_box() const override;
	Rescaling draw(cimg_library::CImg<unsigned char>& img, const rectangle2i box) const override;
};

}

#endif // SHAPE_HPP