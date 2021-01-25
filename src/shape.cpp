/**
 * shape.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "shape.hpp"

#include <cmath>

namespace dragon8
{

PointsState Shape::gen_state(const uint32_t n, RGen& rgen) const
{
	PointsState state;
	for (uint32_t i = 0; i < n; i++)
		state.push_back(gen_point(rgen));
	return state;
}

ShapeCircle::ShapeCircle(const double r)
	: r(r)
{
}

vec2d ShapeCircle::gen_point(RGen& rgen) const
{
	std::uniform_real_distribution<> dist(-r, r);
	double x, y;
	do
	{
		x = dist(rgen);
		y = dist(rgen);
	} while (x * x + y * y > r * r);
	return vec2d(x, y);
}

vec2d ShapeCircle::bound(const vec2d vfrom, const vec2d vto) const
{
	if (vto.len2() <= r * r)
		return vto;
	return vto.get_unit() * r;
	/*const vec2d vmove = vto - vfrom;
	const double vmlen2 = vmove.len2();
	const double b = vfrom.dot(vmove);
	const double t = sqrt(b * b + (r * r - vfrom.len2()) * vmlen2) / vmlen2;
	return vfrom + vmove * t;*/
}

rectangle2d ShapeCircle::get_box() const
{
	return rectangle2d(-r, -r, r, r);
}

Rescaling ShapeCircle::draw(cimg_library::CImg<unsigned char>& img, const rectangle2i box) const
{
	const vec2d center = (vec2d(box.c0) + vec2d(box.c1)) / 2;
	const double rad = std::min(box.get_width(), box.get_height()) / 2;

	const point ccoors = center.get_rounded();
	const unsigned char white[] = { 255, 255, 255 };
	img.draw_circle(ccoors.x, ccoors.y, round(rad), white, 1.0);

	const double sc = rad / r;
	return Rescaling(center, sc, sc);
}

ShapePolygon::ShapePolygon(const std::vector<vec2d>& verts) : verts(verts)
{
}

vec2d ShapePolygon::gen_point(RGen& rgen) const
{
	// TODO: implement
	return vec2d();
}

vec2d ShapePolygon::bound(const vec2d vfrom, const vec2d vto) const
{
	// TODO: implement
	return vec2d();
}

rectangle2d ShapePolygon::get_box() const
{
	// TODO: implement
	return rectangle2d(0, 0, 0, 0);
}

Rescaling ShapePolygon::draw(cimg_library::CImg<unsigned char>& img, const rectangle2i box) const
{
	// TODO: implement
	return Rescaling(vec2d(), 0, 0);
}

}