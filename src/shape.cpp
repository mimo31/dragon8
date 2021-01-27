/**
 * shape.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "shape.hpp"

#include <cmath>
#include <vector>

namespace dragon8
{

// for shape-drawing
const unsigned char white[] = { 255, 255, 255 };

PointsState Shape::gen_state(const uint32_t n, RGen& rgen) const
{
	PointsState state;
	for (uint32_t i = 0; i < n; i++)
		state.push_back(gen_point(rgen));
	return state;
}

bool intersects(const vec2d v0, const vec2d v1, const vec2d w0, const vec2d w1)
{
	const vec2d v = v1 - v0, w = w1 - v0;
	const vec2d w0v = w0 - v0, w1v = w1 - v0;
	if (v.cross(w) == 0 && v.cross(w0 - v0) == 0)
	{
		const double v1d = v.dot(v), w0d = v.dot(w0v), w1d = v.dot(w1v);
		const double wdfrom = std::min(w0d, w1d), wdto = std::max(w0d, w1d);
		return wdfrom <= v1d && 0 <= wdto;
	}
	const double w0c = v.cross(w0v), w1c = v.cross(w1v);
	if ((w0c > 0 && w1c > 0) || (w0c < 0 && w1c < 0))
		return false;
	const vec2d v0w = v0 - w0, v1w = v1 - w0;
	const double v0c = w.cross(v0w), v1c = v.cross(v1w);
	if ((v0c > 0 && v1c > 0) || (v0c < 0 && v1c < 0))
		return false;
	return true;
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
	img.draw_circle(ccoors.x, ccoors.y, round(rad), white, 1.0);

	const double sc = rad / r;
	return Rescaling(center, sc, sc);
}

ShapePolygon::ShapePolygon(const std::vector<vec2d>& verts) : verts(verts)
{
}

int32_t quad(const vec2d v)
{
	if (v.x > 0 && v.y >= 0)
		return 0;
	if (v.x <= 0 && v.y > 0)
		return 1;
	if (v.x < 0 && v.y <= 0)
		return 2;
	return 3;
}

bool ShapePolygon::is_inside(const vec2d p) const
{
	int32_t quadstravelled = 0;
	vec2d lastpos = verts[verts.size() - 1] - p;
	int32_t lastquad = quad(lastpos);
	for (uint32_t i = 0; i < verts.size(); i++)
	{
		const vec2d nv = verts[i] - p;
		if (nv.is_zero())
			return true;
		const int32_t nq = quad(nv);
		if (lastquad == nq)
		{
			lastpos = nv;
		}
		else if (nq == lastquad + 1 || (nq == 0 && lastquad == 3))
		{
			lastquad = nq;
			lastpos = nv;
			quadstravelled++;
		}
		else if (nq == lastquad - 1 || (nq == 3 && lastquad == 0))
		{
			lastquad = nq;
			lastpos = nv;
			quadstravelled--;
		}
		else
		{
			const vec2d vec = nv - lastpos;
			if (vec.get_lrot().dot(-lastpos) > 0)
				quadstravelled += 2;
			else
				quadstravelled -= 2;
			lastquad = nq;
			lastpos = nv;
		}
	}
	return (quadstravelled & 7) >= 2 && (quadstravelled & 7) <= 6;
	/*double mnx = 0, mny = 0;
	for (const vec2d v : verts)
	{
		mnx = std::min(v.x, mnx);
		mny = std::min(v.y, mny);
	}
	const vec2d b(mnx - 1, mny - 1);
	bool inside = intersects(b, p, verts[0], verts[verts.size() - 1]);
	for (uint32_t i = 0; i < verts.size() - 1; i++)
		inside ^= intersects(b, p, verts[i], verts[i + 1]);
	return inside;*/
}

vec2d ShapePolygon::gen_point(RGen& rgen) const
{
	const rectangle2d box = get_box();
	std::uniform_real_distribution<> xdist(box.c0.x, box.c1.x);
	std::uniform_real_distribution<> ydist(box.c0.y, box.c1.y);
	vec2d p;
	do
	{
		p.x = xdist(rgen);
		p.y = ydist(rgen);
	} while (!is_inside(p));
	return p;
}

double line_segment_dist2(const vec2d v, const vec2d s0, const vec2d s1)
{
	const vec2d vs = v - s0;
	const vec2d s = s1 - s0;
	if (s.dot(vs) < 0)
		return vs.len2();
	else if (s.dot(vs) > s.dot(s))
		return (vs - s).len2();
	else
		return s.cross(vs) * s.cross(vs) / s.len2();
}

vec2d ShapePolygon::bound(const vec2d vfrom, const vec2d vto) const
{
	if (is_inside(vto))
		return vto;
	double leastdist2 = line_segment_dist2(vto, verts[0], verts[verts.size() - 1]);
	uint32_t besti = verts.size() - 1;
	for (uint32_t i = 0; i < verts.size() - 1; i++)
	{
		const double dist2 = line_segment_dist2(vto, verts[i], verts[i + 1]);
		if (dist2 < leastdist2)
		{
			leastdist2 = dist2;
			besti = i;
		}
	}
	const vec2d s0 = verts[besti], s1 = verts[besti != verts.size() - 1 ? besti + 1 : 0];
	const vec2d vs = vto - s0;
	const vec2d s = s1 - s0;
	if (s.dot(vs) < 0)
		return s0;
	else if (s.dot(vs) > s.dot(s))
		return s1;
	else
		return s0 + s * s.dot(vs) / s.len2();
}

rectangle2d ShapePolygon::get_box() const
{
	double mnx = verts[0].x, mny = verts[0].y, mxx = verts[0].x, mxy = verts[0].y;
	for (uint32_t i = 1; i < verts.size(); i++)
	{
		const double x = verts[i].x, y = verts[i].y;
		mnx = std::min(mnx, x);
		mny = std::min(mny, y);
		mxx = std::max(mxx, x);
		mxy = std::max(mxy, y);
	}
	return rectangle2d(mnx, mny, mxx, mxy);
}

Rescaling ShapePolygon::draw(cimg_library::CImg<unsigned char>& img, const rectangle2i box) const
{
	const rectangle2d real_box = get_box();
	Rescaling resc;
	if (real_box.get_width() * box.get_height() < real_box.get_height() * box.get_width())
	{
		const double x0at = (box.get_width() - real_box.get_width() / real_box.get_height() * box.get_height()) / 2;
		const double x1at = (box.get_width() + real_box.get_width() / real_box.get_height() * box.get_height()) / 2;
		const double orx = x0at - real_box.c0.x / real_box.get_width() * (x1at - x0at);
		const double ory = -real_box.c0.y / real_box.get_height() * box.get_height();
		const double sx = (x1at - x0at) / real_box.get_width();
		const double sy = box.get_height() / real_box.get_height();
		resc = Rescaling(vec2d(orx, ory), sx, sy);
	}
	else
	{
		const double orx = -real_box.c0.x / real_box.get_width() * box.get_width();
		const double y0at = (box.get_height() - real_box.get_height() / real_box.get_width() * box.get_width()) / 2;
		const double y1at = (box.get_height() + real_box.get_height() / real_box.get_width() * box.get_width()) / 2;
		const double ory = y0at - real_box.c0.y / real_box.get_height() * (y1at - y0at);
		const double sx = (y1at - y0at) / real_box.get_height();
		const double sy = box.get_width() / real_box.get_width();
		resc = Rescaling(vec2d(orx, ory), sx, sy);
	}
	cimg_library::CImg<int> points(verts.size(), 2);
	for (uint32_t i = 0; i < verts.size(); i++)
	{
		const point p = resc.map(verts[i]).get_rounded();
		points(i, 0) = p.x;
		points(i, 1) = p.y;
	}
	img.draw_polygon(points, white);
	return resc;
}

}