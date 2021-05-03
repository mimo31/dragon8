/**
 * graphics.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#include "graphics.hpp"

#include <iostream>

#include "CImg.h"

#include "rescaling.hpp"

namespace dragon8
{

void write_image(const ShapePtr shape, const PointsState& state, const std::string& filename, const double circ_rad)
{
	using namespace cimg_library;

	constexpr int32_t w = 1024;
	constexpr int32_t h = 1024;
	constexpr int32_t xbrd = 16, ybrd = 16;
	CImg<unsigned char> img(w, h, 1, 3, 0);

	const rectangle2i draw_area = rectangle2i(xbrd, ybrd, w - 1 - xbrd, h - 1 - ybrd);

	const Rescaling resc = shape->draw(img, draw_area);

	constexpr int32_t pr = 4;

	const unsigned char pcol[] = { 0, 0, 255 };
	const unsigned char circcol[] = { 247, 255, 247 };
	const unsigned char circline[] = { 0, 191, 0 };

	for (const vec2d v : state)
	{
		const point p = resc.map(v).get_rounded();
		if (circ_rad != 0)
		{
			img.draw_circle(p.x, p.y, circ_rad * resc.sx, circcol, 1.0);
			img.draw_circle(p.x, p.y, circ_rad * resc.sx, circline, 1.0, -1);
		}
		img.draw_circle(p.x, p.y, pr, pcol, 1.0);
	}

	if (circ_rad != 0)
	{
		CImg<unsigned char> mask(w, h, 1, 3, 0);
		shape->draw(mask, draw_area);

		for (uint32_t y = 0; y < h; y++)
		{
			for (uint32_t x = 0; x < w; x++)
			{
				if (mask(x, y) == 0)
				{
					img(x, y, 0) = img(x, y, 1) = img(x, y, 2) = 0;
				}
			}
		}
	}

	img.save(filename.c_str());
}

void write_image(const ShapePtr shape, const PointsState& state, const std::string& filename, const double circ_rad, const std::vector<vec2d>& pots, const vec<std::pair<vec2d, vec2d>> &edges)
{
	using namespace cimg_library;

	constexpr int32_t w = 1024;
	constexpr int32_t h = 1024;
	constexpr int32_t xbrd = 16, ybrd = 16;
	CImg<unsigned char> img(w, h, 1, 3, 0);

	const rectangle2i draw_area = rectangle2i(xbrd, ybrd, w - 1 - xbrd, h - 1 - ybrd);

	const Rescaling resc = shape->draw(img, draw_area);

	constexpr int32_t pr = 4;

	const unsigned char pcol[] = { 0, 0, 255 };
	const unsigned char potcol[] = { 255, 0, 0 };
	const unsigned char edgecol[] = { 0, 255, 0 };
	const unsigned char circcol[] = { 247, 255, 247 };
	const unsigned char circline[] = { 0, 191, 0 };

	for (uint32_t y = ybrd; y < h - ybrd; y++)
	{
		for (uint32_t x = xbrd; x < w - xbrd; x++)
		{
			const vec2d p(x, y);
			double dist2;
			int32_t ind = -1;
			for (uint32_t i = 0; i < state.size(); i++)
			{
				const vec2d v = state[i];
				const double d2 = (resc.map(v) - p).len2();
				if (ind == -1 || d2 < dist2)
				{
					dist2 = d2;
					ind = i;
				}
			}
			unsigned char shade = 255 * ind / state.size();
			img(x, y, 0) = shade;
			img(x, y, 1) = shade;
			img(x, y, 2) = shade;
		}
	}

	for (const vec2d v : state)
	{
		const point p = resc.map(v).get_rounded();
		if (circ_rad != 0)
		{
			img.draw_circle(p.x, p.y, circ_rad * resc.sx, circcol, 1.0);
			img.draw_circle(p.x, p.y, circ_rad * resc.sx, circline, 1.0, -1);
		}
		img.draw_circle(p.x, p.y, pr, pcol, 1.0);
	}

	for (const vec2d v : pots)
	{
		const point p = resc.map(v).get_rounded();
		img.draw_circle(p.x, p.y, pr, potcol, 1.0);
	}

	for (const auto e : edges)
	{
		const point p0 = resc.map(e.first).get_rounded();
		const point p1 = resc.map(e.second).get_rounded();
		img.draw_line(p0.x, p0.y, p1.x, p1.y, edgecol, 1.0);
	}
	

	/*for (uint32_t i = 0; i < state.size(); i++)
	{
		for (uint32_t j = i + 1; j < state.size(); j++)
		{
			const vec2d mid = (state[i] + state[j]) / 2;
			const vec2d v = (state[i] - state[j]).get_lrot();
			const vec2d e0 = resc.map(mid + v);
			const vec2d e1 = resc.map(mid - v);
			const point E0 = e0.get_rounded(), E1 = e1.get_rounded();
			img.draw_line(E0.x, E0.y, E1.x, E1.y, pcol, 1.0, -1);
		}
	}*/

	if (circ_rad != 0)
	{
		CImg<unsigned char> mask(w, h, 1, 3, 0);
		shape->draw(mask, draw_area);

		for (uint32_t y = 0; y < h; y++)
		{
			for (uint32_t x = 0; x < w; x++)
			{
				if (mask(x, y) == 0)
				{
					img(x, y, 0) = img(x, y, 1) = img(x, y, 2) = 0;
				}
			}
		}
	}

	img.save(filename.c_str());
}

}