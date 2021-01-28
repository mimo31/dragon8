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

	const rectangle2i draw_area = rectangle2i(xbrd, ybrd, w - 1 - 2 * xbrd, h - 1 - 2 * ybrd);

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

}