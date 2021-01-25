/**
 * graphics.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#include "graphics.hpp"

#include "CImg.h"

#include "rescaling.hpp"

namespace dragon8
{

void write_image(const ShapePtr shape, const PointsState& state, const std::string& filename)
{
	using namespace cimg_library;

	constexpr int32_t w = 1920;
	constexpr int32_t h = 1080;
	CImg<unsigned char> img(w, h, 1, 3, 0);

	const rectangle2d box = shape->get_box();

	Rescaling resc = shape->draw(img, rectangle2i(0, 0, w - 1, h - 1));

	constexpr int32_t pr = 4;

	const unsigned char pcol[] = { 0, 0, 255 };

	for (const vec2d v : state)
	{
		const point p = resc.map(v).get_rounded();
		img.draw_circle(p.x, p.y, pr, pcol, 1.0);
	}	

	img.save(filename.c_str());
}

}