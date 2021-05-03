/**
 * graphics.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/24
 */
#ifndef GRAPHICS_HPP
#define GRAPHICS_HPP

#include "points-state.hpp"
#include "shape.hpp"
#include "vec.hpp"

namespace dragon8
{

void write_image(const ShapePtr shape, const PointsState& state, const std::string& filename, const double circ_rad);
void write_image(const ShapePtr shape, const PointsState& state, const std::string& filename, const double circ_rad, const vec<vec2d>& pots, const vec<std::pair<vec2d, vec2d>> &edges);

}

#endif // GRAPHICS_HPP