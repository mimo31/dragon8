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

namespace dragon8
{

void write_image(const ShapePtr shape, const PointsState& state, const std::string& filename);

}

#endif // GRAPHICS_HPP