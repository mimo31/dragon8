/**
 * vec2d.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#ifndef VEC2D_HPP
#define VEC2D_HPP

#include <cmath>
#include <ostream>

#include "point.hpp"

namespace dragon8
{

/**
 * Struct containing two doubles representing a 2D vector
 */
struct vec2d
{
	double x;
	double y;

	vec2d() = default;

	constexpr vec2d(const double x, const double y) : x(x), y(y)
	{
	}

	vec2d(const point p) : x(p.x), y(p.y)
	{
	}

	double len2() const
	{
		return x * x + y * y;
	}

	double len() const
	{
		return sqrt(x * x + y * y);
	}

	vec2d get_unit() const
	{
		return *this / this->len();
	}

	vec2d get_lrot() const
	{
		return vec2d(-y, x);
	}

	bool operator==(const vec2d &other) const
	{
		return x == other.x && y == other.y;
	}

	vec2d &operator+=(const vec2d &other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}

	vec2d &operator-=(const vec2d &other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}

	vec2d operator+(const vec2d &other) const
	{
		return vec2d(x + other.x, y + other.y);
	}

	vec2d operator-(const vec2d &other) const
	{
		return vec2d(x - other.x, y - other.y);
	}

	vec2d operator*(const double a) const
	{
		return vec2d(a * x, a * y);
	}

	vec2d operator/(const double a) const
	{
		return vec2d(x / a, y / a);
	}

	vec2d operator-() const
	{
		return vec2d(-x, -y);
	}

	double dist2(const vec2d &other) const
	{
		return (*this - other).len2();
	}

	double dist(const vec2d &other) const
	{
		return (*this - other).len();
	}

	vec2d mid(const vec2d &other) const
	{
		return (*this + other) / 2;
	}

	vec2d operator^(const vec2d &other) const
	{
		return mid(other);
	}

	double dot(const vec2d &other) const
	{
		return x * other.x + y * other.y;
	}

	double cross(const vec2d &other) const
	{
		return x * other.y - y * other.x;
	}

	bool is_zero() const
	{
		return x == 0 && y == 0;
	}

	bool inside(const double x0, const double y0, const double x1, const double y1) const
	{
		return x0 <= x && x <= x1 && y0 <= y && y <= y1;
	}

	point get_rounded() const
	{
		return point(round(x), round(y));
	}
};

vec2d operator*(const double a, const vec2d& v);

std::ostream& operator<<(std::ostream& os, const vec2d& v);

}

#endif // VEC2D_HPP
