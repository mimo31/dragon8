/**
 * vec2.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/04/26
 */
#ifndef VEC2_HPP
#define VEC2_HPP

#include <cmath>
#include <ostream>

#include "point.hpp"

namespace dragon8
{

/**
 * Struct containing two doubles representing a 2D vector
 */
template<typename float_t>
struct vec2
{
	float_t x;
	float_t y;

	vec2() = default;

	constexpr vec2(const float_t x, const float_t y) : x(x), y(y)
	{
	}

	vec2(const point p) : x(p.x), y(p.y)
	{
	}

	float_t len2() const
	{
		return x * x + y * y;
	}

	float_t len() const
	{
		return sqrt(x * x + y * y);
	}

	vec2 get_unit() const
	{
		return *this / this->len();
	}

	vec2 get_lrot() const
	{
		return vec2(-y, x);
	}

	bool operator==(const vec2 &other) const
	{
		return x == other.x && y == other.y;
	}

	vec2 &operator+=(const vec2 &other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}

	vec2 &operator-=(const vec2 &other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}

	vec2 operator+(const vec2 &other) const
	{
		return vec2d(x + other.x, y + other.y);
	}

	vec2 operator-(const vec2 &other) const
	{
		return vec2d(x - other.x, y - other.y);
	}

	vec2 operator*(const float_t a) const
	{
		return vec2d(a * x, a * y);
	}

	vec2 operator/(const float_t a) const
	{
		return vec2d(x / a, y / a);
	}

	vec2 operator-() const
	{
		return vec2d(-x, -y);
	}

	float_t dist2(const vec2 &other) const
	{
		return (*this - other).len2();
	}

	float_t dist(const vec2 &other) const
	{
		return (*this - other).len();
	}

	vec2 mid(const vec2 &other) const
	{
		return (*this + other) / 2;
	}

	vec2 operator^(const vec2 &other) const
	{
		return mid(other);
	}

	float_t dot(const vec2 &other) const
	{
		return x * other.x + y * other.y;
	}

	float_t cross(const vec2 &other) const
	{
		return x * other.y - y * other.x;
	}

	bool is_zero() const
	{
		return x == 0 && y == 0;
	}

	bool inside(const float_t x0, const float_t y0, const float_t x1, const float_t y1) const
	{
		return x0 <= x && x <= x1 && y0 <= y && y <= y1;
	}

	point get_rounded() const
	{
		return point(round(x), round(y));
	}
};

template<typename float_t>
vec2<float_t> operator*(const float_t a, const vec2<float_t>& v);

template<typename float_t>
std::ostream& operator<<(std::ostream& os, const vec2<float_t>& v);

}

#endif // VEC2_HPP