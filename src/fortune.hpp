/**
 * fortune.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/31
 */
#ifndef FORTUNE_HPP
#define FORTUNE_HPP

#include <queue>
#include <set>
#include <vector>

#include "shape.hpp"
#include "vec2d.hpp"

namespace dragon8
{

struct LineSegment
{
	vec2d v0, v1;

	LineSegment(const vec2d v0, const vec2d v1);

	vec2d get_normal() const;
};

template <typename T>
using vec = std::vector<T>;

class Fortune
{
private:
	struct BeachSec
	{
		bool is_point;
		double x;
		const double *t;
		int32_t i, li, ri;
		vec2d p, l, r;

		BeachSec() = default;
		BeachSec(const double x);
		BeachSec(const double *t, const int32_t i, const int32_t li, const int32_t ri, const vec2d p, const vec2d l, const vec2d r);

		bool operator==(const BeachSec &other) const;

		double lx() const;
		double rx() const;
		double midx() const;
		bool operator<(const BeachSec &other) const;

		bool will_be_removed() const;
		double remove_time() const;

		double y_at(const double x) const;
	};

	enum EventType
	{
		NEW_POINT, REMOVE_SEC, BORDER_HIT
	};

	struct Event
	{
		EventType type;
		double time;
		uint32_t pi;
		BeachSec sec;
		vec2d hit_point;

		Event(const double time, const BeachSec sec, const vec2d hit_point);
		Event(const double time, const BeachSec sec);
		Event(const double time, const uint32_t pi);

		bool operator<(const Event &other) const;
	};

	ShapePtr container;
	vec<vec2d> points;
	vec<LineSegment> segs;

	double t;
	std::priority_queue<Event> events;
	std::set<BeachSec> beach;

	void submit_point(const vec2d p, const double d2, const bool border);
	void submit_point_if_border(const vec2d p, const double d2);
	void check_beach_integrity() const;
	void print_beach() const;
	void add_border_line(const BeachSec& left);
	void resolve_top();

	void add_firsts();

public:
	vec2d furthest;
	double dist2 = -1;
	std::vector<vec2d> pots;

	Fortune(ShapePtr container, const vec<vec2d>& points, const vec<LineSegment>& segs);
};

}

#endif // FORTUNE_HPP