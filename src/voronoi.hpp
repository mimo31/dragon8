/**
 * voronoi.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/04/26
 */
#ifndef VORONOI_HPP
#define VORONOI_HPP

#include <queue>
#include <set>

#include "shape.hpp"
#include "vec.hpp"
#include "vec2.hpp"

namespace dragon8
{

template<typename float_t>
struct FurthestPointData
{
	vec2<float_t> point;
	float_t dist;

	FurthestPointData(const vec2<float_t> point, const float_t dist) : point(point), dist(dist)
	{
	}
};

template<typename float_t>
class Voronoi
{
private:
	ShapePtr container;
	vec<vec2<float_t>> sites;

	float_t t;

	struct BeachLinePoint
	{
	private:
		const bool is_dummy;
		const float_t dummy_x;
		const uint32_t point_id;
		const uint32_t left_site, right_site;
		const float_t *const t;
		const vec<vec2<float_t>> *sites;

		float_t get_x() const;

	public:
		BeachLinePoint(float_t dummy_x);

		BeachLinePoint(uint32_t point_id, uint32_t left_site, uint32_t right_site, const float_t *t, const vec<vec2<float_t>> *sites);

		bool operator<(const BeachLinePoint &other) const;
	};

	struct BeachLinePointData
	{
	private:
		/// -1 iff no point on the left / on the right
		uint32_t left_id, right_id;
		const uint32_t left_site, right_site;
	public:
		BeachLinePointData(const uint32_t left_id, const uint32_t right_id, const uint32_t left_site, const uint32_t right_site)
			: left_id(left_id), right_id(right_id), left_site(left_site), right_site(right_site)
		{
		}
	};

	std::set<BeachLinePoint> beach_line;
	vec<BeachLinePointData> beach_line_points;

	struct FortuneEvent
	{
	private:
		const float_t time;
	public:
		const bool is_site;
		const uint32_t site_id;
		const uint32_t left_site, rem_site, right_site;

		FortuneEvent(float_t time, uint32_t site_id);
		FortuneEvent(float_t time, uint32_t left_site, uint32_t rem_site, uint32_t right_site);

		bool operator<(const FortuneEvent &other) const;
	};

	std::priority_queue<FortuneEvent> events;

	void potential_add_vertex_event(uint32_t left_id, uint32_t rem_id, uint32_t right_id);

	void submit_vertex(uint32_t site0_id, uint32_t site1_id, uint32_t site2_id);

public:
	Voronoi(ShapePtr container);
	void init_sites(const vec<vec2<float_t>> &sites);
	void init_compute();
	void add_site(vec2<float_t> p);
	void remove_site(uint32_t i);
	FurthestPointData<float_t> get_furthest() const;
};

}

#endif // VORONOI_HPP