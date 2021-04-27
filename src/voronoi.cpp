/**
 * voronoi.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/04/26
 */
#include "voronoi.hpp"

#include <cassert>

namespace dragon8
{

template<typename float_t>
Voronoi<float_t>::Voronoi(const ShapePtr container)
	: container(container)
{
}

template<typename float_t>
void Voronoi<float_t>::init_sites(const vec<vec2<float_t>> &sites)
{
	sites = sites;
}

template<typename float_t>
void Voronoi<float_t>::init_compute()
{
	for (uint32_t i = 0; i < sites.size(); i++)
		events.push(FortuneEvent(sites[i].y, i));
	const uint32_t site0 = events.top().site_id;
	events.pop();
	while (!events.empty())
	{
		const FortuneEvent event = events.top();
		events.pop();

		if (event.is_site)
		{
			const uint32_t site_id = event.site_id;
			if (!beach_line.empty())
			{
				const float_t dummy_x = sites[site_id].x;
				const BeachLinePoint dummy(dummy_x);
				auto lb = beach_line.lower_bound(dummy);
				const uint32_t next_id = beach_line_points.size();
				if (lb != beach_line.end())
				{
					const BeachLinePoint put_before = *lb;
					const uint32_t put_into_site = put_before.left_site;
					const BeachLinePoint p0(next_id, put_into_site, site_id, &t, &sites);
					const BeachLinePoint p1(next_id + 1, site_id, put_into_site, &t, &sites);
					beach_line.insert(p0);
					beach_line.insert(p1);
					const uint32_t prev_point_id = beach_line_points[put_before.point_id].left_id;
					beach_line_points.push_back(BeachLinePointData(prev_point_id, next_id + 1, put_into_site, site_id));
					beach_line_points.push_back(BeachLinePointData(next_id, put_before.point_id, site_id, put_into_site));
					if (prev_point_id != -1)
					{
						beach_line_points[prev_point_id].right_id = next_id;
						potential_add_vertex_event(beach_line_points[prev_point_id].left_site, put_into_site, site_id);
					}
					beach_line_points[put_before.point_id].left_id = next_id + 1;
					potential_add_vertex_event(site_id, put_into_site, put_before.right_site);
				}
				else
				{
					const BeachLinePoint put_after = *std::prev(beach_line.end());
					const uint32_t put_into_site = put_after.right_site;
					const BeachLinePoint p0(next_id, put_into_site, site_id, &t, &sites);
					const BeachLinePoint p1(next_id + 1, site_id, put_into_site, &t, &sites);
					beach_line.insert(p0);
					beach_line.insert(p1);
					beach_line_points.push_back(BeachLinePointData(put_after.point_id, next_id + 1, put_into_site, site_id));
					const uint32_t next_point_id = beach_line_points[put_after.point_id].right_id;
					beach_line_points.push_back(BeachLinePointData(next_id, next_point_id, site_id, put_into_site));
					if (next_point_id != -1)
					{
						beach_line_points[next_point_id].left_id = next_id + 1;
						potential_add_vertex_event(site_id, put_into_site, beach_line_points[next_point_id].right_site);
					}
					beach_line_points[put_after.point_id].right_id = next_id;
					potential_add_vertex_event(put_after.left_site, put_into_site, site_id);
				}
			}
			else
			{
				const BeachLinePoint p0(0, site0, site_id, &t, &sites);
				const BeachLinePoint p1(1, site_id, site0, &t, &sites);
				beach_line.insert(p0);
				beach_line.insert(p1);
				beach_line_points.push_back(BeachLinePointData(-1, 1, site0, site_id));
				beach_line_points.push_back(BeachLinePointData(0, -1, site_id, site0));
			}
		}
		else
		{
			const BeachLinePoint left_point(-1, event.left_site, event.rem_site, &t, &sites);
			auto it = beach_line.find(left_point);
			if (it == beach_line.end())
				continue;
			const uint32_t left_point_id = it->point_id;
			const uint32_t right_point_id = beach_line_points[left_point_id].right_id;
			assert(right_point_id != -1);
			if (beach_line_points[right_point_id].right_site != event.right_site)
				continue;
			submit_vertex(event.left_site, event.rem_site, event.right_site);
			const uint32_t new_point_id = beach_line_points.size();
			const uint32_t ll_point_id = beach_line_points[left_point_id].left_id;
			const uint32_t rr_point_id = beach_line_points[right_point_id].right_id;
			beach_line_points.push_back(BeachLinePointData(ll_point_id, rr_point_id, event.left_site, event.right_site));
			if (ll_point_id != -1)
				beach_line_points[ll_point_id].right_id = new_point_id;
			if (rr_point_id != -1)
				beach_line_points[rr_point_id].left_id = new_point_id;
			
			beach_line.insert(BeachLinePoint(new_point_id, event.left_site, event.right_site, &t, &sites));
			beach_line.erase(BeachLinePoint(-1, event.left_site, event.rem_site, &t, &sites));
			beach_line.erase(BeachLinePoint(-1, event.rem_site, event.right_site, &t, &sites));
		}
	}
}

template<typename float_t>
vec2<float_t> circumc(const vec2<float_t> a, const vec2<float_t> b, const vec2<float_t> c)
{
	const vec2<float_t> v1 = a - b, v2 = c - b;
	const double s1 = v1.dot(a ^ b), s2 = v2.dot(c ^ b);
	const double det = v1.x * v2.y - v1.y * v2.x;
	return vec2<float_t>(
		(s1 * v2.y - v1.y * s2) / det,
		(v1.x * s2 - s1 * v2.x) / det
	);
}

template<typename float_t>
void Voronoi<float_t>::potential_add_vertex_event(const uint32_t left_id, const uint32_t rem_id, const uint32_t right_id)
{
	const vec2<float_t> left_site = sites[left_id], rem_site = sites[rem_id], right_site = sites[right_id];
	const vec2<float_t> connector = right_site - left_site, u = rem_site - left_site;
	if (connector.get_lrot().dot(u) <= 0)
		return;
	const vec2<float_t> cc = circumc(left_site, rem_site, right_site);
	events.push(FortuneEvent(cc.y + (cc - left_site).len(), left_site, rem_site, right_site));
}

template<typename float_t>
float_t Voronoi<float_t>::BeachLinePoint::get_x() const
{
	if (is_dummy)
		return dummy_x;
	return 0; // TODO implement
}

template<typename float_t>
Voronoi<float_t>::BeachLinePoint::BeachLinePoint(const float_t dummy_x) : is_dummy(true), dummy_x(dummy_x)
{
}

template<typename float_t>
Voronoi<float_t>::BeachLinePoint::BeachLinePoint(const uint32_t point_id, const uint32_t left_site, const uint32_t right_site, const float_t *const t, const vec<vec2<float_t>> *const sites)
	: is_dummy(false), point_id(point_id), left_site(left_site), right_site(right_site), t(t), sites(sites)
{
}

template<typename float_t>
bool Voronoi<float_t>::BeachLinePoint::operator<(const BeachLinePoint &other) const
{
	const float_t this_x = get_x();
	const float_t other_x = other.get_x();
	return this_x < other_x || (this_x == other_x && (left_site < other.left_site || (left_site == other.left_site && right_site < other.right_site)));
}

template<typename float_t>
Voronoi<float_t>::FortuneEvent::FortuneEvent(const float_t time, const uint32_t site_id) : time(time), is_site(true), site_id(site_id)
{
}

template<typename float_t>
Voronoi<float_t>::FortuneEvent::FortuneEvent(const float_t time, const uint32_t left_site, const uint32_t rem_site, const uint32_t right_site)
	: is_site(false), left_site(left_site), rem_site(rem_site), right_site(right_site)
{
}

template<typename float_t>
bool Voronoi<float_t>::FortuneEvent::operator<(const FortuneEvent &other) const
{
	return time < other.time;
}

}