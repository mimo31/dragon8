/**
 * voronoi.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/04/26
 */
#include "voronoi.hpp"

#include <cassert>

#include "print.hpp"

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
	this->sites = sites;
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
		cout << "event" << endl;
		const FortuneEvent event = events.top();
		events.pop();

		y0 = event.time;

		cout << "beach line: ";
		print_beach_line();

		if (event.is_site)
		{
			const uint32_t site_id = event.site_id;
			if (!beach_line_points.empty())
			{
				const float_t search_x = sites[site_id].x;
				const int32_t put_before_id = find_put_before(search_x);
				const int32_t put_after_id = put_before_id != -1 ? beach_line_points[put_before_id].left_neighbor_id : most_right_id;
				const uint32_t put_into_site = put_before_id != -1 ? beach_line_points[put_before_id].left_site : beach_line_points[most_right_id].right_site;
				put_before(put_before_id, put_into_site, site_id);
				put_before(put_before_id, site_id, put_into_site);
				/*const BeachLinePoint dummy(dummy_x);
				auto lb = beach_line.lower_bound(dummy);
				const uint32_t next_id = beach_line_points.size();*/
				// TODO: deal with potential add vertex calls
				if (put_before_id != -1)
					potential_add_vertex_event(site_id, put_into_site, beach_line_points[put_before_id].right_site);
				if (put_after_id != -1)
					potential_add_vertex_event(beach_line_points[put_after_id].left_site, put_into_site, site_id);
				/*if (lb != beach_line.end())
				{
					const BeachLinePoint put_before = *lb;
					const uint32_t put_into_site = put_before.left_site;
					const BeachLinePoint p0(next_id, put_into_site, site_id, &t, &sites);
					const BeachLinePoint p1(next_id + 1, site_id, put_into_site, &t, &sites);
					beach_line.insert(p0);
					beach_line.insert(p1);
					const int32_t prev_point_id = beach_line_points[put_before.point_id].left_id;
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
					const int32_t next_point_id = beach_line_points[put_after.point_id].right_id;
					beach_line_points.push_back(BeachLinePointData(next_id, next_point_id, site_id, put_into_site));
					if (next_point_id != -1)
					{
						beach_line_points[next_point_id].left_id = next_id + 1;
						potential_add_vertex_event(site_id, put_into_site, beach_line_points[next_point_id].right_site);
					}
					beach_line_points[put_after.point_id].right_id = next_id;
					potential_add_vertex_event(put_after.left_site, put_into_site, site_id);
				}*/
			}
			else
			{
				const BeachLinePoint p0(site0, site_id, -1, 1, -1, -1, 1, 1);
				const BeachLinePoint p1(site_id, site_id, -1, -1, 0, 0, -1, 0);
				beach_line_points.push_back(p0);
				beach_line_points.push_back(p1);
				root_id = 0;
				most_left_id = 0;
				most_right_id = 1;
				/*const BeachLinePoint p1(1, site_id, site0, &t, &sites);
				const BeachLinePoint p0(0, site0, site_id, &t, &sites);
				const BeachLinePoint p1(1, site_id, site0, &t, &sites);
				beach_line.insert(p0);
				beach_line.insert(p1);
				beach_line_points.push_back(BeachLinePointData(-1, 1, site0, site_id));
				beach_line_points.push_back(BeachLinePointData(0, -1, site_id, site0));*/
			}
		}
		else
		{
			auto p0it = sites_to_ind.find(PointSiteData(event.left_site, event.rem_site));
			if (p0it == sites_to_ind.end())
				continue;
			
			
			/*const BeachLinePoint left_point(-1, event.left_site, event.rem_site, &t, &sites);
			auto it = beach_line.find(left_point);
			if (it == beach_line.end())
				continue;*/
			const uint32_t left_point_id = p0it->second;
			const uint32_t right_point_id = beach_line_points[left_point_id].right_neighbor_id;
			//assert(right_point_id != -1);
			if (beach_line_points[right_point_id].right_site != event.right_site)
				continue;
			submit_vertex(event.left_site, event.rem_site, event.right_site);
			//const uint32_t new_point_id = beach_line_points.size();
			//const int32_t ll_point_id = beach_line_points[left_point_id].left_neighbor_id;
			const int32_t rr_point_id = beach_line_points[right_point_id].right_neighbor_id;
			erase(left_point_id);
			erase(right_point_id);
			put_before(rr_point_id, event.left_site, event.right_site);
			//beach_line_points.push_back(BeachLinePointData(ll_point_id, rr_point_id, event.left_site, event.right_site));
			/*if (ll_point_id != -1)
				beach_line_points[ll_point_id].right_id = new_point_id;
			if (rr_point_id != -1)
				beach_line_points[rr_point_id].left_id = new_point_id;
			
			beach_line.insert(BeachLinePoint(new_point_id, event.left_site, event.right_site, &t, &sites));
			beach_line.erase(BeachLinePoint(-1, event.left_site, event.rem_site, &t, &sites));
			beach_line.erase(BeachLinePoint(-1, event.rem_site, event.right_site, &t, &sites));*/
		}
		cout << "after add ";
		print_beach_line();
	}
}

template<typename float_t>
int32_t Voronoi<float_t>::find_put_before(const float_t search_x) const
{
	int32_t from_node = -1;
	int32_t next_node = root_id;
	bool from_upper = false;
	while (next_node != -1)
	{
		const float_t cur_x = beach_line_points[next_node].get_x(y0, sites);
		from_node = next_node;
		from_upper = search_x < cur_x;
		next_node = from_upper ? beach_line_points[next_node].left_child_id : beach_line_points[next_node].right_child_id;
	}
	return from_upper ? from_node : beach_line_points[from_node].right_neighbor_id;
}

template<typename float_t>
void Voronoi<float_t>::put_before(const int32_t before_id, const uint32_t left_site, const uint32_t right_site)
{
	int32_t left_neighbor_id;
	int32_t right_neighbor_id;
	int32_t parent_id;
	const uint32_t new_id = beach_line_points.size();
	auto put_as_left_child = [&](const uint32_t par_id)
	{
		parent_id = par_id;
		BeachLinePoint &parent = beach_line_points[par_id];
		parent.left_child_id = new_id;
		left_neighbor_id = parent.left_neighbor_id;
		right_neighbor_id = par_id;
		parent.left_neighbor_id = new_id;
		if (left_neighbor_id == -1)
			most_left_id = new_id;
	};
	auto put_as_right_child = [&](const uint32_t par_id)
	{
		parent_id = par_id;
		BeachLinePoint &parent = beach_line_points[par_id];
		parent.right_child_id = new_id;
		left_neighbor_id = par_id;
		right_neighbor_id = parent.right_neighbor_id;
		parent.right_neighbor_id = new_id;
		if (right_neighbor_id == -1)
			most_right_id = new_id;
	};

	if (before_id != -1)
	{
		if (beach_line_points[before_id].left_child_id == -1)
			put_as_left_child(before_id);
		else
			put_as_right_child(beach_line_points[before_id].left_neighbor_id);
	}
	else
		put_as_right_child(most_right_id);
	
	beach_line_points.push_back(BeachLinePoint(left_site, right_site, -1, -1, parent_id, left_neighbor_id, right_neighbor_id, 0));

	rotate_upward(parent_id);
	sites_to_ind[PointSiteData(left_site, right_site)] = new_id;
}

template<typename float_t>
void Voronoi<float_t>::possibly_rotate(const uint32_t point_id)
{
	BeachLinePoint &point = beach_line_points[point_id];
	if (point.left_child_id == -1 && point.right_child_id == -1)
	{
		point.depth = 0;
		return;
	}
	auto update_depth = [this](const uint32_t p_id)
	{
		BeachLinePoint &p = beach_line_points[p_id];
		const uint32_t d0 = p.left_child_id != -1 ? beach_line_points[p.left_child_id].depth + 1 : 0;
		const uint32_t d1 = p.right_child_id != -1 ? beach_line_points[p.right_child_id].depth + 1 : 0;
		p.depth = std::max(d0, d1);
	};
	auto rotate_left_heavy = [&]
	{
		int32_t *const entry_to_parent = get_holder(point_id);
		const uint32_t old_left_child_id = point.left_child_id;
		BeachLinePoint &old_left_child = beach_line_points[old_left_child_id];
		const bool need_deep_split = old_left_child.right_child_id != -1 && (old_left_child.left_child_id == -1 ||
			beach_line_points[old_left_child.right_child_id].depth > beach_line_points[old_left_child.left_child_id].depth);
		if (need_deep_split)
		{
			const uint32_t old_right_left_child_id = old_left_child.right_child_id;
			BeachLinePoint &old_right_left_child = beach_line_points[old_right_left_child_id];
			point.left_child_id = old_right_left_child.right_child_id;
			old_left_child.right_child_id = old_right_left_child.left_child_id;
			old_right_left_child.left_child_id = old_left_child_id;
			old_right_left_child.right_child_id = point_id;

			old_right_left_child.parent_id = point.parent_id;
			*entry_to_parent = old_right_left_child_id;
			point.parent_id = old_right_left_child_id;
			old_left_child.parent_id = old_right_left_child_id;
			if (point.left_child_id != -1)
				beach_line_points[point.left_child_id].parent_id = point_id;
			if (old_left_child.right_child_id != -1)
				beach_line_points[old_left_child.right_child_id].parent_id = old_left_child_id;
			update_depth(point_id);
			update_depth(old_left_child_id);
			update_depth(old_right_left_child_id);
		}
		else
		{
			point.left_child_id = old_left_child.right_child_id;
			old_left_child.right_child_id = point_id;

			old_left_child.parent_id = point.parent_id;
			*entry_to_parent = old_left_child_id;
			point.parent_id = old_left_child_id;
			if (point.left_child_id != -1)
				beach_line_points[point.left_child_id].parent_id = point_id;
			update_depth(point_id);
			update_depth(old_left_child_id);
		}
	};
	auto rotate_right_heavy = [&]
	{
		int32_t *const entry_to_parent = get_holder(point_id);
		const uint32_t old_right_child_id = point.right_child_id;
		BeachLinePoint &old_right_child = beach_line_points[old_right_child_id];
		const bool need_deep_split = old_right_child.left_child_id != -1 && (old_right_child.right_child_id == -1 || 
			beach_line_points[old_right_child.left_child_id].depth > beach_line_points[old_right_child.right_child_id].depth);
		if (need_deep_split)
		{
			const uint32_t old_left_right_child_id = old_right_child.left_child_id;
			BeachLinePoint &old_left_right_child = beach_line_points[old_left_right_child_id];
			point.right_child_id = old_left_right_child.left_child_id;
			old_right_child.left_child_id = old_left_right_child.right_child_id;
			old_left_right_child.left_child_id = point_id;
			old_left_right_child.right_child_id = old_right_child_id;

			old_left_right_child.parent_id = point.parent_id;
			*entry_to_parent = old_left_right_child_id;
			point.parent_id = old_left_right_child_id;
			old_right_child.parent_id = old_left_right_child_id;
			if (point.right_child_id != -1)
				beach_line_points[point.right_child_id].parent_id = point_id;
			if (old_right_child.left_child_id != -1)
				beach_line_points[old_right_child.left_child_id].parent_id = old_right_child_id;
			update_depth(point_id);
			update_depth(old_right_child_id);
			update_depth(old_left_right_child_id);
		}
		else
		{
			point.right_child_id = old_right_child.left_child_id;
			old_right_child.left_child_id = point_id;

			old_right_child.parent_id = point.parent_id;
			*entry_to_parent = old_right_child_id;
			point.parent_id = old_right_child_id;
			if (point.right_child_id != -1)
				beach_line_points[point.right_child_id].parent_id = point_id;
			update_depth(point_id);
			update_depth(old_right_child_id);
		}
	};
	if (point.left_child_id != -1 && beach_line_points[point.left_child_id].depth > 1
		&& (point.right_child_id == -1 || beach_line_points[point.left_child_id].depth > beach_line_points[point.right_child_id].depth + 1))
		rotate_left_heavy();
	else if (point.right_child_id != -1 && beach_line_points[point.right_child_id].depth > 1
		&& (point.left_child_id == -1 || beach_line_points[point.right_child_id].depth > beach_line_points[point.left_child_id].depth + 1))
		rotate_right_heavy();
	else
		update_depth(point_id);
}

template<typename float_t>
void Voronoi<float_t>::rotate_upward(const uint32_t point_id)
{
	int32_t next_point_id = point_id;
	while (next_point_id != -1)
	{
		const uint32_t original_depth = beach_line_points[next_point_id].depth;
		int32_t *holder = get_holder(next_point_id);
		possibly_rotate(next_point_id);
		next_point_id = *holder;
		if (beach_line_points[next_point_id].depth == original_depth)
			break;
		next_point_id = beach_line_points[next_point_id].parent_id;
	}
}

template<typename float_t>
int32_t *Voronoi<float_t>::get_holder(const uint32_t point_id)
{
	BeachLinePoint &point = beach_line_points[point_id];
	if (point.parent_id != -1)
	{
		BeachLinePoint &parent = beach_line_points[point.parent_id];
		if (parent.left_child_id == int32_t(point_id))
			return &parent.left_child_id;
		return &parent.right_child_id;
	}
	return &root_id;
}

template<typename float_t>
void Voronoi<float_t>::erase(const uint32_t point_id)
{
	BeachLinePoint &point = beach_line_points[point_id];
	if (point.left_neighbor_id != -1)
		beach_line_points[point.left_neighbor_id].right_neighbor_id = point.right_neighbor_id;
	if (point.right_neighbor_id != -1)
		beach_line_points[point.right_neighbor_id].left_neighbor_id = point.left_neighbor_id;
	int32_t *const holder = get_holder(point_id);
	if (point.left_child_id == -1)
	{
		*holder = point.right_child_id;
		if (point.right_child_id != -1)
			beach_line_points[point.right_child_id].parent_id = point.parent_id;
		rotate_upward(point.parent_id);
		return;
	}
	if (point.right_child_id == -1)
	{
		*holder = point.left_child_id;
		if (point.left_child_id != -1)
			beach_line_points[point.left_child_id].parent_id = point.parent_id;
		rotate_upward(point.parent_id);
		return;
	}
	const uint32_t leaf_id = point.right_neighbor_id;
	BeachLinePoint &leaf = beach_line_points[leaf_id];
	*holder = leaf_id;
	const uint32_t leaf_parent_id = leaf.parent_id;
	BeachLinePoint &leaf_parent = beach_line_points[leaf_parent_id];
	leaf_parent.left_child_id = -1;
	leaf.left_child_id = point.left_child_id;
	leaf.right_child_id = point.right_child_id;
	beach_line_points[leaf.left_child_id].parent_id = leaf_id;
	beach_line_points[leaf.right_child_id].parent_id = leaf_id;
	leaf.parent_id = point.parent_id;
	leaf.depth = point.depth;
	rotate_upward(leaf_parent_id);
	sites_to_ind.erase(sites_to_ind.find(PointSiteData(point.left_site, point.right_site)));
}

template<typename float_t>
void Voronoi<float_t>::BeachLinePoint::print(const float_t y0, const vec<vec2<float_t>> &sites) const
{
	cout << "bpoint(lsite = " << left_site << ", rsite = " << right_site << ", x = " << get_x(y0, sites) << ")";
}

template<typename float_t>
void Voronoi<float_t>::print_beach_line() const
{
	int32_t next_id = most_left_id;
	while (next_id != -1)
	{
		const BeachLinePoint &point = beach_line_points[most_left_id];
		point.print(y0, sites);
		cout << ' ';
		next_id = point.right_neighbor_id;
	}
	cout << endl;
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
	cout << "potential" << endl;
	const vec2<float_t> left_site = sites[left_id], rem_site = sites[rem_id], right_site = sites[right_id];
	const vec2<float_t> connector = right_site - left_site, u = rem_site - left_site;
	cout << "dot = " << connector.get_lrot().dot(u) << endl;
	if (connector.get_lrot().dot(u) <= 0)
		return;
	const vec2<float_t> cc = circumc(left_site, rem_site, right_site);
	events.push(FortuneEvent(cc.y + (cc - left_site).len(), left_id, rem_id, right_id));
}

template<typename float_t>
void Voronoi<float_t>::submit_vertex(const uint32_t site0_id, const uint32_t site1_id, const uint32_t site2_id)
{
	std::cout << "submitted vertex between " << site0_id << " " << site1_id << " " << site2_id << std::endl;	
}

template<typename float_t>
float_t Voronoi<float_t>::BeachLinePoint::get_x(const float_t y0, const vec<vec2<float_t>> &sites) const
{
	const vec2<float_t> l = sites[left_site], r = sites[right_site];
	const vec2<float_t> s0 = l ^ r;
	const vec2<float_t> conn = s0 - l;
	const float_t a2 = conn.len2();
	const vec2<float_t> v = -conn.get_lrot();
	const float_t yd = s0.y - y0;
	const float_t lam = (v.y * yd + sqrt(a2 * std::max(yd * yd - v.x * v.x, float_t(0)))) / (v.x * v.x);
	return (s0 + lam * v).x;
}

template<typename float_t>
Voronoi<float_t>::BeachLinePoint::BeachLinePoint(const uint32_t left_site, const uint32_t right_site, const int32_t left_child_id, const int32_t right_child_id, const int32_t parent_id, const int32_t left_neighbor_id, const int32_t right_neighbor_id, const uint32_t depth)
	: left_site(left_site), right_site(right_site), left_child_id(left_child_id), right_child_id(right_child_id), parent_id(parent_id), left_neighbor_id(left_neighbor_id), right_neighbor_id(right_neighbor_id), depth(depth)
{
}

/*
template<typename float_t>
bool Voronoi<float_t>::BeachLinePoint::operator<(const BeachLinePoint &other) const
{
	const float_t this_x = get_x();
	const float_t other_x = other.get_x();
	return this_x < other_x || (this_x == other_x && (left_site < other.left_site || (left_site == other.left_site && right_site < other.right_site)));
}*/

template<typename float_t>
Voronoi<float_t>::FortuneEvent::FortuneEvent(const float_t time, const uint32_t site_id) 
	: time(time), is_site(true), site_id(site_id)
{
}

template<typename float_t>
Voronoi<float_t>::FortuneEvent::FortuneEvent(const float_t time, const uint32_t left_site, const uint32_t rem_site, const uint32_t right_site)
	: time(time), is_site(false), left_site(left_site), rem_site(rem_site), right_site(right_site)
{
}

template<typename float_t>
bool Voronoi<float_t>::FortuneEvent::operator<(const FortuneEvent &other) const
{
	return time < other.time;
}

}