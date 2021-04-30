/**
 * voronoi.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/04/26
 */
#include "voronoi.hpp"

#include <cassert>
#include <functional>

#include "print.hpp"

namespace dragon8
{

constexpr double eps = 1e-8;

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

constexpr bool print_info = true;
constexpr bool do_checks = true;

template<typename float_t>
void Voronoi<float_t>::init_compute()
{
	for (uint32_t i = 0; i < sites.size(); i++)
		events.push(FortuneEvent(sites[i].y, i));
	const uint32_t site0 = events.top().site_id;
	events.pop();
	while (!events.empty())
	{
		if (do_checks)
			check_beach_line();
		const FortuneEvent event = events.top();
		events.pop();

		y0 = event.time;
		if (print_info)
		{
			cout << "\n----\ny0 = " << y0 << endl;
			cout << "beach line: ";
			print_beach_line();
		}

		if (event.is_site)
		{
			const uint32_t site_id = event.site_id;
			if (print_info)
				cout << "adding site_id " << site_id << " (= " << sites[site_id] << ")";
			if (!beach_line_points.empty())
			{
				if (print_info)
					cout << " (beach line non-empty)" << endl;
				const float_t search_x = sites[site_id].x;
				const int32_t put_before_id = find_put_before(search_x);
				const int32_t put_after_id = put_before_id != -1 ? beach_line_points[put_before_id].left_neighbor_id : most_right_id;
				const uint32_t put_into_site = put_before_id != -1 ? beach_line_points[put_before_id].left_site : beach_line_points[most_right_id].right_site;
				put_before(put_before_id, put_into_site, site_id);
				put_before(put_before_id, site_id, put_into_site);
				if (put_before_id != -1)
					potential_add_vertex_event(site_id, put_into_site, beach_line_points[put_before_id].right_site);
				if (put_after_id != -1)
					potential_add_vertex_event(beach_line_points[put_after_id].left_site, put_into_site, site_id);
				if (do_checks)
					check_beach_line();
			}
			else
			{
				if (print_info)
					cout << " (beach line empty)" << endl;
				const BeachLinePoint p0(site0, site_id, -1, 1, -1, -1, 1, 1);
				const BeachLinePoint p1(site_id, site0, -1, -1, 0, 0, -1, 0);
				beach_line_points.push_back(p0);
				beach_line_points.push_back(p1);
				root_id = 0;
				most_left_id = 0;
				most_right_id = 1;
				sites_to_ind[PointSiteData(site0, site_id)] = 0;
				sites_to_ind[PointSiteData(site_id, site0)] = 1;
				if (do_checks)
					check_beach_line();
			}
		}
		else
		{
			if (print_info)
				cout << "vertex event between sites " << event.left_site << ", " << event.rem_site << ", " << event.right_site << endl;
			auto p0it = sites_to_ind.find(PointSiteData(event.left_site, event.rem_site));
			if (p0it == sites_to_ind.end())
				continue;
			
			const uint32_t left_point_id = p0it->second;
			const uint32_t right_point_id = beach_line_points[left_point_id].right_neighbor_id;
			if (beach_line_points[right_point_id].right_site != event.right_site)
				continue;
			submit_vertex(event.left_site, event.rem_site, event.right_site);
			const int32_t rr_point_id = beach_line_points[right_point_id].right_neighbor_id;
			erase(left_point_id);
			erase(right_point_id);
			put_before(rr_point_id, event.left_site, event.right_site);
			const BeachLinePoint &new_point = beach_line_points.back();
			if (new_point.left_neighbor_id != -1)
				potential_add_vertex_event(beach_line_points[new_point.left_neighbor_id].left_site, event.left_site, event.right_site);
			if (new_point.right_neighbor_id != -1)
				potential_add_vertex_event(event.left_site, event.right_site, beach_line_points[new_point.right_neighbor_id].right_site);
			if (do_checks)
				check_beach_line();
		}
		//cout << "after add ";
		//print_beach_line();
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
	int32_t parent_id;
	const uint32_t new_id = beach_line_points.size();
	auto put_as_left_child = [&](const uint32_t par_id)
	{
		parent_id = par_id;
		BeachLinePoint &parent = beach_line_points[par_id];
		parent.left_child_id = new_id;
	};
	auto put_as_right_child = [&](const uint32_t par_id)
	{
		parent_id = par_id;
		BeachLinePoint &parent = beach_line_points[par_id];
		parent.right_child_id = new_id;
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
	
	const int32_t left_neighbor_id = before_id != -1 ? beach_line_points[before_id].left_neighbor_id : most_right_id;
	if (left_neighbor_id != -1)
		beach_line_points[left_neighbor_id].right_neighbor_id = new_id;
	else
		most_left_id = new_id;
	if (before_id != -1)
		beach_line_points[before_id].left_neighbor_id = new_id;
	else
		most_right_id = new_id;
	
	beach_line_points.push_back(BeachLinePoint(left_site, right_site, -1, -1, parent_id, left_neighbor_id, before_id, 0));

	rotate_upward(parent_id);
	sites_to_ind[PointSiteData(left_site, right_site)] = new_id;

	if (do_checks)
		check_beach_line_data_structure();
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
	if (point.left_child_id != -1 && beach_line_points[point.left_child_id].depth >= 1
		&& (point.right_child_id == -1 || beach_line_points[point.left_child_id].depth > beach_line_points[point.right_child_id].depth + 1))
		rotate_left_heavy();
	else if (point.right_child_id != -1 && beach_line_points[point.right_child_id].depth >= 1
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
	if (do_checks)
		check_beach_line_data_structure();
	BeachLinePoint &point = beach_line_points[point_id];
	if (point.left_neighbor_id != -1)
		beach_line_points[point.left_neighbor_id].right_neighbor_id = point.right_neighbor_id;
	else
		most_left_id = point.right_neighbor_id;
	if (point.right_neighbor_id != -1)
		beach_line_points[point.right_neighbor_id].left_neighbor_id = point.left_neighbor_id;
	else
		most_right_id = point.left_neighbor_id;
	int32_t *const holder = get_holder(point_id);
	if (point.left_child_id == -1)
	{
		*holder = point.right_child_id;
		if (point.right_child_id != -1)
			beach_line_points[point.right_child_id].parent_id = point.parent_id;
		rotate_upward(point.parent_id);
	}
	else if (point.right_child_id == -1)
	{
		*holder = point.left_child_id;
		if (point.left_child_id != -1)
			beach_line_points[point.left_child_id].parent_id = point.parent_id;
		rotate_upward(point.parent_id);
	}
	else
	{
		const uint32_t leaf_id = point.right_neighbor_id;
		BeachLinePoint &leaf = beach_line_points[leaf_id];
		*holder = leaf_id;
		const uint32_t leaf_parent_id = leaf.parent_id;
		leaf.parent_id = point.parent_id;
		leaf.depth = point.depth;
		if (leaf_parent_id != point_id)
		{
			BeachLinePoint &leaf_parent = beach_line_points[leaf_parent_id];
			leaf_parent.left_child_id = leaf.right_child_id;
			if (leaf.right_child_id != -1)
				beach_line_points[leaf.right_child_id].parent_id = leaf_parent_id;
			leaf.left_child_id = point.left_child_id;
			leaf.right_child_id = point.right_child_id;
			beach_line_points[leaf.left_child_id].parent_id = leaf_id;
			beach_line_points[leaf.right_child_id].parent_id = leaf_id;
			rotate_upward(leaf_parent_id);
		}
		else
		{
			leaf.left_child_id = point.left_child_id;
			beach_line_points[leaf.left_child_id].parent_id = leaf_id;
			rotate_upward(leaf_id);
		}
	}
	
	const auto it = sites_to_ind.find(PointSiteData(point.left_site, point.right_site));
	sites_to_ind.erase(it);
	
	if (do_checks)
		check_beach_line_data_structure();
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
		const BeachLinePoint &point = beach_line_points[next_id];
		point.print(y0, sites);
		cout << ' ';
		next_id = point.right_neighbor_id;
	}
	cout << endl;
}


template<typename float_t>
void Voronoi<float_t>::check_tree_structure() const
{
	vec<bool> visited(beach_line_points.size());
	std::fill(visited.begin(), visited.end(), false);
	const std::function<void(uint32_t)> visit = [&](const int32_t point_id)
	{
		if (point_id == -1)
			return;
		const BeachLinePoint &point = beach_line_points[point_id];
		assert(!visited[point_id]);
		visited[point_id] = true;
		assert(point.left_child_id == -1 || beach_line_points[point.left_child_id].parent_id == point_id);
		assert(point.right_child_id == -1 || beach_line_points[point.right_child_id].parent_id == point_id);
		visit(point.left_child_id);
		visit(point.right_child_id);
	};
	visit(root_id);
}

template<typename float_t>
void Voronoi<float_t>::check_avl_structure() const
{
	const std::function<void(uint32_t)> visit = [&](const int32_t point_id)
	{
		if (point_id == -1)
			return;
		const BeachLinePoint &point = beach_line_points[point_id];
		const uint32_t d0 = point.left_child_id != -1 ? beach_line_points[point.left_child_id].depth + 1 : 0;
		const uint32_t d1 = point.right_child_id != -1 ? beach_line_points[point.right_child_id].depth + 1 : 0;
		assert(d0 + 1 >= d1 && d1 + 1 >= d0);
		assert(point.depth == std::max(d0, d1));
		visit(point.left_child_id);
		visit(point.right_child_id);
	};
	visit(root_id);
}

template<typename float_t>
vec<uint32_t> Voronoi<float_t>::get_tree_vertices() const
{
	vec<uint32_t> verts;
	const std::function<void(uint32_t)> visit = [&](const int32_t point_id)
	{
		if (point_id == -1)
			return;
		const BeachLinePoint &point = beach_line_points[point_id];
		visit(point.left_child_id);
		verts.push_back(point_id);
		visit(point.right_child_id);
	};
	visit(root_id);
	return verts;
}

template<typename float_t>
void Voronoi<float_t>::check_linked_list_structure() const
{
	vec<bool> visited(beach_line_points.size());
	std::fill(visited.begin(), visited.end(), false);
	int32_t next_id = most_left_id;
	while (next_id != -1)
	{
		assert(!visited[next_id]);
		visited[next_id] = true;
		const BeachLinePoint &point = beach_line_points[next_id];
		assert(point.left_neighbor_id == -1 || beach_line_points[point.left_neighbor_id].right_neighbor_id == next_id);
		assert((point.left_neighbor_id == -1) == (next_id == most_left_id));
		assert((point.right_neighbor_id == -1) == (next_id == most_right_id));
		next_id = point.right_neighbor_id;
	}
}

template<typename float_t>
vec<uint32_t> Voronoi<float_t>::get_linked_list_vertices() const
{
	vec<uint32_t> verts;
	int32_t next_id = most_left_id;
	while (next_id != -1)
	{
		verts.push_back(next_id);
		next_id = beach_line_points[next_id].right_neighbor_id;
	}
	return verts;
}

template<typename float_t>
void Voronoi<float_t>::check_map_structure() const
{
	vec<bool> visited(beach_line_points.size());
	std::fill(visited.begin(), visited.end(), false);
	for (const auto &entry : sites_to_ind)
	{
		const PointSiteData site_data = entry.first;
		const uint32_t ind = entry.second;
		assert(ind < beach_line_points.size());
		const BeachLinePoint &point = beach_line_points[ind];
		assert(point.left_site == site_data.left_site);
		assert(point.right_site == site_data.right_site);
		assert(!visited[ind]);
		visited[ind] = true;
	}
}

template<typename float_t>
vec<uint32_t> Voronoi<float_t>::get_map_vertices() const
{
	vec<uint32_t> verts;
	for (const auto &entry : sites_to_ind)
		verts.push_back(entry.second);
	return verts;
}

template<typename float_t>
void Voronoi<float_t>::check_beach_line_data_structure() const
{
	check_tree_structure();
	check_avl_structure();
	check_linked_list_structure();
	check_map_structure();

	vec<uint32_t> tree_verts = get_tree_vertices();
	vec<uint32_t> ll_verts = get_linked_list_vertices();
	assert(tree_verts == ll_verts);
	
	vec<uint32_t> map_verts = get_map_vertices();
	std::sort(tree_verts.begin(), tree_verts.end());
	std::sort(map_verts.begin(), map_verts.end());
	assert(tree_verts == map_verts);
}

template<typename float_t>
void Voronoi<float_t>::check_beach_line_content() const
{
	int32_t next_id = most_left_id;
	while (next_id != -1)
	{
		const BeachLinePoint &point = beach_line_points[next_id];
		if (point.right_neighbor_id != -1)
			assert(beach_line_points[point.right_neighbor_id].left_site == point.right_site);
		next_id = beach_line_points[next_id].right_neighbor_id;
	}
}

template<typename float_t>
void Voronoi<float_t>::check_beach_line() const
{
	check_beach_line_data_structure();
	check_beach_line_content();
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
	events.push(FortuneEvent(cc.y - (cc - left_site).len(), left_id, rem_id, right_id));
}

template<typename float_t>
void Voronoi<float_t>::submit_vertex(const uint32_t site0_id, const uint32_t site1_id, const uint32_t site2_id)
{
	vertices.push_back(VoronoiVertex(site0_id, site1_id, site2_id));	
}

template<typename float_t>
vec<vec2<float_t>> Voronoi<float_t>::get_vertices() const
{
	vec<vec2<float_t>> verts;
	verts.reserve(vertices.size());
	for (const VoronoiVertex &v : vertices)
		verts.push_back(circumc(sites[v.site0], sites[v.site1], sites[v.site2]));
	return verts;
}

template<typename float_t>
float_t Voronoi<float_t>::BeachLinePoint::get_x(const float_t y0, const vec<vec2<float_t>> &sites) const
{
	const vec2<float_t> l = sites[left_site], r = sites[right_site];
	const vec2<float_t> s0 = l ^ r;
	const vec2<float_t> conn = s0 - l;
	if (abs(conn.y) < eps)
		return s0.x;
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

template class Voronoi<double>;

}