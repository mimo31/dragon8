/**
 * voronoi.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/04/26
 */
#ifndef VORONOI_HPP
#define VORONOI_HPP

#include <unordered_map>
#include <queue>
#include <set>

#include "print.hpp"
#include "shape.hpp"
#include "vec.hpp"
#include "vec2.hpp"

namespace dragon8
{

struct PointSiteData
{
	uint32_t left_site;
	uint32_t right_site;

	PointSiteData() = default;

	PointSiteData(const uint32_t left_site, const uint32_t right_site) : left_site(left_site), right_site(right_site)
	{
	}

	bool operator==(const PointSiteData &other) const
	{
		return left_site == other.left_site && right_site == other.right_site;
	}
};

}

namespace std
{
	template<>
	struct hash<dragon8::PointSiteData>
	{
		size_t operator()(const dragon8::PointSiteData &site_data) const
		{
			return (((71 * size_t(site_data.left_site)) << 32) + 127 * site_data.right_site);
		}
	};
}

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

	float_t y0;

	struct BeachLinePoint
	{
	public:
		uint32_t left_site, right_site;
		int32_t left_child_id, right_child_id, parent_id;
		int32_t left_neighbor_id, right_neighbor_id;
		/// Depth of the subtree below and including this point
		uint32_t depth;

		BeachLinePoint() = default;

		BeachLinePoint(uint32_t left_site, uint32_t right_site, int32_t left_child_id, int32_t right_child_id, int32_t parent_id, int32_t left_neighbor_id, int32_t right_neighbor_id, uint32_t depth);

		float_t get_x(float_t y0, const vec<vec2<float_t>> &sites) const;

		void print(float_t y0, const vec<vec2<float_t>> &sites) const;
	};

	/// Index (id) of the point at the root of the AVL tree
	int32_t root_id = -1;
	int32_t most_left_id = -1;
	int32_t most_right_id = -1;
	vec<BeachLinePoint> beach_line_points;
	std::unordered_map<PointSiteData, uint32_t> sites_to_ind;

	/**
	 * @param search_x value of x to search for in the beach line
	 * @return index (id) of the beach line point before which a point with the given x should be put (or -1 iff it should be put at the end)
	 */
	int32_t find_put_before(float_t search_x) const;
	/**
	 * Inserts a new point into the beach line by updating the AVL tree (updates pointers in the tree and potentially rotates the tree to balance it)
	 * @param before_id index (id) of the point before which the point should be added (or -1 iff it should be put at the end)
	 * @param left_site left site for the point to add
	 * @param right_site right site for the point to add
	 */
	void put_before(int32_t before_id, uint32_t left_site, uint32_t right_site);
	/**
	 * Rotates edges coming out of the specified vertex in the AVL tree above the beach line iff the vertex is unbalanced.
	 * Always updates the depth at the specified point, but may invalidate depths at its ancestors.
	 * Assumes that all descedants of the specified vertex are balanced and have correct depths.
	 * @param point_id index (id) of the point to check and possibly rotate
	 */
	void possibly_rotate(uint32_t point_id);
	void rotate_upward(uint32_t point_id);
	/**
	 * @param point_id index (id) of the point a reference to which we want
	 * @return pointer to the variable that holds the reference to the specified point
	 */
	int32_t *get_holder(uint32_t point_id);

	void erase(uint32_t point_id);

	struct FortuneEvent
	{
	public:
		float_t time;
		bool is_site;
		uint32_t site_id;
		uint32_t left_site, rem_site, right_site;

		FortuneEvent(float_t time, uint32_t site_id);
		FortuneEvent(float_t time, uint32_t left_site, uint32_t rem_site, uint32_t right_site);

		bool operator<(const FortuneEvent &other) const;
	};

	std::priority_queue<FortuneEvent> events;

	void potential_add_vertex_event(uint32_t left_id, uint32_t rem_id, uint32_t right_id);

	void submit_vertex(uint32_t site0_id, uint32_t site1_id, uint32_t site2_id);

	void print_beach_line() const;

	void check_tree_structure() const;
	void check_avl_structure() const;
	vec<uint32_t> get_tree_vertices() const;
	void check_linked_list_structure() const;
	vec<uint32_t> get_linked_list_vertices() const;
	void check_map_structure() const;
	vec<uint32_t> get_map_vertices() const;
	void check_beach_line_data_structure() const;
	void check_beach_line_content() const;
	void check_beach_line() const;

	struct VoronoiVertex
	{
		uint32_t site0, site1, site2;

		VoronoiVertex() = default;

		VoronoiVertex(const uint32_t site0, const uint32_t site1, const uint32_t site2) : site0(site0), site1(site1), site2(site2)
		{
		}
	};

	vec<VoronoiVertex> vertices;

public:
	Voronoi(ShapePtr container);
	void init_sites(const vec<vec2<float_t>> &sites);
	void init_compute();
	void add_site(vec2<float_t> p);
	void remove_site(uint32_t i);
	vec<vec2<float_t>> get_vertices() const;
	FurthestPointData<float_t> get_furthest() const;
};

}

#endif // VORONOI_HPP