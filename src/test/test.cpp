/**
 * test.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/27
 */
#include "test.hpp"

#include <iostream>

#include "application.hpp"
#include "evaluator.hpp"
#include "graphics.hpp"
#include "shape.hpp"
#include "voronoi.hpp"

namespace dragon8
{

#define TEST_ASSERT(x) if (!(x)) { std::cout << __FUNCTION__ << " on line " << __LINE__ << " failed an assert" << std::endl; failed = true; }

void Tester::test_shape_polygon()
{
	const std::vector<vec2d> verts{ vec2d(0, 0), vec2d(1, 0), vec2d(.2, .2), vec2d(0, 1) };
	const ShapePolygon poly(verts);

	TEST_ASSERT(!poly.is_inside(vec2d(.8, .8)));
	TEST_ASSERT(!poly.is_inside(vec2d(.8, .9)));
	TEST_ASSERT(poly.is_inside(vec2d(.1, .1)));
	TEST_ASSERT(poly.is_inside(vec2d(.5, .05)));

	TEST_ASSERT(poly.get_box().c0 == vec2d(0, 0));
	TEST_ASSERT(poly.get_box().c1 == vec2d(1, 1));

	const std::vector<vec2d> verts2{ vec2d(0, 0), vec2d(3, 3), vec2d(0, 3), vec2d(3, 0), vec2d(1, 2.5), vec2d(2, 2.5) };
	const ShapePolygon poly2(verts2);

	TEST_ASSERT(!poly2.is_inside(vec2d(1.5, 1.55)));
	TEST_ASSERT(poly2.is_inside(vec2d(.5, 2.9)));
	TEST_ASSERT(poly2.is_inside(vec2d(2.5, .6)));
	TEST_ASSERT(!poly2.is_inside(vec2d(1.5, 2.45)));

	TEST_ASSERT(poly2.get_box().c0 == vec2d(0, 0));
	TEST_ASSERT(poly2.get_box().c1 == vec2d(3, 3));
}

void Tester::test_no_interior_polygon()
{
	const std::vector<vec2d> verts0{ vec2d(0, 0), vec2d(0, 1), vec2d(0, 0), vec2d(1, 0), vec2d(1, -1) };
	const std::vector<vec2d> verts1{ vec2d(0, 0), vec2d(0, 1), vec2d(0, 2) };
	const std::vector<vec2d> verts2{ vec2d(0, 0), vec2d(0, 1), vec2d(0, 0), vec2d(1, 0) };
	const std::vector<vec2d> verts3{ vec2d(0, 0), vec2d(0, 1), vec2d(1, 1), vec2d(0, 1), vec2d(0, 2), vec2d(-1, 2), vec2d(0, 2), vec2d(0, 1) };
	const std::vector<vec2d> verts4{ vec2d(0, 0), vec2d(0, .1), vec2d(.1, .1), vec2d(0, .1), vec2d(0, .2), vec2d(-.1, .2), vec2d(0, .2), vec2d(0, .1) };
	TEST_ASSERT(has_interior(verts0));
	TEST_ASSERT(!has_interior(verts1));
	TEST_ASSERT(!has_interior(verts2));
	TEST_ASSERT(!has_interior(verts3));
	TEST_ASSERT(!has_interior(verts4));
}

void Tester::test_to_nneg_int()
{
	uint32_t val;
	TEST_ASSERT(to_nneg_int("32", val) == NumberParseResult::OK && val == 32);
	TEST_ASSERT(to_nneg_int("00000000000000", val) == NumberParseResult::OK && val == 0);
	TEST_ASSERT(to_nneg_int("000000000000001", val) == NumberParseResult::OK && val == 1);
	TEST_ASSERT(to_nneg_int(" \t3 \n", val) == NumberParseResult::OK && val == 3);
	TEST_ASSERT(to_nneg_int("3240932  ", val) == NumberParseResult::OK && val == 3240932);

	TEST_ASSERT(to_nneg_int("29435872938476289347698273456", val) == NumberParseResult::TOO_BIG);
	TEST_ASSERT(to_nneg_int("\n\n\n  29435872933338476289347698273456\t", val) == NumberParseResult::TOO_BIG);

	TEST_ASSERT(to_nneg_int("", val) == NumberParseResult::SYNTAX_ERROR);
	TEST_ASSERT(to_nneg_int("\n\n\n  \t", val) == NumberParseResult::SYNTAX_ERROR);
	TEST_ASSERT(to_nneg_int(" a", val) == NumberParseResult::SYNTAX_ERROR);
	TEST_ASSERT(to_nneg_int(" 83 3 ", val) == NumberParseResult::SYNTAX_ERROR);
	TEST_ASSERT(to_nneg_int("-1", val) == NumberParseResult::SYNTAX_ERROR);
}

void Tester::test_all()
{
	test_shape_polygon();
	test_no_interior_polygon();
	test_to_nneg_int();
}

void draw_voronoi(const PointsState &points, const ShapePtr container, const std::string &filename)
{
	Voronoi<double> voron1(container);
	voron1.init_sites(points);
	voron1.init_compute();
	const vec<ProjectiveEdge<double>> raw_edges = voron1.get_edges();
	vec<std::pair<vec2d, vec2d>> realized_edges;
	vec<vec2d> extended_vertices = voron1.get_finite_vertices();
	for (const ProjectiveEdge<double> e : raw_edges)
	{
		if (e.v0.at_infinity && e.v1.at_infinity)
			continue;
		const vec<vec2d> inters = container->get_intersections(e);
		extended_vertices.insert(extended_vertices.end(), inters.begin(), inters.end());
		if (!e.v0.at_infinity && !e.v1.at_infinity)
		{
			realized_edges.push_back(std::make_pair(e.v0.v, e.v1.v));
			continue;
		}
		vec2d p = e.v0.v, v = e.v1.v;
		if (e.v0.at_infinity)
			std::swap(p, v);
		realized_edges.push_back(std::make_pair(p, p + v.get_unit() * 2));
	}
	write_image(container, points, filename, 0, extended_vertices, realized_edges);
}

void Tester::test_voronoi()
{
	Voronoi<double> voronoi(nullptr);
	/*const vec<vec2<double>> sites0({ vec2<double>(0, 0), vec2<double>(1, 1), vec2<double>(4, 3) });
	voronoi.init_sites(sites0);
	voronoi.init_compute();*/
	const std::vector<vec2d> verts{ vec2d(0, 0), vec2d(1, 0), vec2d(1, 1), vec2d(0, 1) };
	//const ShapePtr shape = std::make_shared<ShapeCircle>(1);
	const ShapePtr shape = std::make_shared<ShapePolygon>(verts);
	std::random_device dev;
	std::mt19937 rgen(dev());
	
	constexpr uint32_t point_count = 10;
	/*for (const vec2d p : st)
		cout << p << endl;
	*/
	//const PointsState st({ vec2d(0.460016, 0.802885), vec2d(0.14821, 0.672051), vec2d(0.77762, 0.170662), vec2d(0.958275, 0.894059) });
	//const PointsState st({ vec2d(0.4, 0.4), vec2d(0.4, 0.6), vec2d(0.6, 0.4), vec2d(0.6, 0.6) });
	/*PointsState st;
	constexpr int w = 10, h = 10;
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
			st.push_back(vec2d((i + .5) / (2 * w) + .25, (j + .5) / (2 * h) + .25));
	}*/
	std::uniform_int_distribution<> index_dist(0, point_count - 1);
	constexpr uint32_t tries = UINT32_MAX;
	constexpr uint32_t iterations = 4096;
	PointsState best_state = shape->gen_state(point_count, rgen);
	double best_score = MinDistEvaluator().evalv(best_state);
	for (uint32_t tr = 0; tr < tries; tr++)
	{
		PointsState st = shape->gen_state(point_count, rgen);
		for (uint32_t i = 0; i < iterations; i++)
		{
			const uint32_t move_index = index_dist(rgen);
			std::swap(st[move_index], st.back());
			st.pop_back();
			Voronoi<double> voron1(shape);
			voron1.init_sites(st);
			voron1.init_compute();
			const vec2d to_add = voron1.get_furthest().point;
			//cout << "adding " << to_add << endl;
			st.push_back(to_add);
			/*if (i % 32 == 0)
			{
				cout << "drawing " << i << endl;
				draw_voronoi(st, shape, "alg_test" + std::to_string(i) + ".png");
			}*/
		}
		const double score = MinDistEvaluator().evalv(st);
		cout << "try " << tr << '\n';
		cout << "score = " << score << ", best = " << best_score << endl;
		if (score > best_score)
		{
			best_score = score;
			best_state = st;
			draw_voronoi(best_state, shape, "test_found.png");
		}
	}
}

void Tester::run_tests()
{
	test_voronoi();
	//test_all();
	if (!failed)
		std::cout << "ALL TESTS PASSED" << std::endl;
}

}