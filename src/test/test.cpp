/**
 * test.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/27
 */
#include "test.hpp"

#include <iostream>

#include "application.hpp"
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

void Tester::test_voronoi()
{
	Voronoi<double> voronoi(nullptr);
	/*const vec<vec2<double>> sites0({ vec2<double>(0, 0), vec2<double>(1, 1), vec2<double>(4, 3) });
	voronoi.init_sites(sites0);
	voronoi.init_compute();*/
	const std::vector<vec2d> verts{ vec2d(0, 0), vec2d(1, 0), vec2d(1, 1), vec2d(0, 1) };
	const ShapePtr unit_square = std::make_shared<ShapePolygon>(verts);// unit_square(verts);
	std::random_device dev;
	std::mt19937 rgen(dev());
	
	//const PointsState st = unit_square->gen_state(10000, rgen);
	/*for (const vec2d p : st)
		cout << p << endl;
	*/
	//const PointsState st({ vec2d(0.460016, 0.802885), vec2d(0.14821, 0.672051), vec2d(0.77762, 0.170662), vec2d(0.958275, 0.894059) });
	//const PointsState st({ vec2d(0.4, 0.4), vec2d(0.4, 0.6), vec2d(0.6, 0.4), vec2d(0.6, 0.6) });
	PointsState st;
	constexpr int w = 3, h = 2;
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
			st.push_back(vec2d((i + .5) / (2 * w) + .25, (j + .5) / (2 * h) + .25));
	}
	Voronoi<double> voron1(nullptr);
	voron1.init_sites(st);
	voron1.init_compute();
	write_image(ShapePtr(unit_square), st, "voronoi_test.png", 0, voron1.get_vertices());
}

void Tester::run_tests()
{
	test_voronoi();
	//test_all();
	if (!failed)
		std::cout << "ALL TESTS PASSED" << std::endl;
}

}