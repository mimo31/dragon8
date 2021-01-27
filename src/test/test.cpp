/**
 * test.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/27
 */
#include "test.hpp"

#include <iostream>

#include "../shape.hpp"

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

void Tester::test_all()
{
	test_shape_polygon();
}

void Tester::run_tests()
{
	test_all();
	if (!failed)
		std::cout << "ALL TESTS PASSED" << std::endl;
}

}