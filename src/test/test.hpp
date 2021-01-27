/**
 * test.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/27
 */
#ifndef TEST_HPP
#define TEST_HPP

namespace dragon8
{

class Tester
{
private:
	bool failed = false;
	void test_shape_polygon();
	void test_all();

public:
	void run_tests();
};

}

#endif // TEST_HPP