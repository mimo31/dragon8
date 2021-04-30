/**
 * main.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include <cstdlib>

#include "application.hpp"

#define TESTING

#ifdef TESTING
#include "test/test.hpp"
#endif

int main()
{
	#ifdef TESTING
	dragon8::Tester().run_tests();
	#else
	dragon8::Application app;
	app.run();
	#endif
	return EXIT_SUCCESS;
}
