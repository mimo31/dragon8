cmake_minimum_required(VERSION 3.16)

project(dragon8)

add_executable(${PROJECT_NAME}
	application.cpp
	dist-inv-solver.cpp
	evaluator.cpp
	graphics.cpp
	main.cpp
	min-dist-solver.cpp
	rescaling.cpp
	shape.cpp
	solver.cpp
	test/test.cpp
	vec2d.cpp
)

target_link_libraries(${PROJECT_NAME} -pthread -lX11)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
