cmake_minimum_required(VERSION 3.16)

project(dragon8)

add_executable(${PROJECT_NAME}
	application.cpp
	dist-inv-solver.cpp
	evaluator.cpp
	main.cpp
	shape.cpp
	solver.cpp
)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
