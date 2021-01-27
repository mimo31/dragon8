/**
 * application.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#ifndef APPLICATION_HPP
#define APPLICATION_HPP

#include "optimalization-type.hpp"
#include "shape.hpp"

namespace dragon8
{

enum NumberParseResult
{
	OK, TOO_BIG, SYNTAX_ERROR
};

NumberParseResult to_nneg_int(const std::string& s, uint32_t& val);
bool has_interior(const std::vector<vec2d>& verts);

class Application
{
private:
	ShapePtr shape;
	bool print_coors;
	bool generate_image;
	std::string image_filename;
	OptimalizationType opt_type;

	const uint32_t MAX_N = 1'000'000;
	uint32_t n;

	void read_n(std::istream&, std::ostream&);
	void read_shape(std::istream&, std::ostream&);
	void read_opt_type(std::istream&, std::ostream&);
	void read_output_spec(std::istream&, std::ostream&);
	void read_input(std::istream&, std::ostream&);

public:
	void run();
};

}

#endif // APPLICATION_HPP