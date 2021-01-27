/**
 * application.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/23
 */
#include "application.hpp"

#include <iostream>
#include <memory>
#include <sstream>

#include "dist-inv-solver.hpp"
#include "evaluator.hpp"
#include "graphics.hpp"
#include "min-dist-solver.hpp"
#include "points-state.hpp"

namespace dragon8
{

void Application::run()
{
	read_input(std::cin, std::cout);

	std::cout.precision(15);

	const PointsState ps = opt_type == OptimalizationType::MIN_DIST_INV ? DistInvSolver(shape, n).solve() : MinDistSolver(shape, n).solve();
	const double score = opt_type == OptimalizationType::MIN_DIST_INV ? DistInvEvaluator().evalv(ps) : MinDistEvaluator().evalv(ps);

	std::cout << "score = " << score << std::endl;

	if (print_coors)
	{
		std::cout << "best configuration found:\n";
		for (const vec2d v : ps)
		{
			std::cout << v.x << " " << v.y << '\n';
		}
		std::cout << std::endl;
	}

	if (generate_image)
		write_image(shape, ps, "out.png");
}

NumberParseResult to_nneg_int(const std::string& s, uint32_t& val)
{
	uint32_t frst = 0;
	while (frst < s.length() && isspace(s[frst]))
		frst++;
	uint32_t lst = s.length();
	while (lst >= 1 && isspace(s[lst - 1]))
		lst--;
	if (lst <= frst)
		return NumberParseResult::SYNTAX_ERROR;
	while (frst < s.length() && s[frst] == '0')
		frst++;
	for (uint32_t i = frst; i < lst; i++)
		if (s[i] < '0' || s[i] > '9')
			return NumberParseResult::SYNTAX_ERROR;
	if (lst - frst > 9)
		return NumberParseResult::TOO_BIG;
	val = 0;
	for (uint32_t i = frst; i < lst; i++)
	{
		val *= 10;
		val += s[i] - '0';
	}
	return NumberParseResult::OK;
}

std::string getln(std::istream& is)
{
	std::string ln;
	getline(is, ln);
	return ln;
}

bool has_interior(const std::vector<vec2d>& verts)
{
	if (verts.size() < 3)
		return false;
	const uint32_t n = verts.size();
	double area = (verts[0].x - verts[n - 1].x) * (verts[0].y + verts[n - 1].y);
	for (uint32_t i = 0; i < n - 1; i++)
		area += (verts[i + 1].x - verts[i].x) * (verts[i + 1].y + verts[i].y);
	return area != 0;
}

using std::endl;

typedef std::string str;

bool get_bool(const str& descr, std::istream& is, std::ostream& os)
{
	while (true)
	{
		os << descr << endl;
		const str ln = getln(is);
		if (ln.size() != 0 && (ln[0] == 'y' || ln[0] == 'Y'))
		{
			os << "You chose YES." << endl;
			return true;
		}
		else if (ln.size() != 0 && (ln[0] == 'n' || ln[0] == 'N'))
		{
			os << "You chose NO." << endl;
			return false;
		}
		os << "Enter 'y' (for YES) or 'n' (for NO)." << endl;
	}
}

void Application::read_n(std::istream& is, std::ostream& os)
{
	while (true)
	{
		os << "How many points?" << endl;
		const str ln = getln(is);
		const NumberParseResult pres = to_nneg_int(ln, n);
		if (pres == NumberParseResult::OK)
			break;
		else if (pres == NumberParseResult::TOO_BIG)
			os << "Only " << MAX_N << " points are allowed at most." << endl;
		else if (pres == NumberParseResult::SYNTAX_ERROR)
			os << "Invalid number." << endl;
	}
}

struct ShapePreset
{
	str name;
	ShapePtr shape;

	ShapePreset(const str name, const ShapePtr shape)
		: name(name), shape(shape)
	{
	}
};

std::vector<vec2d> generate_regular_ngon_verts(const uint32_t n)
{
	std::vector<vec2d> verts;
	for (uint32_t i = 0; i < n; i++)
		verts.push_back(vec2d(sin(i * 2 * M_PI / n), cos(i * 2 * M_PI / n)));
	return verts;
}

std::vector<vec2d> generate_comb_verts()
{
	std::vector<vec2d> verts;
	constexpr uint32_t points = 20;
	for (uint32_t i = 0; i < points; i++)
	{
		verts.push_back(vec2d(i / double(points), .1));
		verts.push_back(vec2d((i + .5) / double(points), .9));
	}
	verts.push_back(vec2d(1, .1));
	verts.push_back(vec2d(1, 0));
	verts.push_back(vec2d(0, 0));
	return verts;
}

const uint32_t preset_count = 10;
const ShapePreset shape_presets[preset_count] = 
{
	ShapePreset("unit circle (r = 1)", std::make_shared<ShapeCircle>(1.0)),
	ShapePreset("unit square", std::make_shared<ShapePolygon>(std::vector<vec2d>{ vec2d(0, 0), vec2d(1, 0), vec2d(1, 1), vec2d(0, 1) })),
	ShapePreset("equilateral triangle", std::make_shared<ShapePolygon>(std::vector<vec2d>{ vec2d(0, 0), vec2d(1, 0), vec2d(.5, sqrt(3) / 2) })),
	ShapePreset("right isosceles triangle", std::make_shared<ShapePolygon>(std::vector<vec2d>{ vec2d(0, 0), vec2d(1, 0), vec2d(0, 1) })),
	ShapePreset("1.0 x .1 rectangle", std::make_shared<ShapePolygon>(std::vector<vec2d>{ vec2d(0, 0), vec2d(1, 0), vec2d(1, .1), vec2d(0, .1) })),
	ShapePreset("\"bent\" quadrilateral", std::make_shared<ShapePolygon>(std::vector<vec2d>{ vec2d(0, 0), vec2d(1, 0), vec2d(.2, .2), vec2d(0, 1) })),
	ShapePreset("regular pentagon", std::make_shared<ShapePolygon>(generate_regular_ngon_verts(5))),
	ShapePreset("regular hexagon", std::make_shared<ShapePolygon>(generate_regular_ngon_verts(6))),
	ShapePreset("4-pointed star", std::make_shared<ShapePolygon>(std::vector<vec2d>{
		vec2d(0, 1),
		vec2d(.2, .2),
		vec2d(1, 0),
		vec2d(.2, -.2),
		vec2d(0, -1),
		vec2d(-.2, -.2),
		vec2d(-1, 0),
		vec2d(-.2, .2)
	})),
	ShapePreset("20-pointed comb", std::make_shared<ShapePolygon>(generate_comb_verts()))
};

void Application::read_shape(std::istream& is, std::ostream& os)
{
	bool custom_circle = false;
	bool custom_polygon = false;
	while (true)
	{
		os << "Choose the container." << endl;
		for (uint32_t i = 0; i < preset_count; i++)
			os << "[" << (i + 1) << "] " << shape_presets[i].name << std::endl;
		os << "[" << preset_count + 1 << "] custom circle" << std::endl;
		os << "[" << preset_count + 2 << "] custom polygon" << std::endl;
		const str ln = getln(is);
		uint32_t selected;
		const NumberParseResult res = to_nneg_int(ln, selected);
		if (res == OK && selected != 0 && selected <= preset_count + 2)
		{
			if (selected <= preset_count)
				shape = shape_presets[selected - 1].shape;
			else if (selected == preset_count + 1)
				custom_circle = true;
			else
				custom_polygon = true;
			break;
		}
		os << "Enter an integer between 1 and " << preset_count << "." << endl;
	}

	if (custom_circle)
	{
		os << "You chose custom circle." << endl;
		double r;
		while (true)
		{
			os << "Enter the radius. (Leave blank for default r = 1.0.)" << endl;
			const str ln = getln(is);
			if (ln.length() != 0)
			{
				std::stringstream ss(ln);
				ss >> r;
				ss >> std::ws;
				if (ss.fail() || !ss.eof())
					os << "Invalid number." << endl;
				else
					break;
			}
			else
			{
				r = 1.0;
				break;
			}
		}
		shape = std::make_shared<ShapeCircle>(r);
	}
	else if (custom_polygon)
	{
		os << "You now need to specify the coordinates of the vertices of the custom polygon." << endl;
		while (true)
		{
			std::vector<vec2d> points;
			while (true)
			{
				const bool point_required = points.size() <= 2;
				const str note = point_required ? "" : " or leave blank to finish";
				os << "Enter the coordinates (format \"x y\") of vertex " << points.size() + 1 << note << "." << endl;
				const str ln = getln(is);
				if (ln.length() == 0 && !point_required)
					break;
				std::stringstream ss(ln);
				double x, y;
				ss >> x >> y;
				ss >> std::ws;
				if (ss.fail() || !ss.eof())
					os << "Invalid input." << endl;
				else
				{
					points.push_back(vec2d(x, y));
				}
			}
			if (points.size() != 0)
			{
				if (has_interior(points))
				{
					shape = std::make_shared<ShapePolygon>(points);
					break;
				}
				else
					os << "Your polygon has no interior. Starting over." << endl;
			}
		}
	}
}

void Application::read_opt_type(std::istream& is, std::ostream& os)
{
	while (true)
	{
		os << "What do you want to optimize?" << endl;
		os << "[1] minimize sum 1 / d_{i,j}" << endl;
		os << "[2] maximize min d_{i,j}" << endl;
		const str ln = getln(is);
		if (ln.size() != 0 && ln[0] == '1')
		{
			opt_type = OptimalizationType::MIN_DIST_INV;
			break;
		}
		else if (ln.size() != 0 && ln[0] == '2')
		{
			opt_type = OptimalizationType::MAX_MIN_DIST;
			break;
		}
	}
}

void Application::read_output_spec(std::istream& is, std::ostream& os)
{
	print_coors = get_bool("Do you want to have the result printed to console (y/n)?", is, os);

	generate_image = get_bool("Do you want to have the result saved to an image file (y/n)?", is, os);

	if (!print_coors && !generate_image)
		os << "WARNING: You have no way of getting the result." << endl;

	if (generate_image)
	{
		os << "Enter the output image filename (filename will be \"{what you enter}.png\"). (Leave blank for default \"out.png\".)" << endl;
		str ln = getln(is);
		if (ln.size() == 0)
			ln = "out";
		image_filename = ln + ".png";
	}
}

void Application::read_input(std::istream& is, std::ostream& os)
{
	read_n(is, os);
	read_shape(is, os);
	read_opt_type(is, os);
	read_output_spec(is, os);
}

}