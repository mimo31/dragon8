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
#include "points-state.hpp"

namespace dragon8
{

void Application::run()
{
	//const ShapePtr shape = std::make_shared<ShapeCircle>(1);
	//const uint32_t n = 10000;

	read_input(std::cin, std::cout);

	std::cout.precision(15);

	const PointsState ps = DistInvSolver(shape, n).solve();

	std::cout << "score = " << DistInvEvaluator().evalv(ps) << std::endl;

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

enum NumberParseResult
{
	OK, TOO_BIG, SYNTAX_ERROR
};

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

bool intersects(const vec2d v0, const vec2d v1, const vec2d w0, const vec2d w1)
{
	const vec2d v = v1 - v0, w = w1 - v0;
	const vec2d w0v = w0 - v0, w1v = w1 - v0;
	if (v.cross(w) == 0 && v.cross(w0 - v0) == 0)
	{
		const double v1d = v.dot(v), w0d = v.dot(w0v), w1d = v.dot(w1v);
		const double wdfrom = std::min(w0d, w1d), wdto = std::max(w0d, w1d);
		return wdfrom <= v1d && 0 <= wdto;
	}
	const double w0c = v.cross(w0v), w1c = v.cross(w1v);
	if ((w0c > 0 && w1c > 0) || (w0c < 0 && w1c < 0))
		return false;
	const vec2d v0w = v0 - w0, v1w = v1 - w0;
	const double v0c = w.cross(v0w), v1c = v.cross(v1w);
	if ((v0c > 0 && v1c > 0) || (v0c < 0 && v1c < 0))
		return false;
	return true;
}

bool vertices_distinct(std::vector<vec2d>& verts)
{
	for (uint32_t i = 0; i < verts.size(); i++)
	{
		for (uint32_t j = i; j < verts.size(); j++)
		{
			if (verts[i] == verts[j])
				return false;
		}
	}
	return true;
}

bool has_no_intersections_open(std::vector<vec2d>& verts)
{
	if (!vertices_distinct(verts))
		return false;
	for (uint32_t i = 0; i < verts.size() - 1; i++)
	{
		for (uint32_t j = i + 1; j < verts.size() - 1; j++)
		{
			if (intersects(verts[i], verts[i + 1], verts[j], verts[j + 1]))
				return false;
		}
	}
	return true;
}

bool has_no_intersections_closing(std::vector<vec2d>& verts)
{
	const vec2d v0 = verts[0], v1 = verts[verts.size() - 1];
	for (uint32_t i =  0; i < verts.size() - 1; i++)
	{
		if (intersects(v0, v1, verts[i], verts[i + 1]))
			return false;
	}
	return true;
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

void Application::read_shape(std::istream& is, std::ostream& os)
{
	bool circle;
	while (true)
	{
		os << "Distributing points inside a circle or a polygon? (c/p)" << endl;
		const str ln = getln(is);
		if (ln.length() != 0)
		{
			const char frst = ln[0];
			if (frst == 'c' || frst == 'C')
			{
				circle = true;
				break;
			}
			else if (frst == 'p' || frst == 'P')
			{
				circle = false;
				break;
			}
		}
		os << "Enter 'c' (for circle) or 'p' (for polygon)." << endl;
	}

	os << "You chose " << (circle ? "circle" : "polygon") << "." << endl;

	if (circle)
	{
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
	else
	{
		while (true)
		{
			std::vector<vec2d> points;
			while (true)
			{
				const bool point_required = points.size() <= 2;
				const str note = point_required ? " or leave blank to finish" : "";
				os << "Enter the coorindates (format \"x y\") of vertex " << points.size() + 1 << note << "." << endl;
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
					if (!has_no_intersections_open(points))
					{
						os << "Your polygon was intersecting itself. Starting over." << endl;
						points.clear();
						break;
					}
				}
			}
			if (points.size() != 0)
			{
				if (has_no_intersections_closing(points))
				{
					shape = std::make_shared<ShapePolygon>(points);
					break;
				}
				else
					os << "Your polygon was intersecting itself. Starting over." << endl;
			}
		}
	}
}

void Application::read_opt_type(std::istream& is, std::ostream& os)
{
	while (true)
	{
		os << "What do you want to optimize?" << endl;
		os << "[1] minimize the sum 1 / d_{i,j}" << endl;
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