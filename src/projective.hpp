/**
 * projective.hpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/05/03
 */
#ifndef PROJECTIVE_HPP
#define PROJECTIVE_HPP

#include "vec2.hpp"

namespace dragon8
{

template<typename float_t>
struct ProjectiveVertex
{
	bool at_infinity;
	vec2<float_t> v;

	ProjectiveVertex() = default;
	ProjectiveVertex(const bool at_infinity, const vec2<float_t> v) : at_infinity(at_infinity), v(v)
	{
	}

};

template<typename float_t>
struct ProjectiveEdge
{
	ProjectiveVertex<float_t> v0, v1;

	ProjectiveEdge() = default;
	ProjectiveEdge(const ProjectiveVertex<float_t> &v0, const ProjectiveVertex<float_t> &v1) : v0(v0), v1(v1)
	{
	}
};

}

#endif // PROJECTIVE_HPP