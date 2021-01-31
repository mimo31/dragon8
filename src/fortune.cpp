/**
 * fortune.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/31
 */
#include "fortune.hpp"

#define DBG

#ifdef DBG
#include <cassert>
#endif

#include <iostream>

namespace dragon8
{

using std::cout;
using std::endl;

LineSegment::LineSegment(const vec2d v0, const vec2d v1)
	: v0(v0), v1(v1)
{
}

vec2d LineSegment::get_normal() const
{
	return (v1 - v0).get_lrot();
}

Fortune::BeachSec::BeachSec(const double x)
	: is_point(true), x(x)
{
}

Fortune::BeachSec::BeachSec(const double *t, const int32_t i, const int32_t li, const int32_t ri, const vec2d p, const vec2d l, const vec2d r)
: is_point(false), t(t), i(i), li(li), ri(ri), p(p), l(l), r(r)
{
}

bool Fortune::BeachSec::operator==(const BeachSec &other) const
{
	return is_point == other.is_point && x == other.x && i == other.i && li == other.li && ri == other.ri;
}

double Fortune::BeachSec::lx() const
{
	if (*t == p.y)
		return p.x;
	if (*t == l.y)
		return l.x;
	if (l.x == p.x)
	{
		const vec2d mid = l ^ p;
		const double d2 = (*t - mid.y) * (*t - mid.y);
		const double dx = sqrt(d2 - mid.dist2(l));
		return mid.x + (l.y > p.y ? dx : -dx);
	}
	if (l.y == p.y)
		return (p.x + l.x) / 2;
	const vec2d mid = l ^ p;
	const vec2d vp = (p - l).get_lrot();
	const vec2d v = (l.x < p.x) ^ (vp.y > 0) ? -vp : vp;
	const double mh = *t - mid.y;
	const double bp = mh * v.y;
	const double s = (-bp + sqrt(bp * bp - v.x * v.x * ((mid - p).len2() - mh * mh))) / (v.x * v.x);
	return mid.x + v.x * s;
}

double Fortune::BeachSec::rx() const
{
	if (*t == p.y)
		return p.x;
	if (*t == r.y)
		return r.x;
	if (r.x == p.x)
	{
		const vec2d mid = r ^ p;
		const double d2 = (*t - mid.y) * (*t - mid.y);
		const double dx = sqrt(d2 - mid.dist2(r));
		return mid.x + (r.y > p.y ? -dx : dx);
	}
	if (r.y == p.y)
		return (p.x + r.x) / 2;
	const vec2d mid = r ^ p;
	const vec2d vp = (p - r).get_lrot();
	const vec2d v = (r.x > p.x) ^ (vp.y > 0) ? -vp : vp;
	const double mh = *t - mid.y;
	const double bp = mh * v.y;
	const double s = (-bp + sqrt(bp * bp - v.x * v.x * ((mid - p).len2() - mh * mh))) / (v.x * v.x);
	return mid.x + v.x * s;
}

double Fortune::BeachSec::midx() const
{
	return (lx() + rx()) / 2;
}

bool Fortune::BeachSec::operator<(const BeachSec &other) const
{
	if (is_point && other.is_point)
		return x < other.x;
	if (is_point && !other.is_point)
		return other.li != -1 && x < other.lx();
	if (!is_point && other.is_point)
		return ri != -1 && rx() <= other.x;
	if ((li == -1 && other.li != -1) || (ri != -1 && other.ri == -1))
		return true;
	if (li == -1 || other.li == -1 || ri == -1 || other.ri == -1)
		return false;
	return midx() < other.midx();
}

bool Fortune::BeachSec::will_be_removed() const
{
	if (li == -1 || ri == -1)
		return false;
	if (p.y < l.y && p.y < r.y)
		return true;
	const vec2d v = r - l;
	const vec2d rot = v.get_lrot();
	return (p - l).dot(rot) < 0;
}

vec2d circumc(const vec2d a, const vec2d b, const vec2d c)
{
	const vec2d v1 = a - b, v2 = c - b;
	const double s1 = v1.dot(a ^ b), s2 = v2.dot(c ^ b);
	const double det = v1.x * v2.y - v1.y * v2.x;
	return vec2d(
		(s1 * v2.y - v1.y * s2) / det,
		(v1.x * s2 - s1 * v2.x) / det
	);
}

double Fortune::BeachSec::remove_time() const
{
	const vec2d joint = circumc(l, p, r);
	return joint.y + joint.dist(p);
}

double Fortune::BeachSec::y_at(const double x) const
{
	const vec2d base(x, *t);
	if (*t == p.y)
	{
		const vec2d pp = li != -1 ? l : r;
		const vec2d mid = base ^ pp;
		const vec2d v = (base - pp).get_lrot();
		const vec2d parpoint = mid + v * (x - mid.x) / v.x;
		return parpoint.y;
	}
	const vec2d mid = base ^ p;
	const vec2d v = (base - p).get_lrot();
	const vec2d parpoint = mid + v * (x - mid.x) / v.x;
	return parpoint.y;
}

Fortune::Event::Event(const double time, const BeachSec sec, const vec2d hit_point)
	: type(Fortune::EventType::BORDER_HIT), time(time), sec(sec), hit_point(hit_point)
{
	#ifdef DBG
	assert(!std::isnan(time));
	#endif
}

Fortune::Event::Event(const double time, const BeachSec sec)
	: type(Fortune::EventType::REMOVE_SEC), time(time), sec(sec)
{
	#ifdef DBG
	assert(!std::isnan(time));
	#endif
}

Fortune::Event::Event(const double time, const uint32_t pi)
	: type(Fortune::EventType::NEW_POINT), time(time), pi(pi)
{
	#ifdef DBG
	assert(!std::isnan(time));
	#endif
}

bool Fortune::Event::operator<(const Fortune::Event &other) const
{
	return time > other.time;
}

void Fortune::submit_point(const vec2d p, const double d2, const bool border)
{
	#ifdef DBG
	for (const vec2d q : points)
	{
		assert(!(p == q));
	}
	#endif
	pots.push_back(p);
	std::cout << "pot at " << p << std::endl;
	if (d2 > dist2 || dist2 == -1)
	{
		if (border || container->is_inside(p))
		{
			dist2 = d2;
			furthest = p;
		}
	}
}

/*void Fortune::submit_point_if_border(const vec2d p, const double d2)
{
	constexpr double eps = 1e-13;

	for (const LineSegment &ls : segs)
	{
		if (ls.v0.dist2(p) < eps * eps || ls.v1.dist2(p) < eps * eps)
		{
			submit_point(p, d2, true);
			continue;
		}
		const vec2d dir = ls.v1 - ls.v0;
		const double s0 = dir.dot(ls.v0);
		const double s1 = dir.dot(ls.v1);
		if (dir.dot(p) >= s0 && dir.dot(p) <= s1 && abs(ls.get_normal().dot(ls.v0 - p)) < eps)
			submit_point(p, d2, true);
	}
}*/

void Fortune::check_beach_integrity() const
{
	if (beach.size() == 0)
		return;
	assert(beach.begin()->li == -1);
	assert(std::prev(beach.end())->ri == -1);
	auto curit = beach.begin();
	auto nxtit = std::next(curit);
	for (; nxtit != beach.end(); curit++, nxtit++)
	{
		assert(curit->ri == nxtit->i);
		assert(curit->i == nxtit->li);
	}
}

void Fortune::print_beach() const
{
	std::cout << "(t=" << t << ")";
	for (const BeachSec &bs : beach)
	{
		std::cout << "[" << (bs.li != -1 ? bs.lx() : 0) << ", " << bs.li << ", " << bs.i << ", " << bs.ri << ", " << (bs.ri != -1 ? bs.rx() : 0) << "]";
	}
	std::cout << std::endl;
}

void Fortune::add_border_line(const BeachSec &left)
{
	const vec2d v = left.r - left.p;
	const vec2d dir = v.get_lrot();
	const vec2d startp(left.rx(), left.y_at(left.rx()));
	cout << "r = " << left.r << ", p = " << left.p << endl;
	const double s = v.dot(left.p ^ left.r);
	for (const LineSegment &ls : segs)
	{
		const double d0 = v.dot(ls.v0), d1 = v.dot(ls.v1);
		if ((d0 < s && s < d1) || (d1 < s && s < d0))
		{
			const vec2d n2 = ls.get_normal();
			const double s2 = n2.dot(ls.v0);
			const double det = v.x * n2.y - v.y * n2.x;
			const vec2d inter(
				(s * n2.y - v.y * s2) / det,
				(v.x * s2 - s * n2.x) / det);
			auto bit = beach.lower_bound(BeachSec(inter.x));
			const double yalready = bit->y_at(inter.x);
			if (yalready < inter.y && dir.dot(inter) > dir.dot(startp))
			{
				const double time = inter.y + left.p.dist(inter);
				#ifdef DBG
				assert(time >= t);
				#endif
				/*cout << "inter = " << inter << endl;
				cout << "!!! " << left.p.dist(inter) << " == " << left.r.dist(inter) << endl;
				const double deltat = time - t;
				BeachSec cpy = left;
				double tm = time;
				cpy.t = &tm;
				cout << "??? " << cpy.rx() << " == " << inter.x << " and " << cpy.y_at(cpy.rx()) << " == " << inter.y << endl;
				tm = t;
				for (uint32_t i = 0; i < 5; i++)
				{
					cout << "t = " << tm << endl;
					const double ftx = cpy.rx(), fty = cpy.y_at(ftx);
					cout << "ftx = " << ftx << ", fty = " << fty << endl;
					cout << "cal t = " << fty + cpy.p.dist(vec2d(ftx, fty)) << endl;
					tm += .05;
					//cout << cpy.rx() << ", " << cpy.y_at(cpy.rx()) << endl;
				}*/
				events.push(Event(time, left, inter));
			}
		}
	}
}

void Fortune::resolve_top()
{
	const Event ev = events.top();
	events.pop();
	std::cout << "processing event at t = " << ev.time << std::endl;
	if (ev.type == EventType::NEW_POINT)
	{
		t = ev.time;
		const uint32_t pi = ev.pi;
		const vec2d p = points[pi];

		std::cout << "will add point " << p << "(" << pi << ")" << std::endl;
		std::cout << "before add: ";
		print_beach();

		constexpr double eps = 1e-13;
		auto bit = beach.lower_bound(BeachSec(p.x));
		if (bit->ri != -1 && fabs(bit->rx() - p.x) < eps)
			bit = std::next(bit);
		const BeachSec orig = *bit;
		if (orig.li != -1 && /*orig.lx() == p.x)*/fabs(orig.lx() - p.x) < eps)
		{
			const auto prvit = std::prev(bit);
			const BeachSec prv = *prvit;
			beach.erase(bit);
			beach.erase(prvit);
			const BeachSec newleft(&t, prv.i, prv.li, pi, prv.p, prv.l, p);
			const BeachSec newcenter(&t, pi, prv.i, orig.i, p, prv.p, orig.p);
			const BeachSec newright(&t, orig.i, pi, orig.ri, orig.p, p, orig.r);
			beach.insert(newleft);
			beach.insert(newcenter);
			beach.insert(newright);
			#ifdef DBG
			check_beach_integrity();
			#endif
			const vec2d inter(p.x, newleft.y_at(p.x));
			submit_point(inter, inter.dist2(p), false);
			add_border_line(newleft);
			add_border_line(newcenter);
			if (newright.ri != -1)
				add_border_line(newright);
			if (newleft.will_be_removed())
				events.push(Event(newleft.remove_time(), newleft));
			if (newright.will_be_removed())
				events.push(Event(newright.remove_time(), newright));
		}
		else
		{
			beach.erase(bit);
			const BeachSec newleft(&t, orig.i, orig.li, pi, orig.p, orig.l, p);
			const BeachSec newcenter(&t, pi, orig.i, orig.i, p, orig.p, orig.p);
			const BeachSec newright(&t, orig.i, pi, orig.ri, orig.p, p, orig.r);
			/*const auto lit = */ beach.insert(newleft);
			beach.insert(newcenter);
			beach.insert(newright);
			#ifdef DBG
			check_beach_integrity();
			#endif
			/*if (newleft.li != -1)
				add_border_line(*std::prev(lit.first));*/
			add_border_line(newleft);
			add_border_line(newcenter);
			if (newright.ri != -1)
				add_border_line(newright);
			if (newleft.will_be_removed())
				events.push(Event(newleft.remove_time(), newleft));
			if (newright.will_be_removed())
				events.push(Event(newright.remove_time(), newright));
			//const vec2d inter(p.x, newleft.y_at(p.x));
			//submit_point_if_border(inter, inter.dist2(p));
		}
		std::cout << "after add: ";
		print_beach();
	}
	else if (ev.type == EventType::REMOVE_SEC)
	{
		const vec2d crcc = circumc(ev.sec.l, ev.sec.p, ev.sec.r);
		std::cout << "possible remove arc [" << ev.sec.li << ", " << ev.sec.i << ", " << ev.sec.ri << "] " << crcc << std::endl;
		const auto bit = beach.find(ev.sec);
		if (bit != beach.end())
		{
			std::cout << "arc exists -- removing" << std::endl;
			submit_point(crcc, crcc.dist2(ev.sec.p), false);
			BeachSec prv = *std::prev(bit);
			BeachSec nxt = *std::next(bit);
			beach.erase(std::next(bit));
			beach.erase(std::prev(bit));
			beach.erase(bit);
			prv.ri = nxt.i;
			prv.r = nxt.p;
			nxt.li = prv.i;
			nxt.l = prv.p;
			/*const auto prvit = */ beach.insert(prv);
			beach.insert(nxt);
			std::cout << "after removal: ";
			print_beach();
			#ifdef DBG
			check_beach_integrity();
			#endif

			t = ev.time;

			/*if (prv.li != -1)
					add_border_line(*std::prev(prvit.first));*/
			add_border_line(prv);
			if (nxt.ri != -1)
				add_border_line(nxt);
			if (prv.will_be_removed())
				events.push(Event(prv.remove_time(), prv));
			if (nxt.will_be_removed())
				events.push(Event(nxt.remove_time(), nxt));
		}
		std::cout << "after time update: ";
		print_beach();
	}
	else if (ev.type == EventType::BORDER_HIT)
	{
		#ifdef DBG
		assert(ev.time >= t);
		#endif
		t = ev.time;
		const auto bit = beach.find(ev.sec);
		if (bit != beach.end() && ev.sec == *bit)
		{
			cout << "border hit" << endl;
			submit_point(ev.hit_point, ev.hit_point.dist2(ev.sec.p), true);
			cout << "p = " << ev.hit_point << endl;
			cout << "left i = " << ev.sec.i << ", right i = " << ev.sec.ri << endl;
			const auto psea = beach.lower_bound(BeachSec(ev.hit_point.x));
			cout << "psea i = " << psea->i << ", par y = " << psea->y_at(ev.hit_point.x) << endl;
			cout << "left-right border x = " << ev.sec.rx() << ", y = " << ev.sec.y_at(ev.sec.rx()) << endl;
		}
	}
	std::cout << std::endl;
}

void Fortune::add_firsts()
{
	vec<uint32_t> firsts;
	const double y0 = points[events.top().pi].y;
	firsts.push_back(events.top().pi);
	events.pop();
	while (!firsts.empty() && points[events.top().pi].y == y0)
	{
		firsts.push_back(events.top().pi);
		events.pop();
	}
	std::sort(firsts.begin(), firsts.end(), [this](const uint32_t i0, const uint32_t i1) { return points[i0].x < points[i1].x; });
	// TODO: add border lines
	for (uint32_t i = 0; i < firsts.size(); i++)
	{
		const int32_t prevind = i != 0 ? firsts[i - 1] : -1;
		const int32_t nextind = i < firsts.size() - 1 ? firsts[i + 1] : -1;
		const vec2d prv = prevind != -1 ? points[prevind] : vec2d();
		const vec2d nxt = nextind != -1 ? points[nextind] : vec2d();
		beach.insert(BeachSec(&t, firsts[i], prevind, nextind, points[firsts[i]], prv, nxt));
	}

	#ifdef DBG
	check_beach_integrity();
	#endif
}

Fortune::Fortune(const ShapePtr container, const vec<vec2d> &points, const vec<LineSegment> &segs)
	: container(container), points(points), segs(segs)
{
	for (uint32_t i = 0; i < points.size(); i++)
		events.push(Event(points[i].y, i));

	add_firsts();
	
	while (!events.empty())
		resolve_top();
}

}