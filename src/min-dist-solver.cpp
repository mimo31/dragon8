/**
 * min-dist-solver.cpp
 * 
 * Author: Viktor Fukala
 * Created on 2021/1/26
 */
#include "min-dist-solver.hpp"

#include <iostream>
#include <limits>
#include <queue>
#include <set>

#include "evaluator.hpp"
#include "graphics.hpp"

template <typename T>
using vec = std::vector<T>;

using std::cout;
using std::endl;

namespace dragon8
{

MinDistSolver::MinDistSolver(const std::shared_ptr<Shape> shape, const uint32_t n)
	: Solver(shape, n)
{
}

void MinDistSolver::shift_solve(PointsState& ps, const uint32_t iters, RGen& rgen) const
{
	double bestscore = MinDistEvaluator().evalv(ps);
	PointsState beststate = ps;
	for (uint32_t it = 0; it < iters; it++)
	{
		/*double leastd2 = std::numeric_limits<double>::max();
		uint32_t li, lj;
		for (uint32_t i = 0; i < n; i++)
		{
			for (uint32_t j = i + 1; j < n; j++)
			{
				const double d2 = (ps[i] - ps[j]).len2();
				if (d2 < leastd2)
				{
					leastd2 = d2;
					li = i;
					lj = j;
				}
			}
		}
		std::uniform_int_distribution<> dist(0, 1);
		if (dist(rgen))
			std::swap(li, lj);*/
		
		std::uniform_int_distribution<> dist(0, n - 1);
		const uint32_t lj = dist(rgen);
		double mdist2 = std::numeric_limits<double>::max();
		uint32_t li;
		for (uint32_t i = 0; i < n; i++)
		{
			if (i == lj)
				continue;
			const double dist2 = (ps[lj] - ps[i]).len2();
			if (dist2 < mdist2)
			{
				mdist2 = dist2;
				li = i;
			}
		}

		double mxmove = std::numeric_limits<double>::max();
		int32_t limiting = -1;
		const vec2d s0 = ps[li], s1 = ps[lj], s = s1 - s0;
		for (uint32_t i = 0; i < n; i++)
		{
			const vec2d p = ps[i];
			if (s.dot(p - s1) <= 0)
				continue;
			const double moveavail = .5 * (p - s0).len2() / s.dot(p - s0) - 1;
			if (moveavail < mxmove)
			{
				mxmove = moveavail;
				limiting = i;
			}
		}
		if (limiting == -1)
			mxmove = .5;
		ps[lj] = shape->bound(s1, s1 + s * mxmove);
		/*if (limiting != -1)
		{
			mxmove = std::numeric_limits<double>::max();
			const vec2d sp0 = (ps[li] + ps[limiting]) / 2, sp1 = ps[lj], sp = sp1 - sp0;
			for (uint32_t i = 0; i < n; i++)
			{
				const vec2d p = ps[i];
				if (sp.dot(p - sp1) <= 0)
					continue;
				const double moveavail = ((p + s0) / 2 - sp1).dot(p - s0) / sp.dot(p - s0);
				mxmove = std::min(mxmove, moveavail);
			}
			if (mxmove == std::numeric_limits<double>::max())
				mxmove = .5;
			ps[lj] = shape->bound(sp1, sp1 + sp * mxmove);
		}
		const double curscore = MinDistEvaluator().evalv(ps);
		if (curscore > bestscore)
		{
			bestscore = curscore;
			beststate = ps;
		}*/
	}
	//ps = beststate;
}

void MinDistSolver::jump_solve(PointsState& ps, const uint32_t iters, RGen& rgen) const
{
	for (uint32_t it = 0; it < iters; it++)
	{
		double leastd2 = std::numeric_limits<double>::max();
		uint32_t li, lj;
		for (uint32_t i = 0; i < n; i++)
		{
			for (uint32_t j = i + 1; j < n; j++)
			{
				const double d2 = (ps[i] - ps[j]).len2();
				if (d2 < leastd2)
				{
					leastd2 = d2;
					li = i;
					lj = j;
				}
			}
		}
		std::uniform_int_distribution<> dist(0, 1);
		while (it < iters)
		{
			if (dist(rgen))
				std::swap(li, lj);
			const vec2d np = shape->gen_point(rgen);
			bool better = true;
			for (uint32_t i = 0; i < n; i++)
			{
				if (i == li)
					continue;
				if ((np - ps[i]).len2() <= leastd2)
				{
					better = false;
					break;
				}
			}
			if (better)
			{
				ps[li] = np;
				break;
			}
			it++;
		}
	}
}

bool MinDistSolver::fit_enough(const vec2d p0, const double d) const
{
	const rectangle2d box = shape->get_box();
	const double h = sqrt(3) / 2 * d;
	uint32_t cou = 0;
	int32_t yi0 = 0;
	while ((yi0 - 1) * h + p0.y >= box.c0.y)
		yi0--;
	int32_t xi0 = 0;
	while ((xi0 - 1) * d + p0.x >= box.c0.x)
		xi0--;
	xi0--;
	for (int32_t yi = yi0; yi * h + p0.y <= box.c1.y; yi++)
	{
		const double sh = yi % 2 ? 0 : d / 2;
		for (int32_t xi = xi0; xi * d + sh + p0.x <= box.c1.x; xi++)
		{
			const vec2d p(xi * d + sh + p0.x, yi * h + p0.y);
			if (shape->is_inside(p))
			{
				cou++;
				if (cou >= n)
				{
					return true;
				}
			}
		}
	}
	return false;
}

PointsState MinDistSolver::get_fit(const vec2d p0, const double d) const
{
	const rectangle2d box = shape->get_box();
	const double h = sqrt(3) / 2 * d;
	PointsState points;
	int32_t yi0 = 0;
	while ((yi0 - 1) * h + p0.y >= box.c0.y)
		yi0--;
	int32_t xi0 = 0;
	while ((xi0 - 1) * d + p0.x >= box.c0.x)
		xi0--;
	xi0--;
	for (int32_t yi = yi0; yi * h + p0.y <= box.c1.y; yi++)
	{
		const double sh = yi % 2 ? 0 : d / 2;
		for (int32_t xi = xi0; xi * d + sh + p0.x <= box.c1.x; xi++)
		{
			const vec2d p(xi * d + sh + p0.x, yi * h + p0.y);
			if (shape->is_inside(p))
			{
				points.push_back(p);
				if (points.size() >= n)
					return points;
			}
		}
	}
	return points;
}

PointsState MinDistSolver::regular_fit(const vec2d p0, double& score) const
{
	double d = shape->get_box().get_width();
	double dm, dp;
	if (fit_enough(p0, d))
	{
		while (fit_enough(p0, d * 2))
			d *= 2;
		dp = d * 2;
		dm = d;
	}
	else
	{
		while (!fit_enough(p0, d / 2))
			d /= 2;
		dp = d;
		dm = d / 2;
	}
	for (uint32_t i = 0; i < 50; i++)
	{
		const double mid = (dm + dp) / 2;
		if (fit_enough(p0, mid))
			dm = mid;
		else
			dp = mid;
	}
	score = dm;
	return get_fit(p0, dm);
}

PointsState MinDistSolver::search_regular_fits(const uint32_t iters, RGen& rgen) const
{
	PointsState best;
	double score = 0;
	for (uint32_t i = 0; i < iters; i++)
	{
		const vec2d p0 = shape->gen_point(rgen);
		double sc;
		const PointsState state = regular_fit(p0, sc);
		if (sc > score)
		{
			score = sc;
			best = state;
		}
	}
	return best;
}

vec2d circumc(const vec2d a, const vec2d b, const vec2d c)
{
	const vec2d v1 = a - b, v2 = c - b;
	const double s1 = v1.dot((a + b) / 2), s2 = v2.dot((c + b) / 2);
	const double det = v1.x * v2.y - v1.y * v2.x;
	return vec2d(
		(s1 * v2.y - v1.y * s2) / det,
		(v1.x * s2 - s1 * v2.x) / det
	);
}

struct BeachSec
{
	bool is_point;
	double x;
	const double *t;
	int32_t i, li, ri;
	vec2d p, l, r;

	BeachSec() = default;

	BeachSec(const double x)
		: is_point(true), x(x)
	{
	}

	BeachSec(const double *t, const int32_t i, const int32_t li, const int32_t ri, const vec2d p, const vec2d l, const vec2d r)
		: is_point(false), t(t), i(i), li(li), ri(ri), p(p), l(l), r(r)
	{
	}

	bool operator==(const BeachSec &other) const
	{
		return is_point == other.is_point && x == other.x && i == other.i && li == other.li && ri == other.ri;
	}

	double lx() const
	{
		if (*t == p.y)
			return p.x;
		if (*t == l.y)
			return l.x;
		const vec2d mid = (l + p) / 2;
		const vec2d vp = (p - l).get_lrot();
		const vec2d v = (l.x < p.x) ^ (vp.y > 0) ? -vp : vp;
		const double mh = *t - mid.y;
		const double bp = mh * v.y;
		const double s = (-bp + sqrt(bp * bp - v.x * v.x * ((mid - p).len2() - mh * mh))) / (v.x * v.x);
		return mid.x + v.x * s;
	}

	double rx() const
	{
		if (*t == p.y)
			return p.x;
		if (*t == r.y)
			return r.x;
		const vec2d mid = (r + p) / 2;
		const vec2d vp = (p - r).get_lrot();
		const vec2d v = (r.x > p.x) ^ (vp.y > 0) ? -vp : vp;
		const double mh = *t - mid.y;
		const double bp = mh * v.y;
		const double s = (-bp + sqrt(bp * bp - v.x * v.x * ((mid - p).len2() - mh * mh))) / (v.x * v.x);
		return mid.x + v.x * s;
	}

	double midx() const
	{
		return (lx() + rx()) / 2;
	}

	bool operator<(const BeachSec& other) const
	{
		if (is_point && other.is_point)
			return x < other.x;
		if (is_point && !other.is_point)
			return other.li != -1 && x < other.lx();
		if (!is_point && other.is_point)
			return ri != -1 && rx() < other.x;
		if ((li == -1 && other.li != -1) || (ri != -1 && other.ri == -1))
			return true;
		if (li == -1 || other.li == -1 || ri == -1 || other.ri == -1)
			return false;
		return midx() < other.midx();
	}

	bool will_be_removed() const
	{
		if (li == -1 || ri == -1)
			return false;
		if (p.y < l.y && p.y < r.y)
			return true;
		
		//return (p.y < l.y && p.x <= l.x) || (p.y < r.y && p.x >= r.x);
		const vec2d v = r - l;
		const vec2d rot = v.get_lrot();
		return (p - l).dot(rot) < 0;
	}

	double remove_time() const
	{
		const vec2d joint = circumc(l, p, r);
		return joint.y + (joint - p).len();
	}

	double y_at(const double x) const
	{
		const vec2d base(x, *t);
		if (*t == p.y)
		{
			const vec2d pp = li != -1 ? l : r;
			const vec2d mid = (base + pp) / 2;
			const vec2d v = (base - pp).get_lrot();
			const vec2d parpoint = mid + v * (x - mid.x) / v.x;
			return parpoint.y;
		}
		const vec2d mid = (base + p) / 2;
		const vec2d v = (base - p).get_lrot();
		const vec2d parpoint = mid + v * (x - mid.x) / v.x;
		return parpoint.y;
	}

};

enum FortuneEventType
{
	NEW_POINT, REMOVE_SEC, BORDER_HIT
};

struct FortuneEvent
{
	FortuneEventType type;
	double time;
	uint32_t pi;
	BeachSec sec;
	vec2d hit_point;

	FortuneEvent(const double time, const BeachSec sec, const vec2d hit_point)
		: type(FortuneEventType::BORDER_HIT), time(time), sec(sec), hit_point(hit_point)
	{
	}

	FortuneEvent(const double time, const BeachSec sec)
		: type(FortuneEventType::REMOVE_SEC), time(time), sec(sec)
	{
	}

	FortuneEvent(const double time, const uint32_t pi)
		: type(FortuneEventType::NEW_POINT), time(time), pi(pi)
	{
	}

	bool operator<(const FortuneEvent& other) const
	{
		return time > other.time;
	}
};

struct LineSegment
{
	vec2d v0, v1;

	LineSegment(const vec2d v0, const vec2d v1)
		: v0(v0), v1(v1)
	{
	}

	vec2d get_normal() const
	{
		return (v1 - v0).get_lrot();
	}
};

class Fortune
{
private:
	double t;
	std::priority_queue<FortuneEvent> events;
	vec<vec2d> points;
	vec<LineSegment> segs; 
	std::set<BeachSec> beach;

	void submit_point(const vec2d p, const double dist2)
	{
		pots.push_back(p);
		std::cout << "pot at " << p << std::endl;
	}

	void print_beach()
	{
		std::cout << "(t=" << t << ")";
		for (const BeachSec& bs : beach)
		{
			std::cout << "[" << (bs.li != -1 ? bs.lx() : 0) << ", " << bs.li << ", " << bs.i << ", " << bs.ri << ", " << (bs.ri != -1 ? bs.rx() : 0) << "]";
		}
		std::cout << std::endl;
	}

	void add_border_line(const BeachSec& left)
	{
		const vec2d v = left.r - left.p;
		const vec2d dir = v.get_lrot();
		const vec2d startp(left.rx(), left.y_at(left.rx()));
		cout << "r = " << left.r << ", p = " << left.p << endl;
		const double s = v.dot((left.p + left.r) / 2);
		for (const LineSegment& ls : segs)
		{
			const double d0 = v.dot(ls.v0), d1 = v.dot(ls.v1);
			if ((d0 < s && s < d1) || (d1 < s && s < d0))
			{
				const vec2d n2 = ls.get_normal();
				const double s2 = n2.dot(ls.v0);
				const double det = v.x * n2.y - v.y * n2.x;
				const vec2d inter(
					(s * n2.y - v.y * s2) / det,
					(v.x * s2 - s * n2.x) / det
				);
				auto bit = beach.lower_bound(BeachSec(inter.x));
				const double yalready = bit->y_at(inter.x);
				if (yalready < inter.y && dir.dot(inter) > dir.dot(startp))
				{
					const double time = inter.y + left.p.dist(inter);
					cout << "inter = " << inter << endl;
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
					}
					events.push(FortuneEvent(time, left, inter));
				}
			}
		}
	}

	void resolve_top()
	{
		const FortuneEvent ev = events.top();
		events.pop();
		std::cout << "processing event at t = " << ev.time << std::endl;
		if (ev.type == FortuneEventType::NEW_POINT)
		{
			t = ev.time;
			const uint32_t pi = ev.pi;
			const vec2d p = points[pi];

			std::cout << "will add point " << p << "(" << pi << ")" << std::endl;
			std::cout << "before add: ";
			print_beach();

			if (beach.empty())
			{
				beach.insert(BeachSec(&t, pi, -1, -1, p, vec2d(), vec2d()));
			}
			else
			{
				const auto bit = beach.lower_bound(BeachSec(p.x));
				const BeachSec orig = *bit;
				beach.erase(bit);
				const BeachSec newleft(&t, orig.i, orig.li, pi, orig.p, orig.l, p);
				const BeachSec newcenter(&t, pi, orig.i, orig.i, p, orig.p, orig.p);
				const BeachSec newright(&t, orig.i, pi, orig.ri, orig.p, p, orig.r);
				/*const auto lit = */beach.insert(newleft);
				beach.insert(newcenter);
				beach.insert(newright);
				/*if (newleft.li != -1)
					add_border_line(*std::prev(lit.first));*/
				add_border_line(newleft);
				add_border_line(newcenter);
				if (newright.ri != -1)
					add_border_line(newright);
				if (newleft.will_be_removed())
					events.push(FortuneEvent(newleft.remove_time(), newleft));
				if (newright.will_be_removed())
					events.push(FortuneEvent(newright.remove_time(), newright));
			}
			std::cout << "after add: ";
			print_beach();
		}
		else if (ev.type == FortuneEventType::REMOVE_SEC)
		{
			const vec2d crcc = circumc(ev.sec.l, ev.sec.p, ev.sec.r);
			std::cout << "possible remove arc [" << ev.sec.li << ", " << ev.sec.i << ", " << ev.sec.ri << "] " << crcc << std::endl;
			const auto bit = beach.find(ev.sec);
			if (bit != beach.end())
			{
				std::cout << "arc exists -- removing" << std::endl;
				submit_point(crcc, (crcc - ev.sec.p).len2());
				BeachSec prv = *std::prev(bit);
				BeachSec nxt = *std::next(bit);
				beach.erase(std::next(bit));
				beach.erase(std::prev(bit));
				beach.erase(bit);
				prv.ri = nxt.i;
				prv.r = nxt.p;
				nxt.li = prv.i;
				nxt.l = prv.p;
				/*const auto prvit = */beach.insert(prv);
				beach.insert(nxt);
				std::cout << "after removal: ";
				print_beach();
				
				t = ev.time;

				/*if (prv.li != -1)
					add_border_line(*std::prev(prvit.first));*/
				add_border_line(prv);
				if (nxt.ri != -1)
					add_border_line(nxt);
				if (prv.will_be_removed())
					events.push(FortuneEvent(prv.remove_time(), prv));
				if (nxt.will_be_removed())
					events.push(FortuneEvent(nxt.remove_time(), nxt));
			}
			std::cout << "after time update: ";
			print_beach();
		}
		else if (ev.type == FortuneEventType::BORDER_HIT)
		{
			if (ev.time < t)
				cout << "baaad" << endl;
			t = ev.time;
			const auto bit = beach.find(ev.sec);
			if (bit != beach.end() && ev.sec == *bit)
			{
				cout << "border hit" << endl;
				submit_point(ev.hit_point, ev.hit_point.dist2(ev.sec.p));
				cout << "p = " << ev.hit_point << endl;
				cout << "left i = " << ev.sec.i << ", right i = " << ev.sec.ri << endl;
				const auto psea = beach.lower_bound(BeachSec(ev.hit_point.x));
				cout << "psea i = " << psea->i << ", par y = " << psea->y_at(ev.hit_point.x) << endl;
				cout << "left-right border x = " << ev.sec.rx() << ", y = " << ev.sec.y_at(ev.sec.rx()) << endl;
			}
		}
		std::cout << std::endl;
	}

public:
	vec2d furthest;
	double dist2;
	std::vector<vec2d> pots;

	Fortune(const PointsState& ps, const uint32_t moving_ind, const vec<LineSegment>& segs)
		: t(0), segs(segs)
	{
		points.reserve(ps.size() - (moving_ind < ps.size()));
		for (uint32_t i = 0; i < ps.size(); i++)
		{
			if (i != moving_ind)
			{
				points.push_back(ps[i]);
				const uint32_t ind = points.size() - 1;
				events.push(FortuneEvent(ps[i].y, ind));
			}
		}
		while (!events.empty())
			resolve_top();
	}
};

PointsState MinDistSolver::solve() const
{
	std::random_device rd;
	std::mt19937 rgen(rd());

	PointsState ps = shape->gen_state(n, rgen);
	//PointsState ps{ vec2d(.49, .1), vec2d(.3, .4), vec2d(.8, .35), vec2d(.5, .7) };
	//PointsState ps{ vec2d(.49, .1), vec2d(.3, .4), vec2d(.8, .35) };
	//PointsState ps{ vec2d(.49, .8), vec2d(.3, .4), vec2d(.8, .35) };
	//PointsState ps{ vec2d(.6, .38), /*vec2d(.49, .1), */vec2d(.3, .4), vec2d(.8, .35), vec2d(.5, .7) };
	const vec<LineSegment> unitsquare{
		LineSegment(vec2d(0, 0), vec2d(1, 0)),
		LineSegment(vec2d(1, 0), vec2d(1, 1)),
		LineSegment(vec2d(1, 1), vec2d(0, 1)),
		LineSegment(vec2d(0, 1), vec2d(0, 0))
	};
	Fortune frt(ps, -1, unitsquare);
	write_image(shape, ps, "fortune.png", 0, frt.pots);

	return ps;

	/*

	PointsState ps = shape->gen_state(n, rgen);
	if (n < 200)
	{
		jump_solve(ps, 500'000, rgen);
		shift_solve(ps, 50'000, rgen);
	}
	else if (n < 1000)
	{
		jump_solve(ps, 500'000, rgen);
		shift_solve(ps, 1'000, rgen);
	}
	else if (n < 2000)
	{
		jump_solve(ps, 10'000, rgen);
		shift_solve(ps, 100, rgen);
	}
	else
	{
		jump_solve(ps, 100, rgen);
		shift_solve(ps, 10, rgen);
	}
	uint32_t fit_iters;
	if (n < 50)
	{
		fit_iters = 50'000;
	}
	else if (n < 200)
	{
		fit_iters = 7'000;
	}
	else if (n < 1000)
	{
		fit_iters = 1'000;
	}
	else
	{
		fit_iters = 100;
	}

	PointsState rfit = search_regular_fits(fit_iters, rgen);
	const double fit_score = MinDistEvaluator().evalv(rfit);
	std::cout << "pure fit: " << fit_score << std::endl;
	shift_solve(rfit, 2000, rgen);
	const double shift_score = MinDistEvaluator().evalv(rfit);
	std::cout << "after shift: " << shift_score << std::endl;
	
	return fit_score > MinDistEvaluator().evalv(ps) ? rfit : ps;*/
}

}