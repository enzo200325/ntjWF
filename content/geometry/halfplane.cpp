/**
 * Author: koosaga
 * Date: 
 * License: 
 * Source: 
 * Description: halfplane intersection do koosaga
 * Status: Works fine, used a lot
 */

struct hpi {
	const double eps = 1e-8;
	typedef pair<long double, long double> pi;
	bool z(long double x) const { return fabs(x) < eps; }
	long double ccw(pi a, pi b, pi c) {
		long double dx1 = b.first - a.first;
		long double dy1 = b.second - a.second;
		long double dx2 = c.first - a.first;
		long double dy2 = c.second - a.second;
		return dx1 * dy2 - dy1 * dx2;
	}
	struct line {
		long double a, b, c;
		pi slope() { return pi(a, b); }
	};
	pi cross(line a, line b) {
		long double det = a.a * b.b - b.a * a.b;
		return pi((a.c * b.b - a.b * b.c) / det, (a.a * b.c - a.c * b.a) / det);
	}
	bool bad(line a, line b, line c) {
		if (ccw(pi(0, 0), a.slope(), b.slope()) <= 0)
			return false;
		pi crs = cross(a, b);
		return crs.first * c.a + crs.second * c.b >= c.c;
	}
	bool solve(vector<line> v, vector<pi> &solution) { // ax + by <= c;
		sort(v.begin(), v.end(), [&](line &fuck, line &l) {
			auto [a, b, c] = make_tuple(fuck.a, fuck.b, fuck.c);
			bool flag1 = pi(a, b) > pi(0, 0);
			bool flag2 = pi(l.a, l.b) > pi(0, 0);
			if (flag1 != flag2)
				return flag1 > flag2;
			long double t = ccw(pi(0, 0), pi(a, b), pi(l.a, l.b));
			return z(t) ? c * hypot(l.a, l.b) < l.c * hypot(a, b) : t > 0;
		});
		deque<line> dq;
		for (auto &i : v) {
			if (!dq.empty() && z(ccw(pi(0, 0), dq.back().slope(), i.slope())))
				continue;
			while (dq.size() >= 2 && bad(dq[dq.size() - 2], dq.back(), i))
				dq.pop_back();
			while (dq.size() >= 2 && bad(i, dq[0], dq[1]))
				dq.pop_front();
			dq.push_back(i);
		}
		while (dq.size() > 2 && bad(dq[dq.size() - 2], dq.back(), dq[0]))
			dq.pop_back();
		while (dq.size() > 2 && bad(dq.back(), dq[0], dq[1]))
			dq.pop_front();
		vector<pi> tmp;
		for (int i = 0; i < dq.size(); i++) {
			line cur = dq[i], nxt = dq[(i + 1) % dq.size()];
			if (ccw(pi(0, 0), cur.slope(), nxt.slope()) <= eps)
				return false;
			tmp.push_back(cross(cur, nxt));
		}
		solution = tmp;
		return true;
	}
} hpi;
