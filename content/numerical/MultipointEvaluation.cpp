
/**
 * Author: UFMG
 * Date: 
 * Description: Avalia o polinomio f(x) nos pontos p[0], p[1], ..., p[n-1]
 * Time: O(n log^2(n))
 * Status: 
 */

using poly = vector<mint>;
const int MAGIC = 512;
poly D(poly p) {
	if (p.empty()) return p;
	for (int i = 0; i + 1 < p.size(); i++)
		p[i] = (i + 1) * p[i + 1];
	p.pop_back();
	return p;
}

pair<poly, poly> divslow(const poly& a, const poly& b) {
	poly q, r = a;
	while (r.size() >= b.size()) {
		q.push_back(r.back() / b.back());
		if (q.back() != 0)
			for (int i = 0; i < b.size(); i++)
				r.end()[-i-1] -= q.back() * b.end()[-i-1];
		r.pop_back();
	}
	reverse(q.begin(), q.end());
	return {q, r};
}

// retorna (q, r) : a(x) = b(x) * q(x) + r(x)
pair<poly, poly> divmod(const poly& a, const poly& b) {
	if (a.size() < b.size()) return {{}, a};
	if (max(b.size(), a.size() - b.size()) < MAGIC) return divslow(a, b);
	poly ra = poly(a.rbegin(), a.rend());
	poly rb = poly(b.rbegin(), b.rend());
	int k = a.size() - b.size() + 1;
	rb.resize(k);
	poly irb = inv(move(rb)), q = convolution(ra, irb);
	q = poly(q.rend() - k, q.rend());
	poly r = convolution(move(q), b);
	for (int i = 0; i < r.size(); i++) r[i] = a[i] - r[i];
	while (r.size() > 1 && r.back() == 0) r.pop_back();
	return {q, r};
}


namespace multipoint {
	vector<poly> tree;
	void build(vector<mint>& p) {
		int n = p.size();
		tree.resize(2*n);
		for (int i = 0; i < n; i++) tree[n + i] = {-p[i], 1};
		for (int i = n - 1; i > 0; i--)
			tree[i] = convolution(tree[2*i], tree[2*i + 1]);
	}
	vector<mint> evaluate(poly& f, vector<mint>& p) {
		build(p);
		int n = p.size();
		vector<poly> ans(2 * n);
		ans[1] = divmod(f, tree[1]).second;
		for (int i = 2; i < 2 * n; i++)
			ans[i] = divmod(ans[i/2], tree[i]).second;
		vector<mint> results(n);
		for (int i = 0; i < n; i++) results[i] = ans[n + i][0];
		return results;
	}
	poly prod(vector<mint>& p, int l, int r) {
		if (l == r) return {-p[l], 1};
		int m = (l + r) / 2;
		return convolution(prod(p, l, m), prod(p, m + 1, r));
	}
	poly interpolate(vector<mint>& x, vector<mint>& y) {
		int n = x.size();
		poly p = D(prod(x, 0, n - 1));
		auto d = evaluate(p, x);
		vector<poly> ans(2 * n);
		for (int i = 0; i < n; i++) ans[n + i] = {y[i] / d[i]};
		for (int i = n - 1; i > 0; i--) {
			poly p1 = convolution(tree[2*i], ans[2*i + 1]);
			poly p2 = convolution(tree[2*i + 1], ans[2*i]);
			ans[i] = p1;
			for (int j = 0; j < p1.size(); j++) ans[i][j] += p2[j];
		}
		return ans[1];
	}
}
