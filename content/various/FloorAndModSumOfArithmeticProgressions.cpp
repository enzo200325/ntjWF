/**
 * Author: Simon Lindholm
 * License: CC0
 * Source: Codeforces
 * Description: Given $a[i] = \min_{lo(i) \le k < hi(i)}(f(i, k))$ where the (minimal)
 * optimal $k$ increases with $i$, computes $a[i]$ for $i = L..R-1$.
 * Time: O((N + (hi-lo)) \log N)
 * Status: tested on http://codeforces.com/contest/321/problem/E
 */

ll sumsq(ll n) { return n / 2 * ((n - 1) | 1); }
ll floor_sum(ll a, ll d, ll m, ll n) {
    ll res = d / m * sumsq(n) + a / m * n;
    d %= m;
    a %= m;
    if (!d) return res;
    ll to = (n * d + a) / m;
    return res + (n - 1) * to - floor_sum(m - 1 - a, m, d, to);
}
ll mod_sum(ll a, ll d, ll m, ll n) {
    a = ((a % m) + m) % m;
    d = ((d % m) + m) % m;
    return n * a + d * sumsq(n) - m * floor_sum(a, d, m, n);
}
