/**
 * Author: Brute
 * Date: 
 * License: 
 * Source: 
 * Description: Pollard Rho and Miller Rabin
 * Time: rho is expected to run in O(\sqrt[4]{n}), and miller rabin is O(\log n). Factorize is O(\log n) expected.
*/

namespace MillerRabin {
    inline ll mul_mod(ll a, ll b, ll m) { return (ll)((__int128)a * b % m); }
    inline ll power(ll b, ll e, ll m) {
        ll r = 1;
        b = b % m;
        while (e > 0) {
            if (e & 1) r = mul_mod(r, b, m);
            b = mul_mod(b, b, m), e >>= 1;
        }
        return r;
    }
    inline bool composite(ll n, ll a, ll d, ll s) {
        ll x = power(a, d, n);
        if (x == 1 || x == n - 1 || a % n == 0) return false;
        for (int r = 1; r < s; r++) {
            x = mul_mod(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    }

    // com esses "primos", o teste funciona garantido para n <= 2^64
    int primes[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    // funciona para n <= 3*10^24 com os primos ate 41, mas tem que cuidar com overflow
    // int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41};

    bool prime(ll n) {
        if (n <= 2 || (n % 2 == 0)) return n == 2;
        ll d = n - 1, r = 0;
        while (d % 2 == 0) d /= 2, r++;
        for (int a : primes)
            if (composite(n, a, d, r)) return false;
        return true;
    }
}

namespace PollardRho {
    mt19937 rng((uint32_t)chrono::steady_clock::now().time_since_epoch().count());
    const ll P = 1e6 + 1;
    ll seq[P];
    inline ll add_mod(ll x, ll y, ll m) { return (x += y) < m ? x : x - m; }
    inline ll mul_mod(ll a, ll b, ll m) { return (ll)((__int128)a * b % m); }
    ll rho(ll n) {
        if (n % 2 == 0) return 2;
        if (n % 3 == 0) return 3;
        ll x0 = rng() % n, c = rng() % n;
        while (1) {
            ll x = x0++, y = x, u = 1, v, t = 0;
            ll *px = seq, *py = seq;
            while (1) {
                *py++ = y = add_mod(mul_mod(y, y, n), c, n);
                *py++ = y = add_mod(mul_mod(y, y, n), c, n);
                if ((x = *px++) == y) break;
                v = u;
                u = mul_mod(u, abs(y - x), n);
                if (!u) return gcd(v, n);
                if (++t == 32) {
                    t = 0;
                    if ((u = gcd(u, n)) > 1 && u < n) return u;
                }
            }
            if (t && (u = gcd(u, n)) > 1 && u < n) return u;
        }
    }
}

vector<ll> factorize(ll x) {
    vector<ll> f;
    if (x == 1) return f;
    function<void(ll)> dfs = [&](ll x) {
        if (x == 1) return;
        if (x < Sieve::P) {
            auto fs = Sieve::factorize(x);
            f.insert(f.end(), fs.begin(), fs.end());
        } else if (MillerRabin::prime(x)) {
            f.push_back(x);
        } else {
            ll d = PollardRho::rho(x);
            dfs(d);
            dfs(x / d);
        }
    };
    dfs(x);
    sort(f.begin(), f.end());
    return f;
}
