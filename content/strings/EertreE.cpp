/**
 * Author: UDESC
 * Date: 
 * License: 
 * Source:
 * Description: palindromic tree, palindromic factorization
 * Time: $O(N)$ to build, can use map instead of array on 'to' for memory optmization. palindromic factorization computa o numero de formas de particionar cada
prefixo da string em strings palindromicas
 */

const int N = 2e5 + 15;
const int ALF = 26;

struct eertree {
    int str[N], len[N], lnk[N], cnt[N], first[N], node_cnt, it, last;
    ll palindrome_substring_sum;
    const char norm = 'a';
    array<int, ALF> to[N];

    inline int get(char c) { return c - norm; }

    void init(int n) {
        memset(str, 0, sizeof(int) * (n + 1));
        memset(len, 0, sizeof(int) * (n + 1));
        memset(lnk, 0, sizeof(int) * (n + 1));
        memset(cnt, 0, sizeof(int) * (n + 1));
        for (int i = 0; i <= n; i++)
            for (int j = 0; j < ALF; j++) to[i][j] = 0;
        node_cnt = 2, it = 1, last = 0, str[0] = -1;
        len[0] = 0, len[1] = -1, lnk[0] = 1, lnk[1] = 1;
    }

    void set_string(const string &s) {
        int n = (int)s.size();
        init(n);
        for (int i = 0; i < n; i++) insert(s[i]);
        build_cnt();
    }

    int insert(char ch) {
        int c = get(ch);
        str[it] = c;
        while (str[it - 1 - len[last]] != c) last = lnk[last];
        if (!to[last][c]) {
            int prev = lnk[last];
            while (str[it - 1 - len[prev]] != c) prev = lnk[prev];
            lnk[node_cnt] = to[prev][c];
            len[node_cnt] = len[last] + 2;
            to[last][c] = node_cnt++;
        }
        last = to[last][c];
        first[last] = it;
        cnt[last]++;
        it++;
        return last;
    }

    void build_cnt() {
        ll ans = 0;
        for (int i = it; i > 1; i--) {
            ans += cnt[i];
            cnt[lnk[i]] += cnt[i];
        }
        palindrome_substring_sum = ans;
    }

    inline ll number_of_palindromes() { return palindrome_substring_sum; }
    inline int number_of_distinct_palindromes() { return node_cnt - 2; }
} et;


ll factorization(string s) {
	int n = s.size(), sz = 2;
	eertree PT(n);
	vector<int> diff(n+2), slink(n+2), sans(n+2), dp(n+1);
	dp[0] = 1;
	for (int i = 1; i <= n; i++) {
		PT.insert(s[i-1]);
		if (PT.size()+2 > sz) {
			diff[sz] = PT.len[sz] - PT.len[PT.lnk[sz]];
			if (diff[sz] == diff[PT.lnk[sz]])
				slink[sz] = slink[PT.lnk[sz]];
			else slink[sz] = PT.lnk[sz];
			sz++;
		}
		for (int v = PT.last; PT.len[v] > 0; v = slink[v]) {
			sans[v] = dp[i - (PT.len[slink[v]] + diff[v])];
			if (diff[v] == diff[PT.lnk[v]])
				sans[v] = (sans[v] + sans[PT.lnk[v]]) % MOD;
			dp[i] = (dp[i] + sans[v]) % MOD;
		}
	}
	return dp[n];
}
