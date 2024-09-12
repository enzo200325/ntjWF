/**
 * Author: 
 * Date: 
 * License: 
 * Source: 
 * Description:
 * Time:
 * Status: 
 */

vector<mint> and_convolution(vector<mint> A, vector<mint> B) {
    int n = (int)max(A.size(), B.size());
    int N = 0;
    while ((1 << N) < n) N++;
    A.resize(1 << N);
    B.resize(1 << N);
    vector<mint> C(1 << N);
    for (int j = 0; j < N; j++) {
        for (int i = (1 << N) - 1; i >= 0; i--) {
            if (~i >> j & 1) {
                A[i] += A[i | (1 << j)];
                B[i] += B[i | (1 << j)];
            }
        }
    }
    for (int i = 0; i < 1 << N; i++) C[i] = A[i] * B[i];
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < 1 << N; i++)
            if (~i >> j & 1) C[i] -= C[i | (1 << j)];
    }
    return C;
}

vector<mint> gcd_convolution(vector<mint> A, vector<mint> B) {
    int N = (int)max(A.size(), B.size());
    A.resize(N + 1);
    B.resize(N + 1);
    vector<mint> C(N + 1);
    for (int i = 1; i <= N; i++) {
        mint a = 0;
        mint b = 0;
        for (int j = i; j <= N; j += i) {
            a += A[j];
            b += B[j];
        }
        C[i] = a * b;
    }
    for (int i = N; i >= 1; i--)
        for (int j = 2 * i; j <= N; j += i) C[i] -= C[j];
    return C;
}

vector<mint> lcm_convolution(vector<mint> A, vector<mint> B) {
    int N = (int)max(A.size(), B.size());
    A.resize(N + 1);
    B.resize(N + 1);
    vector<mint> C(N + 1), a(N + 1), b(N + 1);
    for (int i = 1; i <= N; i++) {
        for (int j = i; j <= N; j += i) {
            a[j] += A[i];
            b[j] += B[i];
        }
        C[i] = a[i] * b[i];
    }
    for (int i = 1; i <= N; i++)
        for (int j = 2 * i; j <= N; j += i) C[j] -= C[i];
    return C;
}

vector<mint> or_convolution(vector<mint> A, vector<mint> B) {
    int n = (int)max(A.size(), B.size());
    int N = 0;
    while ((1 << N) < n) N++;
    A.resize(1 << N);
    B.resize(1 << N);
    vector<mint> C(1 << N);
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < 1 << N; i++) {
            if (i >> j & 1) {
                A[i] += A[i ^ (1 << j)];
                B[i] += B[i ^ (1 << j)];
            }
        }
    }
    for (int i = 0; i < 1 << N; i++) C[i] = A[i] * B[i];
    for (int j = N - 1; j >= 0; j--) {
        for (int i = (1 << N) - 1; i >= 0; i--)
            if (i >> j & 1) C[i] -= C[i ^ (1 << j)];
    }
    return C;
}

vector<mint> subset_convolution(vector<mint> A, vector<mint> B) {
    int n = int(max(A.size(), B.size()));
    int N = 0;
    while ((1 << N) < n) N++;
    A.resize(1 << N), B.resize(1 << N);
    vector a(1 << N, vector<mint>(N + 1)), b(1 << N, vector<mint>(N + 1));
    for (int i = 0; i < 1 << N; i++) {
        int popcnt = __builtin_popcount(i);
        a[i][popcnt] = A[i];
        b[i][popcnt] = B[i];
    }
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < 1 << N; i++) {
            if (~i >> j & 1) continue;
            for (int popcnt = 0; popcnt <= N; popcnt++) {
                a[i][popcnt] += a[i ^ (1 << j)][popcnt];
                b[i][popcnt] += b[i ^ (1 << j)][popcnt];
            }
        }
    }
    vector c(1 << N, vector<mint>(N + 1));
    for (int i = 0; i < 1 << N; i++) {
        for (int j = 0; j <= N; j++)
            for (int k = 0; k + j <= N; k++) c[i][j + k] += a[i][j] * b[i][k];
    }
    for (int j = N - 1; j >= 0; j--) {
        for (int i = (1 << N) - 1; i >= 0; i--) {
            if (~i >> j & 1) continue;
            for (int popcnt = 0; popcnt <= N; popcnt++)
                c[i][popcnt] -= c[i ^ (1 << j)][popcnt];
        }
    }
    vector<mint> ans(1 << N);
    for (int i = 0; i < 1 << N; i++) {
        int popcnt = __builtin_popcount(i);
        ans[i] = c[i][popcnt];
    }
    return ans;
}

vector<mint> xor_convolution(vector<mint> A, vector<mint> B) {
    int n = int(A.size());
    for (int rep = 0; rep < 2; rep++) {
        for (int len = n >> 1; len; len >>= 1) {
            for (int i = 0; i < n; i += len << 1) {
                for (int j = 0; j < len; j++) {
                    int id = i + j;
                    mint x = A[id];
                    mint y = A[id + len];
                    A[id] = x + y;
                    A[id + len] = x - y;
                }
            }
        }
        swap(A, B);
    }
    vector<mint> ans(n);
    for (int i = 0; i < n; i++) ans[i] = A[i] * B[i];
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += len << 1) {
            for (int j = 0; j < len; j++) {
                int id = i + j;
                mint x = ans[id];
                mint y = ans[id + len];
                ans[id] = x + y;
                ans[id + len] = x - y;
            }
        }
    }
    return ans;
}

vector<mint> xor_multiply(vector<mint> A, vector<mint> B) {
    int N = 1;
    int n = int(max(A.size(), B.size()));
    while (N < n) N <<= 1;
    A.resize(N);
    B.resize(N);
    auto ans = xor_convolution(A, B);
    for (int i = 0; i < N; i++) ans[i] /= N;

    return ans;
}
