/**
 * Author: dudu
 * Date:
 * License: 
 * Description: Prefix sum of Dirichlet Convolution
*  Dadas duas funções aritméticas $f$ e $g$ com valores computados para $1 \le i \le N^{\frac{2}{3}}$, além de seus prefixos de soma $F$ e $G$ com valores computados para todo $\left\lfloor \frac{N}{i} \right\rfloor$ com $1 \le i \le N^{\frac{1}{3}}$, obtém-se o vetor $H$ tal que
* $H_i = \sum_{j=1}^{\lfloor \frac{n}{i} \rfloor} h(j)$
* para todo $1 \le i \le N^{\frac{1}{3}}$ em $\mathcal{O}(N^{\frac{2}{3}})$, onde $h(n) = \sum_{d \mid n} f(d) \cdot g(\lfloor \tfrac{n}{d} \rfloor)$ é a convolução de Dirichlet de $f$ com $g$. Para atingir essa complexidade, inicialize a estrutura com $T = N^{2/3}$.
* Para obter os demais valores de $H$ (para $1 \le i \le N^{\frac{2}{3}}$) utilize a [convolução de Dirichlet linear](../Dirichlet-Convolution/dirichlet_convolution.cpp).
* funcoes pra mandar pro solve: f(x) -> retorna f(x), g(x) -> retorna g(x), F(x) -> retorna a soma de prefixo de f(x), G(x) -> retorna a soma de prefixo de g(x)
 * Time: $O(N^{\frac{2}{3}}$
 * Status: Works
 */
struct DirichletConvolutionPrefix {
    ll N;
    int T;
    vector<mint> ans;

    DirichletConvolutionPrefix(ll n, int t) : N(n), T(t) { ans.assign(n / T + 1, 0); }

    vector<mint> solve(auto &&f, auto &&F, auto &&g, auto &&G) {
        if (N == 1) return vector<mint>(2, 1);
        ans.assign(N / T + 1, 0);
        for (ll i = 1; i <= N / T; i++) {
            ll now = N / i;
            mint f_sum = 0, g_sum = 0;
            for (int j = 1; (ll)j * j <= now; j++) {
                f_sum += f(j);
                g_sum += g(j);

                ans[i] += G(i * j) * f(j);
                ans[i] += F(i * j) * g(j);
            }
            ans[i] -= f_sum * g_sum;
        }
        return ans;
    }
};
