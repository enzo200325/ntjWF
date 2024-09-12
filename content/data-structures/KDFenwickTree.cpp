/**
 * Author: 
 * Date: 
 * License: 
 * Source: 
 * Description: Fenwick Tree em $k$ dimensões. Faz apenas queries de prefixo e updates pontuais em $O(k\log^k n)$ 
 * Para queries em range, deve-se fazer inclusão-exclusão, porém a complexidade fica exponencial,
 * para $k$ dimensões a query em range é $O(2^k k \log^k n)$. 
 * Time:
 * Status: 
 */

const int MAX = 20;
long long tree[MAX][MAX][MAX][MAX]; // insira o numero de dimensoes aqui

long long query(vector<int> s, int pos = 0) { // s eh a coordenada
    long long sum = 0;
    while (s[pos] >= 0) {
        if (pos < (int)s.size() - 1) {
            sum += query(s, pos + 1);
        } else {
            sum += tree[s[0]][s[1]][s[2]][s[3]];
            // atualizar se mexer no numero de dimensoes
        }
        s[pos] = (s[pos] & (s[pos] + 1)) - 1;
    }
    return sum;
}

void update(vector<int> s, int v, int pos = 0) {
    while (s[pos] < MAX) {
        if (pos < (int)s.size() - 1) {
            update(s, v, pos + 1);
        } else {
            tree[s[0]][s[1]][s[2]][s[3]] += v;
            // atualizar se mexer no numero de dimensoes
        }
        s[pos] |= s[pos] + 1;
    }
}
