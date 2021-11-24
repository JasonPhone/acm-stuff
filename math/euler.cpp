#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const int MAXN = 1e6 + 5;
const int MOD = 1e9 + 7;
// priority_queue<ll, vector<ll>, greater<ll>> minor_que;

int prime[MAXN];
bool vis[MAXN];
int cnt = 0;
ll maxv = -1;
void EulerPrime(int n) {
  for (int i = 2; i <= n; ++i) {
    if (vis[i] == 0) {
      prime[cnt++] = i;
      vis[i] = 1;
    }
    for (int j = 0; i * prime[j] <= n; ++j) {
      vis[i * prime[j]] = 1;
      if (i % prime[j] == 0) break;  // O(n)
    }
  }
}
int main() {
  EulerPrime(100);
  for (int i = 0; i < cnt; ++i) printf("%d ", prime[i]);
  printf("\n");
  return 0;
}
