/**
 * for inverse element from 1 to 3e6
 * O(n) faster than O(nlogn)(quick_pow)
 * inv(x) = (p - p / x) * inv(p % x) % p;
 */
#include <cstdio>
using ll = long long;
const ll MOD = 1e9 + 7;
const ll MAXN = 3e6 + 6;
using namespace std;
ll inv[MAXN], n, p;
int main() {
  scanf("%lld %lld", &n, &p);
  inv[1] = 1;
  printf("1\n");
  for (int i = 2; i <= n; i++) {
    inv[i] = (p - p / i) * inv[p % i] % p;
    printf("%lld\n", inv[i]);
  }
}
