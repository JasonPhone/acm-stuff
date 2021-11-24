#include <cstdio>
// a^(-1) mod p => a^(p - 2) mod p
// n * n * (n + 1) * (n + 1) / 4 = \sum_{1}^{n} i^3
// n * (n + 1) * (2n + 1) / 6 = \sum_{1}^{n} i^2
// a^(b^(c)) % MOD = a^(p) % MOD, where p = b^(c) % (MOD - 1)
// NOTE: gcd(a, MOD) == 1
using ll = long long;
ll MOD = 1e9 + 7;
ll quick_pow_mod(ll bse, ll pwr, ll mod) {
  ll ret = 1;
  bse %= mod;
  while (pwr) {
    if (pwr & 1) ret = ret * bse % mod;
    bse = bse * bse % mod;
    pwr >>= 1;
  }
  return ret;
}
int main() {
  ll x, y;
  ll n, p;
  scanf("%lld %lld", &n, &p);
  for (ll i = 1; i <= n; i++) {
    printf("%lld\n", quick_pow_mod(i, p - 2, p));
  }
  return 0;
}
