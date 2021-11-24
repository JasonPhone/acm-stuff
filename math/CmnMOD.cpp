#include <stdio.h>
using ll = long long;
const ll MN = 2000000;
const ll MOD = 1000000007;

int fac[MN + 5], inv[MN + 5];

ll qpowMod(ll bse, ll pwr) {
  ll ret = 1;
  while (pwr) {
    if (pwr & 1) ret = ret * bse % MOD;
    bse = bse * bse % MOD;
    pwr >>= 1;
  }
  return ret;
}

void init() {
  fac[0] = 1;
  for (int i = 1; i <= MN; i++) fac[i] = 1ll * fac[i - 1] * i % MOD;
  inv[MN] = qpowMod(fac[MN], MOD - 2);
  for (int i = MN - 1; i >= 0; i--) inv[i] = 1ll * inv[i + 1] * (i + 1) % MOD;
}

int C(int n, int m) {
  if (m > n) return 0;
  return 1ll * fac[n] * inv[m] % MOD * inv[n - m] % MOD;
}
int main() {
  init();
  printf("%d\n", C(5, 3));
  return 0;
}
