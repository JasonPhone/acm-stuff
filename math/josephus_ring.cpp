#include <cstdio>
long long josephus(int n, int k) {
  if (n == 1)
    return 0;
  else
    return (josephus(n - 1, k) + k) % n;
}
int main(void) {
  long long n, k;
  scanf("%lld %lld", &n, &k);

  printf("%lld\n", 1 + josephus(n, k));
  return 0;
}
