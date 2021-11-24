#include <cstdio>
#include <algorithm>
using namespace std;
// 较小规模下计算组合数
constexpr int N = 9000, M = 1e9 + 7;
long long C[N][N];
void getCmn() {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
      if (j == 0)
        C[i][j] = 1;
      else
        C[i][j] = (C[i - 1][j - 1] % M + C[i - 1][j] % M) % M;
    }
  }
}
int main() {
  getCmn();
  long long v = -1;
  for (int i = 1; i <= 100; ++i) {
    v = max(v, C[100][i]);
  }
  printf("%lld", v);
  return 0;
}
