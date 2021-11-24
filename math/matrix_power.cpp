#include <bits/stdc++.h>
#define inf 0x3f3f3f3f
using namespace std;
typedef long long ll;
const int N = 205, mod = 998244353, MS = 205;
struct Mat {
  ll a[MS][MS];
  ll n, m;
  Mat(int n = 0, int m = 0) : n(n), m(m) { memset(a, 0, sizeof(a)); }
  Mat operator*(const Mat& B) const {
    Mat C(n, B.m);
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= B.m; j++)
        for (int k = 1; k <= m; k++)
          C.a[i][j] = (C.a[i][j] + a[i][k] * B.a[k][j]) % mod;
    return C;
  }
};
Mat qpow(Mat a, int n) {
  Mat ans(a.n, a.n);
  for (int i = 1; i <= a.n; i++) ans.a[i][i] = 1;
  for (; n; n >>= 1, a = a * a)
    if (n & 1) ans = ans * a;
  return ans;
}
int main() {
  ll n;
  cin >> n;
  string s;
  cin >> s;
  ll now = stol(s);
  Mat A(100, 100);
  A = qpow(A, n);

  Mat B(100, 100);
  B.a[1][1] = 1;
  B = B * A;
  cout << B.a[1][now];
}
