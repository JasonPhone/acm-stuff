#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll MOD = 1e9 + 7;
ll QpowMod(ll bse, ll pwr) {
  ll ret = 1;
  while (pwr) {
    if (pwr & 1) ret = ret * bse % MOD;
    bse = bse * bse % MOD;
    pwr >>= 1;
  }
  return ret;
}
struct Point2 {
  ll x, y;
  Point2() : x(0), y(0) {}
  Point2(ll _x, ll _y) : x(_x), y(_y) {}
  ll Norm2() { return 1ll * x * x + 1ll * y * y; }
  double Norm() { return sqrt(Norm2()); }
  Point2 operator+(const Point2 &po) { return Point2(x + po.x, y + po.y); }
  Point2 operator-(const Point2 &po) {
    // note the direction
    return Point2(x - po.x, y - po.y);
  }
  bool operator==(const Point2 &po) { return x == po.x && y == po.y; }
  bool operator<(const Point2 &po) {
    if (y == po.y) return x < po.x;
    return y < po.y;
  }
};
ll MulCross(const Point2 &p1, const Point2 &p2) {
  return p1.x * p2.y - p1.y * p2.x;
}
ll MulDot(const Point2 &p1, const Point2 &p2) {
  return p1.x * p2.x + p1.y * p2.y;
}
vector<Point2> pts;
bool cmp(int a, int b) {
  if (pts[a].y == pts[b].y) return pts[a].x >= pts[b].x;
  return pts[a].y < pts[b].y;
}
int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);
  cout.tie(0);
  int n;
  cin >> n;
  for (int i = 0; i < n; i++) {
    Point2 t;
    cin >> t.x >> t.y;
    pts.push_back(t);
  }
  sort(pts.begin(), pts.end());
  vector<int> stk;
  stk.push_back(0);
  stk.push_back(1);
  Point2 vec0, vec1;
  for (int i = 2; i < n; i++) {
    if (stk.size() < 2) continue;
    vec0 = pts[stk[stk.size() - 1]] - pts[stk[stk.size() - 2]];
    vec1 = pts[i] - pts[stk.back()];
    while (MulCross(vec0, vec1) < 0) {
      stk.pop_back();
      if (stk.size() < 2) break;
      vec0 = pts[stk[stk.size() - 1]] - pts[stk[stk.size() - 2]];
      vec1 = pts[i] - pts[stk.back()];
    }
    stk.push_back(i);
  }
  stk.push_back(n - 2);
  for (int i = n - 3; i >= 0; i--) {
    if (stk.size() < 2) continue;
    vec0 = pts[stk[stk.size() - 1]] - pts[stk[stk.size() - 2]];
    vec1 = pts[i] - pts[stk.back()];
    while (MulCross(vec0, vec1) < 0) {
      stk.pop_back();
      if (stk.size() < 2) break;
      vec0 = pts[stk[stk.size() - 1]] - pts[stk[stk.size() - 2]];
      vec1 = pts[i] - pts[stk.back()];
    }
    stk.push_back(i);
  }
  stk.pop_back();
  cout << stk.size() << endl;
  for (auto x : stk) {
    cout << pts[x].x << " " << pts[x].y << endl;
  }
  return 0;
}
