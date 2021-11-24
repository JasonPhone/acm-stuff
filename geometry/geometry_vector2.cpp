/**
 * structs of
 * point, vector, segment
 * and some operator overloads
 */
// whether a seg AB intersects with a circle O?
// see the endpoints' tangent point (P, Q) angle
// angles: AOP + BOQ < AOB <==> intersect
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
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
  Point2 operator+(const Point2 &po) {
    return Point2(x + po.x, y + po.y);
  }
  Point2 operator-(const Point2 &po) {
    // note the direction
    return Point2(x - po.x, y - po.y);
  }
  bool operator==(const Point2 &po) {
    return x == po.x && y == po.y;
  }
};
typedef Point2 Vector2;
struct Segment2 {
  Point2 s, e;
  Segment2() {}
  Segment2(Point2 _s, Point2 _e) : s(_s), e(_e) {}
};
ll MulCross(const Point2 &p1, const Point2 &p2) {
  return p1.x * p2.y - p1.y * p2.x;
}
ll MulDot(const Point2 &p1, const Point2 &p2) {
  return p1.x * p2.x + p1.y * p2.y;
}
double DisPointToSeg(Point2 p, Point2 s1, Point2 s2) {
  Point2 v1 = p - s1, v2 = s2 - s1;
  if (MulDot(v2, v1) < 0 || MulDot(v2, v1) > v2.Norm2())
    return min(1.0 * (p - s1).Norm(), 1.0 * (p - s2).Norm());
  return abs(1.0 * MulCross(v2, v1) / v2.Norm());
}
int Dis2PointToSeg_INT(Point2 p, Point2 s1, Point2 s2) {
  // square of distance between two points
  Point2 v = p - s1, u = s2 - s1;
  if (MulDot(u, v) < 0 || MulDot(u, v) > u.Norm2())
    return min((p - s1).Norm2(), (p - s2).Norm2()) % MOD;
  return ((MulCross(v, u) % MOD) * (MulCross(v, u) % MOD)) % MOD *
         QpowMod(u.Norm2() % MOD, MOD - 2) % MOD;
}
int main() { return 0; }
