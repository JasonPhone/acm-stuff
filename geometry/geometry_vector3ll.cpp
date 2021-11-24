#include <bits/stdc++.h>
using namespace std;
using ll = long long;
ll MOD = 1e9 + 7;
struct Point3fra {
  ll x, y, z;
  Point3fra() : x(0), y(0), z(0) {}
  Point3fra(ll _x, ll _y, ll _z) : x(_x), y(_y), z(_z) {}
  ll norm2() { return x * x + y * y + z * z; }
  double norm() { return sqrt(norm2()); }
  Point3fra operator+(const Point3fra &po) {
    return Point3fra(x + po.x, y + po.y, z + po.z);
  }
  Point3fra operator-(const Point3fra &po) {
    return Point3fra(x - po.x, y - po.y, z - po.z);
  }
  bool operator==(const Point3fra &po) {
    return x == po.x && y == po.y && z == po.z;
  }
};
typedef Point3fra Vector3ll;
struct Segment3ll {
  Point3fra s, e;
  Segment3ll() {}
  Segment3ll(Point3fra _s, Point3fra _e): s(_s), e(_e) {}
};
ll mul_dot(const Point3fra &p1, const Point3fra &p2) {
  return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}
Point3fra mul_cross(const Point3fra &p1, const Point3fra &p2) {
  return Point3fra(p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x);
}
int main() {
  Point3fra a{0, 0, 1}, b{1, 1, 1};
  Point3fra c = mul_cross(a, b);
  cout << c.norm() << endl;
  return 0;
}
