#include <bits/stdc++.h>
using namespace std;
using ll = long long;
ll MOD = 1e9 + 7;
struct Fraction {
  ll up, dn;
  Fraction() : up(0), dn(1) {}
  Fraction(ll _up, ll _dn) : up(_up), dn(_dn) {
    ll cd = __gcd(up, dn);
    up /= cd;
    dn /= cd;
    if (dn < 0) {
      dn = abs(dn);
      up = -up;
    }
  }
  void reduce() {
    ll cd = __gcd(up, dn);
    up /= cd;
    dn /= cd;
  }
  Fraction operator+(const Fraction &otr) const {
    ll n_dn = dn / __gcd(otr.dn, dn) * otr.dn;
    ll n_up = n_dn / dn * up + n_dn / otr.dn * otr.up;
    return Fraction(n_up, n_dn);
  }
  Fraction operator-(const Fraction &otr) const {
    ll n_dn = dn / __gcd(otr.dn, dn) * otr.dn;
    ll n_up = n_dn / dn * up - n_dn / otr.dn * otr.up;
    return Fraction(n_up, n_dn);
  }
  Fraction operator*(const Fraction &otr) const {
    ll n_dn = dn * otr.dn;
    ll n_up = up * otr.up;
    // cout << n_up << "/" << n_dn << endl;
    ll cd = __gcd(n_dn, n_up);
    return Fraction(n_up / cd, n_dn / cd);
  }
  Fraction operator/(const Fraction &otr) const {
    Fraction loprd(up, dn), roprd(otr.dn, otr.up);
    return loprd * roprd;
  }
  bool operator==(const Fraction &otr) const {
    ll uup = up, ddn = dn, cd = __gcd(up, dn);
    uup /= up, ddn /= dn;
    ll oup = otr.up, odn = otr.dn;
    cd = __gcd(oup, odn);
    oup /= cd, odn /= cd;
    return up * otr.dn == dn * otr.up;
  }
  bool operator<(const Fraction &otr) const {
    ll uup = up, ddn = dn, cd = __gcd(up, dn);
    uup /= up, ddn /= dn;
    ll oup = otr.up, odn = otr.dn;
    cd = __gcd(oup, odn);
    oup /= cd, odn /= cd;
    return uup * odn < oup * ddn;
  }
  bool operator<=(const Fraction &otr) const {
    Fraction fra{up, dn};
    return fra < otr || fra == otr;
  }
  double real_val() const { return double(up) / double(dn); }
};
struct Point3fra {
  Fraction x, y, z;
  Point3fra() : x(0, 1), y(0, 1), z(0, 1) {}
  Point3fra(Fraction _x, Fraction _y, Fraction _z) : x(_x), y(_y), z(_z) {}
  Fraction norm2() { return (x * x) + (y * y) + (z * z); }
  double norm() { return sqrt(norm2().real_val()); }
  Point3fra operator+(const Point3fra &po) {
    return Point3fra(x + po.x, y + po.y, z + po.z);
  }
  Point3fra operator-(const Point3fra &po) {
    return Point3fra(x - po.x, y - po.y, z - po.z);
  }
  bool operator==(Point3fra &po) {
    return (x == po.x) && (y == po.y) && (z == po.z);
  }
};
typedef Point3fra Vector3fra;
// struct Segment3fra {
//   Point3fra s, e;
//   Segment3fra() {}
//   Segment3fra(Point3fra _s, Point3fra _e) : s(_s), e(_e) {}
// };
Fraction zerofra{0, 1};
Fraction mul_dot(const Point3fra &p1, const Point3fra &p2) {
  return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}
Point3fra mul_cross(const Point3fra &p1, const Point3fra &p2) {
  return Point3fra(p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z,
                   p1.x * p2.y - p1.y * p2.x);
}
Fraction point_to_point2(const Point3fra A, const Point3fra B) {
  return (A - B).norm2();
}
Fraction point_to_seg2(Point3fra P, Point3fra A, Point3fra B) {
  Vector3fra ap = P - A, ab = B - A, bp = P - B, ba = A - B;
  Fraction ret = zerofra;
  if (mul_dot(ap, ab) <= zerofra || mul_dot(bp, ba) <= zerofra) {
    ret = min(ret, point_to_point2(P, A));
    ret = min(ret, point_to_point2(P, B));
  } else {
    Vector3fra pa = A - P, pb = B - P, ab = B - A;
    Fraction up = mul_cross(pa, pb).norm2(), dn = ab.norm2();
    ret = Fraction{up, dn};
  }
  return ret;
}
// TODO seg to seg
int main() {
  Point3fra a{{0, 1}, {0, 1}, {1, 1}}, b{{1, 1}, {1, 1}, {1, 1}};
  Point3fra c = mul_cross(a, b);
  cout << c.norm() << endl;
  return 0;
}
