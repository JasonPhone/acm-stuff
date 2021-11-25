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
  Fraction(Fraction fup, Fraction fdn) {
    Fraction tmp = fup / fdn;
    tmp.reduce();
    up = tmp.up;
    dn = tmp.dn;
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
Fraction frac_zero{0, 1}, frac_one{1, 1};
Fraction mul_dot(const Point3fra &p1, const Point3fra &p2) {
  return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}
Point3fra mul_cross(const Point3fra &p1, const Point3fra &p2) {
  return Point3fra(p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z,
                   p1.x * p2.y - p1.y * p2.x);
}
Point3fra mul_scale(const Point3fra &p1, const Fraction s) {
  return Point3fra(p1.y * s, p1.z * s, p1.x * s);
}
bool is_segs_intersect(Point3fra A, Point3fra B, Point3fra C, Point3fra D) {
  if (A == C || A == D || B == C || B == D) return true;
  Vector3fra nm = mul_cross(B - A, D - C);
  if (mul_dot(C - A, nm) == frac_zero) return true;
  return false;
}
Fraction point_to_point2(Point3fra A, Point3fra B) {
  return (A - B).norm2();
}
Fraction point_to_seg2(Point3fra P, Point3fra A, Point3fra B) {
  Vector3fra ap = P - A, ab = B - A, bp = P - B, ba = A - B;
  Fraction ret = frac_zero;
  if (mul_dot(ap, ab) <= frac_zero || mul_dot(bp, ba) <= frac_zero) {
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
Fraction seg_to_seg2(Point3fra A, Point3fra B, Point3fra C, Point3fra D) {
  Vector3fra ca = A - C, cb = B - C, cd = D - C, ab = B - A, ac = C - A;
  Fraction tmp = mul_dot(mul_cross(ca, cb), cd);
  if (is_segs_intersect(A, B, C, D)) {
    // intersect
    return frac_zero;
  } else if (tmp == frac_zero) {
    // same plane
    Fraction ret = point_to_seg2(A, C, D);
    ret = min(ret, point_to_seg2(B, C, D));
    ret = min(ret, point_to_seg2(C, A, B));
    ret = min(ret, point_to_seg2(D, A, B));
    return ret;
  } else {
    // different plane
    Fraction dn = mul_dot(ab, cd) * mul_dot(ab, cd) - ab.norm2() * cd.norm2();
    Fraction t(ab.norm2() * mul_dot(cd, ac) - mul_dot(ab, cd) * mul_dot(ab, ac),
               dn);
    Fraction s(mul_dot(ab, cd) * mul_dot(cd, ac) - cd.norm2() * mul_dot(ab, ac),
               dn);
    if (frac_zero <= t && t <= frac_one && frac_zero <= s && s <= frac_one) {
      return point_to_seg2(A + mul_scale(ab, t), C, D);
    } else {
      Fraction ret = point_to_seg2(A, C, D);
      ret = min(ret, point_to_seg2(B, C, D));
      ret = min(ret, point_to_seg2(C, A, B));
      ret = min(ret, point_to_seg2(D, A, B));
      return ret;
    }
  }
}

int main() {
  Point3fra a{{0, 1}, {0, 1}, {1, 1}}, b{{1, 1}, {1, 1}, {1, 1}};
  Point3fra c = mul_cross(a, b);
  cout << c.norm() << endl;
  return 0;
}
