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
  }
  Fraction(Fraction fup, Fraction fdn) {
    Fraction tmp = fup / fdn;
    tmp.reduce();
    up = tmp.up;
    dn = tmp.dn;
  }
  void reduce() {
    ll cd = abs(__gcd(up, dn));
    up /= cd;
    dn /= cd;
    neg_sign();
  }
  void neg_sign() {
    if (dn < 0) {
      dn = -dn;
      up = -up;
    }
  }
  Fraction operator+(const Fraction &otr) const {
    ll n_dn = dn / __gcd(otr.dn, dn) * otr.dn;
    ll n_up = n_dn / dn * up + n_dn / otr.dn * otr.up;
    Fraction ret{n_up, n_dn};
    ret.reduce();
    ret.neg_sign();
    return ret;
  }
  Fraction operator-(const Fraction &otr) const {
    ll n_dn = dn / __gcd(otr.dn, dn) * otr.dn;
    ll n_up = n_dn / dn * up - n_dn / otr.dn * otr.up;
    Fraction ret{n_up, n_dn};
    ret.reduce();
    ret.neg_sign();
    return ret;
  }
  Fraction operator*(const Fraction &otr) const {
    ll n_dn = dn * otr.dn;
    ll n_up = up * otr.up;
    // cout << n_up << "/" << n_dn << endl;
    ll cd = abs(__gcd(n_dn, n_up));
    n_up /= cd, n_dn /= cd;
    Fraction ret{n_up, n_dn};
    ret.reduce();
    ret.neg_sign();
    return ret;
  }
  Fraction operator/(const Fraction &otr) const {
    Fraction loprd(up, dn), roprd(otr.dn, otr.up);
    return loprd * roprd;
  }
  bool operator==(const Fraction &otr) const {
    ll uup = up, ddn = dn;
    if (ddn < 0) uup = -uup, ddn = -ddn;
    ll cd = abs(__gcd(up, dn));
    uup /= up, ddn /= dn;
    ll oup = otr.up, odn = otr.dn;
    if (odn < 0) oup = -oup, odn = -odn;
    cd = abs(__gcd(oup, odn));
    oup /= cd, odn /= cd;
    return up * otr.dn == dn * otr.up;
  }
  bool operator<(const Fraction &otr) const {
    ll uup = up, ddn = dn;
    if (ddn < 0) uup = -uup, ddn = -ddn;
    ll cd = abs(__gcd(up, dn));
    ll oup = otr.up, odn = otr.dn;
    if (odn < 0) oup = -oup, odn = -odn;
    cd = abs(__gcd(oup, odn));
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
/******** types done ********/
/******* functions go *******/
Fraction frac_zero{0, 1}, frac_one{1, 1};
Fraction mul_dot(const Point3fra &p1, const Point3fra &p2) {
  return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}
Point3fra mul_cross(const Point3fra &p1, const Point3fra &p2) {
  return Point3fra(p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z,
                   p1.x * p2.y - p1.y * p2.x);
}
Point3fra mul_scale(const Point3fra &p1, const Fraction &s) {
  Fraction sc{s.up, s.dn};
  sc.reduce();
  return Point3fra(p1.x * sc, p1.y * sc, p1.z * sc);
}
bool is_segs_intersect(Point3fra A, Point3fra B, Point3fra C, Point3fra D) {
  Vector3fra ac = C - A, ad = D - A, ca = A - C, cb = B - C;
  Vector3fra nm_abc = mul_cross(B - A, ac);
  Vector3fra nm_abd = mul_cross(B - A, ad);
  Vector3fra nm_acd = mul_cross(D - C, ca);
  Vector3fra nm_bcd = mul_cross(D - C, cb);
  bool flg1 = mul_dot(nm_abc, nm_abd) < frac_zero && mul_cross(nm_abc, nm_abd).norm2() == frac_zero;
  bool flg2 = mul_dot(nm_acd, nm_bcd) < frac_zero && mul_cross(nm_acd, nm_bcd).norm2() == frac_zero;
  return flg1 && flg2;
}
Fraction point_to_point2(Point3fra A, Point3fra B) { return (A - B).norm2(); }
Fraction point_to_seg2(Point3fra P, Point3fra A, Point3fra B) {
  if (A == B) return point_to_point2(P, A);
  Vector3fra ap = P - A, ab = B - A, bp = P - B, ba = A - B;
  if (mul_dot(ap, ab) <= frac_zero || mul_dot(bp, ba) <= frac_zero) {
    Fraction ret = point_to_point2(P, A);
    ret = min(ret, point_to_point2(P, B));
    return ret;
  } else {
    Vector3fra pa = A - P, pb = B - P, ab = B - A;
    Fraction up = mul_cross(pa, pb).norm2(), dn = ab.norm2();
    return Fraction{up, dn};
  }
}
Fraction seg_to_seg2(Point3fra A, Point3fra B, Point3fra C, Point3fra D) {
  Vector3fra ca = A - C, cb = B - C, cd = D - C, ab = B - A, ac = C - A;
  Fraction tmp = mul_dot(mul_cross(ca, cb), cd);
  bool is_intersec = is_segs_intersect(A, B, C, D);
  if (tmp == frac_zero || is_intersec) {
    // same plane or intersect
    if (is_intersec) return frac_zero;
    Fraction ret = point_to_seg2(A, C, D);
    ret = min(ret, point_to_seg2(B, C, D));
    ret = min(ret, point_to_seg2(C, A, B));
    ret = min(ret, point_to_seg2(D, A, B));
    return ret;
  } else {
    // not in same plane, using maxima of two-variable function
    Fraction dn = mul_dot(ab, cd) * mul_dot(ab, cd) - ab.norm2() * cd.norm2();
    Fraction t(ab.norm2() * mul_dot(cd, ac) - mul_dot(ab, cd) * mul_dot(ab, ac),
               dn);
    Fraction s(mul_dot(ab, cd) * mul_dot(cd, ac) - cd.norm2() * mul_dot(ab, ac),
               dn);
    t.reduce();
    s.reduce();
    if (frac_zero < t && t < frac_one && frac_zero < s && s < frac_one) {
      return point_to_point2(A + mul_scale(ab, s), C + mul_scale(cd, t));
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
  ios::sync_with_stdio(false);
  cin.tie(0);
  cout.tie(0);
  int t = 1;
  cin >> t;
  while (t--) {
    ll ax, ay, az, bx, by, bz;
    ll cx, cy, cz, dx, dy, dz;
    cin >> ax >> ay >> az >> bx >> by >> bz;
    cin >> cx >> cy >> cz >> dx >> dy >> dz;
    Point3fra A{{ax, 1}, {ay, 1}, {az, 1}};
    Point3fra B{{bx, 1}, {by, 1}, {bz, 1}};
    Point3fra C{{cx, 1}, {cy, 1}, {cz, 1}};
    Point3fra D{{dx, 1}, {dy, 1}, {dz, 1}};
    Fraction ans = seg_to_seg2(A, B, C, D);
    ans.reduce();
    cout << abs(ans.up) << " " << abs(ans.dn) << endl;
  }
  return 0;
}
