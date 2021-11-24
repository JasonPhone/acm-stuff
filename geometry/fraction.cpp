#include <bits/stdc++.h>
using namespace std;
using ll = long long;
struct Fraction {
  ll up, dn;
  Fraction() : up(0), dn(1) {}
  Fraction(ll _up, ll _dn) : up(_up), dn(_dn) {
    ll cd = __gcd(up, dn);
    up /= cd;
    dn /= cd;
    if (dn < 0) {
      dn = -dn;
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
    ////////////// possible overflow //////////////////
    ll n_up = n_dn / dn * up + n_dn / otr.dn * otr.up;
    return Fraction(n_up, n_dn);
  }
  Fraction operator-(const Fraction &otr) const {
    ll n_dn = dn / __gcd(otr.dn, dn) * otr.dn;
    ////////////// possible overflow //////////////////
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
int main() {
  Fraction a(1, 2), b(3, 6);
  cout << (a * b).real_val() << endl;
  cout << (a - b).real_val() << endl;
  cout << (a == b) << endl;
  return 0;
}
