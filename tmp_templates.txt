/**
 * 单调栈求凸包
 *
 * stk[] 是整型，存的是下标
 * vector<int, int> p[] 存储向量或点
 */
tp = 0;                       // 初始化栈
std::sort(p + 1, p + 1 + n);  // 对点进行排序
stk[++tp] = 1;
//栈内添加第一个元素，且不更新 used，使得 1 在最后封闭凸包时也对单调栈更新
for (int i = 2; i <= n; ++i) {
  while (tp >= 2  // 下一行 * 操作符被重载为叉积
         && (p[stk[tp]] - p[stk[tp - 1]]) * (p[i] - p[stk[tp]]) <= 0)
    used[stk[tp--]] = 0;
  used[i] = 1;  // used 表示在凸壳上
  stk[++tp] = i;
}
int tmp = tp;  // tmp 表示下凸壳大小
for (int i = n - 1; i > 0; --i)
  if (!used[i]) {
    //      ↓求上凸壳时不影响下凸壳
    while (tp > tmp && (p[stk[tp]] - p[stk[tp - 1]]) * (p[i] - p[stk[tp]]) <= 0)
      used[stk[tp--]] = 0;
    used[i] = 1;
    stk[++tp] = i;
  }
for (int i = 1; i <= tp; ++i)  // 复制到新数组中去
  h[i] = p[stk[i]];
int ans = tp - 1;

//////////////////////////////////////////////////
/**
 * query the number of elements equal to val
 * which is previously neighbouring pos (inclusive)
 * call (1, 5, 0) for [1, 1, 1, 0, 0, 0, 0, 0] (start from 1)
 * will get 2
 */
int query_prefix_num(int k, int pos, int val) {
  // val == -1 if seg under this node all not all same
  if (tree[k].val == val) return min(pos, tree[k].r) - tree[k].l + 1;
  if (tree[k].l == tree[k].r) return tree[k].val == val;
  push_down(k);
  int mid = tree[k].l + (tree[k].r - tree[k].l) / 2;
  if (pos > mid) {
    int len = query_prefix_num(k << 1 | 1, pos, val);
    if (len == min(pos, tree[k << 1 | 1].r) - tree[k << 1 | 1].l + 1)
      len += query_prefix_num(k << 1, pos, val);
    return len;
  }
  return query_prefix_num(k << 1, pos, val);
}

////////////////////////////////////////////////////
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
struct Point3 {
  int x, y, z;
  Point3(int x = 0, int y = 0, int z = 0) : x(x), y(y), z(z) {}
  bool operator<(const Point3& u) const {
    return x - u.x < 0 || (x - u.x == 0 && y - u.y < 0) ||
           (x - u.x == 0 && y - u.y == 0 && z - u.z < 0);
  }
  bool operator>(const Point3& u) const { return u < (*this); }
  bool operator==(const Point3& u) const {
    return !(u < (*this) || (*this) < u);
  }
  bool operator!=(const Point3& u) const { return !((*this) == u); }
  bool operator<=(const Point3& u) const { return *this < u || *this == u; }
  bool operator>=(const Point3& u) const { return *this > u || *this == u; }
  Point3 operator+(const Point3& u) const {
    return Point3(x + u.x, y + u.y, z + u.z);
  }
  Point3 operator-(const Point3& u) const {
    return Point3(x - u.x, y - u.y, z - u.z);
  }
  Point3 operator*(const int u) const { return Point3(x * u, y * u, z * u); }
  Point3 operator/(const int u) const { return Point3(x / u, y / u, z / u); }
  void read() { scanf("%d%d%d", &x, &y, &z); }
};

typedef Point3 Vector3;
ll getDot(Vector3 a, Vector3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vector3 getCross(Vector3 a, Vector3 b) {
  return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x);
}
ll getPowerLength(Vector3 u) { return getDot(u, u); }
ll gcd(ll a, ll b) { return b == 0 ? a : gcd(b, a % b); }
struct Rat {
  ll s, m;
  Rat(ll s = 0, ll m = 1) {
    ll d = gcd(s, m);
    s /= d, m /= d;
    if (m < 0) m = -m, s = -s;
    this->s = s;
    this->m = m;
  }
  Rat operator+(const Rat& u) const {
    ll d = gcd(m, u.m);
    return Rat(s * (u.m / d) + u.s * (m / d), m * (u.m / d));
  }
  Rat operator-(const Rat& u) const {
    ll d = gcd(m, u.m);
    return Rat(s * (u.m / d) - u.s * (m / d), m * (u.m / d));
  }
  Rat operator*(const Rat& u) const { return Rat(s * u.s, m * u.m); }
  // Rat operator * (const int& u) const { return Rat(s*u, m); }
  // Rat operator / (const Rat& u) const { return Rat(s*u.m, m*u.s); }
  // Rat operator / (const int& u) const { return Rat(s, m*u); }
  bool operator<(const Rat& u) const { return s * u.m < u.s * m; }
  bool operator>(const Rat& u) const { return u < (*this); }
  bool operator==(const Rat& u) const { return !(u < (*this) || (*this) < u); }
  bool operator!=(const Rat& u) const { return !((*this) == u); }
  bool operator<=(const Rat& u) const { return *this < u || *this == u; }
  bool operator>=(const Rat& u) const { return *this > u || *this == u; }
};
inline int dcmp(Rat u) {
  if (u.s == 0)
    return 0;
  else
    return u.s < 0 ? -1 : 1;
}
Rat getDistancePointToSegment(Point3 p, Point3 a, Point3 b) {
  if (a == b) return getPowerLength(p - a);
  Vector3 v1 = b - a, v2 = p - a, v3 = p - b;
  if (getDot(v1, v2) < 0)
    return getPowerLength(v2);
  else if (getDot(v1, v3) > 0)
    return getPowerLength(v3);
  else
    return Rat(getPowerLength(getCross(v1, v2)), getPowerLength(v1));
}
bool getDistanceLineToLine(Point3 p1, Vector3 u, Point3 p2, Vector3 v, Rat& s) {
  ll b = getDot(u, u) * getDot(v, v) - getDot(u, v) * getDot(u, v);
  if (b == 0) return false;
  ll a = getDot(u, v) * getDot(v, p1 - p2) - getDot(v, v) * getDot(u, p1 - p2);
  s = Rat(a, b);
  return true;
}
const ll inf = 0x3f3f3f3f;
int main() {
  int cas;
  scanf("%d", &cas);
  while (cas--) {
    Point3 a, b, c, d;
    a.read(), b.read(), c.read(), d.read();
    Rat s, t, ans(inf);
    bool flag1 = getDistanceLineToLine(a, b - a, c, d - c, s);
    bool flag2 = getDistanceLineToLine(c, d - c, a, b - a, t);
    if (flag1 && flag2 && s.s > 0 && s.s < s.m && t.s > 0 && t.s < t.m) {
      Vector3 u = b - a, v = d - c;
      Rat x1 = Rat(a.x) + s * u.x, y1 = Rat(a.y) + s * u.y,
          z1 = Rat(a.z) + s * u.z;
      Rat x2 = Rat(c.x) + t * v.x, y2 = Rat(c.y) + t * v.y,
          z2 = Rat(c.z) + t * v.z;
      ans =
          (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1);
    } else {
      ans = min(ans, getDistancePointToSegment(a, c, d));
      ans = min(ans, getDistancePointToSegment(b, c, d));
      ans = min(ans, getDistancePointToSegment(c, a, b));
      ans = min(ans, getDistancePointToSegment(d, a, b));
    }
    printf("%lld %lld\n", ans.s, ans.m);
  }
  return 0;
}
