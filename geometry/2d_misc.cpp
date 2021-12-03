const double EPS = 1e-9;
bool eq(double a, double b) { return abs(a - b) < EPS; } // ==
bool gt(double a, double b) { return a - b > EPS; }      // >
bool lt(double a, double b) { return a - b < -EPS; }     // <
bool ge(double a, double b) { return a - b > -EPS; }     // >=
bool le(double a, double b) { return a - b < EPS; }      // <=
// 直线与直线交点
// DEPENDS eq, d*V, V*V, V+V, V^V
vector<Point> inter(Line a, Line b) {
  double c = a.v ^ b.v;
  if (eq(c, 0)) return {};
  Vec v = 1 / c * Vec{a.P ^ (a.P + a.v), b.P ^ (b.P + b.v)};
  return {{v * Vec{-b.v.x, a.v.x}, v * Vec{-b.v.y, a.v.y}}};
}

// 直线与圆交点
// DEPENDS eq, gt, V+V, V-V, V*V, d*V, len, pedal
vector<Point> inter(Line l, Circle C) {
  Point P = pedal(C.O, l);
  double h = len(P - C.O);
  if (gt(h, C.r)) return {};
  if (eq(h, C.r)) return {P};
  double d = sqrt(C.r * C.r - h * h);
  Vec vec = d / len(l.v) * l.v;
  return {P + vec, P - vec};
}

// 圆与圆的交点 注意内含和相离的情况
// DEPENDS eq, gt, V+V, V-V, d*V, len, r90c
vector<Point> inter(Circle C1, Circle C2) {
  Vec v1 = C2.O - C1.O, v2 = r90c(v1);
  double d = len(v1);
  if (gt(d, C1.r + C2.r) || gt(abs(C1.r - C2.r), d)) return {};
  if (eq(d, C1.r + C2.r) || eq(d, abs(C1.r - C2.r)))
    return {C1.O + C1.r / d * v1};
  double a = ((C1.r * C1.r - C2.r * C2.r) / d + d) / 2;
  double h = sqrt(C1.r * C1.r - a * a);
  Vec av = a / len(v1) * v1, hv = h / len(v2) * v2;
  return {C1.O + av + hv, C1.O + av - hv};
}

// 三角形的重心
Point barycenter(Point A, Point B, Point C) {
  return {(A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3};
}

// 三角形的外心
// DEPENDS r90c, V*V, d*V, V-V, V+V
// NOTE 给定圆上三点求圆，要先判断是否三点共线
Point circumcenter(Point A, Point B, Point C) {
  double a = A * A, b = B * B, c = C * C;
  double d = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
  return 1 / d * r90c(a * (B - C) + b * (C - A) + c * (A - B));
}

// 三角形的内心
// DEPENDS len, d*V, V-V, V+V
Point incenter(Point A, Point B, Point C) {
  double a = len(B - C), b = len(A - C), c = len(A - B);
  double d = a + b + c;
  return 1 / d * (a * A + b * B + c * C);
}

// 三角形的垂心
// DEPENDS V*V, d*V, V-V, V^V, r90c
Point orthocenter(Point A, Point B, Point C) {
  double n = B * (A - C), m = A * (B - C);
  double d = (B - C) ^ (A - C);
  return 1 / d * r90c(n * (C - B) - m * (C - A));
}
