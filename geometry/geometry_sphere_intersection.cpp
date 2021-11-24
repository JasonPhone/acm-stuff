#include <bits/stdc++.h>
using namespace std;
const double PI = acos(-1.0);
struct Sphere {
  double x, y, z, r;
  Sphere() {}
  Sphere(double x, double y, double z, double r) : x(x), y(y), z(z), r(r) {}
};
double IntersectionVolume(Sphere o, Sphere t) {
  // basic formula: V = (3 * r - h) * h * h * PI / 3
  // calculated from spinning surface calculus
  if (o.r < t.r) swap(o, t);
  double dis = sqrt((o.x - t.x) * (o.x - t.x) + (o.y - t.y) * (o.y - t.y) +
                    (o.z - t.z) * (o.z - t.z));
  if (dis <= o.r - t.r) {  // completely in
    return 4.0 / 3 * PI * t.r * t.r * t.r;
  } else if (dis <= o.r) {  // center of the smaller sphere in bigger sphere
    // cosA = (b2 + c2 - a2) / 2bc
    double angleb = acos((t.r * t.r + dis * dis - o.r * o.r) / (2 * t.r * dis));
    double anglea = PI - angleb;
    double l = t.r * cos(anglea);
    double H = o.r - l - dis;
    double h = t.r - l;
    return 4.0 / 3 * PI * t.r * t.r * t.r - PI / 3 * (3 * t.r - h) * h * h +
           PI / 3 * (3 * o.r - H) * H * H;
  } else if (dis < o.r + t.r) {  // normal intersection
    double angler = acos((t.r * t.r + dis * dis - o.r * o.r) / (2 * t.r * dis));
    double angleR = acos((o.r * o.r + dis * dis - t.r * t.r) / (2 * o.r * dis));
    double H = o.r - o.r * cos(angleR);
    double h = t.r - t.r * cos(angler);
    return PI / 3 * (3 * t.r - h) * h * h + PI / 3 * (3 * o.r - H) * H * H;
  } else {
    return 0;
  }
}
double IntersectionSurface(Sphere &o, Sphere &t) {
  // basic formula: S = 2 * PI * r * h
  if (o.r < t.r) swap(o, t);
  double dis = sqrt((o.x - t.x) * (o.x - t.x) + (o.y - t.y) * (o.y - t.y) +
                    (o.z - t.z) * (o.z - t.z));
  if (dis <= o.r - t.r) {  // completely in
    return 4 * PI * t.r * t.r;
  } else if (dis <= o.r) {  // center of the smaller sphere in bigger sphere
    double angleb = acos((t.r * t.r + dis * dis - o.r * o.r) / (2 * t.r * dis));
    double anglea = PI - angleb;
    double l = t.r * cos(anglea);
    double H = o.r - l - dis;
    double h = t.r - l;
    return 4 * PI * t.r * t.r - 2 * PI * t.r * h + 2 * PI * o.r * H;
  } else if (dis < o.r + t.r) {  // normal intersection
    double angler = acos((t.r * t.r + dis * dis - o.r * o.r) / (2 * t.r * dis));
    double angleR = acos((o.r * o.r + dis * dis - t.r * t.r) / (2 * o.r * dis));
    double H = o.r - o.r * cos(angleR);
    double h = t.r - t.r * cos(angler);
    return 2 * PI * t.r * h + 2 * PI * o.r * H;
  } else {
    return 0;
  }
}
int main() {
  Sphere A, B;
  cin >> A.x >> A.y >> A.z >> A.r;
  cin >> B.x >> B.y >> B.z >> B.r;
  cout << fixed << setprecision(10)
       << 4 * PI * (A.r * A.r + B.r * B.r) - IntersectionSurface(A, B) << endl;
  return 0;
}
