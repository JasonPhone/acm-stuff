#include <iostream>
using namespace std;
namespace UNION {
const int MAXN = 100005;
int father[MAXN];
int trank[MAXN];
void init(int n) {
  for (int i = 0; i < n; ++i) {
    father[i] = i;
    trank[i] = 0;
  }
}
int find(int x) {
  if (father[x] == x) {
    return x;
  }
  return father[x] = find(father[x]);
}
void unite(int x, int y) {
  x = find(x);
  y = find(y);
  if (x == y) {
    return;
  }
  if (trank[x] < trank[y]) {
    father[x] = y;
  } else {
    father[y] = x;
    if (trank[x] == trank[y]) {
      trank[x]++;
    }
  }
}
bool in_same(int x, int y) { return find(x) == find(y); }
}  // namespace UNION
