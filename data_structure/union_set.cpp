#include <iostream>
using namespace std;
const int MAXN = 100005;
int father[MAXN];
int trank[MAXN];

void Init(int n) {
  for (int i = 0; i < n; ++i) {
    father[i] = i;
    trank[i] = 0;
  }
}
int Find(int x) {
  if (father[x] == x) {
    return x;
  }
  return father[x] = Find(father[x]);
}
void Unite(int x, int y) {
  x = Find(x);
  y = Find(y);
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
bool inSame(int x, int y) { return Find(x) == Find(y); }
