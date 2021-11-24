/**
 * segment tree for segment maximum
 * single update
 * @note: ALL the range is [l, r]
 * @author jason
 * @date 2021-07-19
 */
#include <algorithm>
#include <cstdio>
using namespace std;
using ll = long long;

const int MAXN = 200005;
const int INF = 0x3f3f3f3f;
/**
 * @param n num size
 * @param m query times
 * "U a b": change num[a] to b
 * "Q a b": get the maximum of range [a, b]
 */
int n, m;
int num[MAXN], seg[MAXN << 2];
/**
 * pick from two children
 * to update node k
 */
void PushUp(int k) {
  seg[k] = max(seg[k * 2], seg[k * 2 + 1]);
}
/**
 * build the tree upon a data array
 * RANGE IS [l, r]
 */
void BuildTree(int k, int l, int r) {
  // [l, r] is empty if l == r
  if (l == r) {
    seg[k] = num[l];
    return;
  }
  int mid = l + (r - l) / 2;
  BuildTree(k * 2, l, mid);
  BuildTree(k * 2 + 1, mid + 1, r);
  PushUp(k);
}
/**
 * update one point value
 * user's call: Update(1, 1, n + 1, pos, val)
 */
void Update(int k, int l, int r, int pos, int val) {
  if (l == r) {
    seg[k] = num[pos] = val;
    return;
  }
  int mid = l + (r - l) / 2;
  if (pos <= mid) {
    Update(k * 2, l, mid, pos, val);
  } else {
    Update(k * 2 + 1, mid + 1, r, pos, val);
  }
  PushUp(k);
}
int Query(int k, int l, int r, int ql, int qr) {
  int ret = -INF;
  if (ql <= l && qr >= r) {
    // query range covers current range
    ret = seg[k];
    return ret;
  }
  int mid = l + (r - l) / 2;
  // avoid unnecessary query when only half range
  if (ql <= mid) ret = max(ret, Query(k * 2, l, mid, ql, qr));
  if (qr > mid) ret = max(ret, Query(k * 2 + 1, mid + 1, r, ql, qr));
  // int vall = Query(k * 2, l, mid, ql, qr);
  // int valr = Query(k * 2 + 1, mid, r, ql, qr);
  // ret = max(vall, valr);
  return ret;
}
int main() {
  while (scanf("%d %d", &n, &m) != EOF) {
    for (int i = 1; i <= n; ++i) scanf("%d", &num[i]);
    /////////////////////////////
    BuildTree(1, 1, n);
    //////////////////////////////
    char oprt;
    int oprd1, oprd2;
    getchar();
    for (int _ = 0; _ < m; ++_) {
      scanf("%c %d %d", &oprt, &oprd1, &oprd2);
      if (oprt == 'U') Update(1, 1, n, oprd1, oprd2);
      if (oprt == 'Q') printf("%d\n", Query(1, 1, n, oprd1, oprd2));
      getchar();
      // printf("done\n");
    }
  }
  return 0;
}
