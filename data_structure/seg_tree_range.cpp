/**
 * https://www.luogu.com.cn/problem/P3373
 * segment tree for segment sum
 * range update (multiply and add)
 * @note all the range is [l, r]
 * @author jason
 * @date 2021-07-19
 */
#include <iostream>
using namespace std;
using ll = long long;
const int MAXN = 200005;

struct Node {
  // TODO modify to fit the need
  ll l, r;
  ll ans, mulv, addv;
  Node() {}
};
Node tree[MAXN << 2];
ll n, m, q, rawValues[MAXN];

void MergeNode(Node &f, const Node &lc, const Node &rc) {
  // TODO VARY based on different problems
  f.ans = (lc.ans + rc.ans) % m;
  f.addv = 0;
  f.mulv = 1;
}
void NodeAdd(int k, ll addv) {

}
void NodeMul(int k, ll mulv) {

}
void SpreadTag(Node &f, Node &sn) {
  // TODO VARY based on different problems
  ll addv = f.addv, mulv = f.mulv;
  sn.ans = (sn.ans * mulv % m + (sn.r - sn.l + 1) % m * addv % m) % m;
  sn.mulv = sn.mulv * mulv % m;
  sn.addv = (sn.addv * mulv % m + addv) % m;
}
void PushUp(int k) {  // up a level
  MergeNode(tree[k], tree[k << 1], tree[k << 1 | 1]);
}
void PushDown(int k) {  // push the lazy tag down a level
  if (!(tree[k].addv == 0 && tree[k].mulv == 1)) {
    SpreadTag(tree[k], tree[k << 1]);
    SpreadTag(tree[k], tree[k << 1 | 1]);
    // TODO reset father's lazy tag
    tree[k].addv = 0;
    tree[k].mulv = 1;
  }
}
void BuildTree(int k, int l, int r) {
  // prepare the nodes
  tree[k].l = l;
  tree[k].r = r;
  if (l == r) {
    // TODO VARY based on different problems
    tree[k].ans = rawValues[l];
    tree[k].addv = 0;
    tree[k].mulv = 1;
  } else {
    int mid = l + (r - l) / 2;
    BuildTree(k << 1, l, mid);
    BuildTree(k << 1 | 1, mid + 1, r);
    PushUp(k);
  }
}
void UpdateSegMul(int k, int l, int r, ll mulv) {
  if (l <= tree[k].l && tree[k].r <= r) {
    // TODO VARY based on problems
    // record the operation for query with smaller range
    tree[k].ans = tree[k].ans * mulv % m;
    tree[k].mulv = tree[k].mulv * mulv % m;
    tree[k].addv = tree[k].addv * mulv % m;
  } else {
    PushDown(k);
    int mid = tree[k].l + (tree[k].r - tree[k].l) / 2;
    if (mid >= l)  // separated update
      UpdateSegMul(k << 1, l, r, mulv);
    if (mid < r) UpdateSegMul(k << 1 | 1, l, r, mulv);
    PushUp(k);
  }
}
void UpdateSegAdd(int k, int l, int r, ll addv) {
  if (l <= tree[k].l && tree[k].r <= r) {
    // TODO VARY based on problems
    tree[k].ans = (tree[k].ans + addv * (tree[k].r - tree[k].l + 1) % m) % m;
    tree[k].addv = (tree[k].addv + addv) % m;
  } else {
    PushDown(k);
    int mid = tree[k].l + (tree[k].r - tree[k].l) / 2;
    if (mid >= l)  // separated update
      UpdateSegAdd(k << 1, l, r, addv);
    if (mid < r) UpdateSegAdd(k << 1 | 1, l, r, addv);
    PushUp(k);
  }
}
void UpdateDot(int k, int pos, ll val) {
  if (tree[k].l == tree[k].r) {
    // TODO VARY based on problems
    // tree[k].sum = val;
  } else {
    PushDown(k);
    int mid = tree[k].l + (tree[k].r - tree[k].l) / 2;
    if (pos <= mid)  // separated update
      UpdateDot(k << 1, pos, val);
    else
      UpdateDot(k << 1 | 1, pos, val);
    PushUp(k);
  }
}
Node Query(int k, int ql, int qr) {
  if (tree[k].l >= ql && tree[k].r <= qr) return tree[k];
  // when not single, push down firstly, then do the query
  PushDown(k);
  int mid = tree[k].l + (tree[k].r - tree[k].l) / 2;
  Node resL, resR, retVal;
  bool hasL = false, hasR = false;
  if (ql <= mid) {
    hasL = true;
    resL = Query(k << 1, ql, qr);
  }
  if (mid < qr) {
    hasR = true;
    resR = Query(k << 1 | 1, ql, qr);
  }
  if (hasL && hasR)
    MergeNode(retVal, resL, resR);
  else if (hasL)
    retVal = resL;
  else if (hasR)
    retVal = resR;
  return retVal;
}
int main() {
  ios::sync_with_stdio(false);
  cin >> n >> q >> m;
  for (int i = 1; i <= n; i++) cin >> rawValues[i];
  ///////////////////////////////
  BuildTree(1, 1, n);
  ///////////////////////////////
  int t, l, r, v;
  while (q--) {
    cin >> t >> l >> r;
    if (t == 3) {
      cout << Query(1, l, r).ans << "\n";
    } else if (t == 1) {
      cin >> v;
      UpdateSegMul(1, l, r, v);
    } else if (t == 2) {
      cin >> v;
      UpdateSegAdd(1, l, r, v);
    }
  }
  return 0;
}
