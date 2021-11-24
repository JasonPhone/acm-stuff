// start from 1
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MAXN = 100005;
ll tree[MAXN];
ll lowbit(int x) { return (x) & (-x); };
void Update(int i, ll x) {
  // increase
  for (int pos = i; pos <= MAXN; pos += lowbit(pos)) {
    tree[pos] += x;
  }
}
ll PrefixQuery(int n) {
  ll ret = 0;
  for (int pos = n; pos; pos -= lowbit(pos)) {
    ret += tree[pos];
  }
  return ret;
}
ll RangeQuery(int ql, int qr) { return PrefixQuery(qr) - PrefixQuery(ql - 1); }
int main() {
  int a[10] = {-1, 4, 2, 1, 5, 6, 7, 2, 1, 4};
  for (int i = 1; i <= 9; i++) {
    Update(i, a[i]);
  }
  for (int i = 1; i <= 9; i++) {
    cout << PrefixQuery(i) << endl;
  }
  return 0;
}
