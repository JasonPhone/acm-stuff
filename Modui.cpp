/**
 * Modui range number of distinct values
 */
#include <bits/stdc++.h>
using namespace std;
#define endl "\n";
#define IOS_ONLY               \
  ios::sync_with_stdio(false); \
  cin.tie(0);                  \
  cout.tie(0);
const int MAXN = 30005, MAXQ = 200005, MAXM = 1000005;
int sq;
struct Query {
  int ql, qr, id;
  bool operator<(const Query &o) const {
    // sqrt(n) partitions, assign sq with sqrt(n) first
    if (ql / sq != o.ql / sq) return ql < o.ql;
    if (ql / sq & 1) return qr < o.qr;  // order by parity
    return qr > o.qr;
  }
} Q[MAXQ];
int A[MAXN], ans[MAXQ], Cnt[MAXM], cur, pl = 1, pr = 0, n;
inline void add(int pos) {
  if (Cnt[A[pos]] == 0) cur++;
  Cnt[A[pos]]++;
}
inline void del(int pos) {
  Cnt[A[pos]]--;
  if (Cnt[A[pos]] == 0) cur--;
}
int main() {
  IOS_ONLY
  cin >> n;
  sq = sqrt(n);
  for (int i = 1; i <= n; ++i) cin >> A[i];
  int q;
  cin >> q;
  for (int i = 0; i < q; ++i) {  // offline query
    cin >> Q[i].ql >> Q[i].qr;
    Q[i].id = i;
  }
  sort(Q, Q + q);  // sort, KEY of modui
  for (int i = 0; i < q; ++i) {
    while (pl > Q[i].ql) add(--pl);
    while (pr < Q[i].qr) add(++pr);
    while (pl < Q[i].ql) del(pl++);
    while (pr > Q[i].qr) del(pr--);
    ans[Q[i].id] = cur;  // store the rasult
  }
  for (int i = 0; i < q; ++i) cout << ans[i] << endl;
  return 0;
}
