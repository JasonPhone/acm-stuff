/**
 * block partition for range modification and sum
 * index of dot starts form 1
 * note:
 *    sum_of_dots = block_mark + dot_value
 *    sum_of_block = block_sum + (block_mark * block_size)
 */
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const int MAXN = 100005, SQ = 1003;
ll n, block_size[SQ], block_st[SQ], block_ed[SQ], in_block[MAXN], mark[SQ], S[SQ];
ll sq;
ll A[MAXN], m;
void init_block() {
  sq = sqrt(n);
  for (int i = 1; i <= sq; i++) {
    // block i is [st[i], ed[i]] in range
    block_st[i] = n / sq * (i - 1) + 1;
    block_ed[i] = n / sq * i;
  }
  // if n is not square
  block_ed[sq] = n;
  // map the dot index to block index
  for (int i = 1; i <= sq; i++)
    for (int j = block_st[i]; j <= block_ed[i]; j++) in_block[j] = i;
  // size of every block
  for (int i = 1; i <= sq; i++) block_size[i] = block_ed[i] - block_st[i] + 1;
}
void range_add(int adl, int adr, ll adv) {
  if (in_block[adl] == in_block[adr]) {
    // same block, enumerate in O(sqrt(n))
    for (int i = adl; i <= adr; i++) {
      A[i] += adv;
      S[in_block[i]] += adv;
    }
  } else {
    // across block boundry, NOTE THE <= and < operators
    // left dots
    for (int i = adl; i <= block_ed[in_block[adl]]; i++) {
      A[i] += adv;
      S[in_block[i]] += adv;
    }
    // right dots
    for (int i = block_st[in_block[adr]]; i <= adr; i++) {
      A[i] += adv;
      S[in_block[i]] += adv;
    }
    // mid block(s)
    for (int i = in_block[adl] + 1; i < in_block[adr]; i++) mark[i] += adv;
  }
}
ll range_query(int ql, int qr) {
  ll ret = 0;
  if (in_block[ql] == in_block[qr]) {
    for (int i = ql; i <= qr; i++) ret += A[i] + mark[in_block[i]];
  } else {
    for (int i = ql; i <= block_ed[in_block[ql]]; i++) ret += A[i] + mark[in_block[i]];
    for (int i = block_st[in_block[qr]]; i <= qr; i++) ret += A[i] + mark[in_block[i]];
    for (int i = in_block[ql] + 1; i < in_block[qr]; i++) ret += S[i] + mark[i] * block_size[i];
  }
  return ret;
}
int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);
  cout.tie(0);
  cin >> n >> m;
  init_block();
  for (int i = 1; i <= n; i++) cin >> A[i];
  // S[i] for sum of block i
  for (int i = 1; i <= sq; i++)
    for (int j = block_st[i]; j <= block_ed[i]; j++) S[i] += A[j];
  while (m--) {
    ll op, opl, opr, val;
    cin >> op;
    if (op == 1) {
      // range add
      cin >> opl >> opr >> val;
      range_add(opl, opr, val);
    } else {
      // range query
      cin >> opl >> opr;
      cout << range_query(opl, opr) << endl;
    }
  }
  return 0;
}
