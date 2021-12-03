/**
 * https://codeforces.com/contest/600/problem/E
 * 树的节点有权，根为 1
 * 一种权占领了一个子树
 * 当且仅当没有其他权在这个子树中出现更多次
 * 求占领每个子树的所有权之和
 * 输入：
 * 节点数
 * 各节点的权
 * 边
 * 输出：
 * 各节点的占领权之和
 *************************
 * 每个节点的答案是其子树的叠加，利用这个性质处理问题
 * 预处理出每个节点子树的 size 和它的重儿子(节点最多子树的儿子)，可以O(n)完成
 * 用 check[i] 表示颜色 i 有没有出现过，ans[i] 表示出现次数
 * 按以下的步骤遍历一个节点：
 * 遍历其非重儿子，获取它的 ans，但不保留遍历后它的 check
 * 遍历它的重儿子，保留它的 check
 * 再次遍历其非重儿子及其父亲，用重儿子的 check
 * 对遍历到的节点进行计算，获取整棵子树的 ans
 */
#include <bits/stdc++.h>
using namespace std;
const int MAXN = 1e5 + 100;
int n, a[MAXN], tot = -1;
int head[MAXN], to[MAXN << 1], nxt[MAXN << 1];
int bson[MAXN], sz[MAXN];
long long ans[MAXN], sum;
int maxc, flag;
int clr[MAXN];
void add(int u, int v) {
  // 链式前向星
  nxt[++tot] = head[u];
  head[u] = tot;
  to[tot] = v;
  nxt[++tot] = head[v];
  head[v] = tot;
  to[tot] = u;
}
void dfs(int u, int f) {
  sz[u] = 1;
  for (int pp = head[u]; pp != -1; pp = nxt[pp]) {
    int nxt_id = to[pp];
    if (nxt_id == f) continue;
    dfs(nxt_id, u);
    sz[u] += sz[nxt_id];
    if (sz[nxt_id] > sz[bson[u]]) bson[u] = nxt_id;
  }
}
void add(int u, int f, int val) {
  clr[a[u]] += val;
  if (clr[a[u]] > maxc) {
    maxc = clr[a[u]];
    /**********  ans  **********/
    sum = a[u];
    /***************************/
  } else if (clr[a[u]] == maxc) {
    /**********  ans  **********/
    sum += a[u];
    /***************************/
  }
  for (int pp = head[u]; pp != -1; pp = nxt[pp]) {
    int nxt_id = to[pp];
    if (nxt_id == flag || nxt_id == f) continue;
    add(nxt_id, u, val);
  }
}
void dfs(int u, int f, bool keep) {
  for (int pp = head[u]; pp != -1; pp = nxt[pp]) {
    int nxt_id = to[pp];
    if (nxt_id == f || nxt_id == bson[u]) continue;
    dfs(nxt_id, u, 0);
  }
  if (bson[u]) {
    dfs(bson[u], u, 1);
    flag = bson[u];
  }
  add(u, f, 1);
  flag = 0;
  /**********  ans  **********/
  ans[u] = sum;
  /***************************/
  if (!keep) {
    add(u, f, -1);
    /**********  ans  **********/
    maxc = sum = 0;
    /***************************/
  }
}
void solve() {
  int u, v;
  // fill(head+1,head+n+2,-1);
  cin >> n;
  fill(head, head + n + 2, -1);
  for (int i = 1; i <= n; ++i) cin >> a[i];
  for (int i = 1; i < n; ++i) {
    cin >> u >> v;
    add(u, v);
  }
  dfs(1, -1);
  dfs(1, -1, 0);
  for (int i = 1; i < n; ++i) cout << ans[i] << " ";
  cout << ans[n] << "\n";
}
int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);
  cout.tie(0);
  solve();
  return 0;
}
