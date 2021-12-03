// 改进版替罪羊树，在另外一些细节上也进行了一些更改，具体看注释
/**
 * 插入一个整数 x。
 * 删除一个整数 x（若有多个相同的数，只删除一个）。
 * 查询整数 x 的排名（排名定义为比当前数小的数的个数 +1）。
 * 查询排名为 x 的数（如果不存在，则认为是排名小于 x 的最大数。保证 x 不会超过当前数据结构中数的总数）。
 * 求 x 的前驱（小于 x，且最大的数）。
 * 求 x 的后继（大于 x，且最小的数）。
 */
#include <bits/stdc++.h>
using namespace std;
#define ls(x) tree[x].ls
#define rs(x) tree[x].rs
#define num(x) tree[x].num
#define val(x) tree[x].val
#define sz(x) tree[x].sz
#define exist(x) !(num(x) == 0 && ls(x) == 0 && rs(x) == 0)
const double ALPHA = 0.7;
const int MAXN = 2e6 + 5;
int n, m;
struct Node {
  int ls, rs, num, val, sz;
} tree[MAXN];            // 改用结构体进行存储
vector<int> FP, FN, FV;  // 存储拉平后的节点编号、数目、值
int cnt = 1;
// 一趟中序遍历，把当前子树拉平并存到 vector 里，返回当前节点的索引
int flatten(int pos) {
  if (exist(ls(pos)))  // 递归地拉平左子树
    flatten(ls(pos));
  int id = FP.size();  // 记下当前节点的索引
  // 如果该节点是已被删除的节点，就略过，否则把相应信息存入 vector
  if (num(pos) != 0) {
    FP.push_back(pos);
    FV.push_back(val(pos));
    FN.push_back(num(pos));
  }
  // 递归地拉平右子树
  if (exist(rs(pos))) flatten(rs(pos));
  return id;
}
// 以 pos 为根节点，以 [l,r] 内的信息重建一棵平衡的树
void rebuild(int pos, int l = 0, int r = FP.size() - 1) {
  int mid = (l + r) / 2, sz1 = 0, sz2 = 0;
  if (l < mid) {
    ls(pos) = FP[(l + mid - 1) / 2];  // 重用节点编号
    rebuild(ls(pos), l, mid - 1);     // 递归地重建
    sz1 = sz(ls(pos));
  } else {
    ls(pos) = 0;
  }
  if (mid < r) {
    rs(pos) = FP[(mid + 1 + r) / 2];
    rebuild(rs(pos), mid + 1, r);
    sz2 = sz(rs(pos));
  } else {
    rs(pos) = 0;
  }
  num(pos) = FN[mid];  // 把存于 vector 中的信息复制过来
  val(pos) = FV[mid];
  sz(pos) = sz1 + sz2 + num(pos);  // 递归确定重建后树的大小
}
// 尝试重构当前子树
void try_restructure(int pos) {
  double k = max(sz(ls(pos)), sz(rs(pos))) / double(sz(pos));
  if (k > ALPHA) {
    FP.clear(), FV.clear(), FN.clear();  // 清空 vector
    int id = flatten(pos);
    // 这里是确保当前节点的编号在重构后不会改变
    swap(FP[id], FP[(FP.size() - 1) / 2]);
    rebuild(pos);
  }
}
// 接下来是普通的二叉查找树
void bst_insert(int v, int pos = 1) {
  if (!exist(pos)) {
    val(pos) = v;
    num(pos) = 1;
  } else if (v < val(pos)) {
    if (!exist(ls(pos))) ls(pos) = ++cnt;
    bst_insert(v, ls(pos));
  } else if (v > val(pos)) {
    if (!exist(rs(pos))) rs(pos) = ++cnt;
    bst_insert(v, rs(pos));
  } else
    num(pos)++;
  sz(pos)++;
  try_restructure(pos);
}
void bst_remove(int v, int pos = 1) {
  sz(pos)--;
  if (v < val(pos))
    bst_remove(v, ls(pos));
  else if (v > val(pos))
    bst_remove(v, rs(pos));
  else
    num(pos)--;
  try_restructure(pos);
}
int bst_countl(int v, int pos = 1) {
  if (v < val(pos))
    return exist(ls(pos)) ? bst_countl(v, ls(pos)) : 0;
  else if (v > val(pos))
    return sz(ls(pos)) + num(pos) + (exist(rs(pos)) ? bst_countl(v, rs(pos)) : 0);
  else
    return sz(ls(pos));
}
int bst_countg(int v, int pos = 1) {
  if (v > val(pos))
    return exist(rs(pos)) ? bst_countg(v, rs(pos)) : 0;
  else if (v < val(pos))
    return sz(rs(pos)) + num(pos) + (exist(ls(pos)) ? bst_countg(v, ls(pos)) : 0);
  else
    return sz(rs(pos));
}
int bst_rank(int v) { return bst_countl(v) + 1; }
int bst_kth(int k, int pos = 1) {
  if (sz(ls(pos)) + 1 > k)
    return bst_kth(k, ls(pos));
  else if (sz(ls(pos)) + num(pos) < k)
    return bst_kth(k - sz(ls(pos)) - num(pos), rs(pos));
  else
    return val(pos);
}
int bst_pre(int v) {
  int r = bst_countl(v);
  return bst_kth(r);
}
int bst_suc(int v) {
  int r = sz(1) - bst_countg(v) + 1;
  return bst_kth(r);
}
int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);
  cout.tie(0);
  cin >> n >> m;
  for (int i = 0; i < n; i++) {
    int a;
    cin >> a;
    bst_insert(a);
  }
  int lasta = 0;
  vector<int> res;
  while (m--) {
    int op, x;
    cin >> op >> x;
    x ^= lasta;
    if (op == 1) // insert
      bst_insert(x);
    else if (op == 2)  // delete
      bst_remove(x);
    else if (op == 3)  // rank
      lasta = bst_rank(x);
    else if (op == 4)  // k-th
      lasta = bst_kth(x);
    else if (op == 5)  // pre
      lasta = bst_pre(x);
    else if (op == 6)  // suc
      lasta = bst_suc(x);
    if (op > 2) {
      res.push_back(lasta);
    }
  }
  int ans = 0;
  for (auto v : res) ans ^= v;
  cout << ans << endl;
  return 0;
}
