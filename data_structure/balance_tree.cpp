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
const int MAXN = 1e5 + 5;
struct Node {
  int ls, rs, num, val, sz;
} tree[MAXN];            // 改用结构体进行存储
vector<int> FP, FN, FV;  // 存储拉平后的节点编号、数目、值
int cnt = 1;
// 一趟中序遍历，把当前子树拉平并存到vector里，返回当前节点的索引
int flatten(int pos) {
  if (exist(ls(pos)))  // 递归地拉平左子树
    flatten(ls(pos));
  int id = FP.size();  // 记下当前节点的索引
  // 如果该节点是已被删除的节点，就略过，否则把相应信息存入vector
  if (num(pos) != 0) {
    FP.push_back(pos);
    FV.push_back(val(pos));
    FN.push_back(num(pos));
  }
  // 递归地拉平右子树
  if (exist(rs(pos))) flatten(rs(pos));
  return id;
}
// 以pos为根节点，以[l,r]内的信息重建一棵平衡的树
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
  num(pos) = FN[mid];  // 把存于vector中的信息复制过来
  val(pos) = FV[mid];
  sz(pos) = sz1 + sz2 + num(pos);  // 递归确定重建后树的大小
}
// 尝试重构当前子树
void try_restructure(int pos) {
  double k = max(sz(ls(pos)), sz(rs(pos))) / double(sz(pos));
  if (k > ALPHA) {
    FP.clear(), FV.clear(), FN.clear();  // 清空vector
    int id = flatten(pos);
    // 这里是确保当前节点的编号在重构后不会改变
    swap(FP[id], FP[(FP.size() - 1) / 2]);
    rebuild(pos);
  }
}
// 接下来是普通的二叉查找树
void insert(int v, int pos = 1) {
  if (!exist(pos)) {
    val(pos) = v;
    num(pos) = 1;
  } else if (v < val(pos)) {
    if (!exist(ls(pos))) ls(pos) = ++cnt;
    insert(v, ls(pos));
  } else if (v > val(pos)) {
    if (!exist(rs(pos))) rs(pos) = ++cnt;
    insert(v, rs(pos));
  } else
    num(pos)++;
  sz(pos)++;
  try_restructure(pos);
}
void remove(int v, int pos = 1) {
  sz(pos)--;
  if (v < val(pos))
    remove(v, ls(pos));
  else if (v > val(pos))
    remove(v, rs(pos));
  else
    num(pos)--;
  try_restructure(pos);
}
int countl(int v, int pos = 1) {
  if (v < val(pos))
    return exist(ls(pos)) ? countl(v, ls(pos)) : 0;
  else if (v > val(pos))
    return sz(ls(pos)) + num(pos) + (exist(rs(pos)) ? countl(v, rs(pos)) : 0);
  else
    return sz(ls(pos));
}
int countg(int v, int pos = 1) {
  if (v > val(pos))
    return exist(rs(pos)) ? countg(v, rs(pos)) : 0;
  else if (v < val(pos))
    return sz(rs(pos)) + num(pos) + (exist(ls(pos)) ? countg(v, ls(pos)) : 0);
  else
    return sz(rs(pos));
}
int rank(int v) { return countl(v) + 1; }
int kth(int k, int pos = 1) {
  if (sz(ls(pos)) + 1 > k)
    return kth(k, ls(pos));
  else if (sz(ls(pos)) + num(pos) < k)
    return kth(k - sz(ls(pos)) - num(pos), rs(pos));
  else
    return val(pos);
}
int pre(int v) {
  int r = countl(v);
  return kth(r);
}
int suc(int v) {
  int r = sz(1) - countg(v) + 1;
  return kth(r);
}
