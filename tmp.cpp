/**
 * 单调栈求凸包
 * stk[] 是整型，存的是下标
 * vector<int, int> p[] 存储向量或点
 */
tp = 0;                       // 初始化栈
std::sort(p + 1, p + 1 + n);  // 对点按坐标排序
stk[++tp] = 1;
//栈内添加第一个元素，且不更新 used，使得 1 在最后封闭凸包时也对单调栈更新
for (int i = 2; i <= n; ++i) {
  // 下一行 * 操作符被重载为叉积
  while (tp >= 2 && (p[stk[tp]] - p[stk[tp - 1]]) * (p[i] - p[stk[tp]]) <= 0)
    used[stk[tp--]] = 0;
  used[i] = 1;  // used 表示在凸壳上
  stk[++tp] = i;
}
int tmp = tp;  // tmp 表示下凸壳大小
for (int i = n - 1; i > 0; --i)
  if (!used[i]) {
    // 求上凸壳时不影响下凸壳
    while (tp > tmp && (p[stk[tp]] - p[stk[tp - 1]]) * (p[i] - p[stk[tp]]) <= 0)
      used[stk[tp--]] = 0;
    used[i] = 1;
    stk[++tp] = i;
  }
for (int i = 1; i <= tp; ++i)  // 复制到新数组中去
  h[i] = p[stk[i]];
int ans = tp - 1;

//////////////////////////////////////////////////
/**
 * query the number of elements equal to val
 * which is previously neighbouring pos (inclusive)
 * call (1, 5, 0) for [1, 1, 1, 0, 0, 0, 0, 0] (start from 1)
 * will get 2
 */
int query_prefix_num(int k, int pos, int val) {
  // val == -1 if seg under this node all not all same
  if (tree[k].val == val) return min(pos, tree[k].r) - tree[k].l + 1;
  if (tree[k].l == tree[k].r) return tree[k].val == val;
  push_down(k);
  int mid = tree[k].l + (tree[k].r - tree[k].l) / 2;
  if (pos > mid) {
    int len = query_prefix_num(k << 1 | 1, pos, val);
    if (len == min(pos, tree[k << 1 | 1].r) - tree[k << 1 | 1].l + 1)
      len += query_prefix_num(k << 1, pos, val);
    return len;
  }
  return query_prefix_num(k << 1, pos, val);
}
