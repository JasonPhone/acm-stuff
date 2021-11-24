using ll = long long;
ll MAXN = 2e5 + 5;
ll n, q[MAXN], tmp[MAXN];
// [l, r]
ll merge_sort(int l, int r) {
  if (l >= r) return 0;
  ll mid = (l + r) >> 1;
  ll res = merge_sort(l, mid) + merge_sort(mid + 1, r);

  ll k = 0, i = l, j = mid + 1;
  while (i <= mid && j <= r) {
    if (q[i] <= q[j])
      tmp[k++] = q[i++];
    else {
      tmp[k++] = q[j++];
      res += mid - i + 1;
    }
  }
  while (i <= mid) tmp[k++] = q[i++];
  while (j <= r) tmp[k++] = q[j++];
  for (ll i = l, j = 0; i <= r; i++, j++) q[i] = tmp[j];
  return res;
}
