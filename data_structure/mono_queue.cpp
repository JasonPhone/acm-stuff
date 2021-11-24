#include <bits/stdc++.h>
// monotonic descending queue, segMax at front
using namespace std;

void getSegMax(vector<int>& v, int k, vector<int>& ans) {
  deque<int> que;
  int n = v.size();
  for (int i = 0; i + 1 < k; ++i) {
    while (!que.empty() && v[que.back()] <= v[i]) que.pop_back();
    que.push_back(i);
  }
  for (int i = k - 1; i < n; ++i) {
    while (!que.empty() && v[que.back()] <= v[i]) que.pop_back();
    que.push_back(i);
    while (que.front() <= i - k) que.pop_front();
    ans.push_back(v[que.front()]);
  }
}
void getSegMin(vector<int>& v, int k, vector<int>& ans) {
  deque<int> que;
  int n = v.size();
  for (int i = 0; i + 1 < k; ++i) {
    while (!que.empty() && v[que.back()] >= v[i]) que.pop_back();
    que.push_back(i);
  }
  for (int i = k - 1; i < n; ++i) {
    while (!que.empty() && v[que.back()] >= v[i]) que.pop_back();
    que.push_back(i);
    while (que.front() <= i - k) que.pop_front();
    ans.push_back(v[que.front()]);
  }
}
int main() {
  vector<int> v = {2, 3, 1, 4, 5, 6, 7, 3};
  vector<int> ans;
  getSegMin(v, 3, ans);
  for (auto itm: ans) {
    cout << itm << " ";
  }
  return 0;
}
