#include <cstdio>
#include <stack>
using namespace std;
stack<int> stk;
// adjacent matrix
int mp[10][10];
// reversed graph
int mpt[10][10];
int vst[10];
int clr[10];
int vn, en;
void dfs1(int s) {
  if (vst[s] == 1) return;
  vst[s] = 1;
  // dfs routine
  for (int i = 1; i <= vn; ++i) {
    if (mp[s][i] < 0x3f3f3f3f) {
      dfs1(i);
    }
  }
  // push
  stk.push(s);
}
void dfs2(int s, int cnt) {
  if (vst[s] == 0) return;
  clr[s] = cnt;
  vst[s] = 0;
  for (int i = 1; i <= vn; ++i) {
    if (mpt[s][i] < 0x3f3f3f3f) {
      dfs2(i, cnt);
    }
  }
}
void init() {
  for (int i = 1; i <= vn; ++i) {
    for (int j = 1; j <= vn; ++j) {
      mp[i][j] = mp[j][i] = 0x3f3f3f3f;
      mpt[i][j] = mpt[j][i] = 0x3f3f3f3f;
    }
    mpt[i][i] = mp[i][i] = 0;
  }
}
void SCC_kor() {
  // init();
  for (int i = 1; i <= vn; ++i) {
    if (vst[i] == 0) dfs1(i);
  }
  int cnt = 1;
  while (!stk.empty()) {
    int s = stk.top();
    stk.pop();
    if (vst[s] == 0) continue;
    dfs2(s, cnt++);
  }
  // vertexes with same value in clr[] is in one SCC
  for (int i = 1; i <= vn; ++i) {
    printf("%d ", clr[i]);
  }
  printf("\n");
}
int main() {
  scanf("%d %d", &vn, &en);
  init();
  for (int i = 1; i <= en; ++i) {
    int fr, to;
    scanf("%d %d", &fr, &to);
    mp[fr][to] = 1;
    mpt[to][fr] = 1;
  }
  SCC_kor();
  return 0;
}
