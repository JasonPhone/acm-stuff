#include <bits/stdc++.h>
using namespace std;
int n, m;
struct node {
  vector<int> nxt;
} g[100000];
int dfn[100000];
int low[100000];
int d[100000];
int col[100000];
int colour;
int cnt[100000];
int stk[100000];
int vis[100000];
int top;
int deep;
void tarjan(int u) {
  dfn[u] = low[u] = ++deep;
  stk[top++] = u;
  vis[u] = 1;
  for (int i = 0; i < g[u].nxt.size(); i++) {
    int v = g[u].nxt[i];
    if (!vis[v]) {
      tarjan(v);
      low[u] = min(low[v], low[u]);
    } else {
      low[u] = min(low[v], low[u]);
    }
  }
  if (dfn[u] == low[u]) {
    int node;
    colour++;
    while (node != u) {
      node = stk[top - 1];
      top--;
      col[node] = colour;
    }
  }
}
