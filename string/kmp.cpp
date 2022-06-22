int nxt[100005], ns, nt;
char t[100005], s[100005];
void get_next() {
  nxt[0] = -1;
  int k = -1, j = 0;
  while (t[j] != '\0') {
    if (k == -1 || t[k] == t[j]) {
      nxt[++j] = ++k;
    } else {
      k = nxt[k];
    }
  }
}
int search() {
  int id = 0;
  for (int i = 0; i < ns; i++) {
    while (id != -1 && s[i] != t[id]) id = nxt[id];
    id++;
    if (id == nt) return i - nt + 1;
  }
}
