int nxt[100005];
char t[100005];
void getNxt() {
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


