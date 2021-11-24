// Manarcher
// find the palindrome in O(n)
#include <bits/stdc++.h>
using namespace std;
char s[100005];
int ps = 0;
int p[100005], ctr, maxr, mirr;
void solve() {
  ctr = maxr = 0;
  for (int i = 0; i < ps; ++i) {
    mirr = 2 * ctr - i;
    if (i < maxr) {
      p[i] = min(maxr - i, p[mirr]);
    } else {
      p[i] = 0;
    }
    while (s[i - 1 - p[i]] == s[i + 1 + p[i]]) {
      p[i]++;
    }
    if (p[i] + i > maxr) {
      ctr = i;
      maxr = p[i] + i;
    }
  }
  int maxi = 0;
  for (int i = 0; i < ps; ++i) {
    maxi = p[maxi] < p[i] ? i : maxi;
  }
  printf("%d\n", p[maxi]);
  for (int i = maxi - p[maxi]; i <= maxi + p[maxi]; ++i) {
    if (s[i] != '#') {
      printf("%c", s[i]);
    }
  }
  printf("\n");
}
int main() {
  int Case = 1;
  while (Case--) {
    char c = getchar();
    s[ps++] = '#';
    while (c != '\n') {
      s[ps++] = c;
      s[ps++] = '#';
      c = getchar();
    }
    solve();
  }
  return 0;
}
