#include <bits/stdc++.h>
using namespace std;
int main() {
  int a[20] = {1, 2, 3, 3, 3, 6, 9, 100000000, 10000002};
  sort(a, a + 9);
  int len = unique(a, a + 9) - a;
  cout << lower_bound(a, a + len, 9) - a << endl;;

  vector<int> v({1, 2, 3, 3, 3, 6, 9, 100000000, 10000002});
  sort(v.begin(), v.end());
  // split operation for potential different evaluate order of parameters
  vector<int>::iterator ed = unique(v.begin(), v.end());
  v.erase(ed, v.end());
  cout << lower_bound(v.begin(), v.end(), 9) - v.begin() << endl;

  return 0;
}
