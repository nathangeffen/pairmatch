#include <iostream>
#include <random>
#include <vector>
#include "ransampl.h"

using namespace std;

mt19937 rng;

int main()
{
  double p[] = {0.1, 0.1, 0.2, 0.4, 0.1, 0.05, 0.05};
  const int n = sizeof(p) / sizeof(double);
  ransampl_ws* ws = ransampl_alloc(n);
  ransampl_set(ws, p);

  vector<int> cumul(n, 0);
  std::uniform_real_distribution<double> uni;
  for (int j = 0; j < 100000000; ++j) {
    int i = ransampl_draw(ws, uni(rng), uni(rng));
    cumul[i] += 1;
  }

  int total = 0;
  for (int i = 0; i < cumul.size(); ++i) {
    cout << cumul[i] << " ";
    total += cumul[i];
  }
  cout << endl;

  for (int i = 0; i < cumul.size(); ++i) {
    cout << (double) cumul[i] / total  << " ";
  }
  cout << endl;


  // Free workspace and terminate:
  ransampl_free( ws );
  return 0;
}
