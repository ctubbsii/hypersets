#include <iostream>
using namespace std;

#define LEVEL 7
#define U_INT64 unsigned long long

int main(void)
{
  U_INT64 cmin, cmax, mask, mid, i;
  mask = (1ull<<LEVEL)-1;
  cmin = 1ull<<(LEVEL*LEVEL-2);
  cmax = (1ull<<(LEVEL*LEVEL))-1;
  
  mid = ((cmax-cmin)/32ull) & (~mask);

  i = 0;
  cout << "RANGE: " << cmin << " - " << cmax << endl;
  do
  {
    cout << i << ": " << cmin << " - " << cmin+mid << endl;
    cmin += mid;
    ++i;
  }
  while (cmin+mid < cmax);

  cout << i << ": " << cmin << " - " << cmax << endl;


  return 0;
}

