#include <vector>
#include <iostream>
#include <sstream>
#include <cctype>
#include <cstring>
#include <cstdlib> // for |exit|
#include "sl_list.h"

typedef atlas::containers::sl_list<int> intlist;

int main()
{
  intlist a;
  a.push_back(3);
  a.push_back(1);
  a.push_back(4);
  a.push_back(13);
  a.reverse();

  std::ostream_iterator<int> lister(std::cout,", ");
  std::copy(a.begin(),a.end(),lister);
  std::cout<<std::endl;

  a.reverse(a.begin(),++++++a.begin());
  std::copy(a.begin(),a.end(),lister);
  std::cout<<std::endl;
}
