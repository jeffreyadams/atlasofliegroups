#include <vector>
#include <iostream>
#include <sstream>
#include <cctype>
#include <cstring>
#include <cstdlib> // for |exit|
#include "mod2_system.h"
#include "bitmap.h"

using atlas::bitmap::BitMap;
using atlas::mod2_system::Mod2_System;

static const int n = 5;
static const int n2 = n*n;

std::vector<unsigned long> neighbours(unsigned int a, unsigned int b)
{
  std::vector<unsigned long> result; result.reserve(5);
  if (a>0)   result.push_back(n*(a-1)+b);
  if (b>0)   result.push_back(n*a+b-1);
  result.push_back(n*a+b);
  if (b<n-1) result.push_back(n*a+b+1);
  if (a<n-1) result.push_back(n*(a+1)+b);
  return result;
}

bool off (const BitMap& lights, bool full)
{
  Mod2_System sys;
  sys.extend(n2);
  for (unsigned int i=0; i<n2; ++i)
  {
    std::vector<unsigned long> lhs = neighbours(i/n,i%n);
    if (not sys.add(lhs.begin(),lhs.end(), lights.isMember(i) ? 1 : 0))
      return false;
  }
  BitMap sol = sys.a_solution();
  std::cout << "Flip switches at";
  for (BitMap::iterator it=sol.begin(); it(); ++it)
    std::cout << " (" << *it/n << ',' << *it%n << ')';
  std::cout << std::endl;

  if (full)
  {
    std::vector<BitMap> kernel = sys.solution_basis();
    if (kernel.empty())
      std::cout << "Solution is unique" << std::endl;
    else
    {
      std::cout << "Neutral sets of switches:" << std::endl;
      for (unsigned int i=0; i<kernel.size(); ++i)
      {
	for (BitMap::iterator it=kernel[i].begin(); it(); ++it)
	  std::cout << " (" << *it/n << ',' << *it%n << ')';
	std::cout << std::endl;
      }
    }
  }
  return true;
}

int main(int argc, char** argv)
{
  --argc, ++argv; // skip program name

  bool full = argv!=NULL and strcmp(*argv,"-f")==0;
  if (full)
    --argc, ++argv;

  BitMap lights(n2);
  while (argc-->0)
  {
    std::istringstream arg (*argv++);
    unsigned int a,b; char c;
    if (not isdigit(arg.peek()))
      arg >> c;
    arg >> a >> c >> b;
    if (arg.fail() or a>=n or b>=n)
    {
      std::cerr << "Improper argument " << arg.str() << std::endl;
      exit(1);
    }
    lights.insert(n*a+b);
  }
  if (not off(lights,full))
    std::cout << "No solutions.\n";
}
