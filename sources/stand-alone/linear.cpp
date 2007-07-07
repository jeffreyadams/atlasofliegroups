#include "filekl.h"

int main(int argc, char** argv)
{
  --argc; ++argv; // skip program name

  if (argc!=1)
  {
    std::cerr <<
      "Usage: linear pol-file\n";
    exit(1);
  }

  std::ifstream pol_file(argv[0]);
  if (not pol_file.is_open())
  {
    std::cerr << "Failure opening file.\n";
    exit(1);
  }

  atlas::filekl::polynomial_info pi(pol_file);

  for (size_t i=0; i<pi.n_polynomials(); ++i)
    if (pi.degree(i)==1)
    {
      std::cout << '#' << i << ": ";
      const std::vector<size_t>& c =pi.coefficients(i);
      bool first=true;
      for (size_t i=c.size(); i-->0;)
	if (c[i]!=0)
	{
	  if (first) first=false; else std::cout << " + ";
	  if (c[i]!=1 or i==0) std::cout << c[i];
	  std::cout << (i==0? "" : "q");
	  if (i>1) std::cout << '^' << i;
	}
      std::cout << '.' << std::endl;
    }
}
