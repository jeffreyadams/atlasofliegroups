#include "filekl_in.h"

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
      for (size_t j=c.size(); j-->0;)
	if (c[j]!=0)
	{
	  if (first) first=false; else std::cout << " + ";
	  if (c[j]!=1 or j==0) std::cout << c[j];
	  std::cout << (j==0? "" : "q");
	  if (j>1) std::cout << '^' << j;
	}
      std::cout << '.' << std::endl;
    }
}
