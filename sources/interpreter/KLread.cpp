#include <string>
#include <fstream>
#include <iostream>

size_t read_bytes(unsigned int n, std::istream& in)
{
  if (n==0) return 0;
  else
    {
      char c; in.get(c); unsigned char low=c;
      return low+(read_bytes(n-1,in)<<8);
    }
}

int main(int argc,char** argv)
{
  std::ifstream in_file;
  if (argc==2) in_file.open(argv[1]);
  else
    {
      std::string file_name;
      std::cout << "File name: " ;
      std::cin >> file_name;
      in_file.open(file_name.c_str());
    }
  if (not in_file.is_open())
    { std::cerr << "Open failed"; exit(1); }

  size_t n_polynomials=read_bytes(4,in_file);

  std::streamoff index_begin=in_file.tellg();
  in_file.seekg(5*n_polynomials,std::ios_base::cur);

  size_t n_coefficients =read_bytes(5,in_file);
  std::streamoff coefficients_begin=in_file.tellg();

  std::cout << n_polynomials << " polynomials, "
	    << n_coefficients << " coefficients.\n";

  while (true)
    {
      std::cout << "index: ";
      size_t i=~0ul; std::cin >> i;
      if (i==size_t(~0ul)) break;
      if (i>=n_polynomials)
	{
	  std::cout << "index too large, limit is " << n_polynomials-1
		    << ".\n";
	  continue;
	}
      in_file.seekg(index_begin+5*i,std::ios_base::beg);
      size_t index=read_bytes(5,in_file);
      size_t next_index=read_bytes(5,in_file);
      if (next_index-index>=33)
	{
	  std::cout << "Degree found too large: " << next_index-index-1
		    << ".\n";
	  continue;
	}
      char* coefficients = new char[next_index-index];

      in_file.seekg(coefficients_begin+index,std::ios_base::beg);
      in_file.read(coefficients, next_index-index);
      if (not in_file.good())
	{
	  std::cout << "Input error reading coefficients.\n";
	  continue;
	}
      for (size_t i=0; i<next_index-index; ++i)
	std::cout << (int)coefficients[i] << "X^" << i
		  << (i+1<next_index-index ? " + " : ".\n");

      delete[] coefficients;
    }

}
