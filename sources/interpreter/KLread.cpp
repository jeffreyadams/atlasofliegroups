#include <string>
#include <fstream>
#include <iostream>

const std::ios_base::openmode binary_out=
			    std::ios_base::out
			  | std::ios_base::trunc
			  | std::ios_base::binary;

const std::ios_base::openmode binary_in=
			    std::ios_base::in
			  | std::ios_base::binary;

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
  if (argc==2) in_file.open(argv[1],binary_in);
  else
    {
      std::string file_name;
      std::cout << "File name: " ;
      std::cin >> file_name;
      in_file.open(file_name.c_str(),binary_in);
    }
  if (not in_file.is_open())
    { std::cerr << "Open failed"; exit(1); }

  size_t n_polynomials=read_bytes(4,in_file);

  std::streamoff index_begin=in_file.tellg();
  read_bytes(10,in_file); // skip initial 2 indices
  unsigned int coefficient_size=read_bytes(5,in_file); // size of One
  std::cout << "Coefficient size " << coefficient_size << ".\n";

  in_file.seekg(index_begin+5*n_polynomials,std::ios_base::beg);
  size_t n_coefficients =read_bytes(5,in_file)/coefficient_size;
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
      size_t length=(next_index-index)/coefficient_size;
      if (length>=33)
	{
	  std::cout << "Degree found too large: " << length-1
		    << ".\n";
	  continue;
	}
      unsigned long* coefficients = new unsigned long[length];

      in_file.seekg(coefficients_begin+index,std::ios_base::beg);
      for (unsigned int i=0; i<length; ++i)
	coefficients[i]=read_bytes(coefficient_size,in_file);

      if (not in_file.good())
	{
	  std::cout << "Input error reading coefficients.\n";
	  continue;
	}
      for (size_t i=0; i<length; ++i)
	std::cout << coefficients[i] << "X^" << i
		  << (i+1<length ? " + " : ".\n");

      delete[] coefficients;
    }

}
