#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <cassert>
#include <iostream>
#include <stdexcept>

#include "filekl.h"
#include "basic_io.h"

namespace { bool verbose=true; }

namespace atlas {
  namespace filekl {

void scan_polynomials
    (const atlas::filekl::block_info& bi
    ,const atlas::filekl::polynomial_info& pi
    ,const atlas::filekl::progress_info& ri
    ,std::ostream& row_out)
{
  using blocks::BlockElt;
  ullong last_l_offset=0;
  ullong prev_offset=0;
  for (size_t l=0; l<=bi.max_length; ++l)
  {
    for (blocks::BlockElt y=bi.start_length[l]; y<bi.start_length[l+1]; ++y)
    {
      ullong offset=pi.coeff_start(ri.first_new_in_row(y+1));
      basic_io::put_int(offset-prev_offset,row_out); prev_offset=offset;
    }
    std::cout<< "Completed l=" << l << " @y=" << bi.start_length[l+1]
	     << ", new coefficients: " << prev_offset-last_l_offset
	     << ", total: " << prev_offset << ".\n";
    last_l_offset=prev_offset;
  }
}


  } // namespace filekl
} // namespace atlas

int main(int argc, char** argv)
{
  --argc; ++argv; // read and skip program name

  if (argc>0 and std::string(*argv)=="-q")
    { verbose=false; --argc; ++argv;}

  if (argc!=4)
  {
    std::cerr <<
      "Usage: polstat [-q] block-file poly-file row-file row-coeff-file\n";
    exit(1);
  }

  std::ifstream block_file(argv[0]);
  std::ifstream pol_file(argv[1]);
  std::ifstream row_file(argv[2]);
  if (not block_file.is_open()
      or not row_file.is_open()
      or not row_file.is_open())
  {
    std::cerr << "Failure opening file(s).\n";
    exit(1);
  }

  atlas::filekl::block_info bi(block_file);
  atlas::filekl::polynomial_info pi(pol_file);
  atlas::filekl::progress_info ri(row_file);

  std::ofstream row_out(argv[3],
		 std::ios_base::out
		 | std::ios_base::trunc
		 | std::ios_base::binary);
  if (not row_out.is_open())
  {
    std::cerr << "Failure opening "<< argv[3] <<" for writing.\n";
    exit(1);
  }

  atlas::filekl::scan_polynomials(bi,pi,ri,row_out);


}
