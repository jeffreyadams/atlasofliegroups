#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <cassert>
#include <iostream>
#include <stdexcept>

#include "filekl.h"
#include "basic_io.h"
#include "tally.h"

namespace { bool verbose=true; }

namespace atlas {
  namespace filekl {

typedef tally::TallyVec<unsigned char> tally_vec;

void scan_polynomials
    (const atlas::filekl::block_info& bi
    ,const atlas::filekl::polynomial_info& pi
    ,const atlas::filekl::progress_info& ri
    ,const std::string file_name_base)
{
  using blocks::BlockElt;

  const size_t deg_limit=bi.max_length/2;

  tally_vec deg_val(deg_limit*(deg_limit+1)/2);
                              // degree-valuation joint distribution
  tally_vec q_is_1(70000000); // distribution of specialisation at 1
  tally_vec q_is_0(9);        // distribution of specialisation at 0
  tally_vec q_is_2(1ul<<31);  // distribution of specialisation at 2
  tally_vec q_min_1_pos(1000000); // non-negative values at q=-1
  tally_vec q_min_1_neg(1000000); // negative values at q=-1
  tally_vec leading(65536);   // distribution of leading coefficients
  tally_vec coeff(70000000);  // distribution of all coefficients

  unsigned int deg_max=0,val_max=0; // maxmial degree, valuation found;
  filekl::KLIndex maxi_deg,maxi_val,maxi_at_1,maxi_at_0,maxi_at_2,
    maxi_at_min1,mini_at_min1,maxi_lead; // indices where maxima first attained

  for (size_t l=0; l<=bi.max_length; ++l)
  {
    for (blocks::BlockElt y=bi.start_length[l]; y<bi.start_length[l+1]; ++y)
    {
      if (verbose)
	std::cerr << y << ", #" << ri.first_new_in_row(y) << '\r'
		  << std::flush;
      for (filekl::KLIndex i=ri.first_new_in_row(y);
	   i<ri.first_new_in_row(y+1); ++i)
      {
	std::vector<size_t> coeffs=pi.coefficients(i);
	if (coeffs.empty()) continue; // skip zero polynomial
	size_t degree=coeffs.size()-1;
	size_t val=0; while (coeffs[val]==0) ++val;
	size_t lead=coeffs[degree];

	if (degree>deg_max) { deg_max=degree; maxi_deg=i; }
	if (val>val_max) { val_max=val; maxi_val=i; }
	deg_val.tally(degree*(degree+1)/2+val);

	tally_vec::Index at_1=lead,at_2=lead,at_0=coeffs[0];
	long int at_min1=lead;
	coeff.tally(lead);
        for (size_t j=degree; j-->0;)
	{
	  size_t c=coeffs[j];
	  coeff.tally(c);
	  at_1+=c; at_2=(at_2<<1)+c; at_min1=c-at_min1;
	}
	if (at_0>=q_is_0.size()) maxi_at_0=i;
	if (at_1>=q_is_1.size()) maxi_at_1=i;
	if (at_2>=q_is_2.size()) maxi_at_2=i;
	if (at_min1>=0)
	{
	  if (unsigned(at_min1)>=q_min_1_pos.size()) maxi_at_min1=i;
	}
	else if (unsigned(~at_min1)>=q_min_1_neg.size()) mini_at_min1=i;
	if (lead>=leading.size()) maxi_lead=i;

	q_is_0.tally(at_0);
	q_is_1.tally(at_1);
	q_is_2.tally(at_2);
	if (at_min1>=0) q_min_1_pos.tally(at_min1);
	else q_min_1_neg.tally(~at_min1);
	leading.tally(lead);

      } // for (i)
    } // for(y)
    if (verbose)
      std::cout<< "Completed l=" << l << " @y=" << bi.start_length[l+1]
	       << ", maximal deg, val, q=1: " << deg_max << ", "
	       << val_max <<  ", " << q_is_1.size()-1 << '.' << std::endl;
  } // for (l)
  std::cout <<
    "Indices for maximal deg, val, q=0, 1, 2, -1, leading, and minimal q=-1:\n#"
    << maxi_deg <<  ", #"
    << maxi_val <<  ", #"
    << maxi_at_0 <<  ", #"
    << maxi_at_1 <<  ", #"
    << maxi_at_2 <<  ", #"
    << maxi_at_min1 <<  ", #"
    << maxi_lead <<  ", #"
    << mini_at_min1 << ".\n";

  std::ofstream stat_out;
  stat_out.open((file_name_base+"-deg_val").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening deg_val file failed");
  for (size_t d=0; d<deg_limit; ++d)
  {
    stat_out << "Degree " << d << ": ";
    for (size_t v=0; v<=d; ++v)
      stat_out << deg_val.multiplicity(d*(d+1)/2+v) << (v<d ? ',' : '\n');
  }
  stat_out.close();

  stat_out.open((file_name_base+"@q=0").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=0 file failed");
  for (tally_vec::Index i=0; i<q_is_0.size(); q_is_0.advance(i))
    stat_out << i << ": " << q_is_0.multiplicity(i) << ".\n";
  stat_out.close();

  stat_out.open((file_name_base+"@q=1").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=1 file failed");
  q_is_1.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"@q=2").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=2 file failed");
  q_is_2.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"@q=-1").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=-1 file failed");

  for (tally_vec::Index i=q_min_1_neg.size(); q_min_1_neg.lower(i);)
    stat_out << ~(long long)i << ": "
	     << q_min_1_neg.multiplicity(i) << ".\n";
  for(tally_vec::Index i=0; i<q_min_1_pos.size(); q_min_1_pos.advance(i))
    stat_out << i << ": " << q_min_1_pos.multiplicity(i) << ".\n";

  stat_out.close();

  stat_out.open((file_name_base+"-leading").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening leading term file failed");
  leading.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-coefficients").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening coefficients file failed");
  coeff.write_to(stat_out);
  stat_out.close();

}


  } // namespace filekl
} // namespace atlas


std::string getFileName(std::string prompt)
{
  std::cerr << prompt;
  std::string name;
  std::cin >> name;
  return name;
}

int main(int argc, char** argv)
{
  --argc; ++argv; // read and skip program name

  if (argc>0 and std::string(*argv)=="-q")
    { verbose=false; --argc; ++argv;}

  if (argc!=3)
  {
    std::cerr <<
      "Usage: polstat [-q] block-file poly-file row-file\n";
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

  std::string file_name_base= getFileName
    ("Base file name for statistics output: ");

  atlas::filekl::scan_polynomials(bi,pi,ri,file_name_base);


}
