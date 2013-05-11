#include <cassert>
#include <cstdlib> // for |exit|
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <iostream>
#include <stdexcept>

#include "filekl_in.h"
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
    ,const tally_vec& prim_mu
    ,const tally_vec& total_mu
    ,const std::string file_name_base)
{
  const size_t deg_limit=(bi.max_length+1)/2; // this is |(max_length-1)/2+1|

  tally_vec deg_val(deg_limit*(deg_limit+1)/2);
                              // degree-valuation joint distribution
  tally_vec q_is_1(65000000); // distribution of specialisation at 1
  tally_vec q_is_0(9);        // distribution of specialisation at 0
  //  tally_vec q_is_2(1ul<<31);  // distribution of specialisation at 2
  //  tally_vec q_min_1_pos(25000); // non-negative values at q=-1
  //  tally_vec q_min_1_neg(25000); // negative values at q=-1
  tally_vec leading(65536);   // distribution of leading coefficients
  tally_vec coeff(11808808);  // distribution of all coefficients

  tally_vec deg_val_t(deg_limit*(deg_limit+1)/2);
                              // degree-valuation joint distribution
  tally_vec q_is_1_t(65000000); // distribution of specialisation at 1
  tally_vec q_is_0_t(9);        // distribution of specialisation at 0
  tally_vec leading_t(65536);   // distribution of leading coefficients
  tally_vec coeff_t(11808808);  // distribution of all coefficients

  unsigned int deg_max=0,val_max=0; // maxmial degree, valuation found;

  for (size_t l=0; l<=bi.max_length; ++l)
  {
    for (BlockElt y=bi.start_length[l]; y<bi.start_length[l+1]; ++y)
    {
      if (verbose)
	std::cerr << y << ", #" << ri.first_new_in_row(y) << '\r'
		  << std::flush;
      for (filekl::KLIndex i=ri.first_new_in_row(y);
	   i<ri.first_new_in_row(y+1); ++i)
      {
        unsigned long long int mu=prim_mu.multiplicity(i);
        unsigned long long int tot_mu=total_mu.multiplicity(i);
	if (mu>tot_mu)
	{
	  std::cerr << i << ": " << mu << '>' << tot_mu << ". \n";
	  throw std::runtime_error("inconsistent multiplicities");
	}

	std::vector<size_t> coeffs=pi.coefficients(i);
	if (coeffs.empty()) continue; // skip zero polynomial
	size_t degree=coeffs.size()-1;
	size_t val=0; while (coeffs[val]==0) ++val;
	size_t lead=coeffs[degree];

	if (degree>deg_max) deg_max=degree;
	if (val>val_max) val_max=val;
	deg_val.tally(degree*(degree+1)/2+val,mu);
	deg_val_t.tally(degree*(degree+1)/2+val,tot_mu);
	if (deg_val.multiplicity(degree*(degree+1)/2)
	    > deg_val_t.multiplicity(degree*(degree+1)/2))
	{
	  std::cerr << "Problem at " << i << ", ("
		    << degree << ',' << val << "), (" << mu << "<=" << tot_mu
		    << "): "
		    << deg_val.multiplicity(degree*(degree+1)/2) << '>'
		    << deg_val_t.multiplicity(degree*(degree+1)/2)
		    << std::endl;
	  throw std::runtime_error("Tally inversion");
	}

	tally_vec::Index at_1=lead,at_0=coeffs[0];
	coeff.tally(lead,mu);
	coeff_t.tally(lead,tot_mu);
        for (size_t j=degree; j-->0;)
	{
	  size_t c=coeffs[j];
	  coeff.tally(c,mu);
	  coeff_t.tally(c,tot_mu);
	  at_1+=c;
	}

	q_is_0.tally(at_0,mu);
	q_is_0_t.tally(at_0,tot_mu);
	q_is_1.tally(at_1,mu);
	q_is_1_t.tally(at_1,tot_mu);
	leading.tally(lead,mu);
	leading_t.tally(lead,tot_mu);

      } // for (i)
    } // for(y)
    if (verbose)
      std::cout<< "Completed l=" << l << " @y=" << bi.start_length[l+1]
	       << ", maximal deg, val, q=1: " << deg_max << ", "
	       << val_max <<  ", " << q_is_1.size()-1 << '.' << std::endl;
  } // for (l)

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

  stat_out.open((file_name_base+"-total-deg_val").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening deg_val file failed");
  for (size_t d=0; d<deg_limit; ++d)
  {
    stat_out << "Degree " << d << ": ";
    for (size_t v=0; v<=d; ++v)
      stat_out << deg_val_t.multiplicity(d*(d+1)/2+v) << (v<d ? ',' : '\n');
  }
  stat_out.close();

  stat_out.open((file_name_base+"@q=0").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=0 file failed");
  for (tally_vec::Index i=0; i<q_is_0.size(); q_is_0.advance(i))
    stat_out << i << ": " << q_is_0.multiplicity(i) << ".\n";
  stat_out.close();

  stat_out.open((file_name_base+"-total@q=0").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=0 file failed");
  for (tally_vec::Index i=0; i<q_is_0_t.size(); q_is_0_t.advance(i))
    stat_out << i << ": " << q_is_0_t.multiplicity(i) << ".\n";
  stat_out.close();

  stat_out.open((file_name_base+"@q=1").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=1 file failed");
  q_is_1.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-total@q=1").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=1 file failed");
  q_is_1_t.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-leading").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening leading term file failed");
  leading.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-total-leading").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening leading term file failed");
  leading_t.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-coefficients").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening coefficients file failed");
  coeff.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-total-coefficients").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening coefficients file failed");
  coeff_t.write_to(stat_out);
  stat_out.close();

}

size_t gcd(size_t a,size_t b)
{
  if (a==1 or b==1) return 1;
  while(a!=0)
  { b%=a;
    if (b==0) return a;
    a%=b;
  }
  return b;

}

void scan_polynomials
    (const atlas::filekl::block_info& bi
    ,const atlas::filekl::polynomial_info& pi
    ,const atlas::filekl::progress_info& ri
    ,const std::string file_name_base)
{
  const size_t deg_limit=(bi.max_length+1)/2; // this is |(max_length-1)/2+1|

  tally_vec deg_val(deg_limit*(deg_limit+1)/2);
                              // degree-valuation joint distribution
  tally_vec q_is_1(65000000); // distribution of specialisation at 1
  tally_vec q_is_0(9);        // distribution of specialisation at 0
  //  tally_vec q_is_2(1ul<<31);  // distribution of specialisation at 2
  //  tally_vec q_min_1_pos(25000); // non-negative values at q=-1
  //  tally_vec q_min_1_neg(25000); // negative values at q=-1
  tally_vec leading(65536);   // distribution of leading coefficients
  tally_vec coeff(11808808);  // distribution of all coefficients
  tally_vec cont(256);        // distribution of contents (=gcd(coefficients))

  unsigned int deg_max=0,val_max=0; // maxmial degree, valuation found;

  for (size_t l=0; l<=bi.max_length; ++l)
  {
    for (BlockElt y=bi.start_length[l]; y<bi.start_length[l+1]; ++y)
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

	if (degree>deg_max) deg_max=degree;
	if (val>val_max) val_max=val;
	deg_val.tally(degree*(degree+1)/2+val);

	tally_vec::Index at_1=lead, contents=lead;
	tally_vec::Index at_0=coeffs[0];
	coeff.tally(lead);
        for (size_t j=degree; j-->0;)
	{
	  size_t c=coeffs[j];
	  coeff.tally(c);
	  at_1+=c;
	  contents=gcd(contents,c);
	}

	q_is_0.tally(at_0);
	q_is_1.tally(at_1);
	leading.tally(lead);
	cont.tally(contents);

      } // for (i)
    } // for(y)
    if (verbose)
      std::cout<< "Completed l=" << l << " @y=" << bi.start_length[l+1]
	       << ", maximal deg, val, q=1: " << deg_max << ", "
	       << val_max <<  ", " << q_is_1.size()-1 << '.' << std::endl;
  } // for (l)

  std::ofstream stat_out;

  stat_out.open((file_name_base+"-uniq-deg_val").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening deg_val file failed");
  for (size_t d=0; d<deg_limit; ++d)
  {
    stat_out << "Degree " << d << ": ";
    for (size_t v=0; v<=d; ++v)
      stat_out << deg_val.multiplicity(d*(d+1)/2+v) << (v<d ? ',' : '\n');
  }
  stat_out.close();

  stat_out.open((file_name_base+"-uniq@q=0").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=0 file failed");
  for (tally_vec::Index i=0; i<q_is_0.size(); q_is_0.advance(i))
    stat_out << i << ": " << q_is_0.multiplicity(i) << ".\n";
  stat_out.close();

  stat_out.open((file_name_base+"-uniq-contents").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening contents file failed");
  for (tally_vec::Index i=1; i<cont.size(); cont.advance(i))
    stat_out << "Contents " << i << ": " << cont.multiplicity(i) << ".\n";
  stat_out.close();

  stat_out.open((file_name_base+"-uniq@q=1").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening q=1 file failed");
  q_is_1.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-uniq-leading").c_str());
  if (not stat_out.is_open())
    throw std::runtime_error("Opening leading term file failed");
  leading.write_to(stat_out);
  stat_out.close();

  stat_out.open((file_name_base+"-uniq-coefficients").c_str());
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

  if (argc!=3 and argc!=5)
  {
    std::cerr <<
      "Usage: polstat [-q] block-file poly-file row-file "
      "[tally-prim-file tally-mult-prim]\n";
    exit(1);
  }

  std::ifstream block_file(argv[0]);
  std::ifstream pol_file(argv[1]);
  std::ifstream row_file(argv[2]);
  if (not block_file.is_open()
      or not pol_file.is_open()
      or not row_file.is_open())
  {
    std::cerr << "Failure opening file(s).\n";
    exit(1);
  }

  std::string file_name_base= getFileName
    ("Base file name for statistics output: ");

  atlas::filekl::block_info bi(block_file);
  atlas::filekl::polynomial_info pi(pol_file);
  atlas::filekl::progress_info ri(row_file);

  if (argc==3)
  {
    atlas::filekl::scan_polynomials(bi,pi,ri,file_name_base);
  }
  else
  {
    std::ifstream tally_file(argv[3]);
    std::ifstream tally_mu_file(argv[4]);
    if (not tally_file.is_open()
	or not tally_mu_file.is_open())
    {
      std::cerr << "Failure opening file(s).\n";
      exit(1);
    }
    atlas::filekl::tally_vec prim_mu(tally_file);
    atlas::filekl::tally_vec tot_mu(tally_mu_file);

    atlas::filekl::scan_polynomials(bi,pi,ri,prim_mu,tot_mu,file_name_base);
  }





}
