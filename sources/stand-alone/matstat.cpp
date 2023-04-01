#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cassert>
#include <stdexcept>

#include "../Atlas.h"
#include "filekl_in.h"
#include "basic_io.h"
#include "tally.h"

#include "blocks.h"

bool verbose=true;
bool with_mu=false;

namespace atlas {
  namespace filekl {

typedef tally::TallyVec<unsigned char> tally_vec;

/* A function that computes the vector saying for any strongly primitive
   element $x$ for $y$ how many $x_0$'s primitivise to $x$ for $y$.
*/
std::vector<unsigned int>
prim_multiplicities(matrix_info& m,BlockElt y)
{
  std::vector<unsigned int> multiplicity(y+1,0);
  for (BlockElt x0=0; x0<=y; ++x0)
  {
    BlockElt x=m.primitivize(x0,y);
    if (x<=y) ++multiplicity[x];
  }
  const strong_prim_list& spy=m.strongly_primitives(y);
  std::vector<unsigned int> result(spy.size());
  for (size_t i=0; i<spy.size(); ++i)
    result[i]=multiplicity[spy[i]];
  return result;
}


void scan_matrix(matrix_info& m, size_t n_pol, bool with_multiplicities,
		 std::ostream& y_out,
		 std::ostream& tally_out,
		 std::ostream& length_out)
{
  tally_vec t(n_pol);
  t.tally(0); // tally the zero polynomial, else it would be found missing
  ullong nr_sp=0,nr_nonzero=0;
  ullong last_nr_sp=0, last_nr_nonzero=0, last_npol=0;
  ullong prev_l_nr_sp=0, prev_l_nr_nonzero=0, prev_l_npol=0;
  size_t l=0;

  for (BlockElt y=0; y<m.block_size(); ++y)
  {
    while (y==m.first_of_length(l+1))
    {
      length_out << "l=" << l << " @y=" << y
		 << " (+" << y-m.first_of_length(l) << "), sp: "
		 << nr_sp-prev_l_nr_sp << ", nonzero: "
		 << nr_nonzero-prev_l_nr_nonzero
		 << ", new polys: " << t.size()-prev_l_npol << std::endl;
      prev_l_nr_sp=nr_sp; prev_l_nr_nonzero=nr_nonzero; prev_l_npol=t.size();
      ++l;
    }

    if (verbose) std::cerr << y << '\r';

    const strong_prim_list& spy=m.strongly_primitives(y);
    std::vector<unsigned int> mu =prim_multiplicities(m,y);

    nr_sp+=spy.size();
    for (size_t i=0; i<spy.size(); ++i)
    {
      nr_nonzero+=mu[i];
      if (with_multiplicities) t.tally(m.find_pol_nr(spy[i],y),mu[i]);
      else t.tally(m.find_pol_nr(spy[i],y));
    }
    basic_io::put_int(nr_sp-last_nr_sp,y_out);    // strongly primitives seen
    basic_io::put_int(nr_nonzero-last_nr_nonzero,y_out); // nonzero pols seen
    basic_io::put_int(t.size()-last_npol,y_out); // distinct polynomials seen
    last_nr_sp=nr_sp; last_nr_nonzero=nr_nonzero; last_npol=t.size();
  } // for(y)

  // summary of final level
  length_out << "l=" << l << " @y=" << m.block_size() << " (+"
	     << m.block_size()-m.first_of_length(l) << "), sp: "
	     << nr_sp-prev_l_nr_sp << ", nonzero: "
	     << nr_nonzero-prev_l_nr_nonzero
	     << ", new polys: " << t.size()-prev_l_npol << std::endl;
  length_out << "Totals. strong prim: " << nr_sp
	     << ", nonzero: " << nr_nonzero
	     << ", polynomials: " << t.size() << '.' << std::endl;

  t.write_to(tally_out);

  tally_vec mu_mu=t.derived<unsigned char>(65536);

  length_out << "frequencies of lowest multiplicities:\n";
  for (ullong i=1; i<=64 and i<mu_mu.size(); mu_mu.advance(i))
    length_out << "multiplicity " << i << " occurred for "
	       << mu_mu.multiplicity(i) << " distinct polynomials.\n";

  length_out << "most frequenly encountered polynomials:\n";

  const unsigned int tier=256;

  std::vector<ullong> high_multiplicities; high_multiplicities.reserve(tier);
  ullong i=mu_mu.size();
  for (size_t j=0; j<tier and mu_mu.lower(i); ++j)
    high_multiplicities.push_back(i);

  std::vector<std::vector<KLIndex> > top_tier
    (high_multiplicities.size(),std::vector<KLIndex>());
  ullong limit=high_multiplicities.back();

  for (KLIndex i=0; i<t.size(); t.advance(i))
  {
    size_t mu=t.multiplicity(i);
    if (mu>=limit)
      for (size_t j=0; j<top_tier.size(); ++j)
	if (mu==high_multiplicities[j])
	{ top_tier[j].push_back(i); break; }
  }

  for (size_t j=0; j<top_tier.size(); ++j)
  {
    if (top_tier[j].size()==1)
      length_out << "Polynomial #" << top_tier[j][0]
		 << " occurs with multiplicity " << high_multiplicities[j]
		 << ".\n" ;
    else
    {
      length_out << "Multiplicity " << high_multiplicities[j]
		 << " occurred for " << top_tier[j].size()
		 << " distinct polynomials";
      for (size_t i=0; i<top_tier[j].size(); ++i)
	length_out << (i==0 ? ':' : ',') << " #" << top_tier[j][i];
      length_out << ".\n";
    }
  }

}


  } // namespace filekl
} // namespace atlas

int main(int argc, char** argv)
{
  --argc; ++argv; // read and skip program name

  if (argc>0 and std::string(*argv)=="-q") { verbose=false; --argc; ++argv;}
  if (argc>0 and std::string(*argv)=="-w") { with_mu=true; --argc; ++argv;}

  //  atlas::constants::initConstants();
  if (argc!=5)
  {
    std::cerr <<
      "Usage: matstat [-q] [-w] block-file matrix-file count"
      " row-file tally-file\n";
    exit(1);
  }

  std::ifstream block_file(argv[0]);
  std::ifstream matrix_file(argv[1]);
  if (not  block_file.is_open() or not matrix_file.is_open())
  {
    std::cerr << "Failure opening file(s).\n";
    exit(1);
  }

  atlas::filekl::matrix_info mi(block_file,matrix_file);

  // test last argument before opening output file
  std::istringstream s(argv[2]);
  unsigned long long int limit=0; s>>limit;

  if (s.fail())
  {
    std::cerr << "Illegal polynomial count: " << argv[2] << ".\n";
    exit(1);
  }

  std::ofstream row_out(argv[3],
		 std::ios_base::out
		 | std::ios_base::trunc
		 | std::ios_base::binary);
  std::ofstream tally_out(argv[4],
		 std::ios_base::out
		 | std::ios_base::trunc
		 | std::ios_base::binary);
  if (not row_out.is_open() or not tally_out.is_open())
  {
    std::cerr << "Failure opening file(s) for writing.\n";
    exit(1);
  }

  atlas::filekl::scan_matrix(mi,limit,with_mu,row_out,tally_out,std::cout);


}
