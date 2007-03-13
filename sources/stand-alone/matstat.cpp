#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <cassert>
#include <iostream>
#include <stdexcept>

#include "filekl.h"
#include "basic_io.h"


bool verbose=true;
bool with_mu=false;

namespace atlas {
  namespace filekl {

class tally_vec
{
  typedef std::map<KLIndex,ullong> map_type;
  std::vector<unsigned char> count; // number of times seen
  map_type overflow; // for multiplicities >=256
  KLIndex max;
  ullong total;

public:
  tally_vec(size_t limit) : count(0), overflow(), max(0), total(0)
  { count.reserve(limit); }

  bool tally (KLIndex i); // increase count for i by 1; tell whether new
  bool tally (KLIndex i,unsigned int multiplicity); // same with multiplicity
  KLIndex size() const { return max+1; } // size of collection now tallied
  ullong multiplicity (KLIndex i) const
  {
    if (i<count.size())
      return count[i]!=255 ? count[i] : overflow.find(i)->second;
    map_type::const_iterator it=overflow.find(i);
    return it==overflow.end() ? 0 : it->second;
  }
  ullong sum() const { return total; }

  void advance(KLIndex& i) const; // like ++ and -- where |i| iterates
  void lower(KLIndex& i) const;   // over indices with nonzero multiplicity

  void write_to(std::ostream& out) const;
};

bool tally_vec::tally(KLIndex i)
{
  ++total;
  if (i<count.size()) // then |i| is already recorded in |count| or |overflow|
  {
    if (count[i]!=255)
      if (++count[i]==255) overflow[i]=255; // create entry upon hitting 255
      else return count[i]==1;  // it just might be the first occurrence of |i|
    else ++overflow[i];
    return false;
  }
  if (i>max) max=i;
  if (i<count.capacity()) // then slot for |i| can be created
  {
    while (count.size()<i) count.push_back(0);
    count.push_back(1); return true;
  }

  // now |i>=count.capacity()|, it must be added to overflow
  std::pair<map_type::iterator,bool> p =overflow.insert(std::make_pair(i,1));
  if (not p.second) ++p.first->second; // if already recorded, increase tally
  return p.second;
}

bool tally_vec::tally(KLIndex i,unsigned int multiplicity)
{
  total+=multiplicity;
  if (i<count.size()) // then |i| is already recorded in |count| or |overflow|
  {
    if (count[i]!=255)
      if (count[i]+multiplicity<255) // then |count[i]| not yet saturated
      {
	bool result= count[i]==0; // this might just be the first occurence
	count[i]+=multiplicity;
	return result;
      }
      else
      {
	overflow[i]=count[i]+multiplicity; // create new entry
	count[i]=255; // and mark |count[i]| as saturated
      }
    else overflow[i]+=multiplicity;
    return false;
  }
  if (i>max) max=i;
  if (i<count.capacity()) // then slot for |i| can be created
  {
    while (count.size()<i) count.push_back(0);
    count.push_back(multiplicity); return true;
  }

  // now |i>=count.capacity()|, it must be added to overflow
  std::pair<map_type::iterator,bool> p
    =overflow.insert(std::make_pair(i,multiplicity));
  if (not p.second) // then it was already recorded
    p.first->second+=multiplicity; // so increase tally for |i| instead
  return p.second; // return whether |i| was previously unrecorded
}

void tally_vec::advance(KLIndex& i) const
{
  if (i>=max)
  { i=size(); return; } // avoid fruitless search
  ++i; // make sure we advance at least by one
  while (i<count.size())
    if (count[i]!=0) return;
    else ++i;

  // now find the first entry |j| in overflow with |j>=i|
  map_type::const_iterator it=overflow.lower_bound(i);
  assert(it!=overflow.end()); // there is at least one |j>=i|
  i=it->first;
}

void tally_vec::lower(KLIndex& i) const
{
  if (i>count.size()) // then find last entry |j| in overflow with |j<i|
  {
    map_type::const_iterator it=overflow.lower_bound(i);
    if (it!=overflow.begin() and (--it)->first>=count.size())
    { i=it->first; return; }

    else i=count.size(); // and fall through to search in |count|
  }

  // now |i<=count.size()|; find last entry |j<i| in count with |count[j]!=0|
  while (i-->0)
    if (count[i]!=0) return;

  // we should not come here; there should be a sentinel |count[0]!=0|
  throw std::runtime_error("lowering hit untallied 0");
}

void tally_vec::write_to(std::ostream& out) const
{
  basic_io::put_int(size(),out);
  basic_io::put_int(overflow.size(),out);
  for (size_t i=0; i<size(); ++i) out.put(count[i]);
  for (map_type::const_iterator it=overflow.begin(); it!=overflow.end(); ++it)
  {
    basic_io::put_int(it->first,out);
    basic_io::put_int(it->second,out);
  }
}

/* A function that computes the vector saying for any strongly primitve
   element $x$ for $y$ how many $x_0$'s primitivise to $x$ for $y$.
*/
std::vector<unsigned int> prim_multiplicities(matrix_info& m,BlockElt y)
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

  tally_vec mu_mu(65536); mu_mu.tally(0); // provide sentinel
  for (KLIndex i=0; i<t.size(); t.advance(i))
    mu_mu.tally(t.multiplicity(i));

  length_out << "frequencies of lowest multiplicities:\n";
  for (ullong i=1; i<=64 and i<mu_mu.size(); mu_mu.advance(i))
    length_out << "multiplicity " << i << " occurred for "
	       << mu_mu.multiplicity(i) << " distinct polynomials.\n";

  length_out << "most frequenly encountered polynomials:\n";

  const unsigned int tier=256;

  std::vector<ullong> high_multiplicities; high_multiplicities.reserve(tier);
  ullong i=mu_mu.size();
  for (size_t j=0; j<tier and i>0; ++j)
  { mu_mu.lower(i); high_multiplicities.push_back(i); }

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
