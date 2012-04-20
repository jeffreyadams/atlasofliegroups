/*
  This is block_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "block_io.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <cassert>
#include <map>

#include "polynomials.h"
#include "blocks.h"
#include "kgb.h"     // |kgb.size()|
#include "complexredgp.h" // |twoRho| in |nu_block::print|
#include "kl.h"

#include "basic_io.h"	// operator |<<|
#include "ioutils.h"	// |digits|
#include "prettyprint.h" // |printWeylElt|


/*****************************************************************************

  Output functions for |Block|, defined in sources/kl/blocks.h

******************************************************************************/

namespace atlas {


// We start with defining a method from |Block_base|

namespace blocks {

// print a derived block object, with virtual |print| after Cayley transforms
std::ostream& Block_base::print_to(std::ostream& strm,
				   bool as_invol_expr) const
{
  // compute maximal width of entry
  int width = ioutils::digits(size()-1,10ul);
  int xwidth = ioutils::digits(xsize()-1,10ul);
  int ywidth = ioutils::digits(ysize()-1,10ul);
  int lwidth = ioutils::digits(length(size()-1),10ul);

  const int pad = 2;

  for (BlockElt z=0; z<size(); ++z)
  {
    // print entry number and corresponding orbit pair
    strm << std::setw(width) << z
	 << '(' << std::setw(xwidth) << x(z)
	 << ',' << std::setw(ywidth) << y(z) << "):";

    // print length
    strm << std::setw(lwidth+pad) << length(z) << std::setw(pad) << "";

    // print descents
    block_io::printDescent(strm,descent(z),rank());

    // print cross actions
    for (weyl::Generator s = 0; s < rank(); ++s)
    {
      strm << std::setw(width+pad);
      if (cross(s,z)==blocks::UndefBlock) strm << '*';
      else strm << cross(s,z);
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < rank(); ++s)
    {
      BlockEltPair p =
	isWeakDescent(s,z) ? inverseCayley(s,z) : cayley(s,z);
      strm << '(' << std::setw(width);
      if (p.first ==blocks::UndefBlock) strm << '*'; else strm << p.first;
      strm << ',' << std::setw(width);
      if (p.second==blocks::UndefBlock) strm << '*'; else strm << p.second;
      strm << ')' << std::setw(pad) << "";
    }

    // finish with derived class specific output
    print(strm,z,as_invol_expr) << std::endl;
  } // |for (z)|

  return strm;
} // |print_on|


std::ostream& Block::print
  (std::ostream& strm, BlockElt z,bool as_invol_expr) const
{
  int cwidth = ioutils::digits(max_Cartan(),10ul);
  strm << std::setw(cwidth) << Cartan_class(z) << std::setw(2) << "";

  // print root datum involution
  if (as_invol_expr)
    prettyprint::printInvolution(strm,involution(z),twistedWeylGroup());
  else
    prettyprint::printWeylElt(strm,involution(z),weylGroup());

  return strm ;
}

std::ostream& gamma_block::print
  (std::ostream& strm, BlockElt z,bool as_invol_expr) const
{
  int xwidth = ioutils::digits(kgb.size()-1,10ul);

  RatWeight ls = local_system(z);
  strm << "(=" << std::setw(xwidth) << parent_x(z)
       << ',' << std::setw(3*ls.size()+3) << ls
       << ')' << std::setw(2) << "";

  // print root datum involution
  if (as_invol_expr)
    prettyprint::printInvolution(strm,involution(z),kgb.twistedWeylGroup());
  else
    prettyprint::printWeylElt(strm,involution(z),kgb.weylGroup());

  return strm ;
}

std::ostream& non_integral_block::print
  (std::ostream& strm, BlockElt z,bool as_invol_expr) const
{
  int xwidth = ioutils::digits(kgb.size()-1,10ul);
  RatWeight ll=y_part(z);

  strm << (survives(z) ? '*' : ' ')
       << "(x=" << std::setw(xwidth) << parent_x(z)
       << ", nu=" << std::setw(2*ll.size()+5) << nu(z);
//strm << ',' << std::setw(2*ll.size()+5) << ll;
  strm << ",lam=rho+" << std::setw(2*ll.size()+3) << lambda_rho(z);
  strm << ')' << std::setw(2) << "";

  const TwistedInvolution& ti = kgb.involution(parent_x(z));
  const TwistedWeylGroup& tW = kgb.twistedWeylGroup();
  // print root datum involution
  if (as_invol_expr) prettyprint::printInvolution(strm,ti,tW);
  else prettyprint::printWeylElt(strm,ti,tW.weylGroup());

  return strm ;
}


} // |namespace blocks|
/*****************************************************************************

        Chapter I -- Functions declared in block_io.h


******************************************************************************/

namespace block_io {



/*
  Synopsis: outputs the unitary elements in the block in the usual format.

  Explanation: as David explained to me, the unitary representations are
  exactly those for which the support of the involution consists entirely
  of descents.

  NOTE: checking that condtion is awkward here, because currently blocks
  do not have direct access to the descents as a bitset!
*/
std::ostream& printBlockU(std::ostream& strm, const Block& block)
{
  using namespace blocks;
  using namespace descents;
  using namespace kgb;
  using namespace prettyprint;

  // compute maximal width of entry
  int width = ioutils::digits(block.size()-1,10ul);
  int xwidth = ioutils::digits(block.xsize()-1,10ul);
  int ywidth = ioutils::digits(block.ysize()-1,10ul);
  int lwidth = ioutils::digits(block.length(block.size()-1),10ul);
  int cwidth = ioutils::digits(block.Cartan_class(block.size()-1),10ul);
  const int pad = 2;

  for (size_t j = 0; j < block.size(); ++j)
  {
    for (size_t s = 0; s < block.rank(); ++s)
      if (block.involutionSupport(j)[s] and not block.isWeakDescent(s,j))
	goto nextj; // representation is not unitary, continue outer loop

    // print entry number and corresponding orbit pair
    strm << std::setw(width) << j;
    strm << '(' << std::setw(xwidth) << block.x(j);
    strm << ',';
    strm << std::setw(ywidth) << block.y(j) << ')';
    strm << ':';

    // print length
    strm << std::setw(lwidth+pad) << block.length(j);

    // print Cartan class
    strm << std::setw(cwidth+pad) << block.Cartan_class(j);
    strm << std::setw(pad) << "";

    // print descents
    printDescent(strm,block.descent(j),block.rank());
    strm << std::setw(pad) << "";

    // print descents filtered
    printDescent(strm,block.descent(j),block.rank(),
		 block.involutionSupport(j));
    strm << std::setw(pad) << "";

#if 0
    // print cross actions
    for (size_t s = 0; s < block.rank(); ++s) {
      BlockElt z = block.cross(s,j);
      if (z == UndefBlock)
	strm << std::setw(width+pad) << '*';
      else
	strm << std::setw(width+pad) << z;
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < block.rank(); ++s)
    {
      BlockEltPair z = block.isWeakDescent(s,j)
                             ? block.inverseCayley(s,j)
	                     : block.cayley(s,j);
      strm << '(' << std::setw(width);
      if (z.first == blocks::UndefBlock)
	strm << '*';
      else
	strm << z.first;
      strm << ',' << std::setw(width);
      if (z.second == blocks::UndefBlock)
	strm << '*';
      else
	strm << z.second;
      strm << ')' << std::setw(pad) << "";
    }
    strm << ' ';
#endif

    // print root datum involution as involution reduced expression
    prettyprint::printInvolution
      (strm,block.involution(j),block.twistedWeylGroup()) << std::endl;

  nextj:
    continue;
  } // |for(j)|

  return strm;
} // |printBlockU|


/*
  Synopsis: outputs the descent status for the various generators
*/
std::ostream& printDescent(std::ostream& strm,
			   const DescentStatus& ds,
			   size_t rank, RankFlags mask)
{
  strm << '[';

  for (size_t s = 0; s < rank; ++s)
  {
    if (s!=0)
      strm << ',';
    if (not mask.test(s))
      strm << "* ";
    else
      switch (ds[s]) {
      case DescentStatus::ComplexDescent:
	strm << "C-";
	break;
      case DescentStatus::ComplexAscent:
	strm << "C+";
	break;
      case DescentStatus::ImaginaryCompact:
	strm << "ic";
	break;
      case DescentStatus::RealNonparity:
	strm << "rn";
	break;
      case DescentStatus::ImaginaryTypeI:
	strm << "i1";
	break;
      case DescentStatus::ImaginaryTypeII:
	strm << "i2";
	break;
      case DescentStatus::RealTypeI:
	strm << "r1";
	break;
      case DescentStatus::RealTypeII:
	strm << "r2";
	break;
      default: // should not happen!
	assert(false);
	break;
      }
  }

  strm << ']';

  return strm;
}

std::ostream& print_KL(std::ostream& f, non_integral_block& block, BlockElt z)
{
  // silently fill the whole KL table
  const kl::KLContext& klc = block.klc(block.size()-1,false);

  typedef Polynomial<int> Poly;
  typedef std::map<BlockElt,Poly> map_type;
  map_type acc; // non-zero $x'\mapsto\sum_{x\downarrow x'}\eps(z/x)P_{x,z}$
  unsigned int parity = block.length(z)%2;
  for (size_t x = 0; x <= z; ++x)
  {
    const kl::KLPol& pol = klc.klPol(x,z);
    if (not pol.isZero())
    {
      Poly p(pol); // convert
      if (block.length(x)%2!=parity)
	p*=-1;
      BlockEltList nb=block.survivors_below(x);
      for (size_t i=0; i<nb.size(); ++i)
      {
	std::pair<map_type::iterator,bool> trial =
	  acc.insert(std::make_pair(nb[i],p));
	if (not trial.second) // failed to create a new entry
	  trial.first->second += p;
      } // |for (i)| in |nb|
    } // |if(pol!=0)|
  } // |for (x<=z)|


  f << (block.singular_simple_roots().any() ? "(cumulated) " : "")
    << "KL polynomials (-1)^{l(" << z << ")-l(x)}*P_{x," << z << "}:\n";
  int width = ioutils::digits(z,10ul);
  for (map_type::const_iterator it=acc.begin(); it!=acc.end(); ++it)
  {
    BlockElt x = it->first;
    const Poly& pol = it->second;
    if (not pol.isZero())
    {
      f << std::setw(width) << x << ": ";
      prettyprint::printPol(f,pol,"q") << std::endl;
    }
  }
  return f;
}

} // namespace block_io|

} // namespace atlas|
