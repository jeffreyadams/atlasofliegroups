/*
  This is block_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "block_io.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <cassert>
#include <map>

#include "polynomials.h"
#include "kgb.h"     // |kgb.size()|
#include "innerclass.h" // |twoRho| in |nu_block::print|
#include "blocks.h"
#include "common_blocks.h"
#include "ext_block.h"
#include "kl.h"
#include "repr.h"

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
  int width = ioutils::digits(size()==0 ? 0 : size()-1,10ul);
  int xwidth = ioutils::digits(max_x(),10ul);
  int ywidth = ioutils::digits(max_y(),10ul);
  int lwidth = ioutils::digits(size()==0 ? 0 : length(size()-1),10ul);

  const int pad = 2;

  bool traditional = dynamic_cast<const Block*>(this)!=nullptr;

  for (BlockElt z=0; z<size(); ++z)
  {
    // print entry number and corresponding orbit pair
    strm << std::setw(width) << z;
    if (traditional) // prining "local" x,y is confusing in other cases
      strm << '(' << std::setw(xwidth) << x(z)
	   << ',' << std::setw(ywidth) << y(z) << "):";
    else
      strm << ':';

    // print length
    strm << std::setw(lwidth+pad) << length(z) << std::setw(pad) << "";

    // print descents
    block_io::printDescent(strm,descent(z),rank());

    // print cross actions
    for (weyl::Generator s = 0; s < rank(); ++s)
    {
      strm << std::setw(width+pad);
      if (cross(s,z)==UndefBlock) strm << '*';
      else strm << cross(s,z);
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < rank(); ++s)
    {
      BlockEltPair p =
	isWeakDescent(s,z) ? inverseCayley(s,z) : cayley(s,z);
      strm << '(' << std::setw(width);
      if (p.first ==UndefBlock) strm << '*'; else strm << p.first;
      strm << ',' << std::setw(width);
      if (p.second==UndefBlock) strm << '*'; else strm << p.second;
      strm << ')' << std::setw(pad) << "";
    }

    // finish with derived class specific output
    print(strm,z,as_invol_expr) << std::endl;
  } // |for (z)|

  return strm;
} // |print_to|


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

std::ostream& param_block::print
  (std::ostream& strm, BlockElt z,bool as_invol_expr) const
{
  const KGB& kgb = rc.kgb();
  unsigned int xwidth = ioutils::digits(highest_x,10ul);
  unsigned int rk = root_datum().semisimpleRank();

  strm << (survives(z) ? '*' : ' ')
       << "(x=" << std::setw(xwidth) << x(z)
       << ",lam_rho=" << std::setw(3*rk+1) << lambda_rho(z)
       << ", nu=" << std::setw(3*rk+3) << nu(z)
       << ')' << std::setw(2) << "";

  const TwistedInvolution& ti = kgb.involution(x(z));
  const TwistedWeylGroup& tW = kgb.twistedWeylGroup();
  // print root datum involution
  if (as_invol_expr) prettyprint::printInvolution(strm,ti,tW);
  else prettyprint::printWeylElt(strm,ti,tW.weylGroup());

  return strm ;
}

std::ostream& common_block::print
  (std::ostream& strm, BlockElt z,bool as_invol_expr) const
{
  const KGB& kgb = rc.kgb();
  unsigned int xwidth = ioutils::digits(highest_x,10ul);
  unsigned int rk = root_datum().semisimpleRank();

  strm << " (x=" << std::setw(xwidth) << x(z)
       << ",gamma-lambda=" << std::setw(5*rk+1) << gamma_lambda(z)
       << ')' << std::setw(2) << "";

  const TwistedInvolution& ti = kgb.involution(x(z));
  const TwistedWeylGroup& tW = kgb.twistedWeylGroup();
  // print root datum involution
  if (as_invol_expr) prettyprint::printInvolution(strm,ti,tW);
  else prettyprint::printWeylElt(strm,ti,tW.weylGroup());

  return strm ;
}


} // |namespace blocks|

namespace ext_block {

std::ostream& ext_block::print_to (std::ostream& strm) const
{
  if (size()==0)
    return strm << "Empty extended block (block is not twist-stable)."
		<< std::endl;
  // compute maximal width of entry
  int width = ioutils::digits(z(size()-1),10ul);
  int lwidth = ioutils::digits(length(size()-1),10ul);

  const int pad = 2;

  for (BlockElt n=0; n<size(); ++n)
  {
    // print parent entry number
    strm << std::setw(width) << z(n);

    // print length
    strm << std::setw(lwidth+pad) << length(n) << std::setw(pad) << "";

    // print descents types
    if (rank()==0)
      strm << '[';
    else
      for (weyl::Generator s=0; s<rank(); ++s)
	strm << (s==0 ? '[' : ',') << descent_code(descent_type(s,n));
    strm << ']';

    // print cross actions
    for (weyl::Generator s = 0; s < rank(); ++s)
    {
      auto csn = cross(s,n);
      if (csn==UndefBlock) strm <<  std::setw(width+pad+1) << '*';
      else strm << std::setw(pad)
		<< (n==csn or epsilon(s,n,csn)>=0 ? '+' : '-' )
		<< std::setw(width) << z(cross(s,n));
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < rank(); ++s)
    {
      strm << '(';
      if (is_complex(descent_type(s,n)) or data[s][n].links.first==UndefBlock)
	strm << std::setw(width+1) << '*';
      else
	strm << (info[n].flips[0][s] ? '-' : '+')
	     << std::setw(width) << z(data[s][n].links.first);
      strm << ',';
      if (has_double_image(descent_type(s,n))
	  and data[s][n].links.second!=UndefBlock)
	strm << (info[n].flips[1][s] ? '-' : '+')
	     << std::setw(width) << z(data[s][n].links.second);
      else
	strm << std::setw(width+1) << '*';
      strm << ')' << std::setw(pad) << "";
    }

    strm << std::endl;
  } // |for (n)|

  return strm;
} // |print_to|

const char* descent_code(DescValue v)
{
  switch(v)
  {
  case one_complex_ascent: return "1C+  ";
  case one_complex_descent: return "1C-  ";
  case one_imaginary_single: return "1i1  ";
  case one_real_pair_fixed: return "1r1f ";
  case one_imaginary_pair_fixed: return "1i2f ";
  case one_real_single: return "1r2  ";
  case one_imaginary_pair_switched: return "1i2s ";
  case one_real_pair_switched: return "1r1s ";
  case one_real_nonparity: return "1rn  ";
  case one_imaginary_compact: return "1ic  ";

  case two_complex_ascent: return "2C+  ";
  case two_complex_descent: return "2C-  ";
  case two_semi_imaginary: return "2Ci  ";
  case two_semi_real: return "2Cr  ";
  case two_imaginary_single_single: return "2i11 ";
  case two_real_double_double: return "2r11 ";
  case two_imaginary_single_double_fixed: return "2i12f";
  case two_real_single_double_fixed: return "2r21f";
  case two_imaginary_double_double: return "2i22 ";
  case two_real_single_single: return "2r22 ";
  case two_imaginary_single_double_switched: return "2i12s";
  case two_real_single_double_switched: return "2r21s";
  case two_real_nonparity: return "2rn  ";
  case two_imaginary_compact: return "2ic  ";

  case three_complex_ascent: return "3C+  ";
  case three_complex_descent: return "3C-  ";
  case three_semi_imaginary: return "3Ci  ";
  case three_real_semi: return "3r   ";
  case three_imaginary_semi: return "3i   ";
  case three_semi_real: return "3Cr  ";
  case three_real_nonparity: return "3rn  ";
  case three_imaginary_compact: return "3ic  ";
  }
  assert(false); return nullptr;
}

} // |namespace ext_block|

/*****************************************************************************

        Chapter I -- Functions declared in block_io.h


******************************************************************************/

namespace block_io {



/*
  Print the unitary elements in the block in the usual block element format.

  Explanation: as David explained to me (Fokko), the unitary representations
  are exactly those for which the support of the involution consists entirely
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
  int width = ioutils::digits(block.size()==0 ? 0 : block.size()-1,10ul);
  int xwidth = ioutils::digits(block.max_x(),10ul);
  int ywidth = ioutils::digits(block.max_y(),10ul);
  int lwidth = ioutils::digits(
        block.size()==0 ? 0 : block.length(block.size()-1),10ul);
  int cwidth = ioutils::digits(
        block.size()==0 ? 0 : block.Cartan_class(block.size()-1),10ul);
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
      if (z.first == UndefBlock)
	strm << '*';
      else
	strm << z.first;
      strm << ',' << std::setw(width);
      if (z.second == UndefBlock)
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


// Output the descent status for the various generators
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


std::ostream& print_twist(std::ostream& strm, const Block_base& block)
{
  if (block.Hermitian_dual(0)==UndefBlock)
    return strm << "Block is not stable under twist" << std::endl;

  std::ostringstream os; BlockElt count=0;

  os << "Elements fixed under twist: ";
  for (BlockElt z=0; z<block.size(); ++z)
    if (block.Hermitian_dual(z)==z)
    {
      if (count>0)
	os << ", ";
      os << z;
      ++count;
    }
  ioutils::foldLine(strm,os.str()) << std::endl;

  os.str(""); // clear string for rewriting
  strm << "Elements interchanged by twist: " << std::endl;
  for (BlockElt z=0; z<block.size(); ++z)
    if (block.Hermitian_dual(z)>z)
      os << '(' << z << ' ' << block.Hermitian_dual(z) << ") ";

  ioutils::foldLine(strm,os.str(),"",") ") << std::endl;

  return strm << "Total " << count << " fixed elements out of " << block.size()
	      << std::endl;
}

std::ostream& print_KL(std::ostream& f, param_block& block, BlockElt z)
{
  // silently fill the whole KL table
  const kl::KLContext& klc = block.klc(z,false);

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
      auto finals=block.finals_for(x);
      for (BlockElt final : finals)
      {
	std::pair<map_type::iterator,bool> trial =
	  acc.insert(std::make_pair(final,p));
	if (not trial.second) // failed to create a new entry
	  trial.first->second += p;
      } // |for (final : finals)|
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
      pol.print(f << std::setw(width) << x << ": ","q") << std::endl;
  }
  return f;
}

} // |namespace block_io||

} // |namespace atlas||
