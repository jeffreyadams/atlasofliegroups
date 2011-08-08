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

#include "blocks.h"

#include "ioutils.h"	// |digits|
#include "prettyprint.h" // |printWeylElt|


/*****************************************************************************

  Output functions for |Block|, defined in sources/kl/blocks.h

******************************************************************************/

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in block_io.h


******************************************************************************/

namespace block_io {

// print a derived block object, with virtual |print| after Cayley transforms
std::ostream& print_block(std::ostream& strm, const Block_base& b)
{
  // compute maximal width of entry
  int width = ioutils::digits(b.size()-1,10ul);
  int xwidth = ioutils::digits(b.xsize()-1,10ul);
  int ywidth = ioutils::digits(b.ysize()-1,10ul);
  int lwidth = ioutils::digits(b.length(b.size()-1),10ul);

  const int pad = 2;

  for (BlockElt z=0; z<b.size(); ++z)
  {
    // print entry number and corresponding orbit pair
    strm << std::setw(width) << z
	 << '(' << std::setw(xwidth) << b.x(z)
	 << ',' << std::setw(ywidth) << b.y(z) << "):";

    // print length
    strm << std::setw(lwidth+pad) << b.length(z) << std::setw(pad) << "";

    // print descents
    printDescent(strm,b.descent(z),b.rank());

    // print cross actions
    for (weyl::Generator s = 0; s < b.rank(); ++s)
    {
      strm << std::setw(width+pad);
      if (b.cross(s,z)==blocks::UndefBlock) strm << '*';
      else strm << b.cross(s,z);
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < b.rank(); ++s)
    {
      BlockEltPair p =
	b.isWeakDescent(s,z) ? b.inverseCayley(s,z) : b.cayley(s,z);
      strm << '(' << std::setw(width);
      if (p.first ==blocks::UndefBlock) strm << '*'; else strm << p.first;
      strm << ',' << std::setw(width);
      if (p.second==blocks::UndefBlock) strm << '*'; else strm << p.second;
      strm << ')' << std::setw(pad) << "";
    }

    // derived class specific output
    b.print(strm,z);

    // print root datum involution
    prettyprint::printWeylElt(strm,b.involution(z),b.weylGroup());

    strm << std::endl;
  }

  return strm;
} // |print_block|

/*
  Print the data from block to strm.

  Explanation: for each parameter, we output the length, Cartan class,
  cross-actions and Cayley or inverse Cayley transform(s) for each generator,
  and the underlying root datum permutation (or rather, the corresponding Weyl
  group element). We use a '*' for absent (inverse) Cayley transforms.
*/
std::ostream& printBlock(std::ostream& strm, const Block& block)
{
  // compute maximal width of entry
  int width = ioutils::digits(block.size()-1,10ul);
  int xwidth = ioutils::digits(block.xsize()-1,10ul);
  int ywidth = ioutils::digits(block.ysize()-1,10ul);
  int lwidth = ioutils::digits(block.length(block.size()-1),10ul);
  int cwidth = ioutils::digits(block.max_Cartan(),10ul);
  const int pad = 2;

  for (BlockElt j = 0; j < block.size(); ++j)
  {
    // print entry number and corresponding orbit pair
    strm << std::setw(width) << j;
    strm << '(' << std::setw(xwidth) << block.x(j);
    strm << ',';
    strm << std::setw(ywidth) << block.y(j) << ')';
    strm << ':';

    // print length
    strm << std::setw(lwidth+pad) << block.length(j) << std::setw(pad) << "";

    // print descents
    printDescent(strm,block.descent(j),block.rank());

    // print cross actions
    for (size_t s = 0; s < block.rank(); ++s) {
      BlockElt z = block.cross(s,j);
      if (z == blocks::UndefBlock)
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

    // print Cartan class
    strm << std::setw(cwidth) << block.Cartan_class(j)
	 << std::setw(pad) << "";

    // print root datum involution
    prettyprint::printWeylElt(strm,block.involution(j),block.weylGroup());

    strm << std::endl;
  }

  return strm;
} // |printBlock|


/*
  Print the data from block to strm.

  Explanation: for each parameter, we output the cross-actions and
  cayley-actions for each generator, the length, and the underlying root
  datum permutation (or rather, the corresponding Weyl group element).
  We use a '*' for undefined cayley actions.

  NOTE: this will print reasonably on 80 columns only for groups that are
  not too large (up to rank 4 or so). We haven't tried to go over to more
  sophisticated formatting for larger groups.

  NOTE: this version outputs involutions in reduced-involution form.
*/
std::ostream& printBlockD(std::ostream& strm, const Block& block)
{
  // compute maximal width of entry
  int width = ioutils::digits(block.size()-1,10ul);
  int xwidth = ioutils::digits(block.xsize()-1,10ul);
  int ywidth = ioutils::digits(block.ysize()-1,10ul);
  int lwidth = ioutils::digits(block.length(block.size()-1),10ul);
  int cwidth = ioutils::digits(block.Cartan_class(block.size()-1),10ul);
  const int pad = 2;

  for (size_t j = 0; j < block.size(); ++j) {
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

    // print cross actions
    for (size_t s = 0; s < block.rank(); ++s) {
      BlockElt z = block.cross(s,j);
      if (z == blocks::UndefBlock)
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

    // print root datum involution as involution reduced expression
    prettyprint::printInvolution
      (strm,block.involution(j),block.twistedWeylGroup()) << std::endl;
  }

  return strm;
}


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
}


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

}

}
