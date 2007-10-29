/*
  This is block_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include <iomanip>
#include <iostream>
#include <sstream>
#include <cassert>

#include "bitmap.h"
#include "block_io.h"
#include "bruhat.h"
#include "set.h"
#include "ioutils.h"
#include "blocks.h"
#include "descents.h"
#include "prettyprint.h"
#include "weyl.h"

/*****************************************************************************

  Output functions for |Block|, defined in sources/kl/blocks.h

******************************************************************************/

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in block_io.h


******************************************************************************/

namespace block_io {


/*
  Print the data from block to strm.

  Explanation: for each parameter, we output the length, Cartan class,
  cross-actions and Cayley or inverse Cayley transform(s) for each generator,
  and the underlying root datum permutation (or rather, the corresponding Weyl
  group element). We use a '*' for absent (inverse) Cayley transforms.

  NOTE: this will print reasonably on 80 columns only for groups that are
  not too large (up to rank 4 or so). We haven't tried to go over to more
  sophisticated formatting for larger groups.
*/
std::ostream& printBlock(std::ostream& strm, const blocks::Block& block)
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
      blocks::BlockElt z = block.cross(s,j);
      if (z == blocks::UndefBlock)
	strm << std::setw(width+pad) << '*';
      else
	strm << std::setw(width+pad) << z;
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < block.rank(); ++s)
    {
      blocks::BlockEltPair z = block.isWeakDescent(s,j)
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

    // print root datum involution
    prettyprint::printWeylElt(strm,block.involution(j),block.weylGroup());

    strm << std::endl;
  }

  return strm;
}


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
std::ostream& printBlockD(std::ostream& strm, const blocks::Block& block)
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
      blocks::BlockElt z = block.cross(s,j);
      if (z == blocks::UndefBlock)
	strm << std::setw(width+pad) << '*';
      else
	strm << std::setw(width+pad) << z;
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < block.rank(); ++s)
    {
      blocks::BlockEltPair z = block.isWeakDescent(s,j)
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
    prettyprint::printInvolution(strm,block.involution(j),block.weylGroup());

    strm << std::endl;
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
std::ostream& printBlockU(std::ostream& strm, const blocks::Block& block)
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

  for (size_t j = 0; j < block.size(); ++j) {
    for (size_t s = 0; s < block.rank(); ++s) {
      if (not block.involutionSupport(j).test(s)) // s is not in the support
	continue; // try next |s|
      if (not block.isWeakDescent(s,j)) // representation is not unitary
	goto nextj;
    }
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
      blocks::BlockElt z = block.cross(s,j);
      if (z == UndefBlock)
	strm << std::setw(width+pad) << '*';
      else
	strm << std::setw(width+pad) << z;
    }
    strm << std::setw(pad+1) << "";

    // print Cayley transforms
    for (size_t s = 0; s < block.rank(); ++s)
    {
      blocks::BlockEltPair z = block.isWeakDescent(s,j)
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
    prettyprint::printInvolution(strm,block.involution(j),block.weylGroup());

    strm << std::endl;

  nextj:
    continue;
  }

  return strm;
}


/*
  Synopsis: outputs the descent status for the various generators
*/
std::ostream& printDescent(std::ostream& strm,
			   const descents::DescentStatus& ds,
			   size_t rank, bitset::RankFlags supp)
{
  using descents::DescentStatus;

  strm << '[';

  for (size_t s = 0; s < rank; ++s) {
    if (s!=0)
      strm << ',';
    if (not supp.test(s))
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
