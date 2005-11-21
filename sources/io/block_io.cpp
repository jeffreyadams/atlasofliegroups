/*
  This is block_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include <iomanip>
#include <iostream>
#include <sstream>

#include "block_io.h"

#include "ioutils.h"
#include "blocks.h"
#include "descents.h"
#include "prettyprint.h"
#include "weyl.h"

/*****************************************************************************

  Input/output functions for the Block data structure, defined in 
  sources/kl/blocks.h

******************************************************************************/

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in block_io.h

  ... explain here when it is stable ...

******************************************************************************/

namespace block_io {

std::ostream& printBlock(std::ostream& strm, const blocks::Block& block)

/*
  Synopsis: outputs the data from block to strm.

  Explanation: for each parameter, we output the cross-actions and 
  cayley-actions for each generator, the length, and the underlying root
  datum permutation (or rather, the corresponding Weyl group element).
  We use a '*' for undefined cayley actions.

  NOTE: this will print reasonably on 80 columns only for groups that are
  not too large (up to rank 4 or so). We haven't tried to go over to more
  sophisticated formatting for larger groups.
*/

{
  using namespace blocks;
  using namespace kgb;
  using namespace prettyprint;

  // compute maximal width of entry
  int width = ioutils::digits(block.size()-1,10ul);
  int xwidth = ioutils::digits(block.xsize()-1,10ul);
  int ywidth = ioutils::digits(block.ysize()-1,10ul);
  int lwidth = ioutils::digits(block.length(block.size()-1),10ul);
  const int pad = 2;

  for (size_t j = 0; j < block.size(); ++j) {
    // print entry number and corresponding orbit pair
    strm << std::setw(width) << j;
    strm << "(" << std::setw(xwidth) << block.x(j);
    strm << ",";
    strm << std::setw(ywidth) << block.y(j) << ")";
    strm << ":";

    // print cross actions
    for (size_t s = 0; s < block.rank(); ++s) {
      BlockElt z = block.cross(s,j);
      if (z == UndefBlock)
	strm << std::setw(width+pad) << '*';
      else
	strm << std::setw(width+pad) << z;
    }
    strm << std::setw(2*pad) << "";

    // print cayley transforms
    for (size_t s = 0; s < block.rank(); ++s) {
      BlockEltPair z = block.cayley(s,j);
      strm << "(";
      if (z.first == UndefBlock)
	strm << std::setw(width) << '*';
      else
	strm << std::setw(width) << z.first;
      strm << ",";
      if (z.second == UndefBlock)
	strm << std::setw(width) << '*';
      else
	strm << std::setw(width) << z.second;
      strm << ")" << std::setw(pad) << "";
    }
    strm << std::setw(pad) << "";

    // print descents
    printDescent(strm,block.descent(j),block.rank());

    // print length
    strm << std::setw(lwidth+pad) << block.length(j);
    strm << "  ";

    // print root datum involution
    printWeylElt(strm,block.tau(j),block.weylGroup());

    strm << std::endl;
  }

  return strm;
}

std::ostream& printDescent(std::ostream& strm, 
			   const descents::DescentStatus& ds, size_t rank)

/*
  Synopsis: outputs the descent status for the various generators
*/

{
  using namespace descents;

  strm << "[";

  for (size_t s = 0; s < rank; ++s) {
    if (s)
      strm << ",";
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
      strm << "*";
      break;
    }
  }

  strm << "]";

  return strm;
}

}

}
