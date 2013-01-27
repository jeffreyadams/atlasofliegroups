/*
  This is block_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef BLOCK_IO_H  /* guard against multiple inclusions */
#define BLOCK_IO_H

#include <iosfwd>

#include "atlas_types.h"

#include "bitset.h"	// for default argument to |printDescent|

namespace atlas {

/******** function declarations *********************************************/

// the main output interface is via methods of |blocks::Block_base|

namespace block_io {

  std::ostream& printBlockU(std::ostream&, const Block&);

  std::ostream& printDescent(std::ostream&, const DescentStatus& ds,
			     size_t rank, RankFlags mask = RankFlags(~0ul));

  std::ostream& print_twist(std::ostream& strm, const Block_base& block);

  std::ostream& print_KL(std::ostream&f, param_block& block, BlockElt z);

}

}

#endif
