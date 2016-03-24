/*
  This is block_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef BLOCK_IO_H  /* guard against multiple inclusions */
#define BLOCK_IO_H

#include <iosfwd>

#include "../Atlas.h"

#include "bitset.h"	// for default argument to |printDescent|
#include "ext_block.h"  // for |ext_block::DescValue|

namespace atlas {

namespace block_io {


/******** function declarations *********************************************/

// the main output interface is via methods of |blocks::Block_base|

  std::ostream& printBlockU(std::ostream&, const Block&);

  std::ostream& printDescent(std::ostream&, const DescentStatus& ds,
			     size_t rank, RankFlags mask = RankFlags(~0ul));

  std::ostream& print_twist(std::ostream& strm, const Block_base& block);

  std::ostream& print_KL(std::ostream&f, param_block& block, BlockElt z);

} // |namespace block_io|

namespace ext_block {

  const char* descent_code(DescValue v);

} // |namespace ext_block|

} // |namespace atlas|


#endif
