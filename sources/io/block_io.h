/*
  This is block_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BLOCK_IO_H  /* guard against multiple inclusions */
#define BLOCK_IO_H

#include <iosfwd>

#include "atlas_types.h"

#include "bitset.h"

namespace atlas {

/******** function declarations *********************************************/

namespace block_io {

  std::ostream& print_block(std::ostream&, const Block_base&);

  std::ostream& printBlockD(std::ostream&, const Block&);

  std::ostream& printBlockU(std::ostream&, const Block&);

  std::ostream& printDescent(std::ostream&, const DescentStatus&,
		     size_t, RankFlags = RankFlags(~0ul));

}

}

#endif
