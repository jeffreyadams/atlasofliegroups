/*!
\file
\brief Type definitions for the class BitSet.
*/
/*
  This is bitset_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BITSET_FWD_H  /* guard against multiple inclusions */
#define BITSET_FWD_H

#include <vector>

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace bitset {

  template<size_t n> class BitSet;

  typedef BitSet<constants::RANK_MAX> RankFlags;
  typedef BitSet<2*constants::RANK_MAX> TwoRankFlags;
  typedef std::vector<RankFlags> RankFlagsList;

}

}

#endif
