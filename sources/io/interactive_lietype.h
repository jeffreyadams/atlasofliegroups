/*
  This is interactive_lietype.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For license information see the LICENSE file
*/

#ifndef INTERACTIVE_LIETYPE_H  /* guard against multiple inclusions */
#define INTERACTIVE_LIETYPE_H

#include <map>
#include <string>

#include "atlas_types.h"
#include "lietype.h"	// |lietype::TypeLetter|

/******** function declarations **********************************************/

namespace atlas {

namespace interactive_lietype {

  bool checkInnerClass(input::InputBuffer&, const LieType&, 
		       bool output=true);

  bool checkLieType(input::InputBuffer&);

  bool checkSimpleLieType(input::InputBuffer&);

  bool checkTotalRank(input::InputBuffer&);

  std::ostream& printRankMessage(std::ostream&, lietype::TypeLetter);

  void readInnerClass(InnerClassType&, input::InputBuffer&,
		      const LieType&);

  void readLieType(LieType&, input::InputBuffer&);

}

}

#endif
