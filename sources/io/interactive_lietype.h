/*
  This is interactive_lietype.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef INTERACTIVE_LIETYPE_H  /* guard against multiple inclusions */
#define INTERACTIVE_LIETYPE_H

#include <map>
#include <string>

#include "input.h"
#include "lietype.h"

/******** function declarations **********************************************/

namespace atlas {

namespace interactive_lietype {

  bool checkInnerClass(input::InputBuffer&, const lietype::LieType&, 
		       bool output=true);

  bool checkLieType(input::InputBuffer&);

  bool checkSimpleLieType(input::InputBuffer&);

  bool checkTotalRank(input::InputBuffer&);

  std::ostream& printRankMessage(std::ostream&, lietype::TypeLetter);

  void readInnerClass(lietype::InnerClassType&, input::InputBuffer&,
		      const lietype::LieType&);

  void readLieType(lietype::LieType&, input::InputBuffer&);

}

}

#endif
