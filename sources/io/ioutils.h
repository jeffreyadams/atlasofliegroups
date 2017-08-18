/*
  This is ioutils.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef IOUTILS_H  /* guard against multiple inclusions */
#define IOUTILS_H

#include <iosfwd>
#include <string>

#include "bigint.h"

namespace atlas {

namespace ioutils {

  /* The class OutputFile was moved from here to interactive.h

     The reason is that its implementation performs user input, which would
     make the rest of this module unusable for atlas, which does not want to
     link in any code performing input
  */

/******** constant declarations **********************************************/

  const size_t LineSize = 79;
  const size_t IndentSize = 4;


/******** function declarations *********************************************/

unsigned int digits(unsigned long n, unsigned int base);
 unsigned int digits(arithmetic::big_int n, unsigned int base);

std::ostream& foldLine(std::ostream&, const std::string& line,
		       const char* preHyphens = " ",
		       const char* postHyphens = "",
		       size_t h = IndentSize,
		       size_t lineSize = LineSize);

std::istream& skipSpaces(std::istream&);

}


}

#endif
