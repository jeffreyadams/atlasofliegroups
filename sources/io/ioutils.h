/*
  This is ioutils.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef IOUTILS_H  /* guard against multiple inclusions */
#define IOUTILS_H

#include <iosfwd>
#include <string>

namespace atlas {

  /* The class OutputFile was moved from here to interactive.h

     The reason is that its implementation performs user input, which would
     make the rest of this module unusable for realex, which does not want to
     link in any code performing input
  */

/******** constant declarations **********************************************/

namespace ioutils {

  const size_t LineSize = 79;
  const size_t IndentSize = 4;

}

/******** function declarations *********************************************/

namespace ioutils {

unsigned long digits(unsigned long, unsigned long);

std::ostream& foldLine(std::ostream&, const std::string&,
		       const char* preHyphens = " ",
		       const char* postHyphens = "",
		       size_t h = IndentSize,
		       size_t lineSize = LineSize);

std::istream& skipSpaces(std::istream&);

}


}

#endif
