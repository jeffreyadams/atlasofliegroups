/*
  This is io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef IO_H  /* guard against multiple inclusions */
#define IO_H

#include <iosfwd>

/******** constants ********************************************************/

namespace atlas {

namespace io {

  const char* const MESSAGE_DIR = "messages/";

}

/******** function declarations ********************************************/

namespace io {

  std::ostream& printFile(std::ostream& strm, const char* a);
  std::ostream& printFile(std::ostream& strm, const char* a, const char* b);

}

}

#endif
