/*
  This is io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "io.h"

#include <fstream>
#include <string>

/****************************************************************************

  This file contains some very simple auxiliary i/o functions.

*****************************************************************************/

namespace atlas {

namespace io {

std::ostream& printFile(std::ostream& strm, const char* a)

/*
  Prints the contents of the file with name a to strm.
*/

{  
  std::ifstream file(a);
  if (file.good())
    strm << file.rdbuf();
  file.close();

  return strm;
}

std::ostream& printFile(std::ostream& strm, const char* a, const char* b)

/*
  Prints the contents of the file with name b+a to strm (b is presumably
  a directory name.)
*/

{  
  std::string aStr = std::string(a);
  std::string bStr = std::string(b);
  std::string str = bStr + aStr;

  return printFile(strm,str.c_str());
}

}

}
