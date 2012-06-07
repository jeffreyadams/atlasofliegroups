/*
  This is io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

#include "io.h"

#include <fstream>
#include <string>

/* the following macro should be set to an absolute path from the atlas
   makefile, but if not we use a relative path from the atlas directory */
#ifndef MESSAGE_DIR_MACRO
#define MESSAGE_DIR_MACRO "messages/"
#endif

/****************************************************************************

  This file contains some very simple auxiliary i/o functions.

*****************************************************************************/

namespace atlas {

namespace io {

const char* MESSAGE_DIR = MESSAGE_DIR_MACRO;

/*
  Prints the contents of the file with name a to strm.
*/
std::ostream& printFile(std::ostream& strm, const char* a)
{
  std::ifstream file(a);
  if (file.good())
    strm << file.rdbuf(); // copy entire file contents to |strm|
  else
    strm << "No help file found at " << a ;
  file.close();

  return strm<< std::endl; // add a newline
}


/*
  Prints the contents of the file with name b+a to strm (b is presumably
  a directory name.)
*/
std::ostream& printFile(std::ostream& strm, const char* a, const char* b)
{
  std::string aStr = std::string(a);
  std::string bStr = std::string(b);
  std::string str = bStr + aStr;

  return printFile(strm,str.c_str());
}

}

}
