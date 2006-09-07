/*
  This is typestring.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "typestring.h"

#include <sstream>

#include "typenumber.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in typestring.h

******************************************************************************/

namespace typestring {

const char* name(size_t t)

/*
  Synopsis: outputs the name of the type number t.

  If t < Unknown, this is an explicitly named type; otherwise, it's name is
  of the form "unknown type # such and such".
*/

{
  using namespace typenumber;

  static std::string str;

  switch (t) {
  case Int:
    str = "int";
    break;
  case Ulong:
    str = "unsigned long";
    break;
  case WeylElt:
    str = "WeylElt";
    break;
  default: // unknown type
    std::ostringstream os;
    os << "unknown type #" << t-Unknown;
    str = os.str();
    break;
  }

  return str.c_str();
}

}

}
