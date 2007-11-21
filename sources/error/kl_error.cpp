/*
  This is kl_error.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "kl_error.h"

#include <iostream>

#include "kl.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- The KLError class

  ... explain here when it is stable ...

******************************************************************************/

namespace kl_error {

void KLError::operator() (const char* mess) const

/*
  Synopsis: outputs error information and exits.
*/

{
  std::cerr << mess << std::endl;
  std::cerr << "line #" << line << " in kl.cpp" << std::endl;
  std::cerr << "x = " << x << "; y = " << y << std::endl;

  exit(0);
}

}

}

