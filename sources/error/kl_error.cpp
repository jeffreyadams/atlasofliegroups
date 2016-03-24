/*
  This is kl_error.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "kl_error.h"

#include <iostream>
#include <cstdlib>

namespace atlas {

/*****************************************************************************

        Chapter I -- The KLError class

******************************************************************************/

namespace kl_error {

void KLError::operator() (const char* mess) const

/*
  Synopsis: outputs error information
*/

{
  std::cerr << mess << std::endl;
  std::cerr << "line #" << line << " in kl.cpp" << std::endl;
  std::cerr << "x = " << x << "; y = " << y << std::endl;
}

}

}

