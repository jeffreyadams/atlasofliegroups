/*
  This is version.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef VERSION_H  /* guard against multiple inclusions */
#define VERSION_H

namespace atlas {

namespace version {

const char* const NAME = "the Atlas of Lie Groups and Representations";
const char* const VERSION = "1.0.8"; // last advanced April 1st, 2020
// version 1.0.6 completes support for cyclotomic polynomials (for Weyl c.f.)
// version 1.0.7 (workshop) has better W_cells support, and fixes some bugs
// version 1.0.8 has sprawling library, to track unipotent representations
const char* const COMPILEDATE = __DATE__ " at " __TIME__;

}

}

#endif
