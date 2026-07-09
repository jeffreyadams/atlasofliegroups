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

  // note that the axis language version is defined in interpreter/main.w

const char* const NAME = "the Atlas of Lie Groups and Representations";
const char* const VERSION = "1.2"; // last advanced December 1, 2024
// version 1.0.6 completes support for cyclotomic polynomials (for Weyl c.f.)
// version 1.0.7 (workshop) has better W_cells support, and fixes some bugs
// version 1.0.8 has sprawling library, to track unipotent representations
// version 1.0.9 has many improvements to deformation computation, alcoves
// version 1.1   was declared (on master) on December 30, 2021
// version 1.1.1 has almost two extra years of fixing; e.g. integer conversions
// version 1.2 has almost two extra years of changes of all kinds of nature
const char* const COMPILEDATE = __DATE__ " at " __TIME__;
constexpr bool debugging =
#ifdef NDEBUG
  false
#else
  true
#endif
  , with_readline =
#ifdef NREADLINE
  false
#else
  true
#endif
  ;
}

}

#endif
