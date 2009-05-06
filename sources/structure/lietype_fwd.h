/*!
\file
\brief Type declarations for namespace lietype.
*/
/*
  This is lietype_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef LIETYPE_FWD_H  /* guard against multiple inclusions */
#define LIETYPE_FWD_H

#include <vector>

namespace atlas {

/******** type declarations **************************************************/

namespace lietype {

  typedef char TypeLetter;
  typedef std::vector<TypeLetter> InnerClassType;
  class SimpleLieType;
  class LieType;

}

}

#endif
