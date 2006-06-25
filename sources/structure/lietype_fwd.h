/*!
\file
\brief Type declarations for namespace lietype.
*/
/*
  This is lietype_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef LIETYPE_FWD_H  /* guard against multiple inclusions */
#define LIETYPE_FWD_H

#include <vector>

namespace atlas {

/******** type declarations **************************************************/

namespace lietype {

  typedef char TypeLetter;
  typedef std::pair<TypeLetter,size_t> SimpleLieType;
  typedef std::vector<SimpleLieType> LieType;
  typedef std::vector<TypeLetter> InnerClassType;

}

}

#endif
