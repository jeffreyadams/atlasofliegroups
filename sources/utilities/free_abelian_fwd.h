/*
  This is free_abelian_fwd.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef FREE_ABELIAN_FWD_H  /* guard against multiple inclusions */
#define FREE_ABELIAN_FWD_H

#include <functional>

namespace atlas {

namespace free_abelian {

template<typename T, typename C=long int, typename Compare=std::less<T> >
  struct Free_Abelian;
template<typename T, typename C=long int, typename Compare=std::less<T> >
  struct Monoid_Ring;
template<typename T, typename C=long int, typename Compare=std::less<T> >
  class Free_Abelian_light;
}

}

#endif
