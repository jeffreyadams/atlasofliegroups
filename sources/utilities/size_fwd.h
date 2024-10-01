/*
  This is size_fwd.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef SIZE_FWD_H  /* guard against multiple inclusions */
#define SIZE_FWD_H

/******** forward type declarations *****************************************/

namespace atlas {

namespace size {

  template<typename C> class SizeType;
  using BaseType = signed char; // safe for RankMax<=128 (for C128 factor 2^255)
  using Size = SizeType<BaseType>;

} // |namespace size|

} // |namespace atlas|

#endif
