/*!
\file
  This is size_fwd.h
*/
/*
  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SIZE_FWD_H  /* guard against multiple inclusions */
#define SIZE_FWD_H

/******** forward type declarations *****************************************/

namespace atlas {

namespace size {

  template<typename C> class SizeType;
  typedef signed char BaseType; // safe for RankMax<=128 (for C128 factor 2^255)
  typedef unsigned char UnsignedBaseType;
  typedef SizeType<BaseType> Size;

}

}

#endif
