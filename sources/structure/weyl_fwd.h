/*
  This is weyl_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef WEYL_FWD_H  /* guard against multiple inclusions */
#define WEYL_FWD_H

#include <vector>

#include "constants.h"
#include "stlvector.h"
#include "typenumber.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace weyl {

  class WeylElt;
  class WeylGroup;

#if 0
  typedef stlvector::Vector<WeylElt>::type WeylEltList;
#endif
  typedef std::vector<WeylElt> WeylEltList;

  typedef unsigned char Generator;
  typedef unsigned char EltPiece;

  typedef std::vector<Generator> WeylWord;
  typedef Generator Twist[constants::RANK_MAX];

}

namespace typenumber {

  template<> struct TypeNumber<weyl::WeylElt> {
    static size_t number() {
      return WeylElt;
    }
  };

}

}

#endif
