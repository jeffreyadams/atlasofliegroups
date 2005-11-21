/*
  This is kgb_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef KGB_FWD_H  /* guard against multiple inclusions */
#define KGB_FWD_H

#include "bitset_fwd.h"

namespace atlas {

/******** type declarations *************************************************/

namespace kgb {

  class KGB;

  typedef size_t KGBElt;
  typedef std::vector<KGBElt> KGBEltList;

  typedef std::pair<KGBElt,KGBElt> KGBEltPair;
  typedef std::vector<KGBEltPair> KGBEltPairList;

  typedef bitset::RankFlags Descent;
  typedef std::vector<Descent> DescentList;

}

}

#endif
