/*
  This is rootdata_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef ROOTDATA_FWD_H  /* guard against multiple inclusions */
#define ROOTDATA_FWD_H

#include "bitmap_fwd.h"
#include "latticetypes_fwd.h"

#include "stlvector.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace rootdata {

  typedef unsigned long RootNbr;

  class RootIterator;
  class RootDatum;

  typedef latticetypes::LatticeElt Root;
#if 0
  typedef stlvector::Vector<RootNbr>::type RootList;
#endif
  typedef std::vector<RootNbr> RootList;

  typedef bitmap::BitMap RootSet;
  typedef std::vector<RootSet> RootSetList;

}

}

#endif
