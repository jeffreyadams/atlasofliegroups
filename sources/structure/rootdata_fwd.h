/*!
\file
  
*/
/*
  This is rootdata_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef ROOTDATA_FWD_H  /* guard against multiple inclusions */
#define ROOTDATA_FWD_H

#include "bitmap_fwd.h"
#include "latticetypes_fwd.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace rootdata {

  /*!
  The roots in a root datum should almost always be referred to by
  their numbers, which refer to their location on the RootDatum list
  of roots.  This type represents such a number.
  */
  typedef unsigned long RootNbr;

  template<typename I> class RootIterator;
  class RootDatum;

  /*!
  Type for a root: an element of the lattice Z^n.
  */
  typedef latticetypes::LatticeElt Root;

  /*!
  List of _numbers_ of roots, referring to the list of roots in
  RootDatum.  This type should therefore have been called RootNbrList.
  */
  typedef std::vector<RootNbr> RootList;
  typedef RootIterator<RootList::const_iterator> WRootIterator;
  /*!
  Type for a bitmap representing a set of roots (such as the set of
  positive roots, or noncompact roots).
  */
  typedef bitmap::BitMap RootSet;
  typedef std::vector<RootSet> RootSetList;

}

}

#endif
