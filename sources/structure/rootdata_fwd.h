/*!
\file
\brief Class declarations and type definitions for RootDatum.   
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
\brief Integer representing the location of a root on the list
maintained by a RootDatum.

  The roots in a root datum should almost always be referred to by
  their numbers.
  */
  typedef unsigned long RootNbr;

  template<typename I> class RootIterator;
  class RootDatum;

  /*!
 \brief Type for a root, coroot, character, or cocharacter: an element
 of the lattice Z^n. 
  */
  typedef latticetypes::LatticeElt Root;

  /*!
\brief List of _numbers_ of roots, referring to the list of roots in
  RootDatum.  

According to Fokko, this type should therefore have been called RootNbrList.
  */
  typedef std::vector<RootNbr> RootList;
  typedef RootIterator<RootList::const_iterator> WRootIterator;
  /*!
\brief BitMap representing a set of roots (such as the set of
  positive roots, or noncompact roots).
  */
  typedef bitmap::BitMap RootSet;
  typedef std::vector<RootSet> RootSetList;

}

}

#endif
