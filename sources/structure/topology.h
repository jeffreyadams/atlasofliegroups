/*!
\file
  This is topology.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef TOPOLOGY_H  /* guard against multiple inclusions */
#define TOPOLOGY_H

#include <vector>

#include "rootdata_fwd.h"
#include "tori_fwd.h"

#include "latticetypes.h"

/******** type declarations *************************************************/

namespace atlas {

namespace topology {

  class Connectivity;

}

/******** function declarations *********************************************/

namespace topology {

  bool isTrivial(const latticetypes::CoeffList&);

}

/******** type definitions **************************************************/

namespace topology {

/* class of which the constructor computes a basis for the dual of the
   component group of a real reductive group (given its most split torus, and
   the root datum that it is related to).
 */

class Connectivity {

 private:
  latticetypes::SmallBitVectorList d_dpi0; // basis of dual component group

 public:
// constructors and destructors
  Connectivity() {}

  Connectivity(const tori::RealTorus& most_split, const rootdata::RootDatum&);

  ~Connectivity() {}

// accessors
  const latticetypes::SmallBitVectorList& dualComponentReps() const {
    return d_dpi0;
  }

// manipulators
  void swap(Connectivity& other) { d_dpi0.swap(other.d_dpi0);}
};

}

}

#endif
