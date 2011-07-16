/*!
\file
  This is topology.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef TOPOLOGY_H  /* guard against multiple inclusions */
#define TOPOLOGY_H

#include <vector>

#include "atlas_types.h"

#include "bitvector.h" // contained in |Connectivity|

/******** type declarations *************************************************/

namespace atlas {


/******** function declarations *********************************************/

namespace topology {

  bool isTrivial(const CoeffList&);

}

/******** type definitions **************************************************/

namespace topology {

/* class of which the constructor computes a basis for the dual of the
   component group of a real reductive group (given its most split torus, and
   the root datum that it is related to).
 */

class Connectivity
{

 private:
  SmallBitVectorList d_dpi0; // basis of dual component group

 public:
// constructors and destructors
  Connectivity() {}

  Connectivity(const tori::RealTorus& most_split, const RootDatum&);

  ~Connectivity() {}

// accessors
  const SmallBitVectorList& dualComponentReps() const
  { return d_dpi0; }

  const size_t component_rank() const { return d_dpi0.size(); }

// manipulators
  void swap(Connectivity& other) { d_dpi0.swap(other.d_dpi0);}

}; // |class Connectivity|

} // |namespace topology|

} // |namespace atlas|

#endif
