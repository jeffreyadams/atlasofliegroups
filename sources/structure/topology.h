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

class Connectivity {

 private:
  size_t d_rank;
  latticetypes::ComponentList d_dpi0;

 public:
// constructors and destructors
  Connectivity() {}

  Connectivity(const tori::RealTorus&, const rootdata::RootDatum&);

  ~Connectivity() {}

// accessors
  const latticetypes::ComponentList& componentReps() const {
    return d_dpi0;
  }

// manipulators
  void swap(Connectivity&);
};

}

}

#endif
