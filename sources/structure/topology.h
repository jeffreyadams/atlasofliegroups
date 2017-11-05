/*
  This is topology.h
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef TOPOLOGY_H  /* guard against multiple inclusions */
#define TOPOLOGY_H

#include "../Atlas.h"

namespace atlas {
namespace topology {

// function definition

SmallBitVectorList
  dual_component_group_basis(const WeightInvolution theta,const RootDatum& rd);

} // |namespace topology|

} // |namespace atlas|

#endif
