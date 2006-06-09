/*!
\file
  This is realredgp.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef REALREDGP_H  /* guard against multiple inclusions */
#define REALREDGP_H

#include "bitmap_fwd.h"
#include "cartanclass_fwd.h"
#include "complexredgp_fwd.h"
#include "poset_fwd.h"
#include "realredgp_fwd.h"
#include "rootdata_fwd.h"
#include "weyl_fwd.h"
#include "tits_fwd.h"

#include "realform.h"
#include "topology.h"

/******** type definitions **************************************************/

namespace atlas {

namespace realredgp {

class RealReductiveGroup {
/*!
 We do not own the complex group; a RealReductiveGroup should be seen
 as a "pointer" to a complex group; when the complex group changes,
 the corresponding pointers are invalidated.
*/
 private:

  enum StatusFlagNames { FullCartan, IsConnected, IsQuasisplit, IsSemisimple, 
			 IsSplit, NumStatusFlags };

  typedef bitset::BitSet<NumStatusFlags> Status;

  // we do not own the complex group; a RealReductiveGroup should be seen
  // as a "pointer" to a complex group; when the complex group changes,
  // the corresponding pointers are invalidated
  complexredgp::ComplexReductiveGroup* d_complexGroup;

  realform::RealForm d_realForm;
  topology::Connectivity d_connectivity;
  Status d_status;

 public:

// constructors and destructors
  RealReductiveGroup() {}

  RealReductiveGroup(complexredgp::ComplexReductiveGroup&, realform::RealForm);

  ~RealReductiveGroup();

// accessors
  const cartanclass::CartanClass& cartan(size_t) const;

  const poset::Poset& cartanOrdering() const;

  const bitmap::BitMap& cartanSet() const;

  const complexredgp::ComplexReductiveGroup& complexGroup() const {
    return *d_complexGroup;
  }

  const latticetypes::ComponentList& componentReps() const {
    return d_connectivity.componentReps();
  }

  const latticetypes::LatticeMatrix& distinguished() const;

  void grading(rootdata::RootSet&) const;

  bool isConnected() const {
    return d_status[IsConnected];
  }

  bool isQuasisplit() const {
    return d_status[IsQuasisplit];
  }

  bool isSemisimple() const {
    return d_status[IsSemisimple];
  }

  bool isSplit() const {
    return d_status[IsSplit];
  }

  size_t kgbSize() const;

  size_t mostSplit() const;

  size_t numCartan() const;

  size_t rank() const;

  realform::RealForm realForm() const {
    return d_realForm;
  }

  const rootdata::RootDatum& rootDatum() const;

  size_t semisimpleRank() const;

  const tits::TitsGroup& titsGroup() const;

  const weyl::WeylGroup& weylGroup() const;

// manipulators
  complexredgp::ComplexReductiveGroup& complexGroup() {
    return *d_complexGroup;
  }

  void fillCartan();

  void swap(RealReductiveGroup&);
};

}

}

#endif
