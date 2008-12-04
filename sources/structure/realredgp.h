/*!
\file
\brief Class definitions and function declarations for RealReductiveGroup.
*/
/*
  This is realredgp.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
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

/*! \brief Represents a real form on a connected reductive complex group,
 determining a real reductive group

 An object of this class is determined by a ComplexReductiveGroup and the
 number of a real form; in addition it stores some data concerning the group
 of real points of the real form

 The complex group is referred to by a pointer that is not owned; we are
 dependent on the owner of the complex group, and once it is destructed, the
 RealReductiveGroup objects referring to it become invalid.
*/
class RealReductiveGroup {

 private:

  enum StatusFlagNames { IsConnected, IsQuasisplit, IsSemisimple, IsSplit,
			 NumStatusFlags };

  typedef bitset::BitSet<NumStatusFlags> Status;

  // we do not own the complex group; a RealReductiveGroup should be seen
  // as dependent on a complex group; when the complex group changes,
  // the dependent RealReductiveGroup objects are invalidated
  complexredgp::ComplexReductiveGroup* d_complexGroup;

  realform::RealForm d_realForm;
  topology::Connectivity d_connectivity;
  Status d_status;

 public:

// constructors and destructors
  RealReductiveGroup() {}

  RealReductiveGroup(complexredgp::ComplexReductiveGroup&, realform::RealForm);

// accessors
  const cartanclass::CartanClass& cartan(size_t) const;

  const bitmap::BitMap& Cartan_set() const;

  const complexredgp::ComplexReductiveGroup& complexGroup() const {
    return *d_complexGroup;
  }

  size_t numInvolutions() const;

  const latticetypes::SmallBitVectorList& dualComponentReps() const {
    return d_connectivity.dualComponentReps();
  }

  const latticetypes::LatticeMatrix& distinguished() const;

  rootdata::RootSet noncompactRoots() const;

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

  size_t KGB_size() const;

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
