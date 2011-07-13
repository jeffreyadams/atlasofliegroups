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
#include "gradings_fwd.h"
#include "poset_fwd.h"
#include "realredgp_fwd.h"

#include "rootdata.h"
#include "cartanclass.h"
#include "weyl.h"
#include "tits.h"
#include "complexredgp.h"
#include "realform.h"
#include "topology.h"
#include "kgb_fwd.h"
#include "bruhat_fwd.h"

/******** type definitions **************************************************/

namespace atlas {

namespace realredgp {

/*! \brief Represents a real form on a connected reductive complex group,
 determining a real reductive group

 An object of this class is determined by a ComplexReductiveGroup and the
 number of a real form; in addition it stores some data concerning the group
 of real points of the real form

 The complex group is referred to by a reference that is not owned; we are
 dependent on the owner of the complex group, and once it is destructed, the
 RealReductiveGroup objects referring to it become invalid.
*/
class RealReductiveGroup
{
  enum StatusFlagNames
  { IsConnected, IsCompact,IsQuasisplit,IsSplit, IsSemisimple,
    NumStatusFlags };

  typedef bitset::BitSet<NumStatusFlags> Status;

  // we do not own the complex group; a RealReductiveGroup should be seen
  // as dependent on a complex group; when the complex group changes,
  // the dependent RealReductiveGroup objects are invalidated
  complexredgp::ComplexReductiveGroup& d_complexGroup;

  realform::RealForm d_realForm; // our identification number
  topology::Connectivity d_connectivity; // characterss of componentn group

  const tits::TitsCoset* d_Tg; // owned pointer; the group is stored here
  kgb::KGB* kgb_ptr; // owned pointer, but initially |NULL|

  Status d_status;

 public:

// constructors and destructors
  RealReductiveGroup(complexredgp::ComplexReductiveGroup&, realform::RealForm);
  ~RealReductiveGroup(); // not inline: type incomplete; deletes d_Tg, kgb_ptr

// accessors
  const complexredgp::ComplexReductiveGroup& complexGroup() const
    { return d_complexGroup; }

  realform::RealForm realForm() const { return d_realForm; }

  const rootdata::RootDatum& rootDatum() const
    { return d_complexGroup.rootDatum(); }

  const tits::TitsCoset& basedTitsGroup() const { return *d_Tg; }
  const tits::TitsGroup& titsGroup() const  { return d_Tg->titsGroup(); }

  const weyl::WeylGroup& weylGroup() const
    { return d_complexGroup.weylGroup(); }

  const weyl::TwistedWeylGroup& twistedWeylGroup() const
    { return d_complexGroup.twistedWeylGroup(); }

  bitmap::BitMap Cartan_set() const
    { return complexGroup().Cartan_set(d_realForm); }

/*!\brief Returns cartan \#cn in the group.

  Precondition: cn belongs to cartanSet().
*/
  const cartanclass::CartanClass& cartan(size_t cn) const
    { return d_complexGroup.cartan(cn); }

  bool isConnected() const { return d_status[IsConnected]; }

  bool isCompact() const { return d_status[IsCompact]; }
  bool isQuasisplit() const { return d_status[IsQuasisplit]; }
  bool isSplit() const { return d_status[IsSplit]; }

  bool isSemisimple() const { return d_status[IsSemisimple]; }

  size_t numCartan() const { return Cartan_set().size(); }

  size_t rank() const { return rootDatum().rank(); };

  size_t semisimpleRank() const { return rootDatum().semisimpleRank(); }

  size_t numInvolutions()
    { return complexGroup().numInvolutions(Cartan_set()); }

/*!\brief Returns the cardinality of K\\G/B.

  Precondition: fillCartan() has been called.
*/
  size_t KGB_size() const { return d_complexGroup.KGB_size(d_realForm); }

  size_t mostSplit() const { return d_complexGroup.mostSplit(d_realForm); }

/*! \brief
  Returns the grading offset (on simple roots) adapted to |G|. This flags the
  simple roots that are noncompact imaginary at the fundamental Cartan in G.

Algorithm: the variable |rset| is first made to flag, among the imaginary
roots of the fundamental Cartan, those that are noncompact for the chosen
representative (in the adjoint fiber) of the real form of |G|. The result is
formed by extracting only the information concerning the presence of the
\emph{simple} roots in |rset|.
*/
  gradings::Grading grading_offset()
    {
      rootdata::RootSet rset= noncompactRoots(); // grading for real form rep
      return cartanclass::restrictGrading(rset,rootDatum().simpleRootList());
    }

  cartanclass::square_class square_class() const
    { return d_complexGroup.fundamental().central_square_class(d_realForm); }

  const size_t component_rank() const
    { return d_connectivity.component_rank(); }
  const bitvector::SmallBitVectorList& dualComponentReps() const
    { return d_connectivity.dualComponentReps(); }

  const latticetypes::LatticeMatrix& distinguished() const
    { return d_complexGroup.distinguished(); }


/*!\brief Returns the set of noncompact imaginary roots for (the
  representative of) the real form.
*/
  rootdata::RootSet noncompactRoots() const
    { return d_complexGroup.noncompactRoots(d_realForm); }

// manipulators
  void swap(RealReductiveGroup&);

  complexredgp::ComplexReductiveGroup& complexGroup()
    { return d_complexGroup; }

  const kgb::KGB& kgb();
  const bruhat::BruhatOrder Bruhat_KGB();

}; // |class RealReductiveGroup|

} // |namespace realredgp|

} // |namespace atlas|

#endif
