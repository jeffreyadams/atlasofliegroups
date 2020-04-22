/*
  This is realredgp.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Class definitions and function declarations for RealReductiveGroup.

#ifndef REALREDGP_H  /* guard against multiple inclusions */
#define REALREDGP_H

#include "bitmap_fwd.h"
#include "poset_fwd.h"
#include "../Atlas.h"

#include "ratvec.h"     // containment of |RatCoweight| field
#include "bitvector.h"  // containment of |SmallBitVectorList|


/******** type definitions **************************************************/

namespace atlas {

namespace realredgp {

/* Represent a real form on a connected reductive complex group,
   determining a real reductive group

 An object of this class is determined by a InnerClass and the
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

  typedef BitSet<NumStatusFlags> Status;

/*
  we do not own the inner class; a |RealReductiveGroup| should be seen as
  dependent on an inner class. In Fokko, when the complex group (and with it the
  |InnerClass|) changes, the dependent |RealReductiveGroups| are invalidated
*/
  InnerClass& d_innerClass;

  RealFormNbr d_realForm; // our identification number
  SmallBitVectorList dual_pi0_gens;

  RatCoweight square_class_cocharacter; // a base coweight for square class
  TorusPart torus_part_x0; // initial |TorusPart| relative to square class base

  const TitsCoset* d_Tg; // owned pointer; the group is stored here
  KGB* kgb_ptr; // owned pointer, but initially |NULL|

  Status d_status;

 public:

// constructors and destructors
  RealReductiveGroup(InnerClass&, RealFormNbr);
  RealReductiveGroup(InnerClass&, RealFormNbr,
		     const RatCoweight& coch, TorusPart x0_torus_part);
  ~RealReductiveGroup(); // not inline: type incomplete; deletes pointers

// accessors
  const InnerClass& innerClass() const { return d_innerClass; }
  // following method forces |const| result, compare with |cbegin| methods
  RealFormNbr realForm() const { return d_realForm; }
  const RootDatum& root_datum() const;
  const TitsCoset& basedTitsGroup() const { return *d_Tg; }
  const TitsGroup& titsGroup() const;
  const WeylGroup& weylGroup() const;
  const TwistedWeylGroup& twistedWeylGroup() const;
  BitMap Cartan_set() const;
  const CartanClass& cartan(size_t cn) const; // Cartan number of parent

  TorusPart x0_torus_part() const { return torus_part_x0; }
  RatCoweight g() const; // |square_class_cocharacter| + $\check\rho$
  const RatCoweight& g_rho_check() const // that is |g()| minus $\check\rho$:
    { return square_class_cocharacter; }
  Grading base_grading() const; // grading (1=noncompact) at square class base

  bool isConnected() const { return d_status[IsConnected]; }

  bool isCompact() const { return d_status[IsCompact]; }
  bool isQuasisplit() const { return d_status[IsQuasisplit]; }
  bool isSplit() const { return d_status[IsSplit]; }

  bool isSemisimple() const;
  size_t numCartan() const;
  size_t rank() const;
  size_t semisimpleRank() const;
  size_t numInvolutions(); // non-|const|, as Cartan classes get generated
  size_t KGB_size() const; // the cardinality of |K\\G/B|.
  size_t mostSplit() const;

/*
  Return the grading offset (on simple roots) adapted to |G|. This flags among
  the simple roots those that are noncompact imaginary at the initial KGB
  element |x0| of G (which lives on the fundamental Cartan).
*/
  Grading grading_offset();

  const size_t component_rank() const;
  const SmallBitVectorList& dualComponentReps() const;
  const WeightInvolution& distinguished() const;

/* Return the set of noncompact imaginary roots for (the representative of)
   the real form.
*/
  RootNbrSet noncompactRoots() const;

// manipulators
  void swap(RealReductiveGroup&);

  InnerClass& innerClass() { return d_innerClass; }

  const KGB& kgb();
  const BruhatOrder& Bruhat_KGB();

// internal methods
 private:
  void construct();  // work common to two constructors
}; // |class RealReductiveGroup|


// Equality test. Requires identical inner class, same base grading and initial KGB
// elements. This should imply the |realForm| numbers of the two are equal.

inline bool operator==
  (const RealReductiveGroup& x, const RealReductiveGroup& y)
{
  return &x.innerClass() == &y.innerClass()
    and x.g_rho_check() == y.g_rho_check()
    and x.x0_torus_part() == y.x0_torus_part()
    and (assert(x.realForm() == y.realForm()),true);
}

inline bool operator!=
  (const RealReductiveGroup& x, const RealReductiveGroup& y)
{ return not(x==y); }

//			   function declarations

TorusPart minimal_torus_part
  (const InnerClass& G, RealFormNbr wrf, const RatCoweight& coch,
   const TwistedInvolution& tw, const RatCoweight& torus_factor
   );

} // |namespace realredgp|

} // |namespace atlas|

#endif
