/*!
\file\brief Class definition for the Subsystem class.
*/
/*
  This is subystem.h

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef SUBSYSTEM_H  /* guard against multiple inclusions */
#define SUBSYSTEM_H

#include "atlas_types.h"

#include "rootdata.h"	// derive from |RootSystem|
#include "weyl.h"	// containment

namespace atlas {

namespace subsystem {

/* The following data type has a purpose specific for representation theory.
   It is associated to a subsystem of the dual root datum (with positive roots
   contained in the positive parent coroots), and is derived from |RootSystem|
   at that dual side (meaning that inherited |rootSystem| methods present that
   subsystem of the dual root system of the parent). It remains however
   attached to the parent root _datum_ (not system), containing and exporting
   a reference to that rootdatum, and has a method to produce a |PreRootDatum|
   for the subsystem of the _parent_ defined by simple coroots for the
   subsystem (a full such root datum is not stored, for efficientcy reasons).
 */

// A subsystem on the dual side of a given root datum
class SubSystem : public RootSystem // new system, subsytem of dual
{
  const RootDatum& rd; // parent root datum
  RootNbrList pos_map; // map positive roots to root number in parent
  RootNbrList inv_map; // partial map back from all parent roots

  struct root_info
  { weyl::Generator simple;
    WeylWord to_simple;
    WeylWord reflection;

  root_info() : simple(~0), to_simple(), reflection() {}
  };
  std::vector<root_info> sub_root;

 public:
  SubSystem(const RootDatum& parent,
	    const RootNbrList& sub_sys // list of simple roots in subsys
           );

  static SubSystem integral // pseudo contructor for integral system
  (const RootDatum& parent, const RatWeight& gamma);

  SubSystem(const SubSystem& s) // copy contructor, not actually used
  : RootSystem(s) // copy base object
  , rd(s.rd) // share this one
  , pos_map(s.pos_map), inv_map(s.inv_map), sub_root(s.sub_root) // copy those
  {
    assert(false); // should never be actually called, but exist nonetheless
  }

  const RootDatum& parent_datum() const { return rd; }

  weyl::Twist twist(const WeightInvolution& theta, WeylWord& ww) const;
  // output value |ww| gives |-^theta| as twisted involution for |sub|

  weyl::Twist parent_twist(const WeightInvolution& theta, WeylWord& ww) const;
  // similar, but |ww| is (for subsystem) on parent side, and anchored at its
  // distinguished involution |delta|: one has |theta=parent(ww).delta|.

  PreRootDatum pre_root_datum() const; // viewed from parent side

  RootNbr parent_nr_simple(weyl::Generator s) const
  { return pos_map[s]; }

  RootNbr to_parent(RootNbr alpha) const; // |pos_map| with some shifting
  RootNbr from_parent(RootNbr alpha) const { return inv_map[alpha]; }

  weyl::Generator simple(weyl::Generator s) const
  { return sub_root[s].simple; } // parent simple root conjugated to |sub.s|

  const WeylWord& to_simple(weyl::Generator s) const
  { return sub_root[s].to_simple; } // parent conjugating word for |simple(s)|

  const WeylWord& reflection(weyl::Generator s) const
  { return sub_root[s].reflection; } // parent reflection corresponding to |s|

  Coweight sub_2rho() const { return rd.dual_twoRho(pos_map); }
  Weight parent_sub_2rho() const { return rd.twoRho(pos_map); }

  // untwisted action of |ww|, as matrix on parent side
  LatticeMatrix action_matrix(const WeylWord& ww) const;

  RootNbrSet positive_roots() const; // for subsystem only
  InvolutionData involution_data (const WeightInvolution& theta) const;

}; // |class SubSystem|

// We have attempted to alleviate |SubSystem| by splitting off the |WeylGroup|
// The following class is for cases where a Weyl group does need to exist
class SubSystemWithGroup : public SubSystem
{
  const WeylGroup sub_W; // Weyl group no reference: built by contructor
 public:
  SubSystemWithGroup(const RootDatum& parent,
		     const RootNbrList& sub_sys // simple roots in subsys
		     );

  static SubSystemWithGroup integral // pseudo contructor for integral system
  (const RootDatum& parent, const RatWeight& gamma);

  SubSystemWithGroup(const SubSystemWithGroup& s) // copy ctor (for pseudo ctor)
  : SubSystem(s) // copy base object
  , sub_W(s.cartanMatrix()) // reconstruct (Weyl group cannot be copied)
  {
    assert(false); // should never be actually called, but exist nonetheless
  }

  const WeylGroup& Weyl_group() const { return sub_W; }
}; // |class SubSystemWithGroup|



} // |namespace subsystem|

} // |namespace atlas|

#endif
