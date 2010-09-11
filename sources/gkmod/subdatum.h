/*!
\file
\brief Class definitions and function declarations for the RootDatum class.
*/
/*
  This is subdatum.h

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SUBDATUM_H  /* guard against multiple inclusions */
#define SUBDATUM_H

#include "rootdata.h" // derive from |rootdata::RootSystem|
#include "prerootdata.h"
#include "weyl.h"
#include "tits.h"
#include "latticetypes_fwd.h"
#include "gradings_fwd.h"
#include "realredgp_fwd.h"
#include "kgb_fwd.h"

namespace atlas {

namespace subdatum {

/* The following data type is specific for representation theory. It is
   associated so a subsystem of the dual root datum (with positive roots
   contained in the positive parent coroots), and is derived from |RootSystem|
   at that dual side. It remains however attached to the parent root datum,
   contains and exports a reference to that rootdatum, and has a method to
   produce a |PreRootDatum| for the subsystem of the parent defined by simple
   coroots for the subsystem (not a full root datum, for efficientcy reasons).
 */

// A subsystem on the dual side of a given root datum
class SubSystem : public rootdata::RootSystem // new system, subsytem of dual
{
  const rootdata::RootDatum& rd; // parent root datum
  const weyl::WeylGroup sub_W; // Weyl group no reference: built in contructor
  rootdata::RootList pos_map; // map positive roots to root number in parent
  rootdata::RootList inv_map; // partial map back from all parent roots

  struct root_info
  { weyl::Generator simple;
    weyl::WeylWord to_simple;
    weyl::WeylWord reflection;

  root_info() : simple(~0), to_simple(), reflection() {}
  };
  std::vector<root_info> sub_root;

 public:
  SubSystem(const rootdata::RootDatum& parent,
	    const rootdata::RootList& sub_sys);

  static SubSystem integral // pseudo contructor for integral system
  (const rootdata::RootDatum& parent, const latticetypes::RatWeight& gamma);

  SubSystem(const SubSystem& s) // copy contructor (used by pseudo contructor)
  : rootdata::RootSystem(s) // copy base object
  , rd(s.rd) // share this one
  , sub_W(s.cartanMatrix()) // reconstruct (Weyl group cannot be copied)
  , pos_map(s.pos_map), inv_map(s.inv_map), sub_root(s.sub_root) // copy those
  {
    assert(false); // should never be actually called, but exist nonetheless
  }

  const rootdata::RootDatum& parent_datum() const { return rd; }
  const weyl::WeylGroup& Weyl_group() const { return sub_W; }

  weyl::Twist twist(const latticetypes::LatticeMatrix& theta,
		    weyl::WeylWord& ww) const;
  // output value |ww| gives |-^theta| as twisted involution for |sub|

  weyl::Twist parent_twist(const latticetypes::LatticeMatrix& theta,
			   weyl::WeylWord& ww) const;
  // now twist |delta| (result) and |ww| are on parent side, |theta=ww.delta|.

  prerootdata::PreRootDatum pre_root_datum() const; // viewed from parent side

  rootdata::RootNbr parent_nr_simple(weyl::Generator s) const
  { return pos_map[s]; }

  rootdata::RootNbr parent_nr(rootdata::RootNbr alpha) const;

  weyl::Generator simple(weyl::Generator s) const
  { return sub_root[s].simple; }

  const weyl::WeylWord& to_simple(weyl::Generator s) const
  { return sub_root[s].to_simple; }

  const weyl::WeylWord& reflection(weyl::Generator s) const
  { return sub_root[s].reflection; }

  latticetypes::Weight sub_2rho() const { return rd.dual_twoRho(pos_map); }
  latticetypes::Weight parent_sub_2rho() const { return rd.twoRho(pos_map); }

  latticetypes::LatticeMatrix action_matrix(const weyl::WeylWord& ww) const;

  gradings::Grading induced(gradings::Grading base_grading) const;

}; // |class SubSystem|

class SubDatum : public SubSystem
{
  weyl::WeylWord base_ww;  // we need this variable mostly in the constructor!
  latticetypes::LatticeMatrix delta; // together with twist: what was missing
  tits::TitsGroup Tg; // twisted Weyl group, plus stuff our base class knows
  weyl::WeylElt ini_tw;

 public:
  SubDatum(realredgp::RealReductiveGroup& GR,
	   const latticetypes::RatWeight& gamma,
	   kgb::KGBElt x);

  const tits::TitsGroup& Tits_group() const {return Tg; }
  weyl::TwistedInvolution init_twisted() const { return ini_tw; }
  const weyl::WeylWord& base_twisted_in_parent() const { return base_ww; }

  latticetypes::LatticeMatrix involution(weyl::TwistedInvolution tw) const;
}; // |class SubDatum|

} // |namespace subdatum|

} // |namespace atlas|

#endif
