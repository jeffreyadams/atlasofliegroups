/*
  Class definition for the Subsystem class.
*/
/*
  This is subystem.h

  Copyright (C) 2010,2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef SUBSYSTEM_H  /* guard against multiple inclusions */
#define SUBSYSTEM_H

#include "../Atlas.h"

#include "rootdata.h"	// derive from |RootSystem|
#include "weyl.h"	// containment

namespace atlas {

namespace subsystem {

/*
  The following data type has a purpose specific for representation theory. It
  is associated to a subdatum of the root datum, whose coroot system typically
  is an "additively closed" part of the set of parent coroots (in the sense that
  any sum of coroots in the part that is a parent coroot must also be in the
  part), for instance the coroots that are integral on a given weight. This
  property is nowhere used however: we take the subsystem to be specified by a
  set of generating positive (co)roots, with the subsystem taken to be their
  closure under the corresponding reflections. We derive from |RootSystem|,
  which base class we use to generate the set of positive (co)roots (using the
  Cartan matrix for the subsystem); we avoid deriving from |RootDatum| to save
  memory, but we provide methods in coordinates using the parent datum.

  The |RootSystem| from which we derive used to be on the dual side to reflect
  the fact that it was mostly defined in terms of coroots, but this had hardly
  any effect (since |RootSystem| methods are mostly duality agnostic), and if
  anything caused terminological confusion. So this side-swapping is no longer
  in effect and roots of the subsystem are also roots of the parent system.
*/

// A subsystem on the dual side of a given root datum
class SubSystem : public RootSystem // new system, subsystem of parent
{
  const RootDatum& rd; // parent root datum
  RootNbrSet which; // subset of parent posroots that are in subsystem
  RootNbrList pos_map; // map positive roots to root number in parent
  RootNbrList inv_map; // map back from parent roots flagged in |which|

  struct root_info
  { weyl::Generator simple; // some simple root $s$ in parent conjugate to root
    WeylWord to_simple; // word $w$ conjugating the root to mentioned simple
    WeylWord reflection; // reflection word for root: $w^{-1}sw$

  root_info() : simple(~0), to_simple(), reflection() {}
  };
  std::vector<root_info> sub_root;

 public:
  SubSystem(const RootDatum& parent,
	    const sl_list<RootNbr>& sub_sys // list of simple roots in subsys
	    // those simple roots must be positive roots of |parent|
           );

  static SubSystem // pseudo contructor for integral system
    integral(const RootDatum& parent, const RatWeight& gamma);

  SubSystem(const SubSystem& s) // copy contructor, used in |common_block|
    : RootSystem(s) // copy base object
    , rd(s.rd) // share this one
    , which(s.which) // copy other fields
    , pos_map(s.pos_map)
    , inv_map(s.inv_map)
    , sub_root(s.sub_root)
  {}

  SubSystem(SubSystem&& s) // move constructor; cannot use |default|
    : RootSystem(std::move(s)) // move base object
    , rd(s.rd) // share this one
    , which(std::move(s.which)) // move other fields
    , pos_map(std::move(s.pos_map))
    , inv_map(std::move(s.inv_map))
    , sub_root(std::move(s.sub_root))
  {}

  const RootDatum& parent_datum() const { return rd; }
  // RootNbr rank() const; // integral rank; inherited from |RootSystem|

  PreRootDatum pre_root_datum() const;

  RootNbr parent_nr_simple(weyl::Generator s) const { return pos_map[s]; }

  RootNbr to_parent(RootNbr alpha) const; // |pos_map| with some shifting
  RootNbr from_parent(RootNbr alpha) const;
  // could be |RootNbr(-1)| if absent, or if |alpha==rd.numPosRoots|

  weyl::Generator simple(unsigned int n) const
  { assert(n<numPosRoots()); // n must be a positive-root index for subsystem
    return sub_root[n].simple; // parent simple root conjugated to |sub.s|
  }

  const WeylWord& to_simple(unsigned int n) const
  { assert(n<numPosRoots()); // n must be a positive-root index for subsystem
    return sub_root[n].to_simple; // parent conjugating word for |simple(s)|
  }

  const WeylWord& reflection(unsigned int n) const
  { assert(n<numPosRoots()); // n must be a positive-root index for subsystem
    return sub_root[n].reflection; // parent reflection corresponding to |n|
  }

  const Weight& simple_root(weyl::Generator s) const
  { return parent_datum().root(parent_nr_simple(s)); }
  const Coweight& simple_coroot(weyl::Generator s) const
  { return parent_datum().coroot(parent_nr_simple(s)); }

  // numbers in parent for the positive (co)roots of the subsystem
  RootNbrSet positive_roots() const; // for subsystem, as |parent| roots
  const RootNbrSet& posroot_subset() const { return which; } // as posroots

  // methods that avoid building full |RootDatum srd(pre_root_datum())|
  template<typename C>
    void simple_reflect(weyl::Generator s,matrix::Vector<C>& v,int d=0) const
  { parent_datum().reflect(parent_nr_simple(s),v,d); }
  template<typename C>
    void simple_coreflect(matrix::Vector<C>& v,weyl::Generator s,int d=0) const
  { parent_datum().coreflect(v,parent_nr_simple(s),d); }

  ext_gens fold_orbits (const WeightInvolution& delta) const; // as for |srd|
  RankFlags singular_generators(const RatWeight& gamma) const; // as for |srd|

  InvolutionData involution_data (const WeightInvolution& theta) const;

}; // |class SubSystem|

// We have attempted to alleviate |SubSystem| by splitting off the |WeylGroup|
// The following class is for cases where a Weyl group does need to exist
class SubSystemWithGroup : public SubSystem
{
  WeylGroup sub_W; // Weyl is no group no reference, but built by contructor
 public:
  SubSystemWithGroup(const RootDatum& parent,
		     const sl_list<RootNbr>& sub_sys // simple roots in subsys
		     );

  static SubSystemWithGroup integral // pseudo contructor for integral system
  (const RootDatum& parent, const RatWeight& gamma);

  // move ctor (for pseudo ctor)
  SubSystemWithGroup(SubSystemWithGroup&& s) = default;

  const WeylGroup& Weyl_group() const { return sub_W; }
}; // |class SubSystemWithGroup|

struct integral_datum_entry // hashable (integral) subset of positive roots
{
  RootNbrSet posroots; // as set of posroot indices, so starting at 0

  integral_datum_entry (const RootNbrSet& p) : posroots(p) {}

  using Pooltype = std::vector<integral_datum_entry>;
  bool operator!=(const integral_datum_entry& other) const
  { return posroots!=other.posroots; }
  size_t hashCode(size_t modulus) const
  {
    size_t h = 2;
    for (unsigned count = 0; count<posroots.capacity(); count+=32)
      h = 243*h + posroots.range(count,32);
    return h&(modulus-1);
  }
 }; // |struct integral_datum_entry|

class integral_datum_item
{
  std::unique_ptr<SubSystem> // pointer level avoids |SubSystem| being moved
    integral; // references full root datum, presents integral datum
  int_Matrix simple_coroots; // convenience, for creating |codec| values

 public:
  struct codec
  {
    const int_Matrix& coroots_matrix;
    const int_Matrix& theta_1_image_basis;
    std::vector<int> diagonal; // inv.factors image $(1-\theta)X^*$ in $X^*/N$
    int_Matrix in, out; // from $X^*/N$ to $-1$ subspace to $(1-\theta)X^*$
    codec (const InnerClass& ic,
	   unsigned int isys, InvolutionNbr inv, const int_Matrix& cmat);
  }; // |struct integral_datum_item::codec|

  integral_datum_item(InnerClass& ic,const RootNbrSet& int_posroots);
  integral_datum_item(integral_datum_item&&)=default; // move, never copy

  const SubSystem& int_system() const { return *integral; }
  codec data(const InnerClass& ic, unsigned int isys, InvolutionNbr inv) const
  { return { ic,isys,inv,simple_coroots }; }
  const int_Matrix& coroots_matrix() const { return simple_coroots; }

}; // |class integral_datum_item|


} // |namespace subsystem|

} // |namespace atlas|

#endif
