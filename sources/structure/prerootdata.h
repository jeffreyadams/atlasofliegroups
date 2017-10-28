/*
  This is prerootdata.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2014-2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  A class |PreRootDatum| to store simple rooots and coroots, and perform some
  tests and operations prior to generating the full root systen and root datum
*/

#ifndef PREROOTDATA_H  /* guard against multiple inclusions */
#define PREROOTDATA_H

#include "../Atlas.h"
#include "matrix.h"

namespace atlas {

/******** type definitions **************************************************/

namespace prerootdata {
  /* Root datum specified by user interaction, having passed validity tests.

  The idea is to construct an abstract root datum (a lattice and a
  subset of "roots," together with the dual lattice and a subset of
  "coroots") specified interactively as a product of simple Lie types,
  then dividing by a specified subgroup of the center of a simply
  connected form.

  The lattice and dual lattice in which the root datum lives are
  always Z^d_rank.  The simple roots, expressed as a list of elements
  of Z^d_rank, are in d_roots, and the simple coroots in d_coroots.

  More serious calculations (like producing the complete list of
  roots) are handled by the RootDatum class, to which PreRootDatum can
  pass its contents.
  */
class PreRootDatum
{
  int_Matrix simple_roots, simple_coroots; // both of same size
  bool prefer_co; // flag indicating preference for coroots in generation
 public:

// constructors and destructors
  PreRootDatum(const WeightList& roots,
               const CoweightList& coroots,
	       size_t rank,
	       bool prefer_co)
    : simple_roots(roots.begin(),roots.end(),rank,tags::IteratorTag())
    , simple_coroots(coroots.begin(),coroots.end(),rank,tags::IteratorTag())
    , prefer_co(prefer_co)
    {}

  PreRootDatum(const LieType& lt, bool prefer_co);
  PreRootDatum(int_Matrix& projector, const PreRootDatum& rd,tags::DerivedTag);
  PreRootDatum(int_Matrix& injector, const PreRootDatum& rd,tags::CoderivedTag);

  PreRootDatum(const PreRootDatum& prd, tags::DualTag)
    : simple_roots(prd.simple_coroots)
    , simple_coroots(prd.simple_roots)
    , prefer_co(not prd.prefer_co)
  {}

// accessors
  bool operator== (const PreRootDatum& prd) const
  { return simple_roots==prd.simple_roots
    and simple_coroots == prd.simple_coroots
    and prefer_co == prd.prefer_co;
  }

  size_t rank() const { return simple_roots.numRows(); }
  size_t semisimple_rank() const { return simple_roots.numColumns(); }

  const int_Matrix& simple_roots_mat() const { return simple_roots; }
  const int_Matrix& simple_coroots_mat() const { return simple_coroots; }

  bool prefer_coroots() const { return prefer_co; }

  Weight simple_root(unsigned int j) const
    { return simple_roots.column(j); }
  Coweight simple_coroot(unsigned int j) const
    { return simple_coroots.column(j); }

  int_Matrix Cartan_matrix() const;

  template<typename C>
    void simple_reflect(weyl::Generator i,matrix::Vector<C>& v) const;
  void simple_reflect(weyl::Generator i, LatticeMatrix& M) const;
  void simple_reflect(LatticeMatrix& M,weyl::Generator i) const;

// manipulators
  // replace by root datum for finite central quotient with weight |sublattice|
  PreRootDatum& quotient(const LatticeMatrix& sublattice);
  PreRootDatum& dualise() { simple_roots.swap(simple_coroots); return *this; }

}; // |class PreRootDatum|

} // |namespace prerootdata|

} // |namespace atlas|

#endif
