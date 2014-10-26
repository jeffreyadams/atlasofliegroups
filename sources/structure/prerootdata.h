/*
  Constructing a root datum from user interaction

  The idea is to construct an abstract root datum (a lattice and a
  subset of "roots," together with the dual lattice and a subset of
  "coroots") specified interactively as a product of simple Lie types,
  then dividing by a specified subgroup of the center of a simply
  connected form.
*/

/*
  This is prerootdata.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations
  Copyright (C) 2014 Marc van Leeuwen

  For license information see the LICENSE file
*/

#ifndef PREROOTDATA_H  /* guard against multiple inclusions */
#define PREROOTDATA_H

#include "atlas_types.h"
#include <stdexcept>

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
  /*!
  \brief List of the simple roots as elements of Z^d_rank, expressed in the
  basis specified by the argument b of the constructor.
  */
  WeightList d_roots;
  /*!
  List of the simple coroots as elements of Z^d_rank, expressed in the
  dual of the basis specified by the argument b of the constructor.
  */
  CoweightList d_coroots;
  /*!
  \brief  Rank of the root datum.
   */
  size_t d_rank;

 public:

// constructors and destructors
  PreRootDatum() {}

  PreRootDatum(const WeightList& roots,
               const CoweightList& coroots,
	       size_t rank)
    : d_roots(roots),d_coroots(coroots), d_rank(rank) {}

  PreRootDatum(const LieType& lt);

  ~PreRootDatum() {}

// accessors
  bool operator== (const PreRootDatum& prd) const {
    return (d_roots == prd.d_roots) and (d_coroots == prd.d_coroots) and
       (d_rank == prd.d_rank);
  }

  /*!
  \brief Rank of the root datum.
  */
  size_t rank() const { return d_rank; }

  // List of the simple roots, in basis that is implicit in constructor.
  const WeightList& simple_roots() const { return d_roots; }
  // List of the simple coroots, in basis that is implicit in constructor.
  const CoweightList& simple_coroots() const { return d_coroots; }

  int_Matrix Cartan_matrix() const;

  template<typename C>
    void simpleReflect(matrix::Vector<C>& v, weyl::Generator i) const;
  void simple_reflect(weyl::Generator i, LatticeMatrix& M) const;
  void simple_reflect(LatticeMatrix& M,weyl::Generator i) const;

// manipulators
  // replace by root datum for finite central quotient with weight |sublattice|
  PreRootDatum& quotient(const LatticeMatrix& sublattice);

};

} // |namespace prerootdata|

} // |namespace atlas|

#endif
