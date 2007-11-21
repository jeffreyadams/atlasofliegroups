/*!
\file
  \brief Constructs a root datum from user interaction: class
  definitions and function declarations.

  The idea is to construct an abstract root datum (a lattice and a
  subset of "roots," together with the dual lattice and a subset of
  "coroots") specified interactively as a product of simple Lie types,
  then dividing by a specified subgroup of the center of a simply
  connected form.
*/

/*
  This is prerootdata.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For license information see the LICENSE file
*/

#ifndef PREROOTDATA_H  /* guard against multiple inclusions */
#define PREROOTDATA_H

#include "latticetypes.h"
#include "lietype.h"

namespace atlas {

/******** type declarations *************************************************/

namespace prerootdata {

  class PreRootDatum;
}

/******** function declarations **********************************************/


namespace prerootdata {

  void cartanMatrix(latticetypes::LatticeMatrix&, const lietype::LieType&);
  void cartanMatrix(latticetypes::LatticeMatrix&, 
		    const lietype::SimpleLieType&);

}

/******** type definitions **************************************************/

namespace prerootdata {
  /*!
  \brief Root datum specified by user interaction.

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
class PreRootDatum{

 private:
  /*!
  \brief List of the simple roots as elements of Z^d_rank, expressed in the
  basis specified by the argument b of the constructor.
  */
  latticetypes::WeightList d_roots;
  /*!
  List of the simple coroots as elements of Z^d_rank, expressed in the
  dual of the basis specified by the argument b of the constructor.
  */
  latticetypes::WeightList d_coroots;
  /*!
  \brief  Rank of the root datum.
   */
  size_t d_rank;

 public:

// constructors and destructors
  PreRootDatum() {}

  PreRootDatum(const lietype::LieType& lt, const latticetypes::WeightList&);

  ~PreRootDatum() {}

// accessors
  bool operator== (const PreRootDatum& prd) const {
    return (d_roots == prd.d_roots) and (d_coroots == prd.d_coroots) and
       (d_rank == prd.d_rank);
  }

  /*!
  \brief List of the simple coroots as elements of Z^d_rank, expressed
  in the dual of the basis specified by argument b of the constructor.
  */
  const latticetypes::WeightList& coroots() const {
    return d_coroots;
  }

  /*!
  \brief Rank of the root datum.
  */
  size_t rank() const {
    return d_rank;
  }

  /*!
\brief  List of the simple roots as elements of Z^d_rank,
expressed in the basis specified by argument b of the constructor.
  */
  const latticetypes::WeightList& roots() const {
    return d_roots;
  }

// manipulators
  void swap(PreRootDatum&);
};

}

}

#endif
