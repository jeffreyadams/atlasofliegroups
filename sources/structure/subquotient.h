/*!
\file
  \brief Class definitions and function declarations for the
  classes NormalSubspace and NormalSubquotient. 

  These classes deal with subspaces and subquotients
  of vector spaces over Z/2Z, elements of which are of type
  BitVector.  The term "normal" refers to the fact that bases for the
  subspaces are row reduced.
*/
/*
  This is subquotient.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef SUBQUOTIENT_H  /* guard against multiple inclusions */
#define SUBQUOTIENT_H

#include "subquotient_fwd.h"

#include "bitset.h"
#include "bitvector.h"

/******** function declarations **********************************************/

namespace atlas {

namespace subquotient {

template<size_t dim>
  void subquotientMap(bitvector::BitMatrix<dim>&, 
		      const NormalSubquotient<dim>&,
		      const NormalSubquotient<dim>&, 
		      const bitvector::BitMatrix<dim>&);

}

/******** type definitions **************************************************/

namespace subquotient {

  /*!
  \brief Subspace of (Z/2Z)^d_rank.

  Elements of the vector space are BitVector's of capacity dim; the
  first d_rank coordinates are significant.  The subspace is specified
  by its (unique) basis d_basis which is in row-reduced form.  (This
  is the significance of "Normal.") This means that the location of
  the leading non-zero bit in each vector of d_basis holds zero in all
  the other basis vectors; and the locations of these leading bits
  increase.

  The locations of the leading bits are flagged by the BitSet
  d_support; the number of set bits in d_support is the dimension of
  NormalSubspace. 
  */
template<size_t dim> class NormalSubspace {

 private:

  /*!
  \brief Dimension of the ambient vector space.
  */
  size_t d_rank;

  /*!
  \brief Basis of NormalSubspace in row-reduced form.
  */
  std::vector<bitvector::BitVector<dim> > d_basis;

  /*!
  \brief BitSet flagging the leading non-zero bits in the basis of
  NormalSubspace.
  */
  bitset::BitSet<dim> d_support;

 public:

// constructors and destructors
  NormalSubspace():d_rank(0) {}

  explicit NormalSubspace(size_t n):d_rank(n) {}
  
  NormalSubspace(const std::vector<bitvector::BitVector<dim> >&, size_t);

  ~NormalSubspace() {}

// accessors
  const bitvector::BitVector<dim>& basis(size_t j) const {
    return d_basis[j];
  }

  const std::vector<bitvector::BitVector<dim> >& basis() const {
    return d_basis;
  }

  size_t dimension() const {
    return d_basis.size();
  }

  const size_t rank() const {
    return d_rank;
  }

  void representative(bitvector::BitVector<dim>&, 
		      const bitvector::BitVector<dim>&) const;

  const bitset::BitSet<dim>& support() const {
    return d_support;
  }

// manipulators
  void apply (const bitvector::BitMatrix<dim>&);

  void swap(NormalSubspace&);

};

  /*!
  \brief Quotient of two subspaces of (Z/2Z)^d_rank.

  Elements of the vector space are BitVector's of capacity dim; the
  first d_rank coordinates are significant.  (The number d_rank is
  owned by d_space and by d_subspace, not directly by
  NormalSubquotient.  The number is accessible by the public member
  function rank().) The larger subspace is specified by the
  NormalSubspace d_space.  The smaller subspace is specified by the
  NormalSubspace d_subspace.

  A consequence of (d_subspace contained in d_space) is that the
  collection of leading bits for d_subspace is a subset of the
  collection of leading bits for d_space.  The difference of these two
  sets is flagged by the BitSet d_support.  The number of set bets in
  d_support is therefore the dimension of the subquotient. 
  */
template<size_t dim> class NormalSubquotient {

 private:

  /*!
  Larger space, expressed as a NormalSubspace.
  */
  NormalSubspace<dim> d_space;

  /*!
  Smaller space, expressed as a NormalSubspace.
  */
  NormalSubspace<dim> d_subspace;
  
  /*!
  \brief Flags the row-reduced basis vectors for the larger
  space spanning the canonical complement to the smaller space.

  These are the basis vectors whose leading bits do _not_ appear as
  leading bits of row-reduced basis vectors for the smaller space.
  The corresponding set of basis vectors for the larger space
  constitute a basis for a cross section of the subquotient in the
  larger space.  In particular, the number of set bits is equal to the
  dimension of the subquotient.

  Notice that what is flagged is positions in the basis of the larger
  space, _not_ the positions of the leading bits in (Z/2Z)^n.
  
  Example: d_space=(1010,0110,0001) d_space.support={0,1,3}
           d_subspace=(0111) d_subspace.support={1}
           d_support={0,2} (corresponding to the basis vectors 1010
           and 0001 for d_space spanning the canonical complement to
           d_subspace).
  */
bitset::BitSet<dim> d_support;

 public:

// constructors and destructors
  NormalSubquotient() {}

  explicit NormalSubquotient(size_t n):d_space(n),d_subspace(n) {}
  
  /*!
  Constructs a subquotient of (Z/2Z)^n.  Here n is the last argument;
  the larger space is generated by the first set of BitVector's; and
  the smaller space is generated by the second set of BitVector's.
  */
  NormalSubquotient(const std::vector<bitvector::BitVector<dim> >&,
		    const std::vector<bitvector::BitVector<dim> >&, size_t);
  
  ~NormalSubquotient() {}

// accessors

  /*!
  \brief Dimension of the subquotient.
  */
  size_t dimension() const {
    return d_space.dimension() - d_subspace.dimension();
  }

  /*!
  \brief Dimension of the ambient (Z/2Z)^rank in which the larger and smaller
  subspaces live.
  */
  size_t rank() const {
    return d_space.rank();
  }

  void representative(bitvector::BitVector<dim>& r, 
		      const bitvector::BitVector<dim>& w) const {
    /*!
    \brief Puts in r the canonical representative of w.
    
    Express the image of w in the subquotient using the image of the
    standard basis of the larger space.  Then r is the corresponding
    combination of those standard basis vectors.  It is assumed that w
    belongs to the larger space.
    */
    d_subspace.representative(r,w);
  }

  /*!
  \brief Cardinality of the subquotient: 2^dimension.
  */
  unsigned long size() const {
    return 1ul << dimension();
  }

  const NormalSubspace<dim>& space() const {
    return d_space;
  }

  const NormalSubspace<dim>& subspace() const {
    return d_subspace;
  }

  bitset::BitSet<dim> significantBits() const {
  // for testing purposes; these are the bits that determine the subquotient
    return d_space.support() & ~d_subspace.support();
  }

  const bitset::BitSet<dim>& support() const {
    return d_support;
  }

  /*!
  \brief Expresses v in the subquotient basis.

  It is assumed that v belongs to d_space, and is equal to its
  canonical representative modulo d_subspace.  The second condition
  can be achieved with a call to representative.
  */
  void toSubquotient(bitvector::BitVector<dim>& v) const {
    v.slice(d_space.support());  //express v in basis of d_space
    v.slice(d_support);
  }

// manipulators
  void apply (const bitvector::BitMatrix<dim>&);

  void swap(NormalSubquotient&);

};

}

}

#include "subquotient_def.h"

#endif
