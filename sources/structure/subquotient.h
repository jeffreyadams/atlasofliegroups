/*!
\file
  \brief Class definitions and function declarations for the class templates
  |Subspace| and |Subquotient|.
*/
/*
  This is subquotient.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SUBQUOTIENT_H  /* guard against multiple inclusions */
#define SUBQUOTIENT_H

#include <cassert>

#include "subquotient_fwd.h"
#include "latticetypes_fwd.h"

#include "bitset.h"
#include "bitvector.h"

/******** function declarations **********************************************/

namespace atlas {

namespace subquotient {

template<size_t dim>
  bitvector::BitMatrix<dim> subquotientMap
   (const Subquotient<dim>&,
    const Subquotient<dim>&,
    const bitvector::BitMatrix<dim>&);

}

/******** type definitions **************************************************/

namespace subquotient {

  /*!
  \brief Subspace of (Z/2Z)^d_rank.

  Elements of the vector space are |BitVector|s of capacity |dim| and actual
  size |d_rank|. The subspace is recorded by its (unique) ordered canonical
  basis |d_basis|. This means its vectors are reduced (in the sense of
  row-reduced matrices); that is (in the binary case) in the position of the
  leading (nonzero) bit of each vector, the other vectors heve a bit $0$, and
  the vectors are ordered by increasing position of their leading bit.

  The locations of the leading bits are flagged by the |BitSet| |d_support|
  (of capacity |d_rank|); the number of set bits in |d_support| is equal to
  the dimension |d_basis.size()| of the |Subspace|.
  */
template<size_t dim> class Subspace
{
  size_t d_rank; ///< Dimension of the ambient vector space.
  bitvector::BitVectorList<dim> d_basis; ///< reduced echelon basis
  bitset::BitSet<dim> d_support; ///< pivot (leading nonzero bit) positions

 public:

// constructors and destructors

  // dummy constructor needed for constructing |Subquotient|, |CartanClass|
  Subspace() : d_rank(0), d_basis(), d_support() {}

  explicit Subspace(size_t n) : d_rank(n), d_basis(), d_support() {}

  Subspace(const bitvector::BitVectorList<dim>&, size_t);

// copy, assignment: the implicitly generated versions will do fine

// accessors
  const bitvector::BitVector<dim>& basis(size_t j) const
  {
    assert(j<d_rank);
    return d_basis[j];
  }

  const bitvector::BitVectorList<dim>& basis() const { return d_basis; }
  size_t dimension() const { return d_basis.size(); }
  const size_t rank() const { return d_rank; }
  const bitset::BitSet<dim>& support() const { return d_support; }

  //! \brief Expresses |v| in the subspace basis.
  bitvector::BitVector<dim> toBasis(bitvector::BitVector<dim> v) const
  {
    assert(contains(v)); // implies |assert(v.size()==rank())|
    v.slice(support());  // express |v| in |d_basis|
    assert(v.size()==dimension());

    return v;
  }

  /*!
  \brief Interprets |v| in the subspace basis and returns external form
  */
  bitvector::BitVector<dim> fromBasis(bitvector::BitVector<dim> v) const
  {
    assert(v.size()==dimension());

    // expand linear combination of basis given by |v|
    return bitvector::combination(d_basis,rank(),v.data());
  }

  bool contains(const bitvector::BitVector<dim>& v) const;
  bool contains(const bitvector::BitVectorList<dim>& m) const;

  // the following methods reduce (finding representative) MODULO the subspace
  bitvector::BitVector<dim>
    representative(const bitvector::BitVector<dim>&) const;

  bitvector::BitVector<dim> mod_image(const bitvector::BitVector<dim>& w)
    const
  {
    return representative(w);
  }
  // destructive version
  void mod_reduce(bitvector::BitVector<dim>& w) const { w=mod_image(w); }

// manipulators
  void apply (const bitvector::BitMatrix<dim>&);

  void swap(Subspace&);

}; // class Subspace

  /*!
  \brief Quotient of two subspaces of (Z/2Z)^d_rank.

  Elements of the vector space are BitVector's of capacity dim; the first
  d_rank coordinates are significant. (The number d_rank is owned both by
  d_space and by d_subspace, not directly by Subquotient. The number is
  accessible by the public member function rank().) The larger subspace is
  specified by the Subspace d_space. The smaller subspace is specified by the
  Subspace d_subspace.

  A consequence of (d_subspace contained in d_space) is that the
  collection of leading bits for d_subspace is a subset of the
  collection of leading bits for d_space.  The difference of these two
  sets is flagged by the BitSet d_support.  The number of set bets in
  d_support is therefore the dimension of the subquotient.
  */
template<size_t dim> class Subquotient
{

 private:

  /*!
  Larger space, expressed as a Subspace.
  */
  Subspace<dim> d_space;

  /*!
  Smaller space, expressed as a Subspace.
  */
  Subspace<dim> d_subspace;

  /*!
  \brief Flags the row-reduced basis vectors for the (larger) |d_space|
  that span the canonical complement to the (smaller) |d_subspace|.

  The bits are set at the positions indexing the basis vectors whose leading
  bits do _not_ appear as leading bits of row-reduced basis vectors for the
  smaller space. The flagged set of basis vectors for |d_space| constitute a
  basis for a cross section of the subquotient in the larger space. In
  particular, the number of set bits is equal to the dimension of the
  subquotient.

  Note that what is flagged is are numbers of basis vectors of the larger
  space, _not_ the positions of their leading bits in $(Z/2Z)^n$. Consequently
  the effective number of bits is |d_space.dimension()| rather than |rank()|,
  although this number is not recorded in the value |d_rel_support| itself.


  Example: d_space=(1010,0110,0001) d_space.support={0,1,3}
           d_subspace=(0111) d_subspace.support={1}
           d_rel_support={0,2} (corresponding to the basis vectors 1010
           and 0001 for d_space spanning the canonical complement to
           d_subspace).
  */
bitset::BitSet<dim> d_rel_support;

 public:

// constructors and destructors
  Subquotient() : d_space(),d_subspace(), d_rel_support() {}

  explicit Subquotient(size_t n)
    : d_space(n), d_subspace(n), d_rel_support()
    {}

  /*! Constructs a subquotient of $(Z/2Z)^n$, formed by the space generated by
  |bsp| modulo the space generated by |bsub|
  */
  Subquotient(const bitvector::BitVectorList<dim>& bsp,
	      const bitvector::BitVectorList<dim>& bsub, size_t n);

// copy, assignment: the implicitly generated versions will do fine

// accessors

  /*!
  \brief Dimension of the subquotient.
  */
  size_t dimension() const
    { return d_space.dimension() - d_subspace.dimension(); }

  /*!
  \brief Dimension of the ambient vector space in which the larger and smaller
  subspaces live.
  */
  size_t rank() const { return d_space.rank(); }

  const Subspace<dim>& space() const { return d_space; }       // numerator
  const Subspace<dim>& subspace() const { return d_subspace; } // denominator

  /* we call this |support| to the outside world, since it flags basis
    representatives for the quotient among the basis for |d_space| */
  const bitset::BitSet<dim>& support() const { return d_rel_support; }

  /*!
  \brief Cardinality of the subquotient: 2^dimension.
  */
  unsigned long size() const { return 1ul << dimension(); }

  /*!
    \brief Puts in |r| the canonical representative of |w| modulo |d_subspace|.

    It is assumed that |w| belongs to the "numerator" subspace |d_space|. Then
    all that needs to be done is reduce modulo the "denominator" |d_subspace|.
    The value remains in $(Z/2Z)^n$; see |toBasis| to express in subquotient.
  */
  void representative(bitvector::BitVector<dim>& r,
		      const bitvector::BitVector<dim>& w) const
    { r=d_subspace.mod_image(w); }

  // destructive version
  void mod_reduce(bitvector::BitVector<dim>& w) const
    { d_subspace.mod_reduce(w); }

  // functional version
  bitvector::BitVector<dim> mod_image(const bitvector::BitVector<dim>& w)
    const
  { return d_subspace.mod_image(w); }

  bitset::BitSet<dim> significantBits() const {
  // for testing purposes; these are the bits that determine the subquotient
    return d_space.support() - d_subspace.support();
  }

  /*!
  \brief Expresses |v| in the subquotient basis.

  It is assumed that |v| belongs to |d_space|.

  We know that the subset of basis vectors of |d_space| selected by
  |d_rel_support| is independent modulo |d_subspace|, so if |v| is any linear
  combination of them, it suffices to extract and pack those coefficients. We
  can reduce to that situation by reducing |v| modulo |d_subspace|, which is
  what the call to |mod_image| below does; this will clear the bits for the
  support of |d_subspace|, while not taking us out of |d_space|.
  */
  bitvector::BitVector<dim> toBasis(const bitvector::BitVector<dim>& v)
    const
  {
    bitvector::BitVector<dim> result=mod_image(v); // mod out subspace

    // implied: |assert(d_space.contains(result))|, |result.size()==rank())|
    result=d_space.toBasis(result);  // express |result| in basis of |d_space|
    assert(result.size()==d_space.dimension());

    result.slice(d_rel_support);  // remove (null) coordinates for |d_subspace|
    assert(result.size()==dimension()); // i.e., dimension of the subquotient

    return result;
  }

  /*!
  \brief Interprets |v| in the subspace basis and returns external form
  */
  bitvector::BitVector<dim> fromBasis(bitvector::BitVector<dim> v) const
  {
    assert(v.size()==dimension());

    v.unslice(d_rel_support,d_space.dimension()); // insert null coordinates
    return d_space.fromBasis(v);
  }

// manipulators
  void apply (const bitvector::BitMatrix<dim>&);

  void swap(Subquotient&);

 }; // |template <size_t dim> class Subquotient|

}

}

#endif
