/*
  This is subquotient.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Classes |Subspace| and |Subquotient| (for vector spaces over $\Z/2\Z$)

#ifndef SUBQUOTIENT_H  /* guard against multiple inclusions */
#define SUBQUOTIENT_H

#include <cassert>

#include "../Atlas.h"

#include "bitvector.h" // containment

/******** function declarations **********************************************/

namespace atlas {

namespace subquotient {

template<size_t dim>
  BitMatrix<dim> subquotientMap
   (const Subquotient<dim>&,
    const Subquotient<dim>&,
    const BitMatrix<dim>&);


/******** type definitions **************************************************/



/*   Object representing a subspace of (Z/2Z)^d_rank in standard form.

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
  BitVectorList<dim> d_basis; ///< reduced echelon basis
  BitSet<dim> d_support; ///< pivot (leading nonzero bit) positions
  unsigned short int d_rank; ///< Dimension of the ambient vector space.

 public:

// constructors and destructors

  // dummy constructor needed for constructing |Subquotient|, |CartanClass|
  Subspace() : d_basis(), d_support(), d_rank(0) {}

  explicit Subspace(size_t n) : d_basis(), d_support(), d_rank(n) {}

  Subspace(const BitVectorList<dim>&, size_t);
  Subspace(const BitMatrix<dim>&); // column span

// copy, assignment: the implicitly generated versions will do fine

// accessors
  const BitVector<dim>& basis(size_t j) const
  {
    assert(j<d_rank);
    return d_basis[j];
  }

  const BitVectorList<dim>& basis() const { return d_basis; }
  size_t dimension() const { return d_basis.size(); }
  const size_t rank() const { return d_rank; }
  const BitSet<dim>& support() const { return d_support; }

  // reverse-canonical basis of perpendicular subspace of dual
  BitVectorList<dim> basis_perp () const;
  BitSet<dim> support_perp () const // complement of |d_support|
  { return BitSet<dim>(d_support).complement(rank()); } // copy and complement

  //! \brief Expresses |v| in the subspace basis.
  BitVector<dim> toBasis(BitVector<dim> v) // by-value
    const
  {
    assert(contains(v)); // implies |assert(v.size()==rank())|
    v.slice(support());  // express |v| in |d_basis|
    assert(v.size()==dimension());

    return v;
  }

  /*!
  \brief Interprets |v| in the subspace basis and returns external form
  */
  BitVector<dim> fromBasis(const BitVector<dim>& v) const
  {
    assert(v.size()==dimension());

    // expand linear combination of basis given by |v|
    return bitvector::combination(d_basis,rank(),v.data());
  }

  bool contains(const BitVector<dim>& v) const;
  bool contains(const BitVectorList<dim>& m) const;

  // the following methods reduce (finding representative) MODULO the subspace
  BitVector<dim>
    representative(const BitVector<dim>&) const;

  BitVector<dim> mod_image(const BitVector<dim>& w)
    const
  {
    return representative(w);
  }
  // destructive version
  void mod_reduce(BitVector<dim>& w) const { w=mod_image(w); }

// manipulators
  void apply (const BitMatrix<dim>&);

  void swap(Subspace&);

}; // |class Subspace|

/*  Object representing a quotient of two subspaces of (Z/2Z)^d_rank.

  Elements of the vector space are |BitVector|'s of capacity |dim|; the first
  |d_rank| coordinates are significant. (The number |d_rank| is contained both
  in |d_space| and in |d_subspace|, not directly in |Subquotient|. The number
  is accessible by the public member function |rank()|.) The larger subspace
  is specified by the |Subspace d_space|. The smaller subspace is specified by
  the |Subspace d_subspace| (not extremely enlightening names).

  A consequence of the fact that |d_subspace| is contained in |d_space| is
  that the collection of leading bits for |d_subspace| is a subset of the
  collection of leading bits for |d_space|. The \emph{difference} between two
  sets is flagged by the |BitSet d_support|. The number of set bets in
  |d_support| is therefore the dimension of the subquotient.
  */
template<size_t dim> class Subquotient
{
  Subspace<dim> d_space; // larger space, expressed as a |Subspace|.
  Subspace<dim> d_subspace; // smaller space, expressed as a Subspace.

/* The field |d_rel_support| flags the row-reduced basis vectors for |d_space|
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
BitSet<dim> d_rel_support;

 public:

// constructors and destructors
  Subquotient() : d_space(),d_subspace(), d_rel_support() {}

  explicit Subquotient(size_t n)
    : d_space(n), d_subspace(n), d_rel_support()
    {}

// Construct a subquotient of $(Z/2Z)^n$
// namely the space generated by |bsp| modulo the space generated by |bsub|
  Subquotient(const BitVectorList<dim>& bsp,
	      const BitVectorList<dim>& bsub, size_t n);

// copy, assignment: the implicitly generated versions will do fine

// accessors

  // Dimension of the subquotient.
  size_t dimension() const
    { return d_space.dimension() - d_subspace.dimension(); }

  // Dimension of the ambient vector space in which subquotient is defined
  size_t rank() const { return d_space.rank(); }

  const Subspace<dim>& space() const { return d_space; }       // numerator
  const Subspace<dim>& denominator() const { return d_subspace; }

  /* we call this |support| to the outside world, since it flags basis
    representatives for the quotient among the basis for |d_space| */
  const BitSet<dim>& support() const { return d_rel_support; }

  const BitVectorList<dim> basis() const // reprentatives for generators
  { BitVectorList<dim> result; result.reserve(d_rel_support.count());
    for (auto it=d_rel_support.begin(); it(); ++it)
      result.push_back(d_space.basis(*it));
    return result;
  }

  // Cardinality of the subquotient: 2^dimension.
  unsigned long size() const { return 1ul << dimension(); }

/* Replace by the canonical representative of |w| modulo |d_subspace|.

  It is assumed that |w| belongs to the "numerator" subspace |d_space|. Then
  all that needs to be done is reduce modulo the "denominator" |d_subspace|.
  The value remains in $(Z/2Z)^n$; see |toBasis| to express in subquotient.
*/
  void mod_reduce(BitVector<dim>& w) const
    { d_subspace.mod_reduce(w); }

  // functional version
  BitVector<dim> mod_image(const BitVector<dim>& w) const
  { return d_subspace.mod_image(w); }

  // for testing purposes;  the bits that determine the subquotient
  BitSet<dim> significantBits() const
  { return d_space.support() - d_subspace.support(); }

/* Expresses |v| in the subquotient basis.

  It is assumed that |v| belongs to |d_space|.

  We know that the subset of basis vectors of |d_space| selected by
  |d_rel_support| is independent modulo |d_subspace|, so if |v| is any linear
  combination of them, it suffices to extract and pack those coefficients. We
  can reduce to that situation by reducing |v| modulo |d_subspace|, which is
  what the call to |mod_image| below does; this will clear the bits for the
  support of |d_subspace|, while not taking us out of |d_space|.
*/
  BitVector<dim> toBasis(const BitVector<dim>& v) const
  {
    BitVector<dim> result=mod_image(v); // mod out subspace

    // implied: |assert(d_space.contains(result))|, |result.size()==rank())|
    result=d_space.toBasis(result);  // express |result| in basis of |d_space|
    assert(result.size()==d_space.dimension());

    result.slice(d_rel_support);  // remove (null) coordinates for |d_subspace|
    assert(result.size()==dimension()); // i.e., dimension of the subquotient

    return result;
  }

  // interpret |v| in the subspace basis and return an external representative
  BitVector<dim> fromBasis(BitVector<dim> v) const // by-value
  {
    assert(v.size()==dimension());

    v.unslice(d_rel_support,d_space.dimension()); // insert null coordinates
    return d_space.fromBasis(v);
  }

// manipulators
  void apply (const BitMatrix<dim>&);

  void swap(Subquotient&);

 }; // |template <size_t dim> class Subquotient|

} // |namespace subquotient|

} // |namespace atlas|

#endif
