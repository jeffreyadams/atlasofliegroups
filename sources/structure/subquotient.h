/*
  This is subquotient.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

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

template<size_t dim> class NormalSubspace {

 private:

  size_t d_rank;
  std::vector<bitvector::BitVector<dim> > d_basis;
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

template<size_t dim> class NormalSubquotient {

 private:

  NormalSubspace<dim> d_space;
  NormalSubspace<dim> d_subspace;
  bitset::BitSet<dim> d_support;

 public:

// constructors and destructors
  NormalSubquotient() {}

  explicit NormalSubquotient(size_t n):d_space(n),d_subspace(n) {}
  
  NormalSubquotient(const std::vector<bitvector::BitVector<dim> >&,
		    const std::vector<bitvector::BitVector<dim> >&, size_t);
  
  ~NormalSubquotient() {}

// accessors
  size_t dimension() const {
    return d_space.dimension() - d_subspace.dimension();
  }

  size_t rank() const {
    return d_space.rank();
  }

  void representative(bitvector::BitVector<dim>& r, 
		      const bitvector::BitVector<dim>& w) const {
    // puts in r the canonical representative of w
    // it is assumed that w belongs to space
    d_subspace.representative(r,w);
  }

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

  void toSubquotient(bitvector::BitVector<dim>& v) const {
    // expresses v in the subquotient basis
    // it is assumed that v belongs to d_subspace
    v.slice(d_space.support());
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
