/*!
  \file
  \brief Definitions of the templates declared in lattice.h.
*/
/*
  This is lattice_def.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <cassert>
#include "matrix.h"
#include "tags.h"

#include "latticetypes.h"

/*****************************************************************************

  This file contains the definitions of the templates declared in lattice.h

******************************************************************************/

namespace atlas {

namespace lattice {

/*!
  Reduces |v| modulo 2. It is assumed that |v.size() <= dim|.
*/
template<size_t dim> void mod2(bitvector::BitVector<dim>& v2,
			       const latticetypes::LatticeElt& v)
{
  v2.resize(v.size());

  for (size_t j = 0; j < v.size(); ++j)
    v2.set_mod2(j,v[j]);
}


/*!
  Reduces all elements in the list modulo 2.
*/
template<size_t dim> void mod2(std::vector<bitvector::BitVector<dim> >& wl2,
			       const latticetypes::WeightList& wl)
{
  wl2.resize(wl.size());

  for (size_t j = 0; j < wl.size(); ++j)
    mod2(wl2[j],wl[j]);
}

/*!
  Reduces the matrix |m| modulo 2 into |m2|; it is assumed that the number
  of rows in |m| is at most |dim|.
*/
template<size_t dim> void mod2(bitvector::BitMatrix<dim>& m2,
			       const latticetypes::LatticeMatrix& m)
{
  assert(m.numRows()<dim);

  m2.resize(m.numRows(),m.numColumns());

  for (size_t i = 0; i < m.numRows(); ++i)
    for (size_t j = 0; j < m.numColumns(); ++j)
      m2.set_mod2(i,j,m(i,j));
}


/*!
  In this template, we assume that |I|, |J| and |O| are respectively input and
  output iterators for type |Weight|, and that |[firstb, lastb[| holds a new
  $\Q$-basis for the current lattice, expressed in terms of the current basis.
  As we range from begin to end, we write the vectors in the new basis and
  output them to |O|.

  Doing the base change amounts to multiplying with the inverse matrix of
  |b|'s matrix.

  NOTE: we don't assume that |[firstb, lastb[| is necessarily a $\Z$-basis of
  the current lattice, only that it is a basis of a full rank sublattice
  containing the vectors in the input range; the new coordinates will then be
  integers. Users should be aware of the "full rank" condition; without it the
  specification still makes sense, but the implementation will fail.
*/
template<typename I, typename J, typename O>
  void baseChange(const I& first, const I& last, O out, const J& firstb,
		  const J& lastb)
{
  latticetypes::LatticeCoeff d;
  latticetypes::LatticeMatrix q =
    latticetypes::LatticeMatrix(firstb,lastb,tags::IteratorTag()).inverse(d);

  for (I i = first; i != last; ++i)
  {
    latticetypes::Weight v = q.apply(*i);
    v/=d;
    *out++ = v;
  }
}


/*! Like baseChange, but we go from a list expressed in terms of |[firstb,
  lastb[| to list expressed in terms of the original basis. This is actually
  easier, as we don't have to invert the base change matrix!
*/
template<typename I, typename J, typename O>
  void inverseBaseChange(const I& first, const I& last, O out, const J& firstb,
			 const J& lastb)
{
  latticetypes::LatticeMatrix q(firstb,lastb,tags::IteratorTag());

  for (I i = first; i != last; ++i, ++out) {
    latticetypes::Weight v(q.numRows());
    q.apply(v,*i);
    *out = v;
  }
}

} // namespace lattice

} // namespace atlas
