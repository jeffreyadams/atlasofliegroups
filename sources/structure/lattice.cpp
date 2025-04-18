/*
  This is lattice.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007--2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  Implementation for namespace lattice.

  This module defines some more general lattice functions.
*/

#include "lattice.h"
#include "bitmap.h" // for (unused) result of |column_echelon|
#include "bigint.h" // using a variable of this type

#include <cassert>

#include "tags.h"
#include "matrix.h"
#include "matreduc.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include <iterator>
#endif

/*****************************************************************************

  This module defines some more general lattice functions.

******************************************************************************/

namespace atlas {

namespace lattice {

 /*
   Functions for working with lattices.

   This namespace defines some general lattice functions. It includes change of
   basis functions, and calculating the orthogonal of a sublattice.
 */

/*****************************************************************************

        Chapter I -- Functions declared in lattice.h

******************************************************************************/


/*
  In this template, we assume that |I|, and |O| are respectively random
  access input and output iterator types for type |Weight|, and that
  |[firstb,lastb[| holds a new $\Q$-basis for the lattice, in particular that
  |lastb-firstb| is equal to the size of the |Weight|s.

  As we iterate from |first| to |last|, we write the vectors in the
  new basis (this is supposed to be possible) and output the result to |O|.

  Doing the base change amounts to applying the inverse of |b|'s matrix.

  We don't assume that |[firstb, lastb[| is necessarily a $\Z$-basis of $\Z^n$,
  but we do assume that it is a basis of a _full rank_ sublattice, and that
  the vectors in |[first,last[| also lie in that range; the new coordinates will
  then be integers. Users should be aware of the "full rank" condition; without
  it the specification would still make sense, but the implementation will fail.
*/
template<typename I, typename O>
  void baseChange(I first, I last, O out, I firstb, I lastb)
{
  arithmetic::big_int d;
  LatticeMatrix q =
    LatticeMatrix(firstb,lastb,lastb-firstb,tags::IteratorTag())
    .inverse(d);

  while (first!=last)
  {
    *out = q*(*first)/d.convert<LatticeCoeff>();
    ++out, ++first;
  }
}

/*
  This (unused) template function is like |baseChange|, but goes from weights
  expressed in terms of |[firstb, lastb[| to ones expressed in terms of the
  original basis. This is easier, as we don't have to invert the matrix!
*/
template<typename I, typename O>
  void inverseBaseChange(I first, I last, O out, I firstb, I lastb)
{
  LatticeMatrix q(firstb,lastb,lastb-firstb,tags::IteratorTag());

  while (first!= last)
  {
    *out = q*(*first);
    ++out, ++first;
  }
}


/*
  Return a basis of the orthogonal of the sublattice generated by the columns
  of |b| in $\Z^r$, where |r| is the number of rows of |b|.

  Algorithm: diagonalise the matrix |b| using row operations |R| (forgetting
  the column operations |C| used) as for Smith form; then image |R*M| is the
  span of the original images of canonical basis vectors (as many as nonzero
  factors), and orthogonal sublattice is spanned by the remaining rows of |R|.

*/
CoweightList perp(const LatticeMatrix& b)
{
  const auto r = b.n_rows();
  LatticeMatrix R,C;
  CoeffList factor = matreduc::diagonalise(b,R,C);

  // drop any final factors 0
  size_t codim=factor.size();

  // collect final rows of |R|, those that annihilate the span of |b|
  CoweightList result; result.reserve(r-codim);
  for (size_t i=codim; i<r; ++i)
    result.push_back(R.row(i));

  return result;
}

LatticeMatrix kernel (LatticeMatrix M)
{
  size_t m= M.n_columns(); // dimension of space to which |M| can be applied
  LatticeMatrix C; bool flip; // flip is dummy
  matreduc::column_echelon(M,C,flip);

  return C.block(0,M.n_columns(),m,m); // final columns span $\ker(M)$
}

LatticeMatrix eigen_lattice (LatticeMatrix M, LatticeCoeff lambda)
{
  return kernel(M-=lambda);
}

LatticeMatrix row_saturate(const LatticeMatrix& M)
{
  size_t n= M.n_columns(); // dimension of space to which |M| can be applied
  CoeffList factor;
  LatticeMatrix b = matreduc::adapted_basis(M.transposed(),factor);

  size_t c=factor.size(); // rank of M (codimension of kernel)

  LatticeMatrix result(c,n); // initial rows of |b.transposed()|
  for (size_t i=0; i<c; ++i)
    result.set_row(i,b.column(i));
  return result;
}


//			    Template instantiation

template
void baseChange
  (WeightList::iterator,
   WeightList::iterator,
   std::back_insert_iterator<WeightList>,
   WeightList::iterator,
   WeightList::iterator);


} // |namespace lattice|

} // |namespace atlas|
