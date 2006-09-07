/*!
\file
\brief Declarations and definitions for templates for BitVector.

This is a vector in the first d_size coordinates of the vector space
(Z/2Z)^dim over the two-element field.  The software envisions dim
between 0 and four times the machine word length (precisely, four
times the constant longBits, which is the number of bits in an
unsigned long integer).

*/
/*
  This is bitvector_def.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include <algorithm>

#include "comparison.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

/*****************************************************************************

        Chapter I -- The BitVector class

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace bitvector {


/*!
\brief Adds b to the bitvector (as the first coordinate), increasing the
  size by one. 

  It is the user's responsibility to make sure that the size does not
  exceed dim.
*/
template<size_t dim> BitVector<dim>& BitVector<dim>::pushBack(bool b)

{
  if (b)
    d_data.set(d_size);

  ++d_size;

  return *this;
}

/*!
\brief Extracts the bits flagged by t, and packs them consecutively
  while preserving their values. The remaining bits are set to zero.

  This is used to express vectors in the basis of a subspace.  A
  subspace of (Z/2Z)^d_size has a unique basis in row-reduced form;
  such a basis is carried by NormalSubspace.  Suppose that t flags the
  leading bits of this row-reduced basis, and that BitVector belongs
  to the subspace.  Then slice(t) will change the coordinates of
  BitVector to those with respect to the row-reduced basis of the subspace.
*/
template<size_t dim> 
void BitVector<dim>::slice(const bitset::BitSet<dim>& t)

{
  size_t c = 0;

  for (size_t j = 0; j < size(); ++j)
    if (t.test(j)) { // insert bit j at position c
      d_data.set(c,d_data[j]);
      ++c;
    }

  for (size_t j = c; j < size(); ++j) // set to zero to be in the clean
    reset(j);

  d_size = c;

  return;
}

}

/*****************************************************************************

        Chapter II -- The BitMatrix class

  ... explain here when it is stable ...

******************************************************************************/

namespace bitvector {

template<size_t dim> 
BitMatrix<dim>::BitMatrix(const std::vector<BitVector<dim> >& b)
  :d_columns(b.size())

/*!
  Constructs the matrix whose columns are given by the vectors in b.

  NOTE : it is assumed that all the vectors in b have the same size.
*/

{
  d_data.reserve(b.size());

  for (size_t j = 0; j < b.size(); ++j)
    d_data.push_back(b[j].data());

  // d_rows is undefined if b.size is zero

  if (b.size())
    d_rows = b[0].size();
}

/******** accessors *********************************************************/

template<size_t dim> 
BitVector<dim>& BitMatrix<dim>::apply(BitVector<dim>& dest, 
				      const BitVector<dim>& source) const

/*!  
  Applies the BitMatrix to the BitVector source and puts the result
  in the BitVector dest. It is assumed that d_columns is equal to
  source.size().
*/

{
  BitVector<dim> v(source); // for the case source = dest

  dest.resize(d_rows);
  combination(dest.d_data,d_data,v.data());

  return dest;
}

template<size_t dim> template<typename I, typename O> 
void BitMatrix<dim>::apply(const I& first, const I& last, O out) const

/*!
  \brief A pipe-version of apply. 
  
  We assume that I is an InputIterator with
  value-type BitVector<dim>, and O an OutputIterator with the same
  value-type. Then we apply our matrix to each vector in [first,last[
  and output it to out.
*/

{
  for (I i = first; i != last; ++i) {
    BitVector<dim> v;
    apply(v,*i);
    *out++ = v;
  }
   
  return;
}

template<size_t dim> void BitMatrix<dim>::column(BitVector<dim>& c, size_t j) 
  const

/*!
  \brief Puts in c the j-th column of the matrix.

  NOTE: we cannot simply return d_data[j] because that would be a BitSet,
  not a BitVector!
*/

{
  c = BitVector<dim>(d_data[j],d_rows);
  return;
}


/*!
  \brief Puts in b the normalized basis of the image of the matrix.
*/
template<size_t dim>
void BitMatrix<dim>::image(std::vector<BitVector<dim> >& b) const

{
  std::vector<BitVector<dim> > b1(d_columns);

  for (size_t j = 0; j < d_columns; ++j)
    column(b1[j],j);  //puts column \#j of the matrix in entry \#j of b1

  bitset::BitSet<dim> t;
  normalize(t,b1);

  b.swap(b1);

  return;
}

template<size_t dim>
void BitMatrix<dim>::kernel(std::vector<BitVector<dim> >& b) const

/*!
  \brief Puts in b the standard basis for the kernel of the matrix. 

  This is essentially normalizing the transpose matrix.
*/

{
  b.clear();

  if (isEmpty())
    return;

  std::vector<BitVector<dim> > eqn(d_rows);

  // "transpose" the matrix in eqn
  
  for (size_t i = 0; i < d_rows; ++i) {
    eqn[i].resize(d_columns);
    for (size_t j = 0; j < d_columns; ++j)
      if (test(i,j))
	eqn[i].set(j);
  }

  // normalize eqn
      
  bitset::BitSet<dim> t; // will flag a subset of [0,d_columns[

  normalize(t,eqn);

  // now we get a relation for each j in [0,d_columns[ not flagged by t
  // the relations are of the form e_j = sum a_{i,j}e_i, j not in t,
  // i in t, where the a_{i,j} are given by the rows in eqn

  for (size_t j = 0; j < d_columns; ++j)
    if (!t.test(j)) {
      BitVector<dim> v(d_columns,j);
      size_t c = 0;
      for (size_t i = 0; i < d_columns; ++i)
	if (t.test(i)) {
	  if (eqn[c].test(j))
	    v.set(i);
	  ++c;
	}
      b.push_back(v);
    }

  return;
}

template<size_t dim> void BitMatrix<dim>::row(BitVector<dim>& r, size_t i) 
  const

/*!
  Puts in r the i-th row of the matrix.
*/

{
  r.reset();
  r.resize(d_columns);

  for (size_t j = 0; j < d_columns; ++j)
    if (test(i,j))
      r.set(j);

  return;
}

/******** manipulators *******************************************************/

template<size_t dim>
BitMatrix<dim>& BitMatrix<dim>::operator+= (const BitMatrix<dim>& m)

/*!
  \brief Increment *this by adding m.

  Precondition: m has the same size as the current matrix.
*/

{
  for (size_t j = 0; j < d_columns; ++j)
    d_data[j] ^= m.d_data[j];

  return *this;
}

template<size_t dim>
BitMatrix<dim>& BitMatrix<dim>::operator*= (const BitMatrix<dim>& m)

/*!
  \brief Increment through right multiplication by m.

  As in apply, this can be done rather efficently by computing a whole column 
  at a time.

  NOTE : of course m.d_rows must be equal to d_columns.
*/

{
  BitMatrix<dim> res(d_rows,m.d_columns);

  for (size_t j = 0; j < m.d_columns; ++j) {
    const bitset::BitSet<dim>& m_j = m.d_data[j];
    bitset::BitSet<dim>& res_j = res.d_data[j];
    combination(res_j,d_data,m_j);
  }

  swap(res);
  
  return *this;
}

template<size_t dim> template<typename I> 
void BitMatrix<dim>::addColumns(const I& first, const I& last)

/*!
  \brief Appends to *this new columns which are the values of I.

  In this template we assume that I is an InputIterator whose value_type
  is BitVector<dim>, and which outputs BitVectors of size d_rows. Then
  we add one column to the matrix for each i in [first,last[.

  NOTE : this can be done very efficiently because our representation
  is in terms of column vectors.
*/

{
  size_t c = 0;

  for (I i = first; i != last; ++i) {
    d_data.push_back((*i).data());
    ++c;
  }

  d_columns += c;

  return;
}

template<size_t dim> template<typename I> 
void BitMatrix<dim>::addRows(const I& first, const I& last)

/*!
  \brief Appends to *this new rows which are the values of I.

  In this template we assume that I is an InputIterator whose value_type
  is BitVector<dim>, and which outputs BitVectors of size d_columns. Then
  we add one row to the matrix for each i in [first,last[.

  NOTE : it is the caller's responsibility to make sure that the number
  of rows does not exceed dim.
*/

{
  size_t c = numRows;

  for (I i = first; i != last; ++i) {
    for (size_t j = 0; j < d_columns; ++j)
      if ((*i).test(j))
	set(c,j);
    ++c;
  }

  numRows = c;

  return;
}

template<size_t dim> void BitMatrix<dim>::cutRows(size_t c)

/*!
  \brief Removes the first c rows of the matrix, by right-shifting the data.
*/

{
  for (size_t j = 0; j < d_columns; ++j)
    d_data[j] >>= c;

  d_rows -= c;
  
  return;
}

template<size_t dim> BitMatrix<dim>& BitMatrix<dim>::invert()

/*!
  \brief Replaces the current matrix by its inverse. 
  
  It is the caller's responsibility to make sure that m is in fact
  invertible (in particular, that it is square).  If necessary, this
  may be done by a call to isInvertible().

  For the algorithm, we use the normalSpanAdd function. We start out with
  a matrix of size (2r,c) (if r,c is the size of our original matrix) 
  containing the identity matrix below the given one. Then we apply
  normalSpanAdd; at the end, the upper part of our matrix is (some permutattion
  of) the identity. Setting it right will finish the job.

  Actually, since this would require working with BitMatrix<2*dim>, and since
  conversion functions don't seem to be forthcoming, we work with a pair
  of matrices, and just copy the code from normalSpanAdd.
*/

{
  BitMatrix<dim> i(d_rows);

  for (size_t j = 0; j < d_rows; ++j)
    i.set(j,j);

  std::vector<size_t> f;

  for (size_t k = 0; k < d_rows; ++k) {// add column k

    for (size_t j = 0; j < k; ++j)
      if (test(k,f[j])) { // set bit f[j] of data[k] to zero
	d_data[k] ^= d_data[j];
	i.d_data[k] ^= i.d_data[j];
      }

    // adjust the basis a

    size_t n = firstBit(d_data[k]);

    for (size_t j = 0; j < k; ++j)
      if (d_data[j].test(n)) {
	d_data[j] ^= d_data[k];
	i.d_data[j] ^= i.d_data[k];
      }

    f.push_back(n);
  }

  // write the appropriate column-permutation of i in the current matrix

  for (size_t j = 0; j < d_rows; ++j)
    d_data[f[j]] = i.d_data[j];

  return *this;
}

template<size_t dim> void BitMatrix<dim>::reset()

/*!
  \brief Resets the matrix to zero.
*/

{
  for (unsigned long j = 0; j < d_data.size(); ++j)
    d_data[j].reset();

  return;
}

template<size_t dim> void BitMatrix<dim>::resize(size_t m, size_t n)

/*!
  \brief Resizes the matrix to m rows, n columns.

  NOTE : it is the caller's responsibility to check that m does not exceed
  dim.
*/

{
  d_data.resize(n);

  d_rows = m;
  d_columns = n;

  return;
}

template<size_t dim> void BitMatrix<dim>::swap(BitMatrix<dim>& m)

/*!
  \brief Swaps contents with m.
*/

{
  d_data.swap(m.d_data);

  size_t tmp = d_rows;
  d_rows = m.d_rows;
  m.d_rows = tmp;

  tmp = d_columns;
  d_columns = m.d_columns;
  m.d_columns = tmp;

  return;
}

template<size_t dim> BitMatrix<dim>& BitMatrix<dim>::transpose()

/*!
  \brief Transposes the matrix. 

  This could have a very simple implementation for square matrices,
  but is quite a bit trickier for rectangular ones!  So we don't try
  to be smart and do the safe thing, which is to make a copy.
*/

{
  BitMatrix<dim> result(d_columns,d_rows);

  for (size_t i = 0; i < d_rows; ++i)
    for (size_t j = 0; j < d_columns; ++j)
      if (test(i,j))
	result.set(j,i);

  d_data.swap(result.d_data);
  d_rows = result.d_rows;
  d_columns = result.d_columns;

  return *this;
}

}

/*****************************************************************************

        Chapter III -- Functions defined in bitvector.h

  ... explain here when it is stable ...

******************************************************************************/

namespace bitvector {

template<size_t dim> 
  void combination(BitVector<dim>& v, const std::vector<BitVector<dim> >& b,
		   const bitset::BitSet<dim>& e)

/*!
  \brief Puts in v the linear combination of the elements of b given by e.

  NOTE : it is the caller's responsibility to check that v has the correct
  size.
*/

{
  v.reset();

  for (size_t i = 0; i < b.size(); ++i)
    if (e.test(i))
      v += b[i];

  return;
}

template<size_t dim> 
  void combination(bitset::BitSet<dim>& v, 
		   const std::vector<bitset::BitSet<dim> >& b,
		   const bitset::BitSet<dim>& e)

/*!
  \brief Puts in v the linear combination of the elements of b given by e.

  NOTE : it is the caller's responsibility to check that v has the correct
  size.
*/

{
  v.reset();

  for (size_t i = 0; i < b.size(); ++i)
    if (e.test(i))
      v ^= b[i];

  return;
}

template<size_t dim>
  void complement(bitset::BitSet<dim>& c, 
		  const std::vector<BitVector<dim> >& b,
		  size_t d)

/*!
  \brief Puts in c a subset of the canonical basis which spans a complementary
  subspace to the subspace spanned by b.

  It is assumed that the vectors in b are all of the same size, but not
  necessarily independent.

  NOTE : we need to pass the dimension in case b is empty; we don't want
  to set bits beyond that (so that s.count(), for instance, yields the
  correct dimension of the complement.)
*/

{
  std::vector<size_t> f;
  std::vector<BitVector<dim> > a;

  for (unsigned long j = 0; j < b.size(); ++j)
    spanAdd(a,f,b[j]);

  // now f contains the indices we _don't_ want

  for (size_t j = 0; j < d; ++j)
    c.set(j);

  for (size_t j = 0; j < f.size(); ++j)
    c.reset(f[j]);

  return;
}

template<size_t dim> 
  bool firstSolution(bitset::BitSet<dim>& c, 
		     const std::vector<BitVector<dim> >& b,
		     const BitVector<dim>& rhs)

/*!
  \brief Puts in c a solution of the system whose columns are b, with the
  given right-hand side.

  Return value is true if there is a solution, false if there is none; in that
  case c is unchanged.
*/

{
  using namespace bitset;

  if (b.size() == 0) {
    if (rhs.isZero()) {
      c.reset();
      return true;
    }
    else
      return false;
  }

  size_t n = b[0].size();
  std::vector<BitSet<dim> > a;
  std::vector<size_t> f;
  BitSet<dim> rh;
  
  for (size_t i = 0; i < n; ++i) {
    BitSet<dim> r;
    // set r to i-th row of transpose matrix
    for (size_t j = 0; j < b.size(); ++j)
      if (b[j].test(i))
	r.set(j);
    bool x = rhs[i];
    // normalize r w.r.t. a
    for (size_t j = 0; j < f.size(); ++j)
      if (r.test(f[j])) {
	r ^= a[j];
	x ^= rh[j];
      }
    // normalize the elements of a if r is new, and add r to a
    if (r.any()) {
      size_t m = r.firstBit();
      for (size_t j = 0; j < a.size(); ++j)
	if (a[j].test(m)) {
	  a[j] ^= r;
	  rh.set(j,rh[j]^x);
	}
      rh.set(a.size(),x);
      a.push_back(r);
      f.push_back(m);
    }
    else if (x)
      return false;
  }

  c.reset();

  for (size_t j = 0; j < f.size(); ++j)
    if (rh[j])
      c.set(f[j]);

  return true;
}

template<size_t dim> 
bool firstSolution(BitVector<dim>& sol,
		   const std::vector<BitVector<dim> >& d_eqn)

/*!
  \brief Solves the system of equations eqn.
  
   eqn holds a system of equations, right-hand side included.
  If the system has at least one solution, we put one in sol, and we return
  true. Otherwise, sol is unchanged, and we return false.
*/

{
  std::vector<BitVector<dim> > eqn(d_eqn); // local copy
  bitset::BitSet<dim> t;

  normalize(t,eqn);

  // now the system is impossible if there is an equation for which the
  // last bit in t is set

  if (eqn.size() and (t.test(eqn[0].size()-1)))
    return false;

  sol.reset();

  if (eqn.size() == 0) // do nothing
    return true;
  
  sol.resize(eqn[0].size()-1); // sol looks only at the lhs

  // we only use the bits flagged by t and set them to the corresponding rhs

  size_t c = 0;
  const size_t rhs = sol.size();

  for (size_t j = 0; j < rhs; ++j)
    if (t.test(j)) {
      if (eqn[c].test(rhs))
	sol.set(j);
      ++c;
    }

  return true;
}

template<size_t dim> void identityMatrix(BitMatrix<dim>& m, size_t n)

/*!
  \brief Puts in m the identity matrix in rank n.

  Precondition: n <= dim;
*/

{
  m.resize(n,n);
  m.reset();

  for (size_t j = 0; j < n; ++j)
    m.set(j,j);

  return;
}

template<size_t dim> void initBasis(std::vector<BitVector<dim> >& b, size_t n)

/*!
  \brief Initializes b to the canonical basis in dimension n.
*/

{
  b.assign(n,BitVector<dim>(n));

  for (size_t j = 0; j < n; ++j)
    b[j].set(j);

  return;
}

template<size_t dim> bool isIndependent(const std::vector<BitVector<dim> >& b)

/*!
  \brief Tells whether the system of bitvectors is independent.
*/

{  
  std::vector<BitVector<dim> > a;
  std::vector<size_t> f;

  for (size_t j = 0; j < b.size(); ++j) {
    spanAdd(a,f,b[j]);
    if (a.size() == j) // first time a dependent vector is found
      return false;
  }

  return true;
}


/*!
  \brief Replaces b by the normal basis of the vector space V it spans. Flags
  in t the leading bits of the normal basis.

  What is flagged in t is the set J in normalSpanAdd.  Those are the
  indices one has to look at to find the coordinates of an element of
  V in the normal basis.)

  NOTE : as we wish to use t to obtain the coordinates in the
  normalized basis by just "slicing", we need to be careful to reorder
  the basis vectors according to their first bit.  Therefore t flags
  the leading bits of the normal basis; every other basis vector has a
  zero in such a leading bit; and the locations of the leading bits
  increase.
*/
template<size_t dim>
  void normalize(bitset::BitSet<dim>& t, std::vector<BitVector<dim> >& b)

{
  using namespace comparison;

  std::vector<BitVector<dim> > a;
  std::vector<size_t> f;

  for (size_t j = 0; j < b.size(); ++j)
    normalSpanAdd(a,f,b[j]);

  t.reset();

  for (size_t j = 0; j < f.size(); ++j)
    t.set(f[j]);

  // reorder the basis elements according to their first bit
  typename std::vector<BitVector<dim> >::iterator first = a.begin();
  typename std::vector<BitVector<dim> >::iterator last = a.end();

  std::sort(first,last,compare(FirstBit<dim>()));

  // commit
  b.swap(a);

  return;
}


/*!  
  \brief Transforms the unordered normal basis a into an unordered
  normal basis for the span of a and v.

  For each subvectorspace V of k^d, and a subset I of {0,...,d-1} such that
  the standard basis vectors flagged by I generate a complementary subspace
  to V, the normal basis of V corresponding to I is the set of vectors in
  V of the form e_j + sum_{i in I}a_{i,j}e_i, for j not in I (in other
  words, we see V as the graph of a map from k^J to k^I, where J is the 
  complement of I.) If moreover we choose I to be as big as possible
  lexicographically (i.e., we choose J to be as small as possible) we get a 
  canonical basis for each subvectorspace of V. 

  This function assumes that a already contains the canonical basis of some 
  subspace, and we add v to the subspace. Then a is not modified if v lies 
  in the subspace generated by a. Otherwise a new element is added, and the
  existing vectors in a are modified accordingly.

  The vector f records the elements of J corresponding to the elements of
  a, in their order of appearance. Note that we refrain from keeping a
  ordered.
*/
template<size_t dim> 
  void normalSpanAdd(std::vector<BitVector<dim> >& a, std::vector<size_t>& f,
		     const BitVector<dim>& v)

{
  if (v.isZero()) // v is the zero vector
    return;

  // reduce v modulo a

  BitVector<dim> w = v;

  for (size_t j = 0; j < a.size(); ++j)
    if (w.test(f[j])) // subtract a[j] from w
      w -= a[j];

  if (w.isZero()) // v is in span of a
    return;

  // adjust the basis a

  size_t n = w.firstBit();

  for (size_t j = 0; j < a.size(); ++j)
    if (a[j].test(n))
      a[j] -= w;

  f.push_back(n);
  a.push_back(w);

  return;
}

template<size_t dim>
  void projection(BitMatrix<dim>& p, const std::vector<BitVector<dim> >& b, 
		  size_t d)

/*!  
  \brief Puts in p the matrix of the projection on the canonical
  complement to the span of b.

  We assume that b holds the normal basis for the
  subspace V that it spans (if not, this can be obtained by a call to
  normalize.) This function then puts in p the matrix of the projection 
  to the canonical complement of V, parallel to V.

  NOTE : we need to pass the dimension in case b is empty.
*/

{
  using namespace bitset;

  // flag the indices for the canonical complement

  BitSet<dim> c;
  std::vector<size_t> f;

  bitset::set(c,d);

  for (size_t j = 0; j < b.size(); ++j) {
    size_t fb = b[j].firstBit();
    f.push_back(fb);
    c.reset(f[j]);
  }

  // the coordinates of the basis vectors corresponding to indices in s
  // give the projection matrix

  p.resize(d-b.size(),d);
  p.reset();

  size_t k = 0;

  for (size_t j = 0; j < d; ++j)
    if (c.test(j)) {
      p.set(k,j);
      ++k;
    }

  for (size_t j = 0; j < b.size(); ++j) {
    BitVector<dim> v = b[j];
    v.slice(c);
    const BitSet<dim>& vd = v.data(); 
    p.setColumn(f[j],vd);
  }

  return;
}

template<size_t dim>
  void reflectionMatrix(BitMatrix<dim>& m, const BitVector<dim>& a, 
			const BitVector<dim>& a_check)

/*!
  \brief Puts in m the matrix of the reflection defined by a and a_check.

  Precondition: a and a_check have same size; <a,a_check> = 0;

  This is the operator x -> x + <x,a_check>a (we can write + because we are
  in characteristic two.)
*/

{
  identityMatrix(m,a.size());

  for (size_t j = 0; j < a.size(); ++j)
    if (a_check.test(j))
      m.addToColumn(j,a);

  return;
}

template<size_t dim>
  void relations(std::vector<BitVector<dim> >& rel, 
		 const std::vector<BitVector<dim> >& b)

/*!
  \brief Writes in r the relations among the elements in b. 

  In other words, it solves the system of equations defined by the _rows_ in 
  the matrix whose columns are given by b.
*/

{
  using namespace bitset;

  rel.resize(0);

  if (b.size() == 0) // do nothing
    return;

  size_t r = b.size();       // dimension of the source space
  size_t d = b[0].size();    // dimension of the target space
  std::vector<BitVector<dim> > eqn(d);

  // "transpose" b in e
  
  for (size_t i = 0; i < d; ++i) {
    eqn[i].resize(b.size());
    for (size_t j = 0; j < b.size(); ++j)
      if (b[j].test(i))
	eqn[i].set(j);
  }

  // normalize e
      
  BitSet<dim> t; // will flag a subset of [0,r[
  normalize(t,eqn);

  // now we get a relation for each j in [0,d[ not flagged by t
  // the relations are of the form e_j = sum a_{i,j}e_i, j not in t,
  // i in t, where the a_{i,j} are given by the columns in eqn

  for (size_t j = 0; j < r; ++j)
    if (!t.test(j)) {
      BitVector<dim> v(r,j);
      size_t c;
      for (size_t i = 0; i < r; ++i)
	if (t.test(i)) {
	  if (eqn[c].test(j))
	    v.set(i);
	  ++c;
	}
      rel.push_back(v);
    }

  return;
}

template<size_t dim> bool scalarProduct(const BitVector<dim>& vd, 
					const BitVector<dim>& v)

/*!
  \brief Applies the dual vector vd to v.
*/

{
  BitVector<dim> w = v;
  w &= vd;

  return w.count()&1ul;
}

template<size_t dim> void spanAdd(std::vector<BitVector<dim> >& a,
				  std::vector<size_t>& f,
				  const BitVector<dim>& v)

/*!
  \brief Enlarges the basis a to span v.

  It is assumed that a contains a list of independent bitvectors of size
  dim; v is added to the list if it is independent, discarded otherwise.

  NOTE : we will assume that the first set bits of the elements in a are
  all distinct, and that the bits in a given vector in a corresponding
  to first bits of previous vectors are all unset; this then makes it
  easy to check linear dependence. The vector f holds the positions of
  these first bits.
*/

{
  if (v.isZero()) // v is the zero vector
    return;

  // reduce v modulo a

  BitVector<dim> w = v;

  for (unsigned long j = 0; j < a.size(); ++j)
    if (w.test(f[j])) // subtract a[j] from w
      w -= a[j];

  if (w.isZero()) // v is in span of a
    return;

  f.push_back(w.firstBit());
  a.push_back(w);

  return;
}

}

}
