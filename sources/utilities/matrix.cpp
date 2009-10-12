/*!
\file
  This is matrix.cpp. This module contains some simple utilities for matrices.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "matrix.h"

#include <cassert>
#include <limits>
#include <algorithm>
#include "matreduc.h"
#include "intutils.h"
#include "arithmetic.h"


/*****************************************************************************

  This module contains some simple utilities for dense matrices.
  Code is formulated as templates, but only a few instances are generated.

******************************************************************************/

namespace atlas {

namespace matrix {

/*****************************************************************************

        Chapter II -- The Vector class

  This class template adds some arithmetic operations to |std::vector|

******************************************************************************/

//! \brief Adds |v|.
template<typename C>
Vector<C>& Vector<C>::operator+= (const Vector<C>& v)
{
  for (size_t i=0; i<base::size(); ++i)
    (*this)[i] += v[i]; // one may write |base::operator[](i)| for |(*this)[i]|
  return *this;
}

//! \brief Subtracts |v|.
template<typename C>
Vector<C>& Vector<C>::operator-= (const Vector<C>& v)
{
  for (size_t i=0; i<base::size(); ++i)
    (*this)[i] -= v[i];
  return *this;
}

//! \brief Scalar multiplies by |c|
template<typename C>
Vector<C>& Vector<C>::operator*= (C c)
{
  for (size_t i=0; i<base::size(); ++i)
    (*this)[i] *= c;
  return *this;
}

/*! \brief Scalar divides by |c|

  All entries must allow exact division, if not a std::runtime_error is thrown
*/
template<typename C>
Vector<C>& Vector<C>::operator/= (C c) throw (std::runtime_error)
{
  for (typename Vector<C>::iterator it=base::begin(); it!=base::end(); ++it)
    if (*it%c==0)
      *it/=c;
    else throw std::runtime_error("Inexact integer division");
  return *this;
}

template<typename C>
Vector<C>& Vector<C>::negate ()
{
  for (size_t i=0; i<base::size(); ++i)
    (*this)[i] = -(*this)[i];
  return *this;
}

template<typename C>
  C Vector<C>::scalarProduct (const Vector<C>& v) const
{
  C result= C(0);
  for (size_t i=0; i<base::size(); ++i)
    result += (*this)[i] * v[i];

  return result;
}

template<typename C>
  bool Vector<C>::isZero() const
{
  for (size_t i=0; i<base::size(); ++i)
    if ((*this)[i]!=0)
      return false;
  return true;
}

template<typename C>
  Matrix<C> Vector<C>::row_matrix() const
{
  Matrix<C> result(1,base::size());
  for (size_t j=0; j<base::size(); ++j)
    result(0,j)=base::operator[](j);

  return result;
}

template<typename C>
  Matrix<C> Vector<C>::column_matrix() const
{
  Matrix<C> result(base::size(),1);
  for (size_t i=0; i<base::size(); ++i)
    result(i,0)=base::operator[](i);

  return result;
}


} // namespace matrix

/*****************************************************************************

        Chapter II -- The Matrix class

  We implement a matrix as a vector of elements, concatenating the rows.

******************************************************************************/

namespace matrix {

/******** constructors *******************************************************/

//! \brief: Construct the identity matrix of size n.
template<typename C> Matrix<C>::Matrix(size_t n)
  : d_rows(n),d_columns(n),d_data(n*n,C(0))
{
  for (size_t i = 0; i<n; ++i)
    operator()(i,i) = C(1);
}

/*! \brief
  This constructor constructs a matrix from a list of vectors, column-wise.
  It is assumed that all elements of |b| (possibly none) have size |n_rows|.
*/
template<typename C>
Matrix<C>::Matrix(const std::vector<Vector<C> >& b, size_t n_rows)
  : d_rows(n_rows), d_columns(b.size()), d_data(d_rows*d_columns)
{
  for (size_t j = 0; j<d_columns; ++j)
    set_column(j,b[j]);
}

/*!
  Here I is a random access iterator type whose value type should be Vector<C>
  or some other type that can be subscripted giving a value of type C.
  We read the columns of the matrix from the iterator.

  All the vectors obtained from the iterators should have size |n_rows|
*/
template<typename C>
  template<typename I> Matrix<C>::Matrix
    (I first, I last, size_t n_rows, tags::IteratorTag)
  : d_rows(n_rows), d_columns(last-first), d_data(d_rows*d_columns)
{
  I p=first;
  for (size_t j=0; j<d_columns; ++j,++p)
    for (size_t i=0; i<d_rows; ++i)
      (*this)(i,j) = (*p)[i];
}

/*!
  Swaps m with the current matrix.
*/
template<typename C>
void Matrix<C>::swap(Matrix<C>& m)
{
  std::swap(d_rows,m.d_rows);
  std::swap(d_columns,m.d_columns);
  d_data.swap(m.d_data);
}

/******** accessors **********************************************************/

/*!
  Puts the i-th row of the matrix in v.
*/
template<typename C>
void Matrix<C>::get_row(Vector<C>& v, size_t i) const
{
  v.resize(d_columns);

  for (size_t j = 0; j<d_columns; ++j)
    v[j] = (*this)(i,j);
}
/*!
  Puts the j-th column of the matrix in v.
*/
template<typename C>
void Matrix<C>::get_column(Vector<C>& v, size_t j) const
{
  v.resize(d_rows);

  for (size_t i = 0; i<d_rows; ++i)
    v[i] = (*this)(i,j);
}

/*!
  Returns the list of row vectors of m.
*/
template<typename C>
  std::vector<Vector<C> > Matrix<C>::rows() const
{
  std::vector<Vector<C> > result(numRows());
  for (size_t i=0; i<numRows(); ++i)
    get_row(result[i],i);
  return result;
}

/*!
  Returns the list of column vectors of m.
*/
template<typename C>
  std::vector<Vector<C> > Matrix<C>::columns() const
{
  std::vector<Vector<C> > result(numColumns());

  for (size_t j = 0; j<result.size(); ++j)
    get_column(result[j],j);

  return result;
}

template<typename C>
bool Matrix<C>::operator== (const Matrix<C>& m) const
{
  return d_rows==m.d_rows and d_columns==m.d_columns and d_data==m.d_data;
}

/*! \brief
Applies the matrix to the vector w, and returns the result. It is assumed that
the size of w is the number of columns; result size is the number of rows.
*/
template<typename C>
Vector<C> Matrix<C>::apply(const Vector<C>& w) const
{
  Vector<C> result(d_rows);

  for (size_t i = 0; i<d_rows; ++i)
  {
    C c = 0;
    for (size_t j = 0; j<d_columns; ++j)
      c += (*this)(i,j) * w[j];
    result[i] = c;
  }

  return result;
}

/*! \brief
Multiplies the matrix to right to the row-vector w, and returns the result.
It is assumed that the size of w is the number of rows; result size is the
number of columns. This is the proper sense of application for dual space.
*/
template<typename C>
Vector<C> Matrix<C>::right_apply(const Vector<C>& w) const
{
  Vector<C> result(d_columns);

    for (size_t j = 0; j<d_columns; ++j)
  {
    C c = 0;
    for (size_t i = 0; i<d_rows; ++i)
      c += w[i] * (*this)(i,j);
    result[j] = c;
  }

  return result;
}

template<typename C>
Matrix<C> Matrix<C>::operator* (const Matrix<C>&  m) const
{
  assert(d_columns==m.d_rows);
  Matrix<C> result(d_rows,m.d_columns);

  for (size_t i = 0; i<d_rows; ++i)
    for (size_t k = 0; k<m.d_columns; ++k)
    {
      C c=0;
      for (size_t j = 0; j<d_columns; ++j)
	c += (*this)(i,j) * m(j,k);

      result(i,k)=c;
    }

  return result;
}

/*!
  This function constructs expresses our square matrix on the basis b.
*/
template<typename C>
  Matrix<C> Matrix<C>::on_basis(const std::vector<Vector<C> >& b) const
{
  assert (numRows()==numColumns());
  assert (b.size()==numRows());

  Matrix<C> p(b,b.size()); // square matrix
  C d;
  Matrix<C> result(p.inverse(d)* *this *p);
  return result /= d;
}

// manipulators


//! Puts v in the i-th row of the matrix
template<typename C>
  void Matrix<C>::set_row(size_t i, const Vector<C>& v)
{
  assert(v.size()==d_columns);

  for (size_t j = 0; j<d_columns; ++j)
    (*this)(i,j)=v[j];
}

//! Puts v in the j-th column of the matrix
template<typename C>
  void Matrix<C>::set_column(size_t j, const Vector<C>& v)
{
  assert(v.size()==d_rows);

  for (size_t i = 0; i<d_rows; ++i)
    (*this)(i,j)=v[i];
}

/*!
  Adds |v| as new row to the matrix. Entries being stored by row, this is easy
*/
template<typename C>
  void Matrix<C>::add_row(const Vector<C>& v)
{
  assert(v.size()==d_columns);

  ++d_rows;
  d_data.resize(d_rows,d_columns);

  std::copy(v.begin(),v.end(),&d_data[(d_rows-1)*d_columns]);
}

/*!
  Adds |v| as new column to the matrix. This is harder than adding a row
*/
template<typename C>
  void Matrix<C>::add_column(const Vector<C>& v)
{
  assert(v.size()==d_rows);

  size_t old_col=d_columns;
  ++d_columns;
  d_data.resize(d_rows,d_columns);

  typename Vector<C>::iterator dst=d_data.end(),src=dst-d_rows;
  for (size_t i=d_rows; i-->0;)
  {
    *--dst = v[i];
    for (size_t j=old_col; j-->0;)
      *--dst = *--src;
  }
}

/*!
  Incrementation by addition with m. It is assumed that m and *this
  have the same size.
*/
template<typename C>
Matrix<C>& Matrix<C>::operator+= (const Matrix<C>&  m)
{
  assert(numRows()==m.numRows());
  assert(numColumns()==m.numColumns());
  d_data += m.d_data;
  return *this;
}

/*!
  Incrementation by subtraction of m. It is assumed that m and *this
  have the same size.
*/
template<typename C>
Matrix<C>& Matrix<C>::operator-= (const Matrix<C>&  m)
{
  assert(numRows()==m.numRows());
  assert(numColumns()==m.numColumns());
  d_data -= m.d_data;
  return *this;
}


template<typename C>
Matrix<C>& Matrix<C>::operator/= (const C& c) throw (std::runtime_error)
{
  if (c != 1)
    d_data /= c;
  return *this;
}

/*!
  Transposes the matrix. We allow ourselves a temporary vector if the matrix
  is not square; it could likely be done in-place, but the permutation is
  rather complicated!
*/
template<typename C> void Matrix<C>::transpose()
{
  if (d_rows == d_columns) {// matrix is square
    for (size_t j = 0; j<d_columns; ++j)
      for (size_t i = j+1; i<d_rows; ++i) {
	C t = (*this)(i,j);
	(*this)(i,j) = (*this)(j,i);
	(*this)(j,i) = t;
      }
    return;
  }

  // now assume matrix is not square

  C zero = 0; // conversion from int
  Vector<C> v(d_data.size(),zero);
  size_t k = 0;

  // write data for transpose matrix in v

  for (size_t j = 0; j<d_columns; ++j)
    for (size_t i = 0; i<d_rows; ++i) {
      v[k] = (*this)(i,j);
      k++;
    }

  // swap data and row- and column- size

  d_data.swap(v);

  size_t t = d_rows;
  d_rows = d_columns;
  d_columns = t;
}


/*!
  Permutes rows and columns of the matrix according to the
  permutation a, resulting in the matrix of the same operator in the basis
  e_{a[0]}, ... , e_{a[n-1]}. This amounts to conjugating by the permutation
  matrix (delta_{i,a[i]})_{i,j} that transforms from coordinates on the
  standard basis to those on the basis e_{a[0]}, ... , e_{a[n-1]}.

  Precondition : m is a square matrix of order n; a holds a permutation
  of the matrix n;

  Method: the new entry at (i,j) is set to the old entry at (a[i],a[j]), in a
  separate copy (without trying to do the permutation of entries in place)
*/
template<typename C>
void Matrix<C>::permute(const setutils::Permutation& a)
{
  assert (a.size()==d_rows);
  assert (a.size()==d_columns);
  Matrix<C> q(d_rows,d_columns);

  for (size_t j = 0; j<d_columns; ++j)
    for (size_t i = 0; i<d_rows; ++i)
      q(i,j) = (*this)(a[i],a[j]);

  swap(q);
}


/*!
  Here we invert the matrix without catching the denominator. The intent
  is that it should be used for invertible matrices only.
*/
template<typename C> void Matrix<C>::invert()
{
  C d; invert(d);
  assert(d==1);
}

/*!
  This function inverts the matrix M. It is assumed that the coefficents
  are in an integral domain. At the conclusion of the process, d will contain
  the denominator for the inverted matrix (so that the true result is M/d).

  If the matrix is not invertible, the resulting matrix is undefined, and
  d is set to 0.

  We have chosen here to avoid divisions as much as possible. So in fact we
  take the matrix to Smith normal form through elementary operations; from
  there we deduce the smallest possible denominator for the inverse, and
  the inverse matrix. The difference with SmithNormal is that we have to
  keep track of both row and column operations.

  NOTE : probably the Smith normal stuff can be streamlined a bit so that
  functions from there can be called, instead of rewriting things slightly
  differently as we do here. In particular, the accounting of the row
  operations seems a bit different here from there.
*/
template<typename C>
void Matrix<C>::invert(C& d)
{
  assert(d_rows==d_columns);
  size_t n=d_rows;
  if (n==0) // do nothing to matrix, but set |d=1|
  { d=1; return; }

  Matrix<C> row,col;    // for recording column operations
  std::vector<C> diagonal = matreduc::diagonalise(*this,row,col);

  if (diagonal[n-1] == C(0)) // zero entries if any come at end
  { d=0; return; }

  d=1;
  for (size_t i=0; i<n; ++i)
    d=arithmetic::lcm(d,diagonal[i]);

  // finally for |D| "inverse" diagonal matrix w.r.t. |d|, compute |col*D*row|
  for (size_t j=0; j<n; ++j)
  {
    C f = d/diagonal[j];
    for (size_t i=0; i<n; ++i)
      (*this)(i,j)=f*col(i,j);
  }

  *this *= (row);
}

/*
 *
 Auxiliaries for |smithnormal|
 *
 */


/*!
  Whether all the entries in rectangle starting at |(i_min,i_max)| are zero.
*/
template<typename C>
bool Matrix<C>::isZero(size_t i_min, size_t j_min) const
{
  if (d_data.size() == 0) // empty matrix
    return true;

  for (size_t i = i_min; i<d_rows; ++i)
    for (size_t j = j_min; j<d_columns; ++j)
      if ((*this)(i,j))
	return false;

  return true;
}


/*!
  Returns the position of the smallest non-zero entry in absolute value,
  in the region i >= i_min, j >= j_min
*/
template<typename C>
typename Matrix<C>::index_pair Matrix<C>::absMinPos(size_t i_min,
						    size_t j_min) const
{
  C minCoeff = std::numeric_limits<C>::max();
  size_t i_m = d_rows;
  size_t j_m = d_columns;

  for (size_t i = i_min; i<d_rows; ++i)
    for (size_t j = j_min; j<d_columns; ++j)
    {
      C c = (*this)(i,j);
      if (c == 0)
	continue;
      if (intutils::abs(c)<minCoeff) // new smallest value
      {
	minCoeff = intutils::abs(c);
	i_m = i;
	j_m = j;
	if (minCoeff == 1)
	  break;
      }
    }

  return std::make_pair(i_m,j_m);
}



/*!
  Tells if all coefficients of the matrix are divisible by c.
*/
template<typename C>
bool Matrix<C>::divisible(C c) const
{
  for (size_t j = 0; j<d_data.size(); ++j)
    if (d_data[j]%c)
      return false;

  return true;
}


/*!
  Synopsis: constructs the matrix corresponding to the block [r_first,r_last[
  x [c_first,c_last[ of source. This uses that storage is by rows.
*/
template<typename C>
  Matrix<C> Matrix<C>::block(size_t i0, size_t j0, size_t i1, size_t j1) const
{
  Matrix<C> result(i1-i0,j1-j0);
  C* p = &result.d_data[0]; // writing pointer
  for (size_t i=i0; i<i1; ++i)
  {
    const C* q = &(*this)(i,j0);
    p = std::copy(q,q+result.numColumns(),p); // copy a row
  }
  return result;
}



/*!
  Carries out the row operation consisting of adding c times row j
  to row i.
*/
template<typename C>
void Matrix<C>::rowOperation(size_t i, size_t j, const C& c)
{
  for (size_t k = 0; k<rowSize(); ++k)
    (*this)(i,k) += c*(*this)(j,k);
}


/*!
  Carries out the column operation consisting of adding c times column k
  to column j.
*/
template<typename C>
void Matrix<C>::columnOperation(size_t j, size_t k, const C& c)
{
  for (size_t i = 0; i<columnSize(); ++i)
    (*this)(i,j) += c*(*this)(i,k);
}



/*!
  Changes the sign of all the entries in row i.
*/
template<typename C>
void Matrix<C>::rowMultiply(size_t i, C f)
{
  for (size_t j=0; j<d_columns; ++j)
    (*this)(i,j) *= f;
}

/*!
  Changes the sign of all the entries in column j.
*/
template<typename C>
void Matrix<C>::columnMultiply(size_t j, C f)
{
  for (size_t i = 0; i<d_rows; ++i)
    (*this)(i,j) *= f;
}

/*!
  Copies source to the rectangle of the current matrix with upper left corner
  at (r,c) and the appropriate size.
*/
template<typename C>
void Matrix<C>::copy(const Matrix<C>& source, size_t r, size_t c)
{
  for (size_t j = 0; j<source.d_columns; ++j)
    for (size_t i = 0; i<source.d_rows; ++i)
      (*this)(r+i,c+j) = source(i,j);
}




//! Interchanges rows i and j
template<typename C>
void Matrix<C>::swapRows(size_t i, size_t j)
{
  for (size_t k = 0; k<rowSize(); ++k)
    std::swap((*this)(i,k),(*this)(j,k));
}

//! Interchange columns i and j
template<typename C>
void Matrix<C>::swapColumns(size_t i, size_t j)
{
  for (size_t k = 0; k<columnSize(); ++k)
    std::swap((*this)(k,i),(*this)(k,j));
}

/*!
  Erases row i from the matrix.
*/
template<typename C>
void Matrix<C>::eraseRow(size_t i)
{
  typename Vector<C>::iterator first = d_data.begin() + i*d_columns;
  d_data.erase(first,first+d_columns);
  --d_rows;
}

/*!
  Removes column |j| altogether, shifting remaining entries to proper place
*/
template<typename C>
void Matrix<C>::eraseColumn(size_t j)
{
  typename Vector<C>::iterator pos =
     d_data.begin() + j; // position of entry |M(0,j)|
  --d_columns; // already adjust number of columns

  // kill individual entries from left to right, with shifts of new |d_colums|
  // (this is clever but not particularly efficient)
  for (size_t k = 0; k<d_rows; ++k, pos += d_columns)
    d_data.erase(pos);
}



} // namespace matrix


/*****************************************************************************

        Chapter III --- Functions declared in matrix.h

******************************************************************************/

namespace matrix {


/*!
  Synopsis: sets b to the canonical basis in dimension r.
*/
template<typename C>
  std::vector<Vector<C> > standard_basis(size_t r)
{
  std::vector<Vector<C> > result(r,Vector<C>(r,C(0)));

  for (size_t i = 0; i<r; ++i)
    result[i][i] = C(1);

  return result;
}


  /*

    Instantiation of templates (only these are generated)

  */

template std::vector<Vector<int> > standard_basis<int>(size_t n);

template class Vector<int>;           // the main instance used
template class Vector<signed char>;   // used inside root data
template class Vector<unsigned long>; // for |abelian::Homomorphism|
template class Matrix<int>;           // the main instance used
template class Matrix<unsigned long>; // for |abelian::Endomorphism|

template Matrix<int>::Matrix
  (std::vector<Vector<int> >::const_iterator,
   std::vector<Vector<int> >::const_iterator,
   size_t,
   tags::IteratorTag);

template Matrix<int>::Matrix
  (std::vector<Vector<int> >::iterator,
   std::vector<Vector<int> >::iterator,
   size_t,
   tags::IteratorTag);

template Matrix<int>::Matrix
  (Vector<int>*,
   Vector<int>*,
   size_t,
   tags::IteratorTag);



} // |namespace matrix|

} // |namespace atlas|
