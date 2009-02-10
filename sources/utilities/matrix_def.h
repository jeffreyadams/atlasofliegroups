/*!
\file
  This is matrix_def.h. This module contains some simple utilities for
  matrices. When we will need to do stuff for large matrices, we will
  need to look elsewhere, or in any case think a lot more.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <cassert>
#include "intutils.h"

/*****************************************************************************

  This module contains some simple utilities for matrices. When we will need
  to do stuff for large matrices, we will need to look elsewhere, or in any
  case think a lot more.

******************************************************************************/

namespace atlas {

namespace matrix {

namespace {

template<typename C> void blockReduce(Matrix<C>&, size_t, Matrix<C>&,
				      Matrix<C>&);
template<typename C> void blockShape(Matrix<C>&, size_t, Matrix<C>&,
				     Matrix<C>&);
template<typename C> void columnReduce(Matrix<C>&, size_t, size_t, Matrix<C>&);
template<typename C> bool hasBlockReduction(const Matrix<C>&, size_t);
template<typename C> bool hasReduction(const Matrix<C>&, size_t);
template<typename C> typename Matrix<C>::index_pair
  findBlockReduction(const Matrix<C>&, size_t);
template<typename C> typename Matrix<C>::index_pair
  findReduction(const Matrix<C>&, size_t);
template<typename C> void rowReduce(Matrix<C>&, size_t, size_t, Matrix<C>&);

} // namespace

} // namespace matrix

/*****************************************************************************

        Chapter II -- The Vector class

  This class template adds some arithmetic operations to |std::vector|

******************************************************************************/

namespace matrix {

//! \brief Increments by |v|.
template<typename C>
Vector<C>& Vector<C>::operator+= (const Vector<C>& v)
{
  for (size_t i=0; i<base::size(); ++i)
    (*this)[i] += v[i];
  return *this;
}

//! \brief Decrements by |v|.
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

  It is the callers responsibility to check divisibility; remainder is lost.
 */
template<typename C>
Vector<C>& Vector<C>::operator/= (C c)
{
  for (size_t i=0; i<base::size(); ++i)
    (*this)[i] /= c;
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

} // namespace matrix

/*****************************************************************************

        Chapter II -- The Matrix class

  We implement a matrix simply as a vector of elements, concatenating the
  rows.

******************************************************************************/

namespace matrix {

/******** constructors *******************************************************/

/*! \brief
  This constructor constructs a matrix from a bunch of vectors, columnwise.
  It is assumed that all elements of b have the same size.
*/
template<typename C>
Matrix<C>::Matrix(const std::vector<Vector<C> >& b)
{
  if (b.empty()) { // empty matrix
    d_rows = 0;
    d_columns = 0;
    return;
  }

  d_rows = b[0].size();
  d_columns = b.size();
  d_data.resize(d_rows*d_columns);

  for (size_t j = 0; j<d_columns; ++j) {
    for (size_t i = 0; i<d_rows; ++i)
      (*this)(i,j) = b[j][i];
  }
}


/*!
  This constructor constructs a square matrix, which is the matrix
  representing the operator defined by m in the canonical basis, in
  the basis b.
*/
template<typename C>
Matrix<C>::Matrix(const Matrix<C>& m, const std::vector<Vector<C> >& b)
  :d_data(m.d_data), d_rows(m.d_rows), d_columns(m.d_columns)
{
  Matrix<C> p(b);
  invConjugate(*this,p);
}


/*!
  Synopsis: constructs the matrix corresponding to the block [r_first,r_last[
  x [c_first,c_last[ of source.
*/
template<typename C>
Matrix<C>::Matrix(const Matrix<C>& source, size_t r_first, size_t c_first,
		  size_t r_last, size_t c_last)
  :d_data((r_last-r_first)*(c_last-c_first)),
   d_rows(r_last-r_first),
   d_columns(c_last-c_first)
{
  for (size_t j = 0; j<d_columns; ++j)
    for (size_t i = 0; i<d_rows; ++i)
      (*this)(i,j) = source(r_first+i,c_first+j);
}


/*!
  Here I is an iterator whose value type is Weight.

  This constructor constructs a square matrix, which is the matrix
  representing the operator defined by m in the canonical basis, in
  the basis supposed to be contained in [first,last[.
*/
template<typename C> template<typename I>
Matrix<C>::Matrix(const Matrix<C>& m, const I& first, const I& last)
  :d_data(m.d_data), d_rows(m.d_rows), d_columns(m.d_columns)
{
  Matrix<C> p(first,last,tags::IteratorTag());
  invConjugate(*this,p);
}


/*!
  Here I is an iterator type whose value type should be Vector<C>.
  The idea is that we read the columns of the matrix from the iterator.
  However, in order to be able to determine the allocation size,
  and since unfortunately we decided to read the matrix in rows, we
  need to catch the vectors first.

  Of course it is assumed that all the vectors in the range have the
  same size.
*/
template<typename C>
  template<typename I>
Matrix<C>::Matrix(const I& first, const I& last, tags::IteratorTag)
{
  std::vector<Vector<C> > b(first,last);

  if (b.empty()) {
    d_rows = 0;
    d_columns = 0;
    return;
  }

  d_rows = b[0].size();
  d_columns = b.size();
  d_data.resize(d_rows*d_columns);

  for (size_t i = 0; i<d_rows; ++i)
    for (size_t j = 0; j<d_columns; ++j)
      (*this)(i,j) = b[j][i];
}

/******** accessors **********************************************************/

template<typename C>
bool Matrix<C>::operator== (const Matrix<C>& m) const
{
  return d_rows==m.d_rows and d_columns==m.d_columns and d_data==m.d_data;
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

/*!
Applies the matrix to the vector w, and returns the result. It is assumed that
the size of w is the number of columns; v.size() is set to the number of rows.
It is safe to call this method with |v| and |w| the same vector.
*/
template<typename C>
void Matrix<C>::apply(Vector<C>& v, const Vector<C>& w) const
{
  v=apply(w);
}


/*!
  A pipe-version of apply. We assume that I is an InputIterator with
  value-type vector<C>, and O an OutputIterator with the same
  value-type. Then we apply our matrix to each vector in [first,last[
  and output it to out.
*/
template<typename C> template<typename I, typename O>
void Matrix<C>::apply(const I& first, const I& last, O out) const
{
  for (I i = first; i != last; ++i) {
    Vector<C> v(d_rows);
    apply(v,*i);
    *out++ = v;
  }
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

/*!
  Puts the j-th column of the matrix in v.
*/
template<typename C>
void Matrix<C>::set_column(Vector<C>& v, size_t j) const
{
  v.resize(d_rows);

  for (size_t i = 0; i<d_rows; ++i)
    v[i] = (*this)(i,j);
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
  Puts the i-th row of the matrix in v.
*/
template<typename C>
void Matrix<C>::set_row(Vector<C>& v, size_t i) const
{
  v.resize(d_columns);

  for (size_t j = 0; j<d_columns; ++j)
    v[j] = (*this)(i,j);
}

/******** manipulators *******************************************************/

/*!
  Incrementation by addition with m. It is assumed that m and *this
  have the same shape.
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
  have the same shape.
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
  \brief right multiplication by |m|.

  It is of course required that the number of rows of |m| is equal to the
  number of columns of |*this|.
*/
template<typename C>
Matrix<C>& Matrix<C>::operator*= (const Matrix<C>&  m)

{
  // I dont think the code below is faster than just |*this= *this * m|, MvL
  assert(d_columns==m.d_rows);

  C zero = 0; // conversion from int
  Matrix<C> prod(d_rows,m.d_columns,zero);

  for (size_t i = 0; i<d_rows; ++i)
    for (size_t j = 0; j<m.d_columns; ++j)
      for (size_t k = 0; k<d_columns; ++k)
      {
	C t_ik = (*this)(i,k);
	C m_kj = m(k,j);
	prod(i,j) += t_ik * m_kj;
      }

  swap(prod); // replace contents of |*this| by newly computed matrix

  return *this;
}


template<typename C>
Matrix<C>& Matrix<C>::operator/= (const C& c)
{
  if (c != 1)
    d_data /= c;
  return *this;
}


/*!
  Changes the sign of all the entries in column j.
*/
template<typename C>
void Matrix<C>::changeColumnSign(size_t j)
{
  for (size_t i = 0; i<d_rows; ++i) {
    (*this)(i,j) *= -1;
  }
}


/*!
  Changes the sign of all the entries in row i.
*/
template<typename C>
void Matrix<C>::changeRowSign(size_t i)
{
  for (size_t j = 0; j<d_columns; ++j) {
    (*this)(i,j) *= -1;
  }
}


/*!
  Carries out the column operation consisting of adding c times column j
  to column i.
*/
template<typename C>
void Matrix<C>::columnOperation(size_t i, size_t j, const C& c)
{
  for (size_t k = 0; k<columnSize(); ++k)
    (*this)(k,i) += c*(*this)(k,j);
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


/*!
  Copies the column c_s of the matrix m to the column c_d of the current
  matrix. It is assumed that the two columns have the same size.
*/
template<typename C>
void Matrix<C>::copyColumn(const Matrix<C>& m, size_t c_d, size_t c_s)
{
  for (size_t i = 0; i<d_rows; ++i)
    (*this)(i,c_d) = m(i,c_s);
}


/*!
  Copies the column r_s of the matrix m to the column r_d of the current
  matrix. It is assumed that the two columns have the same size.
*/
template<typename C>
void Matrix<C>::copyRow(const Matrix<C>& m, size_t r_d, size_t r_s)
{
  for (size_t j = 0; j<d_columns; ++j)
    (*this)(r_d,j) = m(r_s,j);
}


/*!
  Copies the column r_s of the matrix m to the column r_d of the current
  matrix. It is assumed that the two columns have the same size.
*/
template<typename C>
void Matrix<C>::eraseColumn(size_t j)
{
  iterator pos = begin() + j;
  --d_columns;

  for (size_t k = 0; k<d_rows; ++k, pos += d_columns)
    d_data.erase(pos);
}


/*!
  Erases row i from the matrix.
*/
template<typename C>
void Matrix<C>::eraseRow(size_t i)
{
  iterator first = begin() + i*d_columns;
  d_data.erase(first,first+d_columns);
  --d_rows;
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
  if (d_rows == 0) // do nothing to matrix, but set |d=1|
  { d=1; return; }

  C zero = 0; // conversion from int should work!
  Matrix<C> row(d_rows,d_rows,zero);

  // set res to identity

  for (size_t j = 0; j<d_rows; ++j)
    row(j,j) = 1;

  Matrix<C> col(row);

  // take *this to triangular form

  for (size_t r = 0; r<d_rows; ++r) {

    if (isZero(r,r)) { // null matrix left; matrix is not invertible
      d = 0;
      return;
    }

    index_pair k = absMinPos(r,r);

    if (k.first > r) {
      swapRows(r,k.first);
      row.swapRows(r,k.first);
    }

    if(k.second > r) {
      swapColumns(r,k.second);
      col.swapColumns(r,k.second);
    }

    if ((*this)(r,r)<0) {
      changeRowSign(r);
      row.changeRowSign(r);
    }

    // ensure m(0,0) divides row and column

    while (hasReduction(*this,r)) {
      k = findReduction(*this,r);
      if (k.first > r) { // row reduction
	rowReduce(*this,k.first,r,row);
      }
      else { // column reduction
	columnReduce(*this,k.second,r,col);
      }
    }

    // clean up row and column

    blockShape(*this,r,row,col);
    blockReduce(*this,r,row,col);

  } // next |r|

  // now multiply out diagonal

  for (size_t j = 1; j<d_rows; ++j) {
    (*this)(j,j) *= (*this)(j-1,j-1);
  }

  // and write result to |d| and to |*this|

  d = (*this)(d_rows-1,d_rows-1); // minimal denominator

  for (size_t j = 0; j<d_rows; ++j) {
    (*this)(j,j) = d/(*this)(j,j);
  }

  col *= *this;
  col *= row;
  swap(col);
}


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
  Synopsis: multiply the matrix by -1.
*/
template<typename C>
void Matrix<C>::negate()
{
  d_data.negate();
}


/*!
  Synopsis : permutes rows and columns of the matrix according to the
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
  Matrix<C> q(d_columns);

  for (size_t j = 0; j<d_columns; ++j)
    for (size_t i = 0; i<d_rows; ++i) {
      size_t a_i = a[i];
      size_t a_j = a[j];
      q(i,j) = (*this)(a_i,a_j);
    }

  swap(q);
}


/*!
  Synopsis: changes the size of the matrix to m x n.

  NOTE: this is a bad name, and a bit dangerous. One would expect the data
  of the matrix in the upper left corner to be preserved, but this is not
  done; in other words, the contents of the resized matrix are garbage.
*/
template<typename C>
void Matrix<C>::resize(size_t m, size_t n)
{
  d_data.resize(m*n);
  d_rows = m;
  d_columns = n;
}


/*!
  Synopsis: changes the size of the matrix to m x n, and resets _all_ elements
  to c.

  NOTE: this is a bad name, because it does not behave like resize for vectors.
*/
template<typename C>
void Matrix<C>::resize(size_t m, size_t n, const C& c)
{
  d_data.assign(m*n,c);
  d_rows = m;
  d_columns = n;
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
  Swaps m with the current matrix.
*/
template<typename C>
void Matrix<C>::swap(Matrix<C>& m)
{
  // swap data

  d_data.swap(m.d_data);

  // swap row and column numbers

  std::swap(d_rows,m.d_rows);
  std::swap(d_columns,m.d_columns);
}


/*!
  Interchanges columns i and j
*/
template<typename C>
void Matrix<C>::swapColumns(size_t i, size_t j)
{
  for (size_t k = 0; k<columnSize(); ++k) {
    C tmp = (*this)(k,i);
    (*this)(k,i) = (*this)(k,j);
    (*this)(k,j) = tmp;
  }
}


/*!
  Interchanges columns i and j
*/
template<typename C>
void Matrix<C>::swapRows(size_t i, size_t j)
{
  for (size_t k = 0; k<rowSize(); ++k) {
    C tmp = (*this)(i,k);
    (*this)(i,k) = (*this)(j,k);
    (*this)(j,k) = tmp;
  }
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

} // namespace matrix


/*****************************************************************************

        Chapter III --- Functions declared in matrix.h

******************************************************************************/

namespace matrix {

/*!
  Synopsis: writes in b the list of column vectors of m.
*/
template<typename C>
  void columnVectors(std::vector<Vector<C> >& b, const Matrix<C>& m)
{
  b.resize(m.numColumns());

  for (size_t j = 0; j<b.size(); ++j)
    m.set_column(b[j],j);
}

template<typename C> Matrix<C> row_matrix(Vector<C> v)
{
  Matrix<C> result(1,v.size());
  for (size_t j=0; j<v.size(); ++j)
    result(0,j)=v[j];

  return result;
}

template<typename C> Matrix<C> column_matrix(Vector<C> v)
{
  Matrix<C> result(v.size(),1);
  for (size_t i=0; i<v.size(); ++i)
    result(i,0)=v[i];

  return result;
}



/*!
  Conjugates m by p, i.e. transforms m into pmp^{-1}. It is assumed that
  p is invertible (over the quotient field of the coefficients), and that
  denominators cancel out.
*/
template<typename C>
Matrix<C>& conjugate(Matrix<C>& m, const Matrix<C>& p)
{
  C d;
  Matrix<C> tmp(p*m*p.inverse(d));
  tmp /= d;
  m.swap(tmp);

  return m;
}


/*!
  Synopsis: sets dest equal to the block [firstRow,lastRow[ x [firstColumn,
  lastColumn[ in source.
*/
template<typename C>
void extractBlock(Matrix<C>& dest, const Matrix<C>& source, size_t firstRow,
		  size_t lastRow, size_t firstColumn, size_t lastColumn)
{
  dest.resize(lastRow - firstRow, lastColumn - firstColumn);

  for (size_t j = 0; j<dest.numColumns(); ++j)
    for (size_t i = 0; i<dest.numRows(); ++i)
      dest(i,j) = source(firstRow+i,firstColumn+j);
}


/*!
  Synopsis: extracts the sub-matrix corresponding to the subsets r and c of the
  row and column indices respectively.
*/
template<typename C>
void extractMatrix(Matrix<C>& dest, const Matrix<C>& source,
		   const std::vector<size_t>& r, const std::vector<size_t>& c)
{
  dest.resize(r.size(),c.size());

  for (size_t i = 0; i<c.size(); ++i)
    for (size_t j = 0; j<c.size(); ++j)
      dest(i,j) = source(r[i],c[j]);
}


/*!
  Synopsis: puts in q the identity matrix of size n.
*/
template<typename C> void identityMatrix(Matrix<C>& m, size_t n)
{
  m.resize(n,n,C(0));

  for (size_t i = 0; i<n; ++i)
    m(i,i) = C(1);
}

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


/*!
  Conjugates m by p^{-1}, i.e. transforms m into p^{-1}mp. It is assumed that
  p is invertible (over the quotient field of the coefficients), and that
  denominators cancel out.
*/
template<typename C> Matrix<C>& invConjugate(Matrix<C>& m, const Matrix<C>& p)
{
  C d;
  Matrix<C> tmp(p.inverse(d)*m*p);
  tmp /= d;
  m.swap(tmp);

  return m;
}

} // namespace matrix

/*****************************************************************************

        Chapter IV --- Auxiliary functions

******************************************************************************/

namespace matrix {

namespace {

/*!
  Ensures that m(d,d) divides the block from (d+1,d+1) down.
*/
template<typename C>
void blockReduce(Matrix<C>& m, size_t d, Matrix<C>& r, Matrix<C>& c)
{
  if (m(d,d) == 1) // no reduction
    return;

  while(hasBlockReduction(m,d)) {
    typename Matrix<C>::index_pair k = findBlockReduction(m,d);
    C one = 1; // conversion from int
    m.rowOperation(d,k.first,one);
    r.rowOperation(d,k.first,one);
    while (hasReduction(m,d)) {
      k = findReduction(m,d);
      if (k.first > d) { // row reduction
	rowReduce(m,k.first,d,r);
      }
      else { // column reduction
	columnReduce(m,k.second,d,c);
      }
    }
    blockShape(m,d,r,c);
  }

  if (m(d,d) > 1) // divide
    for (size_t j = d+1; j<m.rowSize(); ++j)
      for (size_t i = d+1; i<m.columnSize(); ++i)
	m(i,j) /= m(d,d);
}


/*!
  Does the final reduction of m to block shape, recording row reductions in
  r and column reductions in c.
*/
template<typename C>
void blockShape(Matrix<C>& m, size_t d, Matrix<C>& r,Matrix<C>& c)
{
  C a = m(d,d);

  for (size_t j = d+1; j<m.rowSize(); ++j) {
    if (m(d,j) == 0)
      continue;
    C q = m(d,j)/a;
    q = -q;
    m.columnOperation(j,d,q);
    c.columnOperation(j,d,q);
  }

  for (size_t i = d+1; i<m.columnSize(); ++i) {
    if (m(i,d) == 0)
      continue;
    C q = m(i,d)/a;
    q = -q;
    m.rowOperation(i,d,q);
    r.rowOperation(i,d,q);
  }
}


/*!
  Does the column reduction for m at place j, and does the same operation
  on r. The reduction consists in subtracting from column j the multiple
  of column d which leaves at (d,j) the remainder of the Euclidian division
  of m(d,j) by m(d,d), and then swapping columns j and d.
*/
template<typename C>
void columnReduce(Matrix<C>& m, size_t j,  size_t d, Matrix<C>& c)
{
  C a = m(d,d);
  C q = intutils::divide(m(d,j),a);
  q = -q;
  m.columnOperation(j,d,q);
  c.columnOperation(j,d,q);
  m.swapColumns(d,j);
  c.swapColumns(d,j);
}


/*!
  Returns the reduction point of m. Assumes that hasBlockReduction(m,r) has
  returned true.
*/
template<typename C>
typename Matrix<C>::index_pair findBlockReduction(const Matrix<C>& m,
					     size_t r)
{
  C a = m(r,r);

  for (size_t j = r+1; j<m.rowSize(); ++j) {
    for (size_t i = r+1; i<m.columnSize(); ++i) {
      if (m(i,j)%a)
	return std::make_pair(i,j);
    }
  }

  // this should never be reached
  return std::make_pair(static_cast<size_t>(0ul),static_cast<size_t>(0ul));
}


/*!
  Returns the reduction point of m. Assumes that hasReduction(m,r) has
  returned true.
*/
template<typename C>
typename Matrix<C>::index_pair findReduction(const Matrix<C>& m,
					     size_t r)
{
  C a = m(r,r);

  for (size_t j = r+1; j<m.rowSize(); ++j) {
    if (m(r,j)%a)
      return std::make_pair(r,j);
  }

  for (size_t i = r+1; i<m.columnSize(); ++i) {
    if (m(i,r)%a)
      return std::make_pair(i,r);
  }

  // this should never be reached
  return std::make_pair(static_cast<size_t>(0ul),static_cast<size_t>(0ul));
}


/*!
  Tells if there is an element in the block under (r,r) which is not divisible
  bu m(r,r)
*/
template<typename C>
bool hasBlockReduction(const Matrix<C>& m, size_t r)
{
  C a = m(r,r);

  if (a == 1)
    return false;

  for (size_t j = r+1; j<m.rowSize(); ++j) {
    for (size_t i = r+1; i<m.columnSize(); ++i) {
      if (m(i,j)%a)
	return true;
    }
  }

  return false;
}


/*!
  Tells if there is an element in the first row or column not divisible by
  m(r,r).
*/
template<typename C>
bool hasReduction(const Matrix<C>& m, size_t r)
{
  C a = m(r,r);

  if (a == 1)
    return false;

  for (size_t j = r+1; j<m.rowSize(); ++j) {
    if (m(r,j)%a)
      return true;
  }

  for (size_t i = r+1; i<m.columnSize(); ++i) {
    if (m(i,r)%a)
      return true;
  }

  return false;
}

/*!
  Does the row reduction for m at place j, and does the same operation
  on r. The reduction consists in subtracting from row i the multiple
  of row d which leaves at (i,d) the remainder of the Euclidian division
  of m(d,j) by m(d,d), and then swapping rows i and d.
*/
template<typename C>
void rowReduce(Matrix<C>& m, size_t i, size_t d, Matrix<C>& r)
{
  C a = m(d,d);
  C q = intutils::divide(m(i,d),a);
  q = -q;
  m.rowOperation(i,d,q);
  r.rowOperation(i,d,q);
  m.swapRows(d,i);
  r.swapRows(d,i);
}

} // namespace

} // namespace matrix

} // namespace atlas
