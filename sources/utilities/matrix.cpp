/*
  This is matrix.cpp. This module contains some simple utilities for matrices.

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "matrix.h"

#include <cassert>
#include <limits>
#include <algorithm>
#include <stdexcept>

#include "matreduc.h"
#include "permutations.h"
#include "arithmetic.h"
#include "polynomials.h" // and most importantly polynomials_def.h


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

// Add |v| to |*this|
template<typename C>
Vector<C>& Vector<C>::operator+= (const Vector<C>& v)
{
  assert(base::size()==v.size());
  auto p = v.begin();
  for (auto it=base::begin(); it!=base::end(); ++it,++p)
    *it += *p;
  return *this;
}

// Subtracts |v| from |*this|
template<typename C>
Vector<C>& Vector<C>::operator-= (const Vector<C>& v)
{
  assert(base::size()==v.size());
  auto p = v.begin();
  for (auto it=base::begin(); it!=base::end(); ++it,++p)
    *it -= *p;
  return *this;
}

// Subtracts |*this| from |v|, and sets |*this| to the result
template<typename C>
Vector<C>& Vector<C>::negate_add (const Vector<C>& v)
{
  assert(base::size()==v.size());
  auto p = v.begin();
  for (auto it=base::begin(); it!=base::end(); ++it,++p)
    *it = *p - *it;
  return *this;
}

template<typename C>
template<typename I>
Vector<C>& Vector<C>::add (I b, C c)
{
  for (auto it=base::begin(); it!=base::end(); ++it,++b)
    *it += *b * c;
  return *this;
}

//! \brief Scalar multiplies by |c|
template<typename C>
Vector<C>& Vector<C>::operator*= (C c)
{
  for (auto it=base::begin(); it!=base::end(); ++it)
    *it *= c;
  return *this;
}

/*
  Scalar-divide by |c|
  All entries must allow exact division, if not a std::runtime_error is thrown
*/
template<typename C>
Vector<C>& operator/= (Vector<C>& v,C c)
{
  if (c==C(0))
    throw std::runtime_error("Vector division by 0");
  for (auto it=v.begin(); it!=v.end(); ++it)
    if (*it%c==C(0))
      *it/=c;
    else throw std::runtime_error("Inexact vector integer division");
  return v;
}

template<typename C>
Vector<C>& divide (Vector<C>& v,C c)
{
  if (c==C(0))
    throw std::runtime_error("Vector division by 0");
  if (c>C(0))
    for (auto it=v.begin(); it!=v.end(); ++it)
      *it = arithmetic::divide(*it,c);
  else
  { c = -c; // we must ensure |c>0| for |arithmetic::divide|
    for (auto it=v.begin(); it!=v.end(); ++it)
      *it = -arithmetic::divide(*it,c);
  }
  return v;
}

template<typename C>
Vector<C>& operator%= (Vector<C>& v,C c)
{
  if (c==C(0))
    throw std::runtime_error("Vector taken modulo 0");
  c=std::abs(c); // we must ensure |c>0| for |arithmetic::remainder|
  for (auto it=v.begin(); it!=v.end(); ++it)
    *it = arithmetic::remainder(*it,c);
  return v;
}

template<typename C>
Vector<C>& Vector<C>::negate ()
{
  for (auto it=base::begin(); it!=base::end(); ++it)
    *it = -*it;
  return *this;
}

template<typename C>
template<typename C1>
  C1 Vector<C>::dot (const Vector<C1>& v) const
{
  assert(base::size()==v.size());
  C1 result= C1(0);
  for (size_t i=0; i<base::size(); ++i)
    result += (*this)[i] * v[i];

  return result;
}

template<typename C>
  bool Vector<C>::isZero() const
{
  for (auto it=base::begin(); it!=base::end(); ++it)
    if (*it!=C(0))
      return false;
  return true;
}


template<typename C>
template<typename C1>
Vector<C1> Vector<C>::scaled (C1 c) const
{
  Vector<C1> result(base::size());
  auto p=result.begin();
  for (auto it=base::begin(); it!=base::end(); ++p,++it)
    *p = *it * c;
  return result;
}

template<typename C>
  Matrix<C> Vector<C>::row_matrix() const
{
  Matrix<C> result(1,base::size());
  auto p=base::begin();
  for (size_t j=0; j<base::size(); ++j,++p)
    result(0,j) = *p;

  return result;
}

template<typename C>
  Matrix<C> Vector<C>::column_matrix() const
{
  Matrix<C> result(base::size(),1);
  auto p=base::begin();
  for (size_t i=0; i<base::size(); ++i,++p)
    result(i,0) = *p;

  return result;
}


} // |namespace matrix|

/*****************************************************************************

        Chapter II -- The Matrix class

  We implement a matrix as a vector of elements, concatenating the rows.

******************************************************************************/

namespace matrix {

/******** constructors *******************************************************/

//! \brief: Construct the identity matrix of size n.
template<typename C> Matrix<C>::Matrix(size_t n) : base(n,n,C(0))
{
  for (size_t i=0; i<n; ++i)
    base::operator()(i,i) = C(1);
}

/*! \brief
  This constructor constructs a matrix from a list of vectors, column-wise.
  It is assumed that all elements of |b| (possibly none) have size |n_rows|.
*/
template<typename C>
Matrix_base<C>::Matrix_base(const std::vector<Vector<C> >& b, size_t n_rows)
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
  template<typename I> Matrix_base<C>::Matrix_base
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
void Matrix_base<C>::swap(Matrix_base<C>& m)
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
void Matrix_base<C>::get_row(Vector<C>& v, size_t i) const
{
  assert(i<d_rows);
  v.resize(d_columns);
  std::copy(at(i,0),at(i+1,0),&v[0] );
}

/*
  Puts the j-th column of the matrix in v.
*/
template<typename C>
void Matrix_base<C>::get_column(Vector<C>& v, size_t j) const
{
  assert(j<d_columns);
  v.resize(d_rows);

  for (size_t i = 0; i<d_rows; ++i)
    v[i] = (*this)(i,j);
}

/*!
  Returns the list of row vectors of m.
*/
template<typename C>
  std::vector<Vector<C> > Matrix_base<C>::rows() const
{
  std::vector<Vector<C> > result(d_rows);
  for (size_t i=0; i<d_rows; ++i)
    get_row(result[i],i);
  return result;
}

/*!
  Returns the list of column vectors of m.
*/
template<typename C>
  std::vector<Vector<C> > Matrix_base<C>::columns() const
{
  std::vector<Vector<C> > result(d_columns);

  for (size_t j = 0; j<d_columns; ++j)
    get_column(result[j],j);

  return result;
}

template<typename C>
bool Matrix_base<C>::operator== (const Matrix_base<C>& m) const
{
  return d_rows==m.d_rows and d_columns==m.d_columns and d_data==m.d_data;
}

template<typename C>
bool Matrix_base<C>::is_zero () const
{
  const auto end=d_data.end();
  for (auto it=d_data.begin(); it!=end; ++it)
    if (*it!=C(0))
      return false;
  return true;
}

/*
Apply the matrix to the vector |w|, and returns the result. It is assumed that
the size of |w| is the number of columns; result size is the number of rows.
*/
template<typename C>
template<typename C1>
Vector<C1> Matrix<C>::operator*(const Vector<C1>& w) const
{
  assert(base::numColumns()==w.size());
  Vector<C1> result(base::numRows());

  for (size_t i=0; i<base::numRows(); ++i)
  {
    C1 c(0);
    for (size_t j=0; j<base::numColumns(); ++j)
      c += (*this)(i,j) * w[j];
    result[i] = c;
  }

  return result;
}

/*
  Multiplies the matrix to right to the row-vector w, and returns the result.
  It is assumed that the size of w is the number of rows; result size is the
  number of columns. This is the proper sense of application for dual space.
*/
template<typename C> template<typename C1>
Vector<C1> Matrix<C>::right_prod(const Vector<C1>& w) const
{
  assert(base::numRows()==w.size());
  Vector<C1> result(base::numColumns());

  for (size_t j=0; j<base::numColumns(); ++j)
  {
    C1 c(0);
    for (size_t i=0; i<base::numRows(); ++i)
      c += w[i] * (*this)(i,j);
    result[j] = c;
  }

  return result;
}

template<typename C>
Matrix<C> Matrix<C>::operator* (const Matrix<C>&  m) const
{
  assert(base::numColumns()==m.numRows());
  Matrix<C> result(base::numRows(),m.numColumns());

  for (size_t i=0; i<base::numRows(); ++i)
    for (size_t k=0; k<m.base::numColumns(); ++k)
    {
      C c(0);
      for (size_t j=0; j<base::numColumns(); ++j)
	c += (*this)(i,j) * m(j,k);

      result(i,k)=c;
    }

  return result;
}

/*!
  This function constructs expresses our square matrix on the basis b.
*/
template<typename C>
  PID_Matrix<C> PID_Matrix<C>::on_basis(const std::vector<Vector<C> >& b) const
{
  assert (base::numRows()==base::numColumns());
  assert (b.size()==base::numRows());

  PID_Matrix<C> p(b,b.size()); // square matrix
  C d;
  PID_Matrix<C> result(p.inverse(d)* *this *p);
  return result /= d;
}

// manipulators


//! Puts v in the i-th row of the matrix
template<typename C>
  void Matrix_base<C>::set_row(size_t i, const Vector<C>& v)
{
  assert(v.size()==d_columns);
  std::copy(&v[0],&v[d_columns],at(i,0));
}

//! Puts v in the j-th column of the matrix
template<typename C>
  void Matrix_base<C>::set_column(size_t j, const Vector<C>& v)
{
  assert(v.size()==d_rows);

  for (size_t i=0; i<d_rows; ++i)
    (*this)(i,j)=v[i];
}

// Add |v| as new row to the matrix. Entries being stored by row, this is easy
template<typename C>
  void Matrix_base<C>::add_row(const Vector<C>& v)
{
  assert(v.size()==d_columns);

  ++d_rows;
  d_data.resize(d_rows*d_columns,C(0));

  std::copy(v.begin(),v.end(),at(d_rows-1,0));
}

// Add |v| as new column to the matrix. This is harder than adding a row
template<typename C>
  void Matrix_base<C>::add_column(const Vector<C>& v)
{
  assert(v.size()==d_rows);

  size_t old_col=d_columns;
  ++d_columns;
  d_data.resize(d_rows*d_columns,C(0));

  typename Vector<C>::iterator dst=d_data.end(),src=dst-d_rows;
  for (size_t i=d_rows; i-->0;)
  {
    *--dst = v[i];
    for (size_t j=old_col; j-->0;)
      *--dst = *--src;
  }
}

// fill |dst| with submatrix starting at offsets |k,l|
template<typename C>
  void Matrix_base<C>::get_block
       (Matrix_base<C>& dst, unsigned k, unsigned l) const
{
  const size_t k_end=k+dst.numRows(), l_end=l+dst.numColumns();
  assert(k_end<=numRows() and l_end<=numColumns());
  auto p=dst.at(0,0);
  for (size_t i=k; i<k_end; ++i)
    p=std::copy(this->at(i,l),this->at(i,l_end),p);
}

// fill submatrix starting at offsets |k,l| from |src|
template<typename C>
  void Matrix_base<C>::set_block
       (unsigned k, unsigned l, const Matrix_base<C>& src)
{
  const size_t r=src.numRows();
  assert(k+r<=numRows() and l+src.numColumns()<=numColumns());
  for (size_t i=0; i<r; ++i)
    std::copy(src.at(i,0),src.at(i+1,0),this->at(k+i,l));
}


// Add matrix |m| to |*this|; the matrices must have the same dimensions
template<typename C>
Matrix<C>& Matrix<C>::operator+= (const Matrix<C>& m)
{
  assert(base::numRows()==m.numRows());
  assert(base::numColumns()==m.numColumns());
  base::d_data += m.d_data;
  return *this;
}

// Subtract matrix |m| to |*this|; the matrices must have the same dimensions
template<typename C>
Matrix<C>& Matrix<C>::operator-= (const Matrix<C>&  m)
{
  assert(base::numRows()==m.numRows());
  assert(base::numColumns()==m.numColumns());
  base::d_data -= m.d_data;
  return *this;
}


template<typename C> Matrix<C> Matrix<C>::transposed() const
#if __GNUC__ > 4 || \
  __GNUC__ == 4 && (__GNUC_MINOR__ > 8 || \
		    __GNUC_MINOR__ == 8 && __GNUC_PATCHLEVEL__ >= 1)
// that is, if compiler version is sufficiently new
  & // though it makes no difference, must explicitly write lvalue ref-qualifier
#endif
{
  Matrix<C> result(base::numColumns(),base::numRows());
  auto p = base::d_data.begin();
  for (size_t i=0; i<base::numRows(); ++i)
    for (size_t j=0; j<base::numColumns(); ++j)
      result(j,i) = *p++; // fill by column, read by row
  return result;
}

/*!
  Transposes the matrix. We allow ourselves a temporary vector if the matrix
  is not square; it could likely be done in-place, but the permutation is
  rather complicated!
*/
template<typename C> Matrix<C>& Matrix<C>::transpose()
{
  if (base::numRows() == base::numColumns()) // matrix is square
    for (size_t j=0; j<base::numColumns(); ++j)
      for (size_t i = j+1; i<base::numRows(); ++i)
	std::swap((*this)(i,j),(*this)(j,i));
  else // now matrix is not square; create a transposed copy
    *this = transposed(); // move-assign from transposed copy
  return *this;
}


/*!
  Here we invert the matrix without catching the denominator. The intent
  is that it should be used for invertible matrices only.
*/
template<typename C> PID_Matrix<C>& PID_Matrix<C>::invert()
{
  C d; invert(d);
  assert(d==C(1));
  return *this;
}

/*!
  This function inverts the matrix M. It is assumed that the coefficents
  are in an integral domain. At the conclusion of the process, d will contain
  the denominator for the inverted matrix (so that the true result is M/d).

  If the matrix is not invertible, the resulting matrix is undefined, and
  d is set to 0.

  We have chosen here to avoid divisions as much as possible. So in fact we
  take the matrix to diagonal form through elementary operations; from there
  we deduce the smallest possible denominator for the inverse, and the inverse
  matrix. We have to keep track of both row and column operations.

*/
template<typename C>
PID_Matrix<C>& PID_Matrix<C>::invert(C& d)
{
  assert(base::numRows()==base::numColumns());
  size_t n=base::numRows();
  if (n==0) // do nothing to matrix, but set |d=1|
  { d=C(1); return *this; }

  PID_Matrix<C> row,col;    // for recording column operations
  std::vector<C> diagonal = matreduc::diagonalise(*this,row,col);

  if (diagonal.size()<n) // insufficient rank for inversion
  { d=C(0); return *this; } // record zero determinant, leave |*this| unchanged

  d=std::abs(diagonal[0]); // guaranteed to exist if we get here
  for (size_t i=1; i<n; ++i)
    d=arithmetic::lcm(d,diagonal[i]); // other diagonal entries are positive

  // finally for |D| "inverse" diagonal matrix w.r.t. |d|, compute |col*D*row|
  for (size_t j=0; j<n; ++j)
  {
    C f = d/diagonal[j];
    for (size_t i=0; i<n; ++i)
      (*this)(i,j)=f*col(i,j);
  }

  return *this *= row;
}

/*!
  Tells if all coefficients of the matrix are divisible by c.
*/
template<typename C>
bool PID_Matrix<C>::divisible(C c) const
{
  for (size_t j=0; j<base::d_data.size(); ++j)
    if (base::d_data[j]%c!=0)
      return false;

  return true;
}


/*!
  Synopsis: constructs the matrix corresponding to the block [i0,i1[ x [j0,j1[
  of source. This implementation uses that storage is by rows.
*/
template<typename C>
  PID_Matrix<C>
    PID_Matrix<C>::block(size_t i0, size_t j0, size_t i1, size_t j1) const
{
  assert(i0<=i1 and i1<=base::numRows());
  assert(j0<=j1 and j1<=base::numColumns());

  PID_Matrix<C> result(i1-i0,j1-j0);
  C* p = result.at(0,0); // writing pointer
  for (size_t i=i0; i<i1; ++i)
    p = std::copy(this->at(i,j0),this->at(i,j1),p); // copy row, advance
  return result;
}



// The row operation consisting of adding |c| times row |k| to row |i|.
template<typename C>
void Matrix<C>::rowOperation(size_t i, size_t k, const C& c)
{
  assert(i<base::numRows() and k<base::numRows());
  if (c!=C(0))
    for (size_t j=0; j<base::numColumns(); ++j)
      (*this)(i,j) += c*(*this)(k,j);
}


// The column operation consisting of adding |c| times column |k| to column |j|.
template<typename C>
void Matrix<C>::columnOperation(size_t j, size_t k, const C& c)
{
  assert(j<base::numColumns() and k<base::numColumns());
  if (c!=C(0))
    for (size_t i=0; i<base::numRows(); ++i)
      (*this)(i,j) += c*(*this)(i,k);
}



/*!
  Changes the sign of all the entries in row i.
*/
template<typename C>
void Matrix<C>::rowMultiply(size_t i, C f)
{
  assert(i<base::numRows());
  if (f!=C(1))
    for (size_t j=0; j<base::numColumns(); ++j)
      (*this)(i,j) *= f;
}

/*!
  Changes the sign of all the entries in column j.
*/
template<typename C>
void Matrix<C>::columnMultiply(size_t j, C f)
{
  assert(j<base::numColumns());
  if (f!=C(1))
    for (size_t i=0; i<base::numRows(); ++i)
      (*this)(i,j) *= f;
}


//! Interchanges rows
template<typename C>
void Matrix<C>::swapRows(size_t i0, size_t i1)
{
  assert(i0<base::numRows() and i1<base::numRows());
  for (size_t k=0; k<base::numColumns(); ++k)
    std::swap((*this)(i0,k),(*this)(i1,k));
}

//! Interchange columns
template<typename C>
void Matrix<C>::swapColumns(size_t j0, size_t j1)
{
  assert(j0<base::numColumns() and j1<base::numColumns());
  for (size_t k=0; k<base::numRows(); ++k)
    std::swap((*this)(k,j0),(*this)(k,j1));
}

/*!
  Erases row i from the matrix.
*/
template<typename C>
void Matrix_base<C>::eraseRow(size_t i)
{
  assert(i<d_rows);
  typename Vector<C>::iterator first = d_data.begin() + i*d_columns;
  d_data.erase(first,first+d_columns);
  --d_rows;
}

/*!
  Removes column |j| altogether, shifting remaining entries to proper place
*/
template<typename C>
void Matrix_base<C>::eraseColumn(size_t j)
{
  assert(j<d_columns);
  typename Vector<C>::iterator
    pos = d_data.begin() + j; // position of entry |M(0,j)|
  --d_columns; // already adjust number of columns

  // kill individual entries from left to right, with shifts of new |d_colums|
  // (this is easily seen to be correct but not particularly efficient)
  for (size_t k=0; k<d_rows; ++k, pos += d_columns)
    d_data.erase(pos);
}



} // |namespace matrix|


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

  for (size_t i=0; i<r; ++i)
    result[i][i] = C(1);

  return result;
}

template<typename C>
  PID_Matrix<C>& operator+= (PID_Matrix<C>& A, C c) // |A=A+c|, avoiding copy
{
  unsigned int i=std::min(A.numRows(),A.numColumns());
  while (i-->0)
    A(i,i) += c;
  return A;
}


  /*

    Instantiation of templates (only these are generated)

  */

 // type abreviations used in these instantiations
typedef arithmetic::Numer_t Num;
typedef polynomials::Polynomial<int> Pol;


template std::vector<Vector<int> > standard_basis<int>(size_t n);
template PID_Matrix<int>& operator+=(PID_Matrix<int>&,int);

template class Vector<int>;           // the main instance used
template class Vector<signed char>;   // used inside root data
template class Vector<unsigned long>; // for |abelian::Homomorphism|
template class Vector<Num>;           // numerators of rational vectors
template class Matrix_base<int>;
template class Matrix<int>;           // the main instance used
template class Matrix_base<unsigned long>; // for |abelian::Endomorphism|
template class Matrix<arithmetic::Split_integer>; // KL matrices eval'd at |s|
template class PID_Matrix<int>;

// template member instances
template Vector<int>& Vector<int>::add(Vector<int>::const_iterator b,int c);
template Vector<Num>& Vector<Num>::add(Vector<int>::const_iterator b,Num c);
template Vector<Num>& Vector<Num>::add(Vector<Num>::const_iterator b,Num c);
template int Vector<int>::dot(Vector<int> const&) const;
template Num Vector<int>::dot(Vector<Num> const&) const;
template int Vector<Num>::dot(Vector<int> const&) const;
template signed char
  Vector<signed char>::dot(const Vector<signed char>&) const;

template Vector<int>& operator/=(Vector<int>&,int);
template Vector<int>& divide (Vector<int>&,int);
template Vector<int>& operator%=(Vector<int>&,int);
template Vector<Num>& operator/=(Vector<Num>&,Num);

template Vector<int> Matrix<int>::operator*(Vector<int> const&) const;
template Vector<Num> Matrix<int>::operator*(Vector<Num> const&) const;
template Vector<int> Matrix<int>::right_prod(const Vector<int>&) const;
template Vector<Num> Matrix<int>::right_prod(const Vector<Num>&) const;

template Matrix_base<int>::Matrix_base
  (std::vector<Vector<int> >::const_iterator,
   std::vector<Vector<int> >::const_iterator,
   size_t,
   tags::IteratorTag);

template Matrix_base<int>::Matrix_base
  (std::vector<Vector<int> >::iterator,
   std::vector<Vector<int> >::iterator,
   size_t,
   tags::IteratorTag);

template Matrix_base<int>::Matrix_base
  (Vector<int>*,
   Vector<int>*,
   size_t,
   tags::IteratorTag);


template class Vector<Pol>;
template class Matrix_base<Pol>;
template class Matrix<Pol>;

template Vector<Pol> Matrix<Pol>::operator*(Vector<Pol> const&) const;

} // |namespace matrix|

} // |namespace atlas|
