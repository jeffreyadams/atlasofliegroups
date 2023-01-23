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

#include "sl_list.h"
#include "matreduc.h"
#include "permutations.h"
#include "arithmetic.h"
#include "bigint.h"
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

// Subtract |v| from |*this|
template<typename C>
Vector<C>& Vector<C>::operator-= (const Vector<C>& v)
{
  assert(base::size()==v.size());
  auto p = v.begin();
  for (auto it=base::begin(); it!=base::end(); ++it,++p)
    *it -= *p;
  return *this;
}

// Subtract |*this| from |v|, and sets |*this| to the result
template<typename C>
Vector<C>& Vector<C>::negate_add (const Vector<C>& v)
{
  assert(base::size()==v.size());
  auto p = v.begin();
  for (auto it=base::begin(); it!=base::end(); ++it,++p)
    *it = *p - *it;
  return *this;
}

// Scalar multiply our vector by |c|
template<typename C>
Vector<C>& Vector<C>::operator*= (C c)
{
  for (auto it=base::begin(); it!=base::end(); ++it)
    *it *= c;
  return *this;
}

/*
  Scalar-divide our vector by |c|
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
  bool Vector<C>::isZero() const
{
  for (auto it=base::begin(); it!=base::end(); ++it)
    if (*it!=C(0))
      return false;
  return true;
}

template<typename C>
  Matrix<C> Vector<C>::row_matrix() const
{
  Matrix<C> result(1,base::size());
  auto p=base::begin();
  for (index_t j=0; j<base::size(); ++j,++p)
    result(0,j) = *p;

  return result;
}

template<typename C>
  Matrix<C> Vector<C>::column_matrix() const
{
  Matrix<C> result(base::size(),1);
  auto p=base::begin();
  for (index_t i=0; i<base::size(); ++i,++p)
    result(i,0) = *p;

  return result;
}


/*****************************************************************************

        Chapter II -- The Matrix class

  We implement a matrix as a vector of elements, concatenating the rows.

******************************************************************************/

/******** constructors *******************************************************/

// Construct the identity matrix of size n.
template<typename C> Matrix<C>::Matrix(index_t n) : base(n,n,C(0))
{
  for (index_t i=0; i<n; ++i)
    base::operator()(i,i) = C(1);
}

/*
  This constructor constructs a matrix from a list of vectors, column-wise.
  It is assumed that all elements of |b| (possibly none) have size |n_rows|.
*/
template<typename C>
  Matrix_base<C>::Matrix_base
    (const std::vector<Vector<C> >& b, index_t n_rows)
  : d_rows(n_rows), d_columns(b.size())
  , d_data(d_rows*b.size()) // computed as |std::size_t| value
{
  for (index_t j = 0; j<d_columns; ++j)
    set_column(j,b[j]);
}

template<typename C>
void Matrix_base<C>::swap(Matrix_base<C>& m)
{
  std::swap(d_rows,m.d_rows);
  std::swap(d_columns,m.d_columns);
  d_data.swap(m.d_data);
}

/******** accessors **********************************************************/

// Assign row |i| of the matrix to |v|.
template<typename C>
void Matrix_base<C>::get_row(Vector<C>& v, index_t i) const
{
  assert(i<d_rows);
  v.resize(d_columns);
  std::copy(at(i,0),at(i+1,0),&v[0] );
}

// Assign column |j| of the matrix to |v|.
template<typename C>
void Matrix_base<C>::get_column(Vector<C>& v, index_t j) const
{
  assert(j<d_columns);
  v.resize(d_rows);

  for (index_t i = 0; i<d_rows; ++i)
    v[i] = (*this)(i,j);
}

// Return the list of row vectors of m.
template<typename C>
  std::vector<Vector<C> > Matrix_base<C>::rows() const
{
  std::vector<Vector<C> > result(d_rows);
  for (index_t i=0; i<d_rows; ++i)
    get_row(result[i],i);
  return result;
}

// Return the list of column vectors of m.
template<typename C>
  std::vector<Vector<C> > Matrix_base<C>::columns() const
{
  std::vector<Vector<C> > result(d_columns);

  for (index_t j = 0; j<d_columns; ++j)
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

template<typename C>
Matrix<C> Matrix<C>::operator* (const Matrix<C>&  m) const
{
  assert(base::n_columns()==m.n_rows());
  Matrix<C> result(base::n_rows(),m.n_columns());

  for (index_t i=0; i<base::n_rows(); ++i)
    for (index_t k=0; k<m.base::n_columns(); ++k)
    {
      C c(0);
      for (index_t j=0; j<base::n_columns(); ++j)
	c += (*this)(i,j) * m(j,k);

      result(i,k)=c;
    }

  return result;
}

// Express our square matrix on the basis |b| rather than the standard basis
template<typename C>
  PID_Matrix<C> PID_Matrix<C>::on_basis(const PID_Matrix<C>& basis) const
{ // both |*this| and |basis| are square, and of same size; no need to |assert|
  arithmetic::big_int d;
  PID_Matrix<C> result(basis.inverse(d)* *this *basis);
  return result /= d.convert<C>();
}

// manipulators


// Puts |v| in the |i|-th row of the matrix
template<typename C>
  void Matrix_base<C>::set_row(index_t i, const Vector<C>& v)
{
  assert(v.size()==d_columns);
  std::copy(&v[0],&v[d_columns],at(i,0));
}

// Put |v| in the |j|-th column of the matrix
template<typename C>
  void Matrix_base<C>::set_column(index_t j, const Vector<C>& v)
{
  assert(v.size()==d_rows);

  for (index_t i=0; i<d_rows; ++i)
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

  index_t old_col=d_columns;
  ++d_columns;
  d_data.resize(d_rows*d_columns,C(0));

  typename Vector<C>::iterator dst=d_data.end(),src=dst-d_rows;
  for (index_t i=d_rows; i-->0;)
  {
    *--dst = v[i];
    for (index_t j=old_col; j-->0;)
      *--dst = *--src;
  }
}

// fill |dst| with submatrix starting at offsets |k,l|
template<typename C>
  void Matrix_base<C>::get_block
       (Matrix_base<C>& dst, index_t k, index_t l) const
{
  const index_t k_end=k+dst.n_rows(), l_end=l+dst.n_columns();
  assert(k_end<=n_rows() and l_end<=n_columns());
  auto p=dst.at(0,0);
  for (index_t i=k; i<k_end; ++i)
    p=std::copy(this->at(i,l),this->at(i,l_end),p);
}

// fill submatrix starting at offsets |k,l| from |src|
template<typename C>
  void Matrix_base<C>::set_block
       (index_t k, index_t l, const Matrix_base<C>& src)
{
  const index_t r=src.n_rows();
  assert(k+r<=n_rows() and l+src.n_columns()<=n_columns());
  for (index_t i=0; i<r; ++i)
    std::copy(src.at(i,0),src.at(i+1,0),this->at(k+i,l));
}


// Add matrix |m| to |*this|; the matrices must have the same dimensions
template<typename C>
Matrix<C>& Matrix<C>::operator+= (const Matrix<C>& m)
{
  assert(base::n_rows()==m.n_rows());
  assert(base::n_columns()==m.n_columns());
  base::d_data += m.d_data;
  return *this;
}

// Subtract matrix |m| to |*this|; the matrices must have the same dimensions
template<typename C>
Matrix<C>& Matrix<C>::operator-= (const Matrix<C>&  m)
{
  assert(base::n_rows()==m.n_rows());
  assert(base::n_columns()==m.n_columns());
  base::d_data -= m.d_data;
  return *this;
}


template<typename C> Matrix<C> Matrix<C>::transposed() const &
{
  Matrix<C> result(base::n_columns(),base::n_rows());
  auto p = base::d_data.begin();
  for (index_t i=0; i<base::n_rows(); ++i)
    for (index_t j=0; j<base::n_columns(); ++j)
      result(j,i) = *p++; // fill by column, read by row
  return result;
}

/*
  Transpose the matrix. We allow ourselves a temporary vector if the matrix
  is not square; it could likely be done in-place, but the permutation is
  rather complicated!
*/
template<typename C> Matrix<C>& Matrix<C>::transpose()
{
  if (base::n_rows() == base::n_columns()) // matrix is square
    for (index_t j=0; j<base::n_columns(); ++j)
      for (index_t i = j+1; i<base::n_rows(); ++i)
	std::swap((*this)(i,j),(*this)(j,i));
  else // now matrix is not square; create a transposed copy
    *this = transposed(); // move-assign from transposed copy
  return *this;
}



template<typename C>
PID_Matrix<C> inverse (PID_Matrix<C> A, arithmetic::big_int& d)
{
  assert(A.n_rows()==A.n_columns());
  index_t n=A.n_rows();
  PID_Matrix<C> result(n,n);
  if (n==0) // do nothing to matrix, but set |d=1|
  { d=arithmetic::big_int(1); return result; }

  bool flip=false;
  PID_Matrix<C> col; // for recording column operations
  auto pivots=matreduc::column_echelon(A,col,flip);
  if (not pivots.full())
  { d=arithmetic::big_int(0); // record zero determinant,
    return result; // return useless value
  }

  std::vector<arithmetic::big_int> denoms(n); // next loop initialises entries
  auto basis = matrix::standard_basis<C>(n);
  for (index_t j=0; j<n; ++j) // from easy to harder, but order isn't vital
    basis[j] = matreduc::echelon_solve(A,pivots,std::move(basis[j]),denoms[j]);

  d=denoms[0];
  for (index_t j=1; j<n; ++j)
    d=arithmetic::lcm(d,denoms[j]);

  // finally for |D| "inverse" diagonal matrix w.r.t. |d|, compute |col*D*row|
  for (index_t j=0; j<n; ++j)
    result.set_column(j,col*basis[j]*(d/denoms[j]).int_val());

  return result;
}

// Whether all coefficients of the matrix are divisible by c.
template<typename C>
bool PID_Matrix<C>::divisible(C c) const
{
  for (index_t j=0; j<base::d_data.size(); ++j)
    if (base::d_data[j]%c!=0)
      return false;

  return true;
}


// Extract the block at indices [i0,i1[ x [j0,j1[ as a |PID_Matrix|
template<typename C>
  PID_Matrix<C>
    PID_Matrix<C>::block(index_t i0, index_t j0, index_t i1, index_t j1) const
{
  assert(i0<=i1 and i1<=base::n_rows());
  assert(j0<=j1 and j1<=base::n_columns());

// The implementation uses that storage is by rows.
  PID_Matrix<C> result(i1-i0,j1-j0);
  C* p = result.at(0,0); // writing pointer
  for (index_t i=i0; i<i1; ++i)
    p = std::copy(this->at(i,j0),this->at(i,j1),p); // copy row, advance
  return result;
}

// Extract the block at indices [j0,j1[ x [i0,i1[, transposed, as a |Matrix|
template<typename C>
  PID_Matrix<C>
    PID_Matrix<C>::transposed_block
     (index_t i0, index_t j0, index_t i1, index_t j1) const
{
  assert(i0<=i1 and i1<=base::n_columns());
  assert(j0<=j1 and j1<=base::n_rows());

  PID_Matrix<C> result(i1-i0,j1-j0);
  for (index_t i=0; i<result.n_rows(); ++i)
    for (index_t j=0; j<result.n_columns(); ++j)
      result(i,j)=(*this)(j0+j,i0+i);
  return result;
}

// The row operation consisting of adding |c| times row |k| to row |i|.
template<typename C>
void Matrix<C>::rowOperation(index_t i, index_t k, const C& c)
{
  assert(i<base::n_rows() and k<base::n_rows());
  if (c!=C(0))
    for (index_t j=0; j<base::n_columns(); ++j)
      (*this)(i,j) += c*(*this)(k,j);
}


// The column operation consisting of adding |c| times column |k| to column |j|.
template<typename C>
void Matrix<C>::columnOperation(index_t j, index_t k, const C& c)
{
  assert(j<base::n_columns() and k<base::n_columns());
  if (c!=C(0))
    for (index_t i=0; i<base::n_rows(); ++i)
      (*this)(i,j) += c*(*this)(i,k);
}



template<typename C>
void Matrix<C>::rowMultiply(index_t i, C f)
{
  assert(i<base::n_rows());
  if (f!=C(1))
    for (index_t j=0; j<base::n_columns(); ++j)
      (*this)(i,j) *= f;
}

template<typename C>
void Matrix<C>::columnMultiply(index_t j, C f)
{
  assert(j<base::n_columns());
  if (f!=C(1))
    for (index_t i=0; i<base::n_rows(); ++i)
      (*this)(i,j) *= f;
}


template<typename C>
void Matrix<C>::swapRows(index_t i0, index_t i1)
{
  assert(i0<base::n_rows() and i1<base::n_rows());
  for (index_t k=0; k<base::n_columns(); ++k)
    std::swap((*this)(i0,k),(*this)(i1,k));
}

template<typename C>
void Matrix<C>::swapColumns(index_t j0, index_t j1)
{
  assert(j0<base::n_columns() and j1<base::n_columns());
  for (index_t k=0; k<base::n_rows(); ++k)
    std::swap((*this)(k,j0),(*this)(k,j1));
}

template<typename C>
void Matrix_base<C>::eraseRow(index_t i)
{
  assert(i<d_rows);
  typename Vector<C>::iterator first = d_data.begin() + i*d_columns;
  d_data.erase(first,first+d_columns);
  --d_rows;
}

template<typename C>
void Matrix_base<C>::eraseColumn(index_t j)
{
  assert(j<d_columns);
  typename Vector<C>::iterator
    pos = d_data.begin() + j; // position of entry |M(0,j)|
  --d_columns; // already adjust number of columns

  // kill individual entries from left to right, with shifts of new |d_colums|
  // (this is easily seen to be correct but not particularly efficient)
  for (index_t k=0; k<d_rows; ++k, pos += d_columns)
    d_data.erase(pos);
}



} // |namespace matrix|


/*****************************************************************************

        Chapter III --- Functions declared in matrix.h

******************************************************************************/

namespace matrix {

  /*

    Instantiation of class templates (only these class instances are generated)

  */

 // type abbreviations used in these instantiations
using Num = long long; // maybe |arithmetic::Numer_t Num|
typedef polynomials::Polynomial<int> Pol;
typedef arithmetic::big_int bigint;

template class Vector<int>;           // the main instance used
template class Vector<signed char>;   // used inside root data
template class Vector<unsigned long>; // for |abelian::Homomorphism|
template class Vector<Num>;           // numerators of rational vectors
template class Matrix_base<signed char>;   // used inside root data
template class Matrix_base<int>;
template class Matrix_base<Num>;
template class Matrix_base<bigint>;
template class Matrix_base<unsigned long>; // for |abelian::Endomorphism|
template class Matrix<int>;           // the main instance used
template class Matrix<short int>;     // for inverse Cartan matrix
template class Matrix<Num>;
template class Matrix<bigint>;
template class Matrix<arithmetic::Split_integer>; // KL matrices eval'd at |s|
template class PID_Matrix<int>;
template class PID_Matrix<short>; // for inverting Cartan matrix
template class PID_Matrix<Num>; // for precision system solving (alcoves)

// template member instances

template Vector<int>& operator/=(Vector<int>&,int);
template Vector<int>& divide (Vector<int>&,int);
template Vector<int>& operator%=(Vector<int>&,int);
template Vector<Num>& operator/=(Vector<Num>&,Num);

template PID_Matrix<int> inverse (PID_Matrix<int> A, arithmetic::big_int& d);
template PID_Matrix<short> inverse (PID_Matrix<short> A, arithmetic::big_int& d);

template class Vector<Pol>;
template class Matrix_base<Pol>;
template class Matrix<Pol>;

template Vector<Pol> Matrix<Pol>::operator*(Vector<Pol> const&) const;

} // |namespace matrix|

} // |namespace atlas|
