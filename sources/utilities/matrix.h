/*
  This is matrix.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef MATRIX_H  /* guard against multiple inclusions */
#define MATRIX_H

#include <cstddef> // for |std::size_t|
#include <vector>
#include <functional> // for |std::reference_wrapper|
#include <stdexcept>
#include <cassert>

#include "matrix_fwd.h"

#include "tags.h"
#include "bigint.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include "constants.h"
#endif

namespace atlas {
namespace matrix {

/******** function definitions ***********************************************/


template<typename C>
  std::vector<Vector<C> > standard_basis(unsigned int n);

template<typename C>
  void initBasis(std::vector<Vector<C> >& v, unsigned int n)
  { v=standard_basis<C>(n); }

template<typename C>
  void row_apply(Matrix<C>& A, const Matrix<C>& ops, // an $r\times r$ matrix
		 unsigned int i); // apply |ops| to rows |i| up to |i+r| of |A|

template<typename C>
  void column_apply(Matrix<C>& A, const Matrix<C>& ops, // an $r\times r$ matrix
		    unsigned int j); // apply |ops| to columns |j| up to |j+r|

template<typename C>
  PID_Matrix<C>& operator+= (PID_Matrix<C>& A, C c); // |A=A+c|, avoiding copy

// non-destructive form takes value parameter, which maybe selects rvalue copy
template<typename C>
  PID_Matrix<C> operator+ (PID_Matrix<C> A, C c) // add scalar matrix
  { return A += c; }

template<typename C>
  PID_Matrix<C>& operator-= (PID_Matrix<C>& A, C c) { return A += -c; }

template<typename C>
  PID_Matrix<C> operator- (PID_Matrix<C> A, C c) { return A += -c; }

template<typename C>
  PID_Matrix<C> operator- (C c, PID_Matrix<C> A) { return A.negate() += c; }

template<typename C>
  PID_Matrix<C> inverse (PID_Matrix<C> A, arithmetic::big_int& d);

template<typename C> PID_Matrix<C> inverse (PID_Matrix<C> A)
 { arithmetic::big_int d; PID_Matrix<C> result=inverse(std::move(A),d);
    if (not d.is_one())
      throw std::runtime_error("Matrix not invertible over the integers");
    return result;
  }

template<typename C>
  void swap(Matrix_base<C>&,Matrix_base<C>&);

/******** type definitions ***************************************************/

template <typename C> class Matrix;

template<typename C>
  class Vector : public std::vector<C>
{
  typedef std::vector<C> base;
public:
  Vector () : base() {}
  explicit Vector (std::size_t n) : base(n) {} // entries remain undefined!
  explicit Vector (const base& b) : base(b) {} // std::vector -> Vector
  Vector (std::size_t n,C c) : base(n,c) {}
  template<typename I>
    Vector (I b, I e) : base(b,e) {} // construction from iterators

  Vector& operator+= (const Vector&);
  Vector& operator-= (const Vector&);
  Vector& operator*= (C);
  Vector& negate (); // negates argument in place
  Vector& negate_add (const Vector& y); // reversed -=, so *this = y - *this

  // the following two methods do not take an end iterator; count is |size()|
  template<typename I> Vector& add(I b,C c) // add |Vector(b,b+size())*c|
  {
    if (c!=C(0)) // we may be able to avoid doing anything at all
      for (auto it=base::begin(); it!=base::end(); ++it,++b)
	*it += *b * c;
    return *this;
  }
  template<typename I> Vector& subtract( I b,C c) { return add(b,-c); }


  template<typename C1> C1 dot (const Vector<C1>& v) const
  {
    assert(base::size()==v.size());
    C1 result = 0;
    for (size_t i=0; i<base::size(); ++i)
      result += (*this)[i] * v[i];
    return result;
  }

  bool isZero() const;

  Vector operator+ (const Vector& v) const &
    { Vector result(*this); return result +=v; }
  Vector operator- (const Vector&v) const &
    { Vector result(*this); return result -=v; }
  Vector operator* (C c) const &
    { Vector result(*this); return result *=c; }
  Vector operator- () const &
    { Vector result(*this); return result.negate(); }

  // when left operand is rvalue reference, use destructive operators
  Vector operator+ (const Vector& v) && { return *this +=v; }
  Vector operator- (const Vector& v) && { return *this -=v; }
  Vector operator* (C c) &&  { return *this *=c; }
  Vector operator- () &&     { return negate(); }

  // when right operand is rvalue reference, use destructive operators too
  Vector operator+ (Vector&& v) const & { return v+=*this; }
  Vector operator- (Vector&& v) const & { return v.negate_add(*this); }

  // when both operands are rvalue references, give priority to the left
  Vector operator+ (Vector&& v) && { return *this +=v; }
  Vector operator- (Vector&& v) && { return *this -=v; }

  template<typename C1>
    Vector<C1> scaled (C1 c) const // like operator*, but forces type to C1
  {
    Vector<C1> result(base::size());
    auto p=result.begin();
    for (auto it=base::begin(); it!=base::end(); ++p,++it)
      *p = *it * c;
    return result;
  }

  Matrix<C> row_matrix() const;    // vector as 1-row matrix
  Matrix<C> column_matrix() const; // vector as 1-column matrix
}; // Vector

// these are external functions, to allow instantiating |Vector<Pol>|
template<typename C>
  Vector<C>& operator/= (Vector<C>& v, C c); // this division must be exact
template<typename C>
  Vector<C>& divide (Vector<C>& v, C c); // this integer division rounds down
template<typename C>
  Vector<C>& operator%= (Vector<C>& v, C c);

// a by-value operand |v| will be copy-constructed only if given by an lvalue
template<typename C>
  Vector<C> operator/ (Vector<C> v,C c) { return v /= c; }


template<typename C> class Matrix_base
{
  unsigned int d_rows;
  unsigned int d_columns;
 protected: // derived classes may sometimes need to acces this
  Vector<C> d_data;   // vector of elements, concatenated by rows

 public:

// constructors
  Matrix_base(): d_rows(0),d_columns(0),d_data() {}

  Matrix_base(unsigned int m, unsigned int n)
    : d_rows(m),d_columns(n),d_data(static_cast<std::size_t>(m)*n) {}
  Matrix_base(unsigned int m, unsigned int n, const C& c)
    : d_rows(m), d_columns(n), d_data(static_cast<std::size_t>(m)*n,c) {}

  Matrix_base (const std::vector<Vector<C> >& b,
	       unsigned int n_rows); // with explicit #rows in case |b| empty

  template<typename I> // from sequence of columns obtained via iterator
    Matrix_base(I begin, I end, unsigned int n_rows, tags::IteratorTag)
      : d_rows(n_rows), d_columns(std::distance(begin,end))
    , d_data(static_cast<std::size_t>(d_rows)*d_columns)
  {
    I p=begin;
    for (unsigned int j=0; j<d_columns; ++j,++p)
      for (unsigned int i=0; i<d_rows; ++i)
	(*this)(i,j) = (*p)[i];
  }


  void swap(Matrix_base&);

  // accessors
  unsigned int numRows() const { return d_rows; }
  unsigned int numColumns() const { return d_columns; }
  unsigned int rowSize() const { return d_columns;  }
  unsigned int columnSize() const { return d_rows; }

  const C& operator() (unsigned int i,unsigned int j)
    const { return d_data[static_cast<std::size_t>(i)*d_columns+j]; }

  void get_row(Vector<C>&, unsigned int) const;
  void get_column(Vector<C>&, unsigned int) const;

  Vector<C> row(unsigned int i) const { return Vector<C>(at(i,0),at(i+1,0)); }
  std::vector<Vector<C> > rows() const;
  Vector<C> column(unsigned int j) const
    { Vector<C> c; get_column(c,j); return c; }
  std::vector<Vector<C> > columns() const;

  Vector<C> partial_row(unsigned int i, unsigned int j, unsigned int l) const
    { return Vector<C>(at(i,j),at(i,l)); }
  Vector<C> partial_column(unsigned int j, unsigned int i, unsigned int k) const
  { Vector<C> result(k-i);
    for (auto it=result.begin(); it!=result.end(); ++it)
      *it = (*this)(i++,j);
    return result;
  }

  bool operator== (const Matrix_base<C>&) const;
  bool operator!= (const Matrix_base<C>& m) const {return not(operator==(m)); }

  bool is_zero() const; // whether all entries are zero
  bool isEmpty() const { return d_data.size() == 0; }
// manipulators
  C& operator() (unsigned int i, unsigned int j)
   { return d_data[static_cast<std::size_t>(i)*d_columns+j]; }

  void set_row(unsigned int,const Vector<C>&);
  void set_column(unsigned int,const Vector<C>&);
  void add_row(const Vector<C>&);
  void add_column(const Vector<C>&);

  // in the following two methods |k,l| are top-left indices of submatrix
  void get_block(Matrix_base<C>& dst, unsigned k, unsigned l) const;
  void set_block(unsigned k, unsigned l, const Matrix_base<C>& src);

// resize undefined; use |Matrix<C>(m,n).swap(M)| instead of |M.resize(m,n)|

  void eraseColumn(unsigned int);
  void eraseRow(unsigned int);
  void clear() { d_rows=d_columns=0; d_data.clear(); }

 protected: // not |private| because |PID_Matrix<C>::block| uses them
  const C* at (unsigned int i,unsigned int j) const { return &operator()(i,j); }
  C* at (unsigned int i,unsigned int j)             { return &operator()(i,j); }
}; // |template<typename C> class Matrix_base|

template<typename C> class Matrix : public Matrix_base<C>
{
 protected:
   typedef Matrix_base<C> base;

 public:
// constructors
  Matrix() : base() {}
  Matrix(const base& M) : base(M) {}
  Matrix(base&& M) : base(std::move(M)) {}

  using base::base; // inherit all constructors
  explicit Matrix(unsigned int n); // square identity matrix

// accessors
  // no non-destructive additive matrix operations; no need for them (yet?)

  Matrix transposed() const &;
  Matrix transposed() &&
    { return Matrix(std::move(*this)).transpose(); }
  Matrix negative_transposed() const & { return transposed().negate(); }
  Matrix negative_transposed() &&
    { return Matrix(std::move(*this)).negate().transpose(); }

  // there is no point in trying to do these matrix multiplications in-place
  template<typename C1> Vector<C1> operator* (const Vector<C1>&) const;
  template<typename C1> Vector<C1> right_prod(const Vector<C1>&) const;

  template<typename C1> void apply_to(Vector<C1>& v) const
    { v= operator*(v); }
  template<typename C1> void right_mult(Vector<C1>& v) const
    { v= right_prod(v); }

  Matrix operator* (const Matrix&) const;

// manipulators

  Matrix& operator+= (const Matrix&);
  Matrix& operator-= (const Matrix&);
  Matrix& operator*= (const Matrix& Q) { return *this = *this * Q; }
  Matrix& leftMult (const Matrix& P)   { return *this = P * *this; }

  Matrix& operator*= (C c) { base::d_data *= c; return *this; }
  Matrix& negate() { base::d_data.negate(); return *this; }
  Matrix& transpose();

  // secondary manipulators

  // |row(i) += c+row(k)|:
  void rowOperation(unsigned int i, unsigned int k, const C& c);
  // |col(j) += c*col(k)|:
  void columnOperation(unsigned int j, unsigned int k, const C& c);

  void rowMultiply(unsigned int i, C f);
  void columnMultiply(unsigned int j,C f);

  void swapColumns(unsigned int, unsigned int);
  void swapRows(unsigned int, unsigned int);

}; // |template<typename C> class Matrix|

// The following derived class is for integer types only
template<typename C> class PID_Matrix : public Matrix<C>
{
  typedef Matrix<C> base;

 public:
  PID_Matrix() : base() {}
  PID_Matrix(const base& M) : base(M) {}
  PID_Matrix(base&& M) : base(std::move(M)) {}

  // forward constructors to Matrix
  using base::base; // inherit all constructors

// manipulators
  PID_Matrix& operator+= (const Matrix<C>& M)
    { base::operator+=(M); return *this; }
  PID_Matrix& operator-= (const Matrix<C>& M)
    { base::operator-=(M); return *this; }
  PID_Matrix& operator*= (const Matrix<C>& Q)
    { base::operator*=(Q); return *this; }
  PID_Matrix& leftMult (const Matrix<C>& P)
    { base::leftMult(P); return *this; }

  PID_Matrix& operator/= (const C& c)
  { if (c != C(1)) base::base::d_data /= c; return *this; }

  PID_Matrix& negate() { base::negate(); return *this; }
  PID_Matrix& transpose() { base::transpose(); return *this; }

// accessors
  PID_Matrix transposed() const  & { return PID_Matrix(base::transposed()); }
  PID_Matrix transposed() &&
    { return PID_Matrix(std::move(*this)).transpose(); }
  PID_Matrix negative_transposed() const &
    { return PID_Matrix(base::negative_transposed()); }
  PID_Matrix negative_transposed() &&
    { return PID_Matrix(std::move(*this)).negate().transpose(); }
  PID_Matrix inverse() const & { return matrix::inverse(*this); }
  PID_Matrix inverse() &&      { return matrix::inverse(std::move(*this)); }
  PID_Matrix inverse(arithmetic::big_int& d) const &
    { return matrix::inverse(*this,d); }
  PID_Matrix inverse(arithmetic::big_int& d) &&
    { return matrix::inverse(std::move(*this),d); }

  using base::operator*;

  PID_Matrix operator* (const Matrix<C>& Q) const
    { return PID_Matrix(base::operator*(Q)); }


  PID_Matrix on_basis(const PID_Matrix<C>& basis) const; // change of basis

  // secondary accessors (for inversion algorithms)
  bool divisible(C) const;
  PID_Matrix block(unsigned int i0, unsigned int j0, // upper left
		   unsigned int i1, unsigned int j1  // lower right
		   ) const;
  PID_Matrix transposed_block // same as |transposed().block(i0,j0,i1,j1)|
    (unsigned int i0, unsigned int j0, // upper left
     unsigned int i1, unsigned int j1  // lower right
     ) const;

}; // |class PID_Matrix|

// auxiliary class to enable making containers of vector references
template<typename C> class Vector_cref
  : public std::reference_wrapper<Vector<C> const>
{
  typedef std::reference_wrapper<Vector<C> const>  base;
public:
  using base::base; // inherit all constructors
  Vector_cref (const Vector_cref& v) = default;

  C operator[] (std::size_t i) const { return base::get()[i]; }
}; // |class Vector_cref|

// instantiations of templated free functions and methods

// Set |b| to the canonical basis (as in identity matrix) in dimension |r|
template<typename C>
  std::vector<Vector<C> > standard_basis(unsigned int r)
{
  std::vector<Vector<C> > result(r,Vector<C>(r,C(0)));

  for (unsigned int i=0; i<r; ++i)
    result[i][i] = C(1);

  return result;
}

template<typename C>
  void swap(Matrix_base<C>& A,Matrix_base<C>& B) { A.swap(B); }

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

  for (unsigned int i=0; i<base::numRows(); ++i)
  {
    C1 c(0);
    for (unsigned int j=0; j<base::numColumns(); ++j)
      c += (*this)(i,j) * w[j];
    result[i] = c;
  }

  return result;
}

/*
  Multiply the matrix to right to the row-vector w, and returns the result.
  It is assumed that the size of w is the number of rows; result size is the
  number of columns. This is the proper sense of application for dual space.
*/
template<typename C> template<typename C1>
Vector<C1> Matrix<C>::right_prod(const Vector<C1>& w) const
{
  assert(base::numRows()==w.size());
  Vector<C1> result(base::numColumns());

  for (unsigned int j=0; j<base::numColumns(); ++j)
  {
    C1 c(0);
    for (unsigned int i=0; i<base::numRows(); ++i)
      c += w[i] * (*this)(i,j);
    result[j] = c;
  }

  return result;
}

template<typename C>
  void row_apply(Matrix<C>& A, const Matrix<C>& ops,
		 unsigned int i) // initial row
{ const auto r=ops.numRows();
  assert(r==ops.numColumns());
  assert(i+r <= A.numRows());
  Vector<C> tmp(r);
  for (unsigned int j=0; j<A.numColumns(); ++j)
  { // |tmp = partial_column(j,i,i+r)|; save values before overwriting
    for (unsigned int k=0; k<r; ++k)
      tmp[k]=A(i+k,j);
    ops.apply_to(tmp); // do row operations on column vector |tmp|
    for (unsigned int k=0; k<r; ++k)
      A(i+k,j)=tmp[k];
  }
}


template<typename C>
  void column_apply(Matrix<C>& A, const Matrix<C>& ops,
		    unsigned int j) // initial column of |A| to act on
{ const auto r=ops.numRows();
  assert(r==ops.numColumns());
  assert(j+r <= A.numColumns());
  Vector<C> tmp(r);
  for (unsigned int i=0; i<A.numRows(); ++i)
  { // |tmp = partial_row(i,j,j+r)|; save values before overwriting
    for (unsigned int l=0; l<r; ++l)
      tmp[l]=A(i,j+l);
    ops.right_mult(tmp); // do column operations on row vector |tmp|
    for (unsigned int l=0; l<r; ++l)
      A(i,j+l)=tmp[l];
  }
} // |column_apply|

template<typename C>
  PID_Matrix<C>& operator+= (PID_Matrix<C>& A, C c) // |A=A+c|, avoiding copy
{
  unsigned int i=std::min(A.numRows(),A.numColumns());
  while (i-->0)
    A(i,i) += c;
  return A;
}


} // |namespace matrix|

} // |namespace atlas|

#endif
