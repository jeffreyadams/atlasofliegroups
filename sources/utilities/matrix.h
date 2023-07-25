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
#include <cstdint> // for |uint32_t|
#include <vector>
#include <functional> // for |std::reference_wrapper|
#include <limits>
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

/*				Type definition
   A deliberate chioce (for now) is to use only 32-bit row/column indices in
   matrices; this saves some space but makes extremely deep or wide matrices
   (for storing huge collections of weights) impossible.
*/
using index_t = std::uint32_t;
constexpr auto index_max = std::numeric_limits<matrix::index_t>::max();


/******** function definitions ***********************************************/


template<typename C>
  std::vector<Vector<C> > standard_basis(index_t n);

template<typename C>
  void initBasis(std::vector<Vector<C> >& v, index_t n)
  { v=standard_basis<C>(n); }

template<typename C>
  void row_apply(Matrix<C>& A, const Matrix<C>& ops, // an $r\times r$ matrix
		 index_t i); // apply |ops| to rows |i| up to |i+r| of |A|

template<typename C>
  void column_apply(Matrix<C>& A, const Matrix<C>& ops, // an $r\times r$ matrix
		    index_t j); // apply |ops| to columns |j| up to |j+r|

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
    if (d!=1)
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

  bool is_zero() const;

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
template<typename C>
  Vector<C> operator% (Vector<C> v,C c) { return v %= c; }

template<typename C> class Matrix_base
{
  index_t d_rows;
  index_t d_columns;
protected: // derived classes may sometimes need to access this
  Vector<C> d_data;   // vector of elements, concatenated by rows

public:

// constructors
  Matrix_base(): d_rows(0),d_columns(0),d_data() {}

  Matrix_base(index_t m, index_t n)
    : d_rows(m),d_columns(n),d_data(static_cast<std::size_t>(m)*n) {}
  Matrix_base(index_t m, index_t n, const C& c)
    : d_rows(m), d_columns(n), d_data(static_cast<std::size_t>(m)*n,c) {}

  Matrix_base (const std::vector<Vector<C> >& b,
	       index_t n_rows); // with explicit #rows in case |b| empty

  template<typename I> // from sequence of columns obtained via iterator
    Matrix_base(I begin, I end, index_t n_rows, tags::IteratorTag)
      : d_rows(n_rows), d_columns(std::distance(begin,end))
    , d_data(static_cast<std::size_t>(d_rows)*d_columns)
  {
    I p=begin;
    for (index_t j=0; j<d_columns; ++j,++p)
      for (index_t i=0; i<d_rows; ++i)
	(*this)(i,j) = (*p)[i];
  }


  void swap(Matrix_base&);

  // accessors
  index_t n_rows() const { return d_rows; }
  index_t n_columns() const { return d_columns; }

  const C& operator() (index_t i,index_t j) const
   { return d_data[static_cast<std::size_t>(i)*d_columns+j]; }

  void get_row(Vector<C>&, index_t) const;
  void get_column(Vector<C>&, index_t) const;

  Vector<C> row(index_t i) const { return Vector<C>(at(i,0),at(i+1,0)); }
  std::vector<Vector<C> > rows() const;
  Vector<C> column(index_t j) const
    { Vector<C> c; get_column(c,j); return c; }
  std::vector<Vector<C> > columns() const;

  Vector<C> partial_row(index_t i, index_t j, index_t l) const
    { return Vector<C>(at(i,j),at(i,l)); }
  Vector<C> partial_column(index_t j, index_t i, index_t k) const
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
  C& operator() (index_t i, index_t j)
   { return d_data[static_cast<std::size_t>(i)*d_columns+j]; }

  void set_row(index_t,const Vector<C>&);
  void set_column(index_t,const Vector<C>&);
  void add_row(const Vector<C>&);
  void add_column(const Vector<C>&);

  // in the following two methods |k,l| are top-left indices of submatrix
  void get_block(Matrix_base<C>& dst, index_t k, index_t l) const;
  void set_block(index_t k, index_t l, const Matrix_base<C>& src);

// resize undefined; use |Matrix<C>(m,n).swap(M)| instead of |M.resize(m,n)|

  void eraseColumn(index_t);
  void eraseRow(index_t);
  void clear() { d_rows=d_columns=0; d_data.clear(); }

 protected: // not |private| because |PID_Matrix<C>::block| uses them
  const C* at (size_t i,size_t j) const { return &operator()(i,j); }
  C* at (size_t i,size_t j)             { return &operator()(i,j); }
}; // |template<typename C> class Matrix_base|

template<typename C> class Matrix : public Matrix_base<C>
{
protected:
  using base = Matrix_base<C>;

public:
// constructors
  Matrix() : base() {}
  Matrix(const base& M) : base(M) {}
  Matrix(base&& M) : base(std::move(M)) {}

  using base::base; // inherit all constructors
  explicit Matrix(index_t n); // square identity matrix

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
  void rowOperation(index_t i, index_t k, const C& c);
  // |col(j) += c*col(k)|:
  void columnOperation(index_t j, index_t k, const C& c);

  void rowMultiply(index_t i, C f);
  void columnMultiply(index_t j,C f);

  void swapColumns(index_t, index_t);
  void swapRows(index_t, index_t);

}; // |template<typename C> class Matrix|

// The following derived class is for integer types only
template<typename C> class PID_Matrix : public Matrix<C>
{
  using base = Matrix<C>;

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
  template<typename D> PID_Matrix<D> entry_type_convert() const;
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
  PID_Matrix block(index_t i0, index_t j0, // upper left
		   index_t i1, index_t j1  // lower right
		   ) const;
  // return a block of the transposed matrix, without calling |transposed|
  PID_Matrix transposed_block // same as |transposed().block(i0,j0,i1,j1)|
    (index_t i0, index_t j0, // upper left
     index_t i1, index_t j1  // lower right
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
  std::vector<Vector<C> > standard_basis(index_t r)
{
  std::vector<Vector<C> > result(r,Vector<C>(r,C(0)));

  for (index_t i=0; i<r; ++i)
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
  assert(base::n_columns()==w.size());
  Vector<C1> result(base::n_rows());

  for (index_t i=0; i<base::n_rows(); ++i)
  {
    C1 c(0);
    for (index_t j=0; j<base::n_columns(); ++j)
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
  assert(base::n_rows()==w.size());
  Vector<C1> result(base::n_columns());

  for (index_t j=0; j<base::n_columns(); ++j)
  {
    C1 c(0);
    for (index_t i=0; i<base::n_rows(); ++i)
      c += w[i] * (*this)(i,j);
    result[j] = c;
  }

  return result;
}

template<typename C>
  void row_apply(Matrix<C>& A, const Matrix<C>& ops, index_t i) // initial row
{ const auto r=ops.n_rows();
  assert(r==ops.n_columns());
  assert(i+r <= A.n_rows());
  Vector<C> tmp(r);
  for (index_t j=0; j<A.n_columns(); ++j)
  { // |tmp = partial_column(j,i,i+r)|; save values before overwriting
    for (index_t k=0; k<r; ++k)
      tmp[k]=A(i+k,j);
    ops.apply_to(tmp); // do row operations on column vector |tmp|
    for (index_t k=0; k<r; ++k)
      A(i+k,j)=tmp[k];
  }
}


template<typename C>
  void column_apply(Matrix<C>& A, const Matrix<C>& ops,
		    index_t j) // initial column of |A| to act on
{ const auto r=ops.n_rows();
  assert(r==ops.n_columns());
  assert(j+r <= A.n_columns());
  Vector<C> tmp(r);
  for (index_t i=0; i<A.n_rows(); ++i)
  { // |tmp = partial_row(i,j,j+r)|; save values before overwriting
    for (index_t l=0; l<r; ++l)
      tmp[l]=A(i,j+l);
    ops.right_mult(tmp); // do column operations on row vector |tmp|
    for (index_t l=0; l<r; ++l)
      A(i,j+l)=tmp[l];
  }
} // |column_apply|

template<typename C> template<typename D>
  PID_Matrix<D> PID_Matrix<C>::entry_type_convert() const
{
  PID_Matrix<D> result(this->n_rows(),this->n_columns());
  // we cannot access |result.d_data| directly, but |&result(0,0)| points there
  // though we access underlying vector, no assumption of layout is made here
  std::copy(this->d_data.begin(),this->d_data.end(),&result(0,0));
  return result;
}

template<typename C>
  PID_Matrix<C>& operator+= (PID_Matrix<C>& A, C c) // |A=A+c|, avoiding copy
{
  index_t i=std::min(A.n_rows(),A.n_columns());
  while (i-->0)
    A(i,i) += c;
  return A;
}


} // |namespace matrix|

} // |namespace atlas|

#endif
