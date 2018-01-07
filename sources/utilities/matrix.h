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
  template<typename I> Vector& add(I b,C c); // add |Vector(b,b+size())*c|
  template<typename I> Vector& subtract( I b,C c) { return add(b,-c); }


  template<typename C1> C1 dot (const Vector<C1>& v) const;
  bool isZero() const;

#if __GNUC__ < 4 || \
  __GNUC__ == 4 && ( __GNUC_MINOR__ < 8 || \
                     __GNUC_MINOR__ == 8 && __GNUC_PATCHLEVEL__ < 1)
  Vector operator+ (const Vector& v) const
    { Vector result(*this); return result +=v; }
  Vector operator- (const Vector&v) const
    { Vector result(*this); return result -=v; }
  Vector operator* (C c) const
    { Vector result(*this); return result *=c; }
  Vector operator- () const
    { Vector result(*this); return result.negate(); }
#else
  Vector operator+ (const Vector& v) const &
    { Vector result(*this); return result +=v; }
  Vector operator- (const Vector&v) const &
    { Vector result(*this); return result -=v; }
  Vector operator* (C c) const &
    { Vector result(*this); return result *=c; }
  Vector operator- () const &
    { Vector result(*this); return result.negate(); }

  Vector operator+ (const Vector& v) &&
  { Vector result(std::move(*this)); return result +=v; }
  Vector operator- (const Vector&v)  &&
    { Vector result(std::move(*this)); return result -=v; }
  Vector operator* (C c) &&
    { Vector result(std::move(*this)); return result *=c; }
  Vector operator- () &&
    { Vector result(std::move(*this)); return result.negate(); }
#endif

  template<typename C1>
    Vector<C1> scaled (C1 c) const; // like operator*, but forces type to C1

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
    : d_rows(m),d_columns(n),d_data(m*n) {}
  Matrix_base(unsigned int m, unsigned int n, const C& c)
    : d_rows(m), d_columns(n), d_data(m*n,c) {}

  Matrix_base (const std::vector<Vector<C> >& b,
	       unsigned int n_rows); // with explicit #rows in case |b| empty

  template<typename I> // from sequence of columns obtained via iterator
    Matrix_base(I begin, I end, unsigned int n_rows, tags::IteratorTag);

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
  Vector<C> column(unsigned int j) const { Vector<C> c; get_column(c,j); return c; }
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

#if __GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__ >= 8
  // that is, if compiler version is sufficiently new
  using base::base; // inherit all constructors
#else
  Matrix(unsigned int m, unsigned int n) : base(m,n) {}
  Matrix(unsigned int m, unsigned int n, const C& c) : base(m,n,c) {}

  Matrix(const std::vector<Vector<C> >&cols, unsigned int n_rows)
    : base(cols,n_rows) {}

  template<typename I> // from sequence of columns obtained via iterator
    Matrix(I begin, I end, unsigned int n_rows, tags::IteratorTag)
    : base(begin,end,n_rows,tags::IteratorTag()) {}
#endif
  explicit Matrix(unsigned int n); // square identity matrix

// accessors
  // no non-destructive additive matrix operations; no need for them (yet?)

#if __GNUC__ < 4 || \
  __GNUC__ == 4 && ( __GNUC_MINOR__ < 8 || \
                     __GNUC_MINOR__ == 8 && __GNUC_PATCHLEVEL__ < 1)
  // that is, if compiler version is too old
  Matrix transposed() const;
  Matrix negative_transposed() const { return transposed().negate(); }
#else
  Matrix transposed() const &;
  Matrix transposed() &&
    { return Matrix(std::move(*this)).transpose(); }
  Matrix negative_transposed() const & { return transposed().negate(); }
  Matrix negative_transposed() &&
    { return Matrix(std::move(*this)).negate().transpose(); }
#endif

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

  Matrix& negate() { base::d_data.negate(); return *this; }
  Matrix& transpose();

  // secondary manipulators

  void rowOperation(unsigned int i, unsigned int k, const C& c); // |row(i) += c+row(k)|
  void columnOperation(unsigned int j, unsigned int k, const C& c); // |col(j) += c*col(k)|

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
#if __GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__ >= 8
  // that is, if compiler version is sufficiently new
  using base::base; // inherit all constructors
#else
  PID_Matrix(unsigned int m, unsigned int n) : base(m,n) {}
  PID_Matrix(unsigned int m, unsigned int n, const C& c) : base(m,n,c) {}

  explicit PID_Matrix(unsigned int n) : base(n) {} // square identity matrix
  PID_Matrix(const std::vector<Vector<C> >&cols, unsigned int n_rows)
    : base(cols,n_rows) {}

  template<typename I> // from sequence of columns obtained via iterator
    PID_Matrix(I begin, I end, unsigned int n_rows, tags::IteratorTag)
    : base(begin,end,n_rows,tags::IteratorTag()) {}
#endif

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
#if __GNUC__ < 4 || \
  __GNUC__ == 4 && ( __GNUC_MINOR__ < 8 || \
                     __GNUC_MINOR__ == 8 && __GNUC_PATCHLEVEL__ < 1)
  // that is, if compiler version is too old
  PID_Matrix transposed() const  { return PID_Matrix(base::transposed()); }
  PID_Matrix negative_transposed() const
    { return PID_Matrix(base::negative_transposed()); }
  PID_Matrix inverse() const     { return matrix::inverse(*this); }
  PID_Matrix inverse(arithmetic::big_int& d) const
    { return matrix::inverse(*this,d); }
#else
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
#endif

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
#if __GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__ >= 8
  // that is, if compiler version is sufficiently new
  using base::base; // inherit all constructors
#else
  expicit Vector_cref (const Vector<C>& v) noexcept : base(v) {}
  Vector_cref (const base& v) noexcept : base(v) {}
#endif
  Vector_cref (const Vector_cref& v) = default;

  C operator[] (std::size_t i) const { return base::get()[i]; }
}; // |class Vector_cref|

} // |namespace matrix|

} // |namespace atlas|

#endif
