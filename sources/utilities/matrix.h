/*!
\file
  This is matrix.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef MATRIX_H  /* guard against multiple inclusions */
#define MATRIX_H

#include <vector>
#include <stdexcept>

#include "matrix_fwd.h"

#include "tags.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include "constants.h"
#endif

namespace atlas {

/******** function definitions ***********************************************/

namespace matrix {

template<typename C>
  std::vector<Vector<C> > standard_basis(size_t n);

template<typename C>
  void initBasis(std::vector<Vector<C> >& v, size_t n)
  { v=standard_basis<C>(n); }

}

/******** type definitions ***************************************************/

namespace matrix {

template <typename C> class Matrix;

template<typename C>
  class Vector : public std::vector<C>
{
  typedef std::vector<C> base;
public:
  Vector () : base() {}
  explicit Vector (size_t n) : base(n) {} // entries remain undefined!
  explicit Vector (const base& b) : base(b) {} // std::vector -> Vector
  Vector (size_t n,C c) : base(n,c) {}
  template<typename I>
    Vector (I b, I e) : base(b,e) {} // construction from iterators

  Vector<C>& operator+= (const Vector<C>&);
  Vector<C>& operator-= (const Vector<C>&);
  Vector<C>& operator*= (C);
  Vector<C>& operator/= (C) throw (std::runtime_error);
  Vector<C>& negate (); // negates argument in place

  template<typename C1> C1 dot (const Vector<C1>& v) const;
  bool isZero() const;

  Vector<C> operator+ (const Vector<C>& v) const
    { Vector<C> result(*this); return result +=v; }
  Vector<C> operator- (const Vector<C>&v) const
    { Vector<C> result(*this); return result -=v; }
  Vector<C> operator* (C c) const
    { Vector<C> result(*this); return result *=c; }
  Vector<C> operator/ (C c) const // throws |std::runtime_error| if not exact
    { Vector<C> result(*this); return result /=c; }
  Vector<C> operator- () const
    { Vector<C> result(*this); return result.negate(); }

  template<typename C1>
    Vector<C1> scaled (C1 c) const // like operator*, but forces tye to C1
  { Vector<C1> result(base::begin(),base::end()); return result *=c; }


  Matrix<C> row_matrix() const;    // vector as 1-row matrix
  Matrix<C> column_matrix() const; // vector as 1-column matrix
}; // Vector


template<typename C> class Matrix_base
{
  size_t d_rows;
  size_t d_columns;
 protected: // derived classes may sometimes need to acces this
  Vector<C> d_data;   // vector of elements, concatenated by rows

 public:

// constructors
  Matrix_base(): d_rows(0),d_columns(0),d_data() {}

  Matrix_base(size_t m, size_t n) : d_rows(m),d_columns(n),d_data(m*n) {}
  Matrix_base(size_t m, size_t n, const C& c)
    : d_rows(m), d_columns(n), d_data(m*n,c) {}

  Matrix_base
    (const std::vector<Vector<C> >&, size_t n_rows); // with explicit #rows

  template<typename I> // from sequence of columns obtained via iterator
    Matrix_base(I begin, I end, size_t n_rows, tags::IteratorTag);

  void swap(Matrix_base&);

  // accessors
  size_t numRows() const { return d_rows; }
  size_t numColumns() const { return d_columns; }
  size_t rowSize() const { return d_columns;  }
  size_t columnSize() const { return d_rows; }

  const C& operator() (size_t i,size_t j)
    const { return d_data[i*d_columns+j]; }

  void get_row(Vector<C>&, size_t) const;
  void get_column(Vector<C>&, size_t) const;

  Vector<C> row(size_t i) const { Vector<C> r; get_row(r,i); return r; }
  std::vector<Vector<C> > rows() const;
  Vector<C> column(size_t j) const { Vector<C> c; get_column(c,j); return c; }
  std::vector<Vector<C> > columns() const;

  bool operator== (const Matrix_base<C>&) const;
  bool operator!= (const Matrix_base<C>& m) const {return not(operator==(m)); }

  bool isEmpty() const { return d_data.size() == 0; }
// manipulators
  C& operator() (size_t i, size_t j) { return d_data[i*d_columns+j]; }

  void set_row(size_t,const Vector<C>&);
  void set_column(size_t,const Vector<C>&);
  void add_row(const Vector<C>&);
  void add_column(const Vector<C>&);

// resize undefined; use |Matrix<C>(m,n).swap(M)| instead of |M.resize(m,n)|

  void eraseColumn(size_t);
  void eraseRow(size_t);
  void reset() { d_data.assign(d_data.size(),C(0)); }
}; // |template<typename C> class Matrix_base|

template<typename C> class Matrix : public Matrix_base<C>
{
  typedef Matrix_base<C> base;

 public:
// constructors
  Matrix() : base() {}

  Matrix(size_t m, size_t n) : base(m,n) {}
  Matrix(size_t m, size_t n, const C& c) : base(m,n,c) {}

  explicit Matrix(size_t n); // square identity matrix
  Matrix(const std::vector<Vector<C> >&cols, size_t n_rows)
    : base(cols,n_rows) {}

  template<typename I> // from sequence of columns obtained via iterator
    Matrix(I begin, I end, size_t n_rows, tags::IteratorTag)
    : base(begin,end,n_rows,tags::IteratorTag()) {}

// accessors
  Matrix<C> transposed() const
    { Matrix<C> result(*this); result.transpose(); return result; }
  Matrix<C> negative_transposed() const
    { Matrix<C> result(*this); result.negate(); result.transpose();
      return result; }

  template<typename C1> Vector<C1> operator* (const Vector<C1>&) const;

  template<typename C1> Vector<C1> right_mult(const Vector<C1>&) const;
  template<typename C1> void apply_to(Vector<C1>& v) const { v= operator*(v); }

  Matrix<C> operator* (const Matrix<C>&) const;

  Matrix<C> inverse() const
    { Matrix<C> result(*this); result.invert(); return result; }
  Matrix<C> inverse(C& d) const
    { Matrix<C> result(*this); result.invert(d); return result; }

  Matrix<C> on_basis(const std::vector<Vector<C> >& basis) const;

// manipulators

  Matrix<C>& operator+= (const Matrix<C>&);
  Matrix<C>& operator-= (const Matrix<C>&);
  Matrix<C>& operator*= (const Matrix<C>& Q)
    { (*this*Q).swap(*this); return *this; }
  Matrix<C>& leftMult (const Matrix<C>& P)
    { (P**this).swap(*this); return *this; }
  Matrix<C>& operator/= (const C& c) throw (std::runtime_error);

  void negate(){ base::d_data.negate(); }
  void transpose();

  void invert();
  void invert(C& d);

  // secondary accessors (mainly for inversion algorithms)

  bool divisible(C) const;
  Matrix<C> block(size_t i0, size_t j0, size_t i1, size_t j1) const;

  // secondary manipulators

  void rowOperation(size_t, size_t, const C&);
  void columnOperation(size_t j, size_t k, const C& c); // |col(j) += c*col(k)|

  void rowMultiply(size_t i, C f);
  void columnMultiply(size_t j,C f);

  void swapColumns(size_t, size_t);
  void swapRows(size_t, size_t);

}; // |template<typename C> class Matrix|

} // |namespace matrix|

} // |namespace atlas|

#endif
