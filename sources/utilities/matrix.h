/*!
\file
  This is matrix.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef MATRIX_H  /* guard against multiple inclusions */
#define MATRIX_H

#include <vector>
#include <stdexcept>

#include "matrix_fwd.h"

#include "setutils.h"
#include "tags.h"

namespace atlas {

/******** function definitions ***********************************************/

namespace matrix {

template<typename C>
  void columnVectors(std::vector<Vector<C> >& b,
		     const Matrix<C>& m);

template<typename C> Matrix<C> row_matrix(Vector<C> v);
template<typename C> Matrix<C> column_matrix(Vector<C> v);

template<typename C>
Matrix<C>& conjugate(Matrix<C>&, const Matrix<C>&);

template<typename C>
void extractBlock(Matrix<C>&, const Matrix<C>&, size_t, size_t, size_t,
		  size_t);

template<typename C>
void extractMatrix(Matrix<C>&, const Matrix<C>&,
		   const std::vector<size_t>&, const std::vector<size_t>&);

template<typename C>
  void identityMatrix(Matrix<C>&, size_t);

template<typename C>
  std::vector<Vector<C> > standard_basis(size_t n);

template<typename C>
  void initBasis(std::vector<Vector<C> >& v, size_t n)
 { v=standard_basis<C>(n); }

template<typename C>
Matrix<C> inv_conjugated(const Matrix<C>& M, const Matrix<C>& g);

template<typename C>
Matrix<C>& invConjugate(Matrix<C>&, const Matrix<C>&);

}

/******** type definitions ***************************************************/

namespace matrix {

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

    C scalarProduct (const Vector<C>&) const;
    C dot(const Vector<C>& v) const { return scalarProduct(v); } // alias
    bool isZero() const;

    Vector<C> operator+ (const Vector<C>& v) const
      { Vector<C> result(*this); return result +=v; }
    Vector<C> operator- (const Vector<C>&v) const
      { Vector<C> result(*this); return result -=v; }
    Vector<C> operator* (C c) const
      { Vector<C> result(*this); return result *=c; }
    Vector<C> operator/ (C c) const
      { Vector<C> result(*this); return result /=c; }
    Vector<C> operator- () const
      { Vector<C> result(*this); return result.negate(); }


}; // Vector


template<typename C> class Matrix {
  /*!
  We implement a matrix simply as a vector of elements, concatenating the
  rows.
  */
 public:

  typedef std::pair<size_t,size_t> index_pair;

 private:

  Vector<C> d_data;
  size_t d_rows;
  size_t d_columns;

 public:

// iterators
  typedef typename Vector<C>::iterator iterator;
  typedef typename Vector<C>::const_iterator const_iterator;

  iterator begin() {
    return d_data.begin();
  }

  iterator end() {
    return d_data.end();
  }

  const_iterator begin() const {
    return d_data.begin();
  }

  const_iterator end() const {
    return d_data.end();
  }

// constructors and destructors
  Matrix()
    {}

  Matrix(size_t m, size_t n)
    :d_data(m*n),d_rows(m),d_columns(n)
    {}

  Matrix(size_t m, size_t n, const C& c)
    :d_data(m*n,c),d_rows(m),d_columns(n)
    {}

  explicit Matrix(size_t n)
    :d_data(n*n),d_rows(n),d_columns(n)
    {}

  explicit Matrix(const std::vector<Vector<C> >&); // ctor from column vectors

  Matrix(const Matrix<C> &, const std::vector<Vector<C> >&); // on other basis

  Matrix(const Matrix<C> &, size_t, size_t, size_t, size_t);

  template<typename I> Matrix(const Matrix<C>&, const I&, const I&);

  template<typename I> Matrix(const I&, const I&, tags::IteratorTag);

  // accessors
  size_t numRows() const { return d_rows; }
  size_t numColumns() const { return d_columns; }
  size_t rowSize() const { return d_columns;  }
  size_t columnSize() const { return d_rows; }

  const C& operator() (size_t i, size_t j) const
  { return d_data[i*d_columns+j]; }

  Vector<C> column(size_t j) const { Vector<C> c; set_column(c,j); return c; }
  std::vector<Vector<C> > columns() const
  {std::vector<Vector<C> > result; columnVectors(result,*this); return result; }
  Vector<C> row(size_t i) const { Vector<C> r; set_row(r,i); return r; }
  std::vector<Vector<C> > rows() const
  {
    std::vector<Vector<C> > result(numRows());
    for (size_t i=0; i<numRows(); ++i)
      set_row(result[i],i);
    return result;
  }

  bool operator== (const Matrix<C>&) const;

  index_pair absMinPos(size_t i_min = 0, size_t j_min = 0) const;

  void apply(Vector<C>&, const Vector<C>&) const;
  Vector<C> apply(const Vector<C>&) const; //functional version
  Vector<C> right_apply(const Vector<C>&) const;

  template<typename I, typename O> void apply(const I&, const I&, O) const;

  void set_column(Vector<C>&, size_t) const;
  void set_row(Vector<C>&, size_t) const;


  bool divisible(C) const;

  Matrix<C> inverse() const
  {
    Matrix<C> result(*this); result.invert(); return result;
  }

  Matrix<C> inverse(C& d) const
  {
    Matrix<C> result(*this); result.invert(d); return result;
  }

  bool isEmpty() const {
    return d_data.size() == 0;
  }

  bool isZero(size_t i_min = 0, size_t j_min = 0) const;


  Matrix<C> transposed() const
  {
    Matrix<C> result(*this); result.transpose(); return result;
  }

  Matrix<C> negative_transposed() const
  {
    Matrix<C> result(*this); result.negate(); result.transpose();
    return result;
  }

  Matrix<C> on_basis(const std::vector<Vector<C> >& basis) const
  {
    return Matrix<C>(*this,basis);
  }

// manipulators
  C& operator() (size_t i, size_t j) {
    return d_data[i*d_columns+j];
  }

  Matrix<C>& operator+= (const Matrix<C>&);

  Matrix<C>& operator-= (const Matrix<C>&);

  Matrix<C>& operator*= (const Matrix<C>&);

  Matrix<C> operator* (const Matrix<C>&) const;

  Matrix<C>& leftMult (const Matrix<C>& p) { return *this=p * *this; }

  Matrix<C>& operator/= (const C& c) throw (std::runtime_error);

  void changeColumnSign(size_t);

  void changeRowSign(size_t);

  void columnOperation(size_t, size_t, const C&);

  void copy(const Matrix<C>&, size_t r = 0, size_t c = 0);

  void copyColumn(const Vector<C>&, size_t);

  void copyRow(const Vector<C>&, size_t);

  void eraseColumn(size_t);

  void eraseRow(size_t);

  void invert();

  void invert(C& d);

  void permute(const setutils::Permutation& a);

  void negate();

  void reset() {
    d_data.assign(d_data.size(),0);
  }

  void resize(size_t, size_t);

  void resize(size_t, size_t, const C&);

  void rowOperation(size_t, size_t, const C&);

  void swap(Matrix&);

  void swapColumns(size_t, size_t);

  void swapRows(size_t, size_t);

  void transpose();
};

}

}

#include "matrix_def.h"

#endif
