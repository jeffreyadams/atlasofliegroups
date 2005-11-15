/*
  This is matrix.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef MATRIX_H  /* guard against multiple inclusions */
#define MATRIX_H

#include <vector>

#include "matrix_fwd.h"

#include "setutils.h"
#include "tags.h"

namespace atlas {

/******** function definitions ***********************************************/

namespace matrix {

template<typename C>
  void columnVectors(std::vector<std::vector<C> >& b, 
		     const Matrix<C>& m);

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
  void initBasis(std::vector<std::vector<C> >&, size_t);

template<typename C>
Matrix<C>& invConjugate(Matrix<C>&, const Matrix<C>&);

template<typename C>
Matrix<C>& leftProd(Matrix<C>&, const Matrix<C>&);

}

/******** type definitions ***************************************************/

namespace matrix {

template<typename C> class Matrix {

 public:

  typedef std::pair<size_t,size_t> index_pair;

 private:

  std::vector<C> d_data;
  size_t d_rows;
  size_t d_columns;

 public:

// iterators
  typedef typename std::vector<C>::iterator iterator;
  typedef typename std::vector<C>::const_iterator const_iterator;

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

  explicit Matrix(const std::vector<std::vector<C> >&);
    
  Matrix(const Matrix<C> &, const std::vector<std::vector<C> >&);
    
  Matrix(const Matrix<C> &, size_t, size_t, size_t, size_t);
    
  template<typename I> Matrix(const Matrix<C>&, const I&, const I&);
    
  template<typename I> Matrix(const I&, const I&, tags::IteratorTag);
    
  virtual ~Matrix() 
    {}

// accessors
  const C& operator() (size_t i, size_t j) const {
    return d_data[i*d_columns+j];
  }

  bool operator== (const Matrix<C>&) const;

  index_pair absMinPos(size_t i_min = 0, size_t j_min = 0) const;

  void apply(std::vector<C>&, const std::vector<C>&) const;

  template<typename I, typename O> void apply(const I&, const I&, O) const;

  void column(std::vector<C>&, size_t) const;

  size_t columnSize() const {
    return d_rows;
  }

  bool divisible(C) const;

  bool isEmpty() const {
    return d_data.size() == 0;
  }

  bool isZero(size_t i_min = 0, size_t j_min = 0) const;

  size_t numRows() const {
    return d_rows;
  }

  size_t numColumns() const {
    return d_columns;
  }

  void row(std::vector<C>&, size_t) const;

  size_t rowSize() const {
    return d_columns;
  }

// manipulators
  C& operator() (size_t i, size_t j) {
    return d_data[i*d_columns+j];
  }

  Matrix<C>& operator+= (const Matrix<C>&);

  Matrix<C>& operator-= (const Matrix<C>&);

  Matrix<C>& operator*= (const Matrix<C>&);

  Matrix<C>& operator/= (const C& c);

  void changeColumnSign(size_t);

  void changeRowSign(size_t);

  void columnOperation(size_t, size_t, const C&);

  void copy(const Matrix<C>&, size_t r = 0, size_t c = 0);

  void copyColumn(const Matrix<C>&, size_t, size_t);

  void copyRow(const Matrix<C>&, size_t, size_t);

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
