/*
  This is bitvector.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef BITVECTOR_H  /* guard against multiple inclusions */
#define BITVECTOR_H

#include "bitvector_fwd.h"

#include <vector>

#include "bitset.h"

/******** function declarations **********************************************/

namespace atlas {

namespace bitvector {

template<size_t dim> 
  void combination(BitVector<dim>&, const std::vector<BitVector<dim> >&,
		   const bitset::BitSet<dim>&);

template<size_t dim> 
  void combination(bitset::BitSet<dim>&, 
		   const std::vector<bitset::BitSet<dim> >&,
		   const bitset::BitSet<dim>&);

template<size_t dim>
  void complement(bitset::BitSet<dim>&, const std::vector<BitVector<dim> >&,
		  size_t);

template<size_t dim> 
  bool firstSolution(bitset::BitSet<dim>&, const std::vector<BitVector<dim> >&,
		     const BitVector<dim>&);

template<size_t dim> bool firstSolution(BitVector<dim>&,
					const std::vector<BitVector<dim> >&);

template<size_t dim> void identityMatrix(BitMatrix<dim>&, size_t);

template<size_t dim> void initBasis(std::vector<BitVector<dim> >&, size_t);

template<size_t dim> bool isIndependent(const std::vector<BitVector<dim> >&);

template<size_t dim> 
  void normalize(bitset::BitSet<dim>&, std::vector<BitVector<dim> >&);

template<size_t dim> 
  void normalSpanAdd(std::vector<BitVector<dim> >&, std::vector<size_t>&,
		     const BitVector<dim>&);

template<size_t dim>
  void projection(BitMatrix<dim>& p, const std::vector<BitVector<dim> >& b, 
		  size_t d);

template<size_t dim>
  void reflectionMatrix(BitMatrix<dim>&, const BitVector<dim>&, 
			const BitVector<dim>&);

template<size_t dim>
  void relations(std::vector<BitVector<dim> >&, 
		 const std::vector<BitVector<dim> >&);

template<size_t dim> 
  bool scalarProduct(const BitVector<dim>&, const BitVector<dim>&);

template<size_t dim> 
  void spanAdd(std::vector<BitVector<dim> >&, std::vector<size_t>&,
	       const BitVector<dim>&);

}

/******** type definitions **************************************************/

namespace bitvector {

template<size_t dim> class BitVector{

  friend BitVector<dim>& BitMatrix<dim>::apply(BitVector<dim>&, 
					       const BitVector<dim>&) const;

 private:

  bitset::BitSet<dim> d_data;
  size_t d_size;

 public:

  // constructors and destructors
  BitVector() 
    {}

  explicit BitVector(size_t n)
    :d_size(n) 
    {}

  BitVector(size_t n, size_t j)
    :d_size(n) {
    set(j);
  }

  BitVector(bitset::BitSet<dim> data, size_t n)
    :d_data(data),
     d_size(n) 
    {}

  ~BitVector() 
    {}

// copy and assignment
  BitVector(const BitVector& v)
    :d_data(v.d_data),
     d_size(v.d_size) 
    {}

  BitVector& operator= (const BitVector& v) {
    d_data = v.d_data;
    d_size = v.d_size;
    return *this;
  }

// accessors

  bool operator< (const BitVector& v) const {
    return d_data < v.d_data;
  }

  bool operator[] (size_t i) const {
    return d_data[i];
  }

  size_t count() {
    return d_data.count();
  }

  const bitset::BitSet<dim>& data() const {
    return d_data;
  }

  size_t firstBit() const {
    return d_data.firstBit();
  }

  bool isZero() const {
    return d_data.none();
  }

  bool nonZero() const {
    return d_data.any();
  }

  size_t size() const {
    return d_size;
  }

  bool test(size_t i) const {
    return d_data.test(i);
  }

// manipulators

  BitVector& operator+= (const BitVector& v) {
    d_data ^= v.d_data; 
    return *this;
  }

  BitVector& operator-= (const BitVector& v) { // same thing as +=
    d_data ^= v.d_data; 
    return *this;
  }

  BitVector& operator&= (const BitVector& v) {
    d_data &= v.d_data; 
    return *this;
  }

  BitVector& operator>>= (size_t pos) {
    d_data >>= pos; 
    return *this;
  }

  BitVector& operator<<= (size_t pos) {
    d_data <<= pos; 
    return *this;
  }

  BitVector& flip(size_t i) {
    d_data.flip(i); 
    return *this;
  }

  BitVector& pushBack(bool);

  BitVector& reset() {
    d_data.reset(); 
    return *this;
  }

  BitVector& reset(size_t i) {
    d_data.reset(i); 
    return *this;
  }

  void resize(size_t n) {
    d_size = n;
  } 

  BitVector& set() {
    set(d_data,d_size); 
    return *this;
  }

  BitVector& set(size_t i) {
    d_data.set(i); 
    return *this;
  }

  BitVector& set(size_t i, bool b) {
    d_data.set(i,b); 
    return *this;
  }

  void slice(const bitset::BitSet<dim>&);
};

template<size_t dim> class FirstBit {

 public:

  typedef const BitVector<dim>& argument_type;
  typedef size_t result_type;

  result_type operator() (argument_type v) const {
    return v.firstBit();
  }
};

// note that the elements in d_data are the _column_ vectors of the
// matrix

template<size_t dim> class BitMatrix {

 private:

  std::vector<bitset::BitSet<dim> > d_data;
  size_t d_rows;
  size_t d_columns;

 public:

// constructors and destructors
  BitMatrix() 
    {}

  explicit BitMatrix(size_t n)
    :d_data(n),
     d_rows(n),
     d_columns(n) 
    {}

  BitMatrix(size_t m, size_t n)
    :d_data(n),
     d_rows(m),
     d_columns(n) 
    {}

  explicit BitMatrix(const std::vector<BitVector<dim> >&);

  ~BitMatrix() 
    {}

// copy and assignment
  BitMatrix(const BitMatrix& m)
    :d_data(m.d_data),
     d_rows(m.d_rows),
     d_columns(m.d_columns) 
    {}

  BitMatrix& operator=(const BitMatrix& m) {
    d_data = m.d_data;
    d_rows = m.d_rows;
    d_columns = m.d_columns;

    return *this;
  }


// accessors

  BitVector<dim>& apply(BitVector<dim>&, const BitVector<dim>&) const;

  template<typename I, typename O> void apply(const I&, const I&, O) const;

  const bitset::BitSet<dim>& column(size_t j) const {
    return d_data[j];
  }

  void column(BitVector<dim>&, size_t) const;

  void image(std::vector<BitVector<dim> >&) const;

  bool isEmpty() const {
    return d_data.size() == 0;
  }

  void kernel(std::vector<BitVector<dim> >&) const;

  size_t numColumns() const {
    return d_columns;
  }

  size_t numRows() const {
    return d_rows;
  }

  void row(BitVector<dim>&, size_t) const;

  bool test(size_t i, size_t j) const {
    return d_data[j].test(i);
  }

// manipulators

  BitMatrix& operator+= (const BitMatrix&);

  BitMatrix& operator*= (const BitMatrix&);

  void addColumn(const bitset::BitSet<dim>& f) {
    d_data.push_back(f); 
    d_columns++;
  }

  void addColumn(const BitVector<dim>& c) {
    addColumn(c.data());
  }

  void addToColumn(size_t j, const BitVector<dim>& v) {
    d_data[j] ^= v.data();
  }

  template<typename I> void addColumns(const I&, const I&);

  template<typename I> void addRows(const I&, const I&);

  void cutRows(size_t);

  BitMatrix& flip(size_t i, size_t j) {
    d_data[j].flip(i); 
    return *this;
  }

  BitMatrix& invert();

  void reset();

  BitMatrix& reset(size_t i, size_t j) {
    d_data[j].reset(i); 
    return *this;
  }

  void resize(size_t n) {
    resize(n,n);
  }

  void resize(size_t m, size_t n);

  BitMatrix& set(size_t i, size_t j) {
    d_data[j].set(i); 
    return *this;
  }

  void setColumn(size_t j, const bitset::BitSet<dim>& data) {
    d_data[j] = data;
  }

  void swap(BitMatrix& m);

  BitMatrix& transpose();
};

}

}

#include "bitvector_def.h"

#endif
