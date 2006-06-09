/*!
\file
\brief Class definitions and function declarations for BitVector.
*/
/*
  This is bitvector.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

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
  /*!  
  \brief A vector d_data in (Z/2Z)^d_size, with d_size a non-negative
  integer less than or equal to dim.

  The class BitVector<dim> represents a vector in the first d_size
  coordinates of the vector space (Z/2Z)^dim over the two-element
  field.  The software envisions dim between 0 and four times the
  machine word length (precisely, four times the constant longBits,
  which is the number of bits in an unsigned long integer).  What is
  now fully implemented allows dim to be one or two times the word
  length (see the discussion in the description of the class
  BitSet<n>).  What is now instantiated seems to be BitVectors of dim
  equal to RANK_MAX or 2*RANK_MAX, (that is 16 or 32).  

  Let m be the smallest integer so that dim is at most m*longBits.  The
  vector is stored as the BitSet d_data, which is a (fixed size) array
  of m words of memory.  (On a 32 bit machine with RANK_MAX=16, this
  means that d_data is a single word of memory.)  We look only at the
  first d_size bits of d_data; but d_size can be changed by
  manipulators (like the member functions resize and pushBack).

  A BitVector should be thought of as a column vector.  A Bitmatrix will
  in general act on it on the left.  [I thought that Fokko always
  prefers matrices acting on the right on row vectors; but if I am
  reading this correctly, that preference didn't make it into this
  software.] 
  */
   
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

/*!
\brief A matrix of d_rows rows and d_columns columns, with entries in Z/2Z.

The number of rows d_rows should be less than or equal to dim, which
in turn is envisioned to be at most four times the machine word size.
The present implementation allows dim at most twice the machine word
size, and what is used is dim equal to RANK_MAX (now 16).  

The at least when d_rows is less than or equal to dim, the matrix can
act on the left on a BitVector of size d_column; in this setting each
column of the matrix is the image of one of the standard basis vectors
in the domain.

What is stored in d_data, as BitSet's, are the d_column column vectors
of the matrix.  Construction of a matrix is therefore most efficient
when columns are added to it.

Notice that the columns are BitSets and not BitVectors.  A BitVector
is a larger data structure than a BitSet, including also an integer
d_size saying how many of the available bits are significant.  In a
BitMatrix this integer must be the same for all of the columns, so it
is easier and safer to store it once for the whole BitMatrix, and also
to modify it just once when the matrix is resized.
*/
template<size_t dim> class BitMatrix {

 private:

  /*!
  A vector of d_rows BitSet's (each of size d_columns), the columns of
  the BitMatrix.
  */
  std::vector<bitset::BitSet<dim> > d_data; 

  /*!
  Number of rows of the BitMatrix.
  */
  size_t d_rows;

  /*!
  Number of columns of the BitMatrix.
  */
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

  /*!
  Returns column j of the BitMatrix, as a BitSet. 
  */
  const bitset::BitSet<dim>& column(size_t j) const {
    return d_data[j];
  }

  void column(BitVector<dim>&, size_t) const;

  void image(std::vector<BitVector<dim> >&) const;

  /*!
  Tests whether the BitMatrix is empty (has zero columns).
  */
  bool isEmpty() const {
    return d_data.size() == 0;
  }

  void kernel(std::vector<BitVector<dim> >&) const;

  /*!
  Returns the number of columns of the BitMatrix.
  */
  size_t numColumns() const {
    return d_columns;
  }

  /*!
  Returns the number of rows of the BitMatrix.
  */
  size_t numRows() const {
    return d_rows;
  }

  void row(BitVector<dim>&, size_t) const;

  /*!
  Tests the (i,j) entry of the BitMatrix.
  */
  bool test(size_t i, size_t j) const {
    return d_data[j].test(i);
  }

// manipulators

  BitMatrix& operator+= (const BitMatrix&);

  BitMatrix& operator*= (const BitMatrix&);

  /*!
 Adds the BitSet f as a new column (the first one, pushing the others
 back) to the BitMatrix.
  */
  void addColumn(const bitset::BitSet<dim>& f) {
    d_data.push_back(f); 
    d_columns++;
  }

  void addColumn(const BitVector<dim>& c) {
    addColumn(c.data());
  }

  /*!
 Adds the BitVector v to column j of the BitMatrix.
  */
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

  /*!
  Resizes the BitMatrix to an n by n square.

  NOTE: it is the caller's responsibility to check that n does not
  exceed dim.
  */
  void resize(size_t n) {
    resize(n,n);
  }

  void resize(size_t m, size_t n);

  /*!
  Puts a 1 in row i and column j of the BitMatrix.
  */
  BitMatrix& set(size_t i, size_t j) {
    d_data[j].set(i); 
    return *this;
  }

  /*!
  Puts the BitSet data in column j of the BitMatrix.
  */
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
