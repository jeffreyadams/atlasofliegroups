/*!
\file
\brief Class definitions and function declarations for BitVector.
*/
/*
  This is bitvector.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BITVECTOR_H  /* guard against multiple inclusions */
#define BITVECTOR_H

#include "bitvector_fwd.h"
#include "latticetypes_fwd.h" // not "latticetypes.h", which includes us!

#include <vector>
#include <cassert>

#include "bitset.h" // so users will see complete types; also needed for inline

/******** function declarations **********************************************/

namespace atlas {

namespace bitvector {

/*!
  \brief Puts in |v| the $Z/2Z$-linear combination of the |BitVector|s of |b|
  (each of size |n|) given by the bits of |e|.
*/
template<size_t dim>
  BitVector<dim> combination
  (const std::vector<BitVector<dim> >& b,
   size_t n,
   const bitset::BitSet<dim>& e);

// version with |BitSet|s instead of |BitVector|s; functional (size no issue)
template<size_t dim>
  bitset::BitSet<dim> combination(const std::vector<bitset::BitSet<dim> >&,
				  const bitset::BitSet<dim>&);


/*!
  \brief Put into |c| a solution of the system with as left hand sides the
  rows of a matrix whose columns are given by |b|, and as right hand sides the
  bits of |rhs|, and return |true|; if no solution exists just return |false|.
*/
template<size_t dim>
  bool firstSolution(bitset::BitSet<dim>& c,
		     const std::vector<BitVector<dim> >& b,
		     const BitVector<dim>& rhs);

/*!
  \brief Either find a solution of the system of equations |eqn|, putting it
  into |sol| and returning |true|, or return |false| if no solition exists.

  Here |eqns| holds a system of equations, the last bit of each being
  interpreted as the right hand side.

  If there is a solution, |sol| will be resized to to number of
  indeterminates, which is one less than the size of each equation; however,
  if the set of equations is empty, |sol| is left unchanged.
*/
template<size_t dimsol, size_t dimeq>
  bool firstSolution(BitVector<dimsol>& sol,
		     const std::vector<BitVector<dimeq> >& eqns);

template<size_t dim> void identityMatrix(BitMatrix<dim>&, size_t);

template<size_t dim> void initBasis(std::vector<BitVector<dim> >&, size_t);

template<size_t dim>
  void normalize(bitset::BitSet<dim>&, std::vector<BitVector<dim> >&);

template<size_t dim>
  void normalSpanAdd(std::vector<BitVector<dim> >&, std::vector<size_t>&,
		     const BitVector<dim>&);
/* unused functions
template<size_t dim>
  void complement(bitset::BitSet<dim>&, const std::vector<BitVector<dim> >&,
		  size_t);

template<size_t dim> bool isIndependent(const std::vector<BitVector<dim> >&);

template<size_t dim>
  void projection(BitMatrix<dim>& p, const std::vector<BitVector<dim> >& b,
		  size_t d);

template<size_t dim>
  void reflectionMatrix(BitMatrix<dim>&, const BitVector<dim>&,
			const BitVector<dim>&);

template<size_t dim>
  void relations(std::vector<BitVector<dim> >&,
		 const std::vector<BitVector<dim> >&);
*/

template<size_t dim>
  void spanAdd(std::vector<BitVector<dim> >&, std::vector<size_t>&,
	       const BitVector<dim>&);

} // namespace bitvector

/******** type definitions **************************************************/

namespace bitvector {
  /*!
  \brief The template class |BitVector<dim>| represents a number |size| with
  |0<=size<=dim|, and a vector in the (Z/2Z)-vector space (Z/2Z)^size.

  The software envisions dim between 0 and four times the machine word length
  (precisely, four times the constant |longBits|, which is the number of bits
  in an unsigned long integer). What is now fully implemented allows |dim| to
  be one or two times the word length (see the discussion in the description
  of the class BitSet<n>). It seems that only BitVector<RANK_MAX> and
  BitVector<2*RANK_MAX> are now instantiated, (that is <16> or <32>), so |dim|
  is not very relevant; one could imagine |dim==longBits| throughout.

  Let the integer |m| be \f$\lceil dim/longBits \rceil\f$ (quotient rounded
  upwards) . The vector is stored as the BitSet<dim> |d_data|, which is a
  (fixed size) array of |m| words of memory. (On a 32 bit machine with
  RANK_MAX=16 one always has |m==1|, so that |d_data| is a single word of
  memory.) We look only at the first |d_size| bits of |d_data|; but |d_size|
  can be changed by manipulators (like the member functions resize and
  pushBack). [Maybe |d_size| never exceeds |RANK_MAX+1|, even for variables
  declared as |BitVector<2*RANK_MAX>| MvL]

  Given the number of methods that are passed on to the |BitSet<dim>| field
  |d_data|, one might wonder if it would not have been better to publicly
  derive from |BitSet<dim>|, MvL.

  A |BitVector| should be thought of as a column vector. A |Bitmatrix| will in
  general act on it on the left, just like a |LatticeMatrix| on a |Weight|.
  This does not stop |BitVector|s from being used most for coweights modulo 2.
  */

template<size_t dim> class BitVector
{
  bitset::BitSet<dim> d_data;
  unsigned short int d_size; // this is more than large enough for all uses

 public:

  // constructors and destructors
  BitVector()
    : d_data(), d_size(0)
    {}

  explicit BitVector(size_t n) // initialised to 0 in $Z/2Z$
    : d_data(), d_size(n)
    {}

  BitVector(size_t n, size_t j) // canonical basis vector $e_j$ in $(Z/2Z)^n$
    : d_data(), d_size(n)
    {
      set(j);
    }

  BitVector(bitset::BitSet<dim> data, size_t n) // view |data| as |BitVector|
    : d_data(data), d_size(n)
    {}

  BitVector(const latticetypes::LatticeElt& weight); // reduce weight mod 2

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

  bool operator< (const BitVector& v) const
  {
    assert(d_size==v.d_size);
    return d_data<v.d_data;
  }

  bool operator== (const BitVector& v) const
  {
    assert(d_size==v.d_size);
    return d_data == v.d_data;
  }

  bool operator!= (const BitVector& v) const
  {
    assert(d_size==v.d_size);
    return d_data != v.d_data;
  }

  bool operator[] (size_t i) const
  {
    assert(i<d_size);
    return d_data[i];
  }

  size_t size() const { return d_size; }

  const bitset::BitSet<dim>& data() const { return d_data; }

  size_t firstBit() const { return d_data.firstBit(); }

  size_t count() { return d_data.count(); }

  bool isZero() const { return d_data.none(); }

  bool nonZero() const { return d_data.any(); }

  //! scalar product with value in $\Z/2$
  bool dot(const BitVector& v) const
  { return ((d_data & v.d_data).count()&1u)!=0; }

  BitVector operator+ (const BitVector& v) const
  { BitVector<dim> result(*this); result+=v; return result; }

// manipulators

  BitVector& operator+= (const BitVector& v)
  {
    assert(d_size==v.d_size);
    d_data ^= v.d_data;
    return *this;
  }

  BitVector& operator-= (const BitVector& v) // same thing as +=
  {
    assert(d_size==v.d_size);
    d_data ^= v.d_data;
    return *this;
  }

  BitVector& operator&= (const BitVector& v) // bitwise multiply
  {
    assert(d_size==v.d_size);
    d_data &= v.d_data;
    return *this;
  }

  BitVector& operator>>= (size_t pos) { d_data >>= pos; return *this; }
  BitVector& operator<<= (size_t pos) { d_data <<= pos; return *this; }

  BitVector& flip(size_t i)
  {
    assert(i<d_size);
    d_data.flip(i);
    return *this;
  }

  BitVector& pushBack(bool);

  BitVector& set(size_t i)
  {
    assert(i<d_size);
    d_data.set(i);
    return *this;
  }

  void set(size_t i, bool b) { assert(i<d_size); d_data.set(i,b); }

  void set_mod2(size_t i, unsigned long v) { set(i,(v&1)!=0); }

  BitVector& reset() { d_data.reset(); return *this; }

  BitVector& reset(size_t i)
  {
    assert(i<d_size);
    d_data.reset(i);
    return *this;
  }

  void resize(size_t n) { assert(n<=dim); d_size = n; }

  void slice(const bitset::BitSet<dim>& mask);
  void unslice(bitset::BitSet<dim> mask, size_t new_size);
}; // class BitVector

/* the following template inherits everything from |std::vector| but after
   some constructors that mimick those of |BitVector|, we also provide a
   constructor that converts from |WeightList|, reducing coefficients mod 2
 */
template<size_t dim> class BitVectorList
 : public std::vector<BitVector<dim> >
{
 public:
  // default constructor
  BitVectorList() : std::vector<BitVector<dim> >() {}

  // dimension-only constructor
  BitVectorList(size_t n) : std::vector<BitVector<dim> >(n) {}

  // fixed element constructor
  BitVectorList(size_t n,BitVector<dim> model)
    : std::vector<BitVector<dim> >(n,model)
    {}

  // also allow explicit consversion (implicit would be too dangerous)
  explicit BitVectorList(const std::vector<BitVector<dim> >& v)
    : std::vector<BitVector<dim> >(v) // copy
    {}

  /* reduction mod 2 is done via range-constructor of vector, which on its
     turn calls |BitVector<dim> (const LatticeElt&)| on the elements */
  BitVectorList(const latticetypes::WeightList& l)
    : std::vector<BitVector<dim> >(l.begin(),l.end())
    {}

  template<typename I> // input iterator
    BitVectorList(I begin, I end)
    : std::vector<BitVector<dim> >(begin,end)
    {}
};


// note that the elements in d_data are the _column_ vectors of the
// matrix

/*!
\brief A matrix of d_rows rows and d_columns columns, with entries in Z/2Z.

The number of rows d_rows should be less than or equal to dim, which
in turn is envisioned to be at most four times the machine word size.
The present implementation allows dim at most twice the machine word
size, and what is used is dim equal to RANK_MAX (now 16).

At least when d_columns is less than or equal to dim, the matrix can
act on the left on a BitVector of size d_columns; in this setting each
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
template<size_t dim> class BitMatrix
{
  /*!
  A vector of d_columns BitSet's (each of size d_rows), the columns of
  the BitMatrix. Thus d_data.size()==d_columns at all times.
  */
  std::vector<bitset::BitSet<dim> > d_data;

  /*!
  Number of rows of the BitMatrix. Cannot exceed template argument dim
  */
  unsigned short int d_rows;

  /*!
  Number of columns of the BitMatrix. This field is redundant, see d_data.
  */
  unsigned short int d_columns;

 public:

// constructors and destructors
  BitMatrix()
    : d_data(), d_rows(0), d_columns(0)
    {}

  explicit BitMatrix(size_t n) // square matrix
    : d_data(n) // make |n| default constructed |BitSet|s as columns
    , d_rows(n)
    , d_columns(n)
    {
      assert(n<=dim);
    }

  BitMatrix(size_t m, size_t n)
    : d_data(n) // make |n| default constructed |BitSet|s as columns
    , d_rows(m)
    , d_columns(n)
    {
      assert(m<=dim); // having |n>dim| is not immediately fatal
    }

  explicit BitMatrix(const std::vector<BitVector<dim> >&,  // set by columns
		     unsigned short int num_rows);  // all of size |num_rows|
  explicit BitMatrix(const latticetypes::LatticeMatrix& m); // set modulo 2

// copy and assignment (implicitly generated ones will do)

// accessors

  //! The (i,j) entry of the BitMatrix.
  bool test(size_t i, size_t j) const
  {
    assert(i<d_rows);
    assert(j<d_columns);
    return d_data[j].test(i);
  }

  BitVector<dim> operator*(const BitVector<dim>& src) const; // matrix * vector


  template<typename I, typename O> void apply(const I&, const I&, O) const;


  /*!
  Tests whether the BitMatrix is empty (has zero columns).
  */
  bool isEmpty() const {
    return d_data.size() == 0;
  }

  // independent vectors spanning the image of the matrix
  BitVectorList<dim> image() const;

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


  void get_row(BitVector<dim>&, size_t) const;
  BitVector<dim> row(size_t i) const
  { BitVector<dim> result; get_row(result,i); return result; }

  //! Column $j$ of the BitMatrix, as a BitVector.
  BitVector<dim> column(size_t j) const
  {
    assert(j<d_columns);
    return BitVector<dim>(d_data[j],d_rows);
  }
  void get_column(BitVector<dim>& c, size_t j) const { c=column(j); }


// manipulators

  BitMatrix& operator+= (const BitMatrix&);

  BitMatrix& operator*= (const BitMatrix&);

  /*!
 Adds the BitSet f as a new column at the end to the BitMatrix.
  */
  void addColumn(const bitset::BitSet<dim>& f) {
    d_data.push_back(f);
    d_columns++;
  }

  void addColumn(const BitVector<dim>& c) {
    assert(c.size()==d_rows);
    addColumn(c.data()); // call previous method, which does |d_columns++|
  }

  /*!
 Adds the BitVector v to column j of the BitMatrix.
  */
  void addToColumn(size_t j, const BitVector<dim>& v) {
    assert(v.size()==d_rows);
    d_data[j] ^= v.data();
  }

  BitMatrix& invert();  // currently unused, but defined anyway

  /*!
  Puts a 1 in row i and column j of the BitMatrix.
  */
  BitMatrix& set(size_t i, size_t j) {
    assert(i<d_rows);
    assert(j<d_columns);
    d_data[j].set(i);
    return *this;
  }

  BitMatrix& reset(size_t i, size_t j) {
    assert(i<d_rows);
    assert(j<d_columns);
    d_data[j].reset(i);
    return *this;
  }

  void set(size_t i, size_t j, bool b) {
    if (b) set(i,j); else reset(i,j);
  }

  void set_mod2(size_t i, size_t j, unsigned long v) {
    set(i,j, (v&1)!=0 );
  }

  void reset();


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
  Puts the BitSet data in column j of the BitMatrix.
  */
  void setColumn(size_t j, const bitset::BitSet<dim>& data)
  {
    assert(j<d_columns);
    d_data[j] = data;
  }

  void swap(BitMatrix& m);

  /*!
  Transposes the BitMatrix

  NOTE: it is the caller's responsibility to check that d_columns does
  not exceed dim.
  */
  BitMatrix& transpose();
}; // |class BitMatrix|

//			  Inlined function definition

inline
  bitvector::BinaryEquation make_equation(const SmallBitVector& lhs, bool rhs)
  {
    BinaryEquation eqn(lhs.data(),lhs.size());
    eqn.pushBack(rhs);
    return eqn;
  }


} // |namespace bitvector|

} // |namespace  atlas|


#endif
