/*
  This is bitvector.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef BITVECTOR_H  /* guard against multiple inclusions */
#define BITVECTOR_H

#include "atlas_types.h"

#include <vector>
#include <cassert>

#include "bitset.h" // so users will see complete types; also needed for inline
#include "matrix.h" // for |matrix::Vector|

/******** function declarations **********************************************/

namespace atlas {

namespace bitvector {

/* Put in |v| the $Z/2Z$-linear combination of the |BitVector|s of |b|
  (each of size |n|) given by the bits of |e|.
*/
template<size_t dim>
  BitVector<dim> combination
  (const std::vector<BitVector<dim> >& b,
   size_t n,
   const BitSet<dim>& e);

// version with |BitSet|s instead of |BitVector|s; functional (size no issue)
template<size_t dim>
  BitSet<dim> combination(const std::vector<BitSet<dim> >&,
				  const BitSet<dim>&);


/*
  Find out whether any combination of the vectos in |b| adds to |rhs|, and if
  so flag such a combination in the bits of |c|. Nothing changes when |false|.
*/
template<size_t dim>
  bool combination_exists(const std::vector<BitVector<dim> >& b,
			  const BitVector<dim>& rhs,
			  BitSet<dim>& c);

/*
  Either find a solution of the system of equations |eqn|, putting it
  into |sol| and returning |true|, or return |false| if no solition exists.

  Here |eqns| holds a system of equations, the last bit of each being
  interpreted as the right hand side.

  If there is a solution, |sol| will be resized to to number of
  indeterminates, which is one less than the size of each equation; however,
  if the set of equations is empty, |sol| is left unchanged.
*/
template<size_t dimsol, size_t dimeq>
  bool solvable(const std::vector<BitVector<dimeq> >& eqns,
		BitVector<dimsol>& sol);

template<size_t dim> void identityMatrix(BitMatrix<dim>&, size_t);

template<size_t dim> void initBasis(std::vector<BitVector<dim> >&, size_t);

template<size_t dim>
  void Gauss_Jordan(BitSet<dim>&, std::vector<BitVector<dim> >&);

template<size_t dim>
  void normalSpanAdd(std::vector<BitVector<dim> >&, std::vector<size_t>&,
		     const BitVector<dim>&);
/* unused functions
template<size_t dim>
  void complement(BitSet<dim>&, const std::vector<BitVector<dim> >&,
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


/******** type definitions **************************************************/

/*
  The template class |BitVector<dim>| represents a number |size| with
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
 public:
  typedef BitSet<dim> base_set;
 private:
  base_set d_data;
  unsigned short int d_size; // this is more than large enough for all uses

 public:

  // constructors and destructors
  explicit BitVector(size_t n) // initialised to 0 in $Z/2Z$
    : d_data(), d_size(n)
    {}

  BitVector(size_t n, size_t j) // canonical basis vector $e_j$ in $(Z/2Z)^n$
    : d_data(), d_size(n)
    {
      set(j);
    }

  BitVector(BitSet<dim> data, size_t n) // view |data| as |BitVector|
    : d_data(data), d_size(n)
    {}

  template<typename C>
    explicit BitVector(const matrix::Vector<C>& weight); // reduce weight mod 2

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

  const BitSet<dim>& data() const { return d_data; }

  size_t firstBit() const { return d_data.firstBit(); }

  size_t count() { return d_data.count(); }

  bool isZero() const { return d_data.none(); }

  bool nonZero() const { return d_data.any(); }

  // scalar product with value in $\Z/2$
  bool dot(const BitVector& v) const
  { return ((d_data & v.d_data).count()&1u)!=0; }

  BitVector operator+ (const BitVector& v) const
  { BitVector<dim> result(*this); result+=v; return result; }

  // the same operation as the previous method, under different name
  // the destinction may sometimes be useful for documentation purposes
  BitVector operator- (const BitVector& v) const
  { BitVector<dim> result(*this); result-=v; return result; }

  template<typename C>
    explicit operator matrix::Vector<C>() const
  { matrix::Vector<C> result(size());
    for (unsigned int i=size(); i-->0;)
      result[i]=d_data[i];
    return result;
  }

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

  void slice(const BitSet<dim>& mask);
  void unslice(BitSet<dim> mask, size_t new_size);
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
  BitVectorList(const std::vector<matrix::Vector<int> >& l)
    : std::vector<BitVector<dim> >(l.begin(),l.end())
    {}

  template<typename I> // input iterator
    BitVectorList(I begin, I end)
    : std::vector<BitVector<dim> >(begin,end)
    {}
};


// note that the elements in d_data are the _column_ vectors of the matrix

/*
  A rectangular matrix with entries in Z/2Z.

  The number of rows |d_rows| should be less than or equal to the template
  parameter |dim|, which in turn is envisioned to be at most four times the
  machine word size. The present |BitSet| implementation allows |dim| at most
  twice the machine word size, and what is used is |dim| equal to |RANK_MAX|.

  At least when |d_columns<=dim|, the matrix can act on the left on a
  |BitVector| of size |d_columns|; in this setting each column of the matrix
  is the image of one of the standard basis vectors in the domain.

  What is stored in |d_data|, as |BitSet|'s, are the column vectors of the
  matrix. Construction of a matrix is therefore most efficient when columns
  are added to it.

  Notice that the columns are stored as |BitSets| and not |BitVectors|. A
  |BitVector| is a larger data structure than a |BitSet|, including also an
  integer |d_size| saying how many of the available bits are significant. In a
  |BitMatrix| this integer must be the same for all of the columns, so it is
  easier and safer to store it once for the whole |BitMatrix|, and also to
  modify it just once when the matrix is resized.
*/
template<size_t dim> class BitMatrix
{
/*
  A vector of |d_columns| |BitSet| values (each to be thought of as a vector
  of size |d_rows|), the columns of the |BitMatrix|. Thus
  |d_data.size()==d_columns| at all times.
*/
  std::vector<BitSet<dim> > d_data;

  // Number of rows of the BitMatrix. Cannot exceed template argument |dim|
  unsigned short int d_rows;

  // Number of columns of the BitMatrix. This field is redundant, see |d_data|.
  unsigned short int d_columns;

 public:

// constructors and destructors
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
  explicit BitMatrix(const matrix::Matrix<int>& m); // set modulo 2

  static BitMatrix identity(unsigned int n)
  { BitMatrix M(n);
    for (size_t i=0; i<n; ++i)
      M.set(i,i);
    return M;
  }

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

  BitVectorList<dim> image() const; // free generators of image of matrix
  BitVectorList<dim> kernel() const; // free generators of kernel of matrix

  size_t numColumns() const { return d_columns; } // the number of columns
  size_t numRows() const { return d_rows; } // the number of rows

  BitVector<dim> row(size_t i) const;

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

  // Add |f| as a new column at the end to the |BitMatrix|.
  void addColumn(const BitSet<dim>& f) {
    d_data.push_back(f);
    ++d_columns;
  }

  void addColumn(const BitVector<dim>& c) {
    assert(c.size()==d_rows);
    addColumn(c.data()); // call previous method, which does |++d_columns|
  }

  // Add the BitVector |v| to column |j| of the BitMatrix.
  void addToColumn(size_t j, const BitVector<dim>& v) {
    assert(v.size()==d_rows);
    d_data[j] ^= v.data();
  }

  // matrix $B$ such that $ABA=A$, will be inverse if $A$ invertible
  BitMatrix section() const;

  // Put a 1 in row |i| and column |j| of the |BitMatrix|.
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


  /* Resize the BitMatrix to an |n| by |n| square.

  NOTE: it is the caller's responsibility to check that |n<=dim|.
  */
  void resize(size_t n) { resize(n,n); }

  void resize(size_t m, size_t n);

  // Puts |data| in column j of the BitMatrix.
  void setColumn(size_t j, const BitSet<dim>& data)
  {
    assert(j<d_columns);
    d_data[j] = data;
  }

  void swap(BitMatrix& m);

/*
  Transpose the BitMatrix

  NOTE: it is the caller's responsibility to check that |d_columns<=dim|.
*/
  BitMatrix& transpose();
}; // |class BitMatrix|

//			  Inlined function definition

inline
  BinaryEquation make_equation(const SmallBitVector& lhs, bool rhs)
  {
    BinaryEquation eqn(lhs.data(),lhs.size());
    eqn.pushBack(rhs);
    return eqn;
  }


} // |namespace bitvector|

} // |namespace  atlas|


#endif
