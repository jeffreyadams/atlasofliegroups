/*!
\file
\brief Class definitions and function declarations for WeylGroup.
*/
/*
  This is weyl.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef WEYL_H  /* guard against multiple inclusions */
#define WEYL_H

#include "weyl_fwd.h"

#include "constants.h"
#include "error.h"
#include "latticetypes.h"
#include "size.h"
#include "tags.h"

/******** type declarations *************************************************/

namespace atlas {

namespace weyl {

  typedef Generator WeylInterface[constants::RANK_MAX];
  typedef std::pair<EltPiece,Generator> ShiftValue;

  class RowBase;
  typedef RowBase ShiftRow;
  typedef RowBase OutRow;

  class Transducer;

}

/******** constant declarations *********************************************/

namespace weyl {

  const unsigned char UndefValue = constants::ucharMax;
  const EltPiece UndefEltPiece = UndefValue;
  const Generator UndefGenerator = UndefValue;
  const unsigned long UndefOrder = 0; // is never a valid value

  const Twist& identityTwist();

}

/******** function declarations *********************************************/

namespace weyl {

  void copy(Twist&, const Twist&);

  void makeReflections(WeylEltList&, const WeylGroup&);

}

/******** type definitions **************************************************/

namespace weyl {

  /*!
\brief Element of a Weyl group.

The representation is described in detail in the description of the
class WeylGroup.  An array of RANK_MAX unsigned char, the ith
representing a shortest length coset representative of a parabolic
subquotient W_{i-1}\\W_i.
  */
class WeylElt {

 private:

  /*!

\brief Represents factorization of Weyl group element w as a product
of shortest length coset representatives for parabolic subquotients
W_{i-1}\\W_i. 

Entry \#i is an unsigned char parametrizing the ith coset
representative w_i for an element of W_{i-1}\\W_i.  Then w =
w_1.w_2...w_n.  
  */
  EltPiece d_data[constants::RANK_MAX];

 public:

// constructors and destructors

/*!
\brief Constructs the identity element of W.
*/
  WeylElt() {
    memset(d_data,0,constants::RANK_MAX);
  }

/*!
\brief Constructs the element of W with all coset representatives
equal to x.

[Hard to see how this could possibly be used except with x=0.  DV 7/24/06.]
*/
  explicit WeylElt(EltPiece x) {
    memset(d_data,x,constants::RANK_MAX);
  }

  ~WeylElt() 
    {}

// copy and assignment
  WeylElt(const WeylElt& w) {
    memcpy(d_data,w.d_data,constants::RANK_MAX);
  }

  WeylElt& operator=(const WeylElt& w) {
    memcpy(d_data,w.d_data,constants::RANK_MAX); 
    return *this;
  }

// accessors


/*!
\brief Returns the jth factor of the Weyl group element.
*/
  EltPiece operator[] (size_t j) const {
    return d_data[j];
  }

  /*!
\brief Tests whether this Weyl group element is lexicographically
strictly less than the Weyl group element following the < sign.
  */
  bool operator< (const WeylElt& w) const {
    return memcmp(d_data,w.d_data,constants::RANK_MAX) < 0;
  }

  /*!
\brief Tests whether this Weyl group element is equal to the Weyl
group element following the == sign.
  */
  bool operator== (const WeylElt& w) const {
    return !memcmp(d_data,w.d_data,constants::RANK_MAX);
  }

  /*!
\brief Tests whether this Weyl group element is not equal to the Weyl
group element following the != sign.
  */
  bool operator!= (const WeylElt& w) const {
    return memcmp(d_data,w.d_data,constants::RANK_MAX);
  }

// manipulators
  EltPiece& operator[] (size_t j) {
    return d_data[j];
  }

  /*!
\brief Sets the Weyl group element equal to 1.
  */
  void clear() {
    memset(d_data,0,constants::RANK_MAX);
  }
};

const WeylElt Identity;

/*!
\brief Represents one row of a transducer table for a Weyl group.
*/
class RowBase {

 private:

  unsigned char d_data[constants::RANK_MAX];

 public:

// constructors and destructors
  RowBase() {
    memset(d_data,UndefValue,constants::RANK_MAX);
  }

  ~RowBase() {}

// copy and assignment
  RowBase(const RowBase& r) {
    memcpy(d_data,r.d_data,constants::RANK_MAX);
  }

  RowBase& operator=(const RowBase& r) {
    memcpy(d_data,r.d_data,constants::RANK_MAX); return *this;
  }

// accessors
  Generator operator[] (size_t j) const {
    return d_data[j];
  }

// manipulators
  Generator& operator[] (size_t j) {
    return d_data[j];
  }
};

/*!
\brief Right multiplication action of simple reflections on a Weyl
group modulo a maximal parabolic subgroup.

In the notation from the description of the class WeylGroup, there
will be one transducer class for each parabolic subquotient
W_{r-1}\\W_r. List the shortest length coset representatives for this
subquotient as x_0,...,x_{N_r-1}.  Recall that the simple roots were
ordered to guarantee that N_r-1 fits in an unsigned char, so each
coset representative can be indexed by an unsigned char.  We wish to
compute the product x_i.s_j for j between 1 and r.  The key
theoretical fact about multiplication is that there are two mutually
exclusive possibilities:

x_i.s_j = x_{i'}  (some i' ne i) 

OR

x_i.s_j = s_k.x_i  (some k < r).

The first possibility is called _transduction_ and the second
_transition_.  Entry \#j of row i of the transduction table is equal
to i' (in the first case) or to i (in the second).  This is stored in
d_shift, and accessed by the member function shift.  The reflection
s_k (which is defined only in the transition case -- that is, when
shift(i,j) = i) is stored in d_out, and accessed by the member
function out.
*/
class Transducer {

 private:

  /*!
\brief Column \#j is the permutation action of s_j on W_{r-1}\\W_r.

Entry i lists the (in its jth entry, for j less than the rank
r) the unsigned char i' so that x_i.s_j has minimal coset
representative x_i'.
  */
  std::vector<ShiftRow> d_shift;

  /*!
\brief Column \#j gives the simple reflections s_k (if any) with x.s_j
= s_k.x, and k<r.

Entry i lists (in its jth entry, for j less than the rank
r) the generator k so that x_i.s_j = s_k.x_i.   
  */
  std::vector<OutRow> d_out;

  /*!
\brief Lengths of the minimal coset representatives.
  */
  std::vector<unsigned long> d_length;

  /*!
\brief Reduced expressions of the minimal coset representatives.
  */
  std::vector<WeylWord> d_piece;

 public:

// constructors and destructors
  Transducer() 
    {}

  Transducer(const latticetypes::LatticeMatrix&, size_t);

  ~Transducer() 
    {}

// accessors


  /*!
\brief Length of minimal coset representative x. 
  */
  unsigned long length(EltPiece x) const {
    return d_length[x];
  }


  /*!
\brief Maximal length of minimal coset representatives.

This is the number of positive roots for the Levi subgroup L_r, minus
the number of positive roots for L_{r-1}. 
  */
  unsigned long maxlength() const {
    return d_length.back();
  }


  /*!
\brief Simple reflection t (strictly preceding s) so that xs = tx.

[I don't know what it returns in the case of transduction, when no
such t exists. DV 7/24/06.]
  */
  Generator out(EltPiece x, Generator s) const {
    return d_out[x][s];
  }

  /*!
\brief Right coset x' defined by x' = xs.

When x' is not equal to s, this is an equality of minimal coset
representatives. When x'=x, the equation for minimal coset representatives
is out(x,s).x = x.s. 
  */
  EltPiece shift(EltPiece x, Generator s) const {
    return d_shift[x][s];
  }

  /*!
\brief Number of cosets W_{r-1}\\W_r. 
  */
  unsigned long size() const {
    return d_shift.size();
  }


  /*!
\brief Reduced decomposition in W (or W_r) of minimal coset representative x.
  */
  const WeylWord& wordPiece(EltPiece x) const {
    return d_piece[x];
  }

// this class should have no manipulators!
};


  /*!
\brief  Represents a  Weyl group.

  The WeylGroup class is a variation on my favourite implementation in terms
  of transducers.

  I have tried to make a careful choice of datatype for the group
  elements in order to strike the right balance between efficiency and
  generality. [Since this is Fokko's mathematics as well as his code,
  I have tried to keep his voice. Occasional amplifications are in
  brackets.  DV 7/23/06.] This has led me to the choice of _fixed
  size_ arrays of unsigned chars, representing the "parabolic
  subquotient" representation of a given element; in other words, in a
  group of rank n, the first n elements of the array are used (but
  since it is fixed size, we have to allocate it to some constant
  value, namely RANK_MAX.)

  [Meaning of "parabolic subquotient representation: one orders the
  simple generators in an appropriate fashion s_1,...,s_n (henceforth
  called the "internal representation). Define W_i = <s_1,...,s_i>,
  so that W_0 = {1} subset W_1 subset W_2 subset ... subset W_n =
  W. Each W_i is a parabolic subgroup of W, so its length function is
  the restriction of that of W. Each right coset in W_{i-1}\\W_i has
  a unique minimal length coset representative. Write N_i for the
  number of cosets (or representatives). Then the cardinality of W is
  N_1.N_2...N_n. More precisely, and element of W has a unique
  factorization as x_1.x_2...x_n, with x_i one of the N_i minimal
  length coset representatives of W_{i-1}\\W_i.

  The key to the implementation is to order the generators in such a
  way that every parabolic subquotient W_{i-1}\\W_i has cardinality
  fitting into an unsigned char: that is, cardinality at most 256.
  (This is possible for any Weyl group of rank at most 128.) A
  WeylElt is an array of unsigned char of size RANK_MAX; the unsigned
  char in entry i is the number of the factor x_i.

  There is a little more about how the group structure is computed in
  this representation in the description of the class Transducer.  The
  mathematical reference is Fokko du Cloux, "A transducer approach to
  Coxeter groups," J. Symbolic Computation 27 (1999), 311-324.]

  This is not so wasteful as it may sound : of course the
  representation as a single number is more compact, but will overflow
  even on 64-bit machines for some groups of rank <= 16, and much
  sooner of course on 32-bit machines; also it imposes some
  computational overhead for packing and unpacking. 
  
  [Meaning of "representation as a single number:" if the cardinality
  of W fits in an unsigned long, then an element of W may be
  represented as the unsigned long

  x_1 + x_2(N_1) + x_3(N_1N_2) + ... +x_n(N_1...N_n).
  
  In this formula I have changed the meaning of x_i: now it means the
  integer between 0 and N_i enumerating the coset representative
  previously called x_i. This representation is called the "packed
  form" in the software.  On a 32 bit machine, this representation
  works for W(E8), but not for W(A12) (which has order 13!, which is
  about 1.5 x 2^{32}).]

  Any variable-size structure like the STL vector uses already three
  unsigned longs for its control structure (the address of the data,
  the size and the capacity), and then it still has to allocate. This
  could perhaps be simplified to just a pointer (after all the size of
  the allocation is known to the group) but you still have the big
  overhead of allocating and deallocating memory from the heap, and
  remembering to delete the pointers when they go out of scope, or
  else use autopointers ...

  And if worst comes to worst and one really is worried about a factor
  2 waste for type E8 (which might be significant if one deals with
  huge containers of group elements), then one can still recompile
  with RANK_MAX=8, which will then give a datatype that in 64 bits is
  the same size as an unsigned long.

  Notice that the unsigned char type miraculously suffices for all
  subquotients of all groups up to rank 128 (indeed, the biggest
  subquotient for B128 is of order 256), _provided_ the generators are
  enumerated in an appropriate order. This forces us to do quite a bit
  of type recognition, which is relegated to the dynkin
  namespace. Because of this reordering, the group carries a little
  interface that will translate back and forth from the external
  ordering and the internal one.
*/
class WeylGroup {

 private:

  size_t d_rank;
  size::Size d_order;
  unsigned long d_maxlength;
  WeylElt d_longest;
  latticetypes::LatticeMatrix d_coxeterMatrix;
  std::vector<Transducer> d_transducer;
  Twist d_twist;
  WeylInterface d_in;
  WeylInterface d_out;

// private member functions
  void leftProdIn(WeylElt&, Generator) const;
  int prodIn(WeylElt&, Generator) const;
  void prodIn(WeylElt&, const WeylWord&) const;

  const WeylWord& wordPiece(const WeylElt& w, size_t j) const {
    const Transducer& tr = d_transducer[j];
    return tr.wordPiece(w[j]);
  }

 public:

// constructors and destructors
  WeylGroup()
    :d_rank(0UL),d_order(1UL),d_maxlength(0UL) {}

  WeylGroup(const latticetypes::LatticeMatrix&, const Twist* = 0);

  WeylGroup(const WeylGroup&, tags::DualTag);

  ~WeylGroup() {}

// accessors
  void conjugacyClass(WeylEltList&, const WeylElt&, bool twisted = true)
    const;

  /*!
\brief Conjugates w by the generator s.
  */
  void conjugate(WeylElt& w, Generator s) const {
    leftProd(w,s);
    prod(w,s);
  }

  bool hasDescent(Generator, const WeylElt&) const;

  bool hasTwistedCommutation(Generator, const WeylElt&) const;

  void invert(WeylElt&) const;

  unsigned long involutionLength(const WeylElt&) const;

  void involutionOut(WeylWord&, const WeylElt&) const;

  void invProd(WeylElt&, const WeylWord&) const;

  Generator leftDescent(const WeylElt&) const;

  void leftProd(WeylElt& w, Generator s) const {
    leftProdIn(w,d_in[s]);
  }

  unsigned long length(const WeylElt&) const;

  const WeylElt& longest() const {
    return d_longest;
  }

  unsigned long maxlength() const {
    return d_maxlength;
  }

  const size::Size& order() const {
    return d_order;
  }

  void out(WeylWord&, const WeylElt&) const;

  void outerTwist(Twist&) const;

  int prod(WeylElt& w, Generator s) const {
    return prodIn(w,d_in[s]);
  }

  void prod(WeylElt&, const WeylElt&) const;

  void prod(WeylElt&, const WeylWord&) const;

  size_t rank() const {
    return d_rank;
  }

  unsigned long toUlong(const WeylElt&) const;

  WeylElt toWeylElt(unsigned long) const;

  void translate(WeylElt&, const WeylInterface&) const;

  void twist(WeylElt&) const;

  Generator twistGenerator(Generator s) const {
    return d_twist[s];
  }

  void twistedConjugate(WeylElt& w, Generator s) const {
    leftProd(w,s);
    prodIn(w,d_twist[d_in[s]]);
  }

// manipulators
  void swap(WeylGroup&);
};

}

}

#endif
