/*!
\file
\brief Class definitions and function declarations for WeylGroup.
*/
/*
  This is weyl.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef WEYL_H  /* guard against multiple inclusions */
#define WEYL_H

#include <cstring>

#include "weyl_fwd.h"
#include "rootdata_fwd.h"

#include "constants.h"
#include "error.h"
#include "latticetypes.h"
#include "size.h"
#include "tags.h"

/******** type declarations *************************************************/

namespace atlas {

namespace weyl {

  class RowBase; // information at one coset element in one transducer table
  typedef RowBase ShiftRow; // the case of a transition entry
  typedef RowBase OutRow;   // the case of a transduction entry

  class Transducer;

}

/******** constant declarations *********************************************/

namespace weyl {

  const unsigned char UndefValue = constants::ucharMax;
  const unsigned long UndefOrder = 0; // is never a valid value

}

/******** type definitions **************************************************/

namespace weyl {

  /*!
\brief A mapping between one interpretation of Generators and another
  */
  class WeylInterface
  {
    Generator d[constants::RANK_MAX];
  public:
    Generator& operator[] (size_t i) { return d[i]; }
    const Generator& operator[] (size_t i) const { return d[i]; }
  };

  /*!
\brief A permutation of the set of Generators, giving a diagram automorphism
  */
  typedef WeylInterface Twist; // use the same implementation


  /*!
\brief Element of a Weyl group.

The representation is described in detail in the description of the
class WeylGroup.  An array of RANK_MAX unsigned char, the ith
representing a shortest length coset representative of a parabolic
subquotient W_{i-1}\\W_i.
  */
class WeylElt {

  friend class WeylGroup; // so we can shield implementation from all others

public:
  /*!
\brief Represents a minimal length coset representative for one of the
  parabolic subquotients W_{i-1}\\W_i.
  */
  typedef unsigned char EltPiece;

 private:

  /*!

\brief Represents factorization of Weyl group element w as a product
of shortest length coset representatives for parabolic subquotients
W_{i-1}\\W_i.

Entry \#i-1 is an unsigned char parametrizing the ith coset
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
    std::memset(d_data,0,sizeof(d_data));
  }

  /*! \brief interpret |w| in weyl group |W| */
  WeylElt(const WeylWord& w, const WeylGroup& W);

// copy and assignment
// use standard definitions (raw copy of array)

// accessors


  /*!
\brief Tests whether this Weyl group element is lexicographically
strictly less than the Weyl group element following the < sign.
  */
  bool operator< (const WeylElt& w) const {
    return std::memcmp(d_data,w.d_data,sizeof(d_data)) < 0;
  }

  /*!
\brief Tests whether this Weyl group element is equal to the Weyl
group element following the == sign.
  */
  bool operator== (const WeylElt& w) const {
    return std::memcmp(d_data,w.d_data,sizeof(d_data))==0;
  }

  /*!
\brief Tests whether this Weyl group element is not equal to the Weyl
group element following the != sign.
  */
  bool operator!= (const WeylElt& w) const {
    return std::memcmp(d_data,w.d_data,sizeof(d_data))!=0;
  }

protected: // these are for |WeylGroup| and |TI_Entry|'s eyes only

/*!
\brief Returns the jth factor of the Weyl group element.
*/
  EltPiece operator[] (size_t j) const {
    return d_data[j];
  }

// manipulators
  EltPiece& operator[] (size_t j) {
    return d_data[j];
  }

public:
// dummy methods that mark transition of interpretation

  // from WeylElt interpretation to TwistedInvolution
  // nothing: cast to TwistedInvolution is allowed, would do construction if
  // TwistedInvolution were a derived class, and is a no-op in reality

  // from TwistedInvolution interpretation to WeylElt
  // calls of these mark referral to the "base" class of |TwistedInvolution|
  const WeylElt& w() const { return *this; }
  WeylElt& contents() { return *this; }

}; // class WeylElt

const WeylElt Identity; // default constructor initialises to identity

/*!
\brief Represents one row of a transducer table for a Weyl group.
*/
class RowBase {

 private:

  unsigned char d_data[constants::RANK_MAX];

 public:

// constructors and destructors
  RowBase() {
    std::memset(d_data,UndefValue,sizeof(d_data));
  }

  ~RowBase() {}

// copy and assignment
  RowBase(const RowBase& r) {
    std::memcpy(d_data,r.d_data,sizeof(d_data));
  }

  RowBase& operator=(const RowBase& r) {
    std::memcpy(d_data,r.d_data,sizeof(d_data)); return *this;
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
will be one Transducer object for each parabolic subquotient
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

The first possibility is called _transition_ and the second _transduction_.
(Confusingly Fokko's 1999 paper interchanges these terms at their definition,
but their usual meaning and the sequel makes clear that this was an error).

The Transducer has tables to describe the two cases. the first table |d_shift|
describes the transistions, namely |d_shift[i][j]==i'| in the first case; the
cases that are tranductions can be distinguished from these by the fact that
|d_shift[i][j]==i|. In these cases, the value |k| emitted by the transduction
is stored in |d_out[i][j]|, which otherwise contains the value |UndefGenerator|
*/

class Transducer {
  // there is one such object for each $r\in\{1,2,\ldots,n\}$
  // but $r$ is not explicitly stored in the Tranducer object

 private:

  /*!
\brief Right multiplication by $s_j$ gives transition |i -> d_shift[i][j]|
  */
  std::vector<ShiftRow> d_shift;

  /*!
\brief If |d_shift[i][j]==i| then $s_j$ transduces in state $i$ to $s_k$
 with $k=d_out[i][j]$ (otherwise |d_out[i][j]==UndefGenerator|).

  In this case $x_i.s_j = s_k.x_i$, so the state $i$ remains unchanged.
  */
  std::vector<OutRow> d_out;

  /*!
\brief Lengths of the minimal coset representatives $x_i$.
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
  unsigned long length(WeylElt::EltPiece x) const { return d_length[x]; }


  /*!
\brief Maximal length of minimal coset representatives.

This is the number of positive roots for the Levi subgroup L_r, minus
the number of positive roots for L_{r-1}.
  */
  unsigned long maxlength() const {
    return d_length.back();
  }


  /*!
\brief Simple reflection t (strictly preceding s) so that xs = tx, if any

In case of a transition, this returns UndefGenerator.
  */
  Generator out(WeylElt::EltPiece x, Generator s) const { return d_out[x][s]; }

  /*!
\brief Right coset x' defined by x' = xs.

When x' is not equal to s, this is an equality of minimal coset
representatives. When x'=x, the equation for minimal coset representatives
is out(x,s).x = x.s.
  */
  WeylElt::EltPiece shift(WeylElt::EltPiece x, Generator s) const
   { return d_shift[x][s]; }

  /*!
\brief Number of cosets W_{r-1}\\W_r.
  */
  unsigned long size() const {
    return d_shift.size();
  }


  /*!
\brief Reduced decomposition in W (or W_r) of minimal coset representative x.
  */
  const WeylWord& wordPiece(WeylElt::EltPiece x) const { return d_piece[x]; }

// this class should have no manipulators!
}; // class Transducer


  /*!
\brief Represents a Weyl group for the purpose of manipulating its elements

  The WeylGroup class is a variation on Fokko's own and favourite
  implementation in terms of transducers.

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
  subquotients of all groups up to rank 128 (indeed, the biggest subquotient
  for B128 is of order 256), _provided_ the generators are enumerated in an
  appropriate order. This forces us to do quite a bit of type recognition,
  which is relegated to the dynkin namespace. Because of this reordering, the
  group carries a little interface that will translate back and forth from the
  external ordering and the internal one. To speed up some computations we
  precompute in |d_min_star| for each generator the smallest other generator
  with which it does not commute (or the generator itself if there are none).
*/
class WeylGroup {

 private:

  size_t d_rank;
  size::Size d_order;
  unsigned long d_maxlength;
  WeylElt d_longest;
  latticetypes::LatticeMatrix d_coxeterMatrix;
  std::vector<Transducer> d_transducer;
  WeylInterface d_in;
  WeylInterface d_out;
  std::vector<Generator> d_min_star;

// private member functions
// these interpret Generators and WeylWords internally, so not for public use!

  // set |w=ws|, return value is $l(ws)-l(w)\in\{+1,-1}$
  int multIn(WeylElt& w, Generator s) const;

  // set |w=wv|, return value is $l(wv)-l(w)-l(v)\in -2\N$
  int multIn(WeylElt& w, const WeylWord& v) const;

  // set |w=sw|, return value is $l(sw)-l(w)\in\{+1,-1}$
  int leftMultIn(WeylElt& w, Generator s) const;

  // get WeylWord for piece |j| of |w|
  const WeylWord& wordPiece(const WeylElt& w, Generator j) const
    { return d_transducer[j].wordPiece(w[j]); }

  /*!
    \brief first generator $<s$ not commuting with |s|, or |s| if none exist
  */
  Generator min_neighbor (Generator s) const { return d_min_star[s]; }

  WeylElt genIn (Generator i) const; // $s_i$ as Weyl group element

// From this point on each Generator or WeylWord uses external numbering
public:

// constructors and destructors
  WeylGroup(const latticetypes::LatticeMatrix& cartan); // from Cartan matrix

// accessors

  WeylElt generator (Generator i) const // $s_i$ as Weyl group element
    { return genIn(d_in[i]); }

  // set |w=ws|, return length change
  int mult(WeylElt& w, Generator s) const { return multIn(w,d_in[s]); }

  // set |w=wv|
  void mult(WeylElt& w, const WeylElt& v) const;

  // set |w=wv|
  void mult(WeylElt&, const WeylWord&) const;

  // set |w=sw| (argument order motivated by modification effect on |w|)
  int leftMult(WeylElt& w, Generator s) const { return leftMultIn(w,d_in[s]); }
  void leftMult(WeylElt& w, const WeylWord& ww) const
  {
    for (size_t i=ww.size(); i-->0; ) // use letters from right to left
      leftMult(w,ww[i]);
  }
  void leftMult(WeylElt& w, const WeylElt& x) const;

  WeylElt prod(const WeylElt& w, Generator s) const
    { WeylElt result=w; mult(result,s); return result; }
  WeylElt prod(Generator s, const WeylElt& w) const
    { WeylElt result=w; leftMult(result,s); return result; }
  WeylElt prod(const WeylElt& w, const WeylElt& v) const
    { WeylElt result=w; mult(result,v); return result; }
  WeylElt prod(const WeylElt& w, const WeylWord& ww) const
    { WeylElt result=w; mult(result,ww); return result; }
  WeylElt prod(const WeylWord& ww,const WeylElt& w) const
    { WeylElt result=w; leftMult(result,ww); return result; }


  /* These additional definitions would be needed if TwistedInvolutions were a
     type distinct from WeylElt (but they are not allowed as it is):

  void leftMult(TwistedInvolution& w, Generator s) const {
    leftMultIn(w.contents(),d_in[s]);
  }
  void leftMult(TwistedInvolution& w, const WeylWord& ww) const {
    leftMult(w.contents(),ww);
  }
  void leftMult(TwistedInvolution& w, const WeylElt& x) const {
    leftMult(w.contents(),x);
  }
  */


  unsigned int length(const WeylElt&) const;

  // return (only) $l(w.s)-l(w)$
  int length_change(WeylElt w, Generator s) const
  {
    return multIn(w,d_in[s]); // discard copied |w|
  }

  // return $l(s.w)-l(w)$
  int length_change(Generator s,const WeylElt& w) const
  {
    return hasDescent(s,w) ? -1 : 1;
  }

  // letter |i| of Weyl word for |w|
  Generator letter(const WeylElt& w, unsigned int i) const;

  void conjugacyClass(WeylEltList&, const WeylElt&) const;

  /*! \brief Conjugates |w| by the generator |s|: |w=sws|. */
  void conjugate(WeylElt& w, Generator s) const
    { leftMult(w,s); mult(w,s); }

  /*! \brief $w^(-1)$ */
  WeylElt inverse(const WeylElt& w) const;

  /*! \brief set |w=w^(-1)| */
  void invert(WeylElt& w) const { w=inverse(w); }

  /*! \brief return a left descent for |w|, or |UndefGenerator| if |w==e|

   In fact the descent whose internal number is lowest is returned, but
   converted to external numbering.
 */
  Generator leftDescent(const WeylElt& w) const;

  const WeylElt& longest() const { return d_longest; }

  // correspondence with dual Weyl group; is coherent with twist on right
  WeylElt opposite (const WeylElt& w) const { return prod(w,d_longest); }

  unsigned long maxlength() const { return d_maxlength; }

  /*! \brief the order of the Weyl group */
  const size::Size& order() const { return d_order; }

  size_t rank() const { return d_rank; }

  WeylWord word(const WeylElt& w) const;
  WeylElt element(const WeylWord& ww) const { return prod(WeylElt(),ww); }

  WeylEltList reflections() const; // list all reflections

  /* give representation of |w| as integral number, supposing this fits. */
  unsigned long toUlong(const WeylElt& w) const;

  /* inverse operation of |toUlong| */
  WeylElt toWeylElt(unsigned long) const;

  bool hasDescent(Generator, const WeylElt&) const;

  // apply automorphism of $(W,S)$ given by |f| in terms of outer numbering
  WeylElt translation(const WeylElt& w, const WeylInterface& f) const;

  void translate(WeylElt& w, const WeylInterface& i) const
    { w=translation(w,i); }

  // standard reflection action of Weyl group using a root datum
  void act(const rootdata::RootDatum& rd,
	   const WeylElt& w,
	   latticetypes::LatticeElt& v) const;
/*!
  \brief Nondestructive version of |act| method
*/
  latticetypes::LatticeElt
    imageBy(const rootdata::RootDatum& rd,
	    const WeylElt& w,
	    latticetypes::LatticeElt v) const
    { act(rd,w,v); return v; }

  void inverseAct(const rootdata::RootDatum& rd,
		  const WeylElt& w,
		  latticetypes::LatticeElt& v) const;

/*!
  \brief Nondestructive version of |inverseAct| method
*/
  latticetypes::LatticeElt
    imageByInverse(const rootdata::RootDatum& rd,
		   const WeylElt& w,
		   latticetypes::LatticeElt v) const
    { inverseAct(rd,w,v); return v; }


// manipulators

}; // |class WeylGroup|

/*
  We have split off in the following class functionality that involves an
  involutive diagram automorphism |delta|. While WeylElt values are assumed to
  live just in the Weyl group W, this derived class implements operations that
  are more naturally interprests them in $W semidirect Z/2Z$ (the second
  factor acting on the first via |delta|), namely in the non-identity coset
  for $W$. This derived class is totally ignorant of the internal numbering
*/
class TwistedWeylGroup
{
  const WeylGroup& W;  // non-owned reference
  const Twist d_twist; // cannot be reference, if dual is to be constructible

public:
  TwistedWeylGroup(const WeylGroup&, const Twist&);

  // construct the "dual" twisted Weyl group: differs by a dual twist
  TwistedWeylGroup(const TwistedWeylGroup&, tags::DualTag);

  const WeylGroup& weylGroup() const { return W; }
  size_t rank() const { return W.rank(); }

  int mult(WeylElt& w, Generator s) const { return W.mult(w,s); }
  void mult(WeylElt& w, const WeylElt& v) const { W.mult(w,v); }
  void mult(WeylElt& w, const WeylWord& ww) const { W.mult(w,ww); }

  int leftMult(WeylElt& w, Generator s) const { return W.leftMult(w,s); }
  void leftMult(WeylElt& w, const WeylWord& ww) const { W.leftMult(w,ww); }
  WeylWord word(const WeylElt& w) const { return W.word(w); }

  WeylElt prod(const WeylElt& w, Generator s) const { return W.prod(w,s); }
  WeylElt prod(Generator s, const WeylElt& w) const { return W.prod(s,w); }
  WeylElt prod(const WeylElt& w, const WeylElt& v) const { return W.prod(w,v); }
  WeylElt prod(const WeylElt& w, const WeylWord& ww) const
   { return W.prod(w,ww); }
  WeylElt prod(const WeylWord& ww,const WeylElt& w) const
   { return W.prod(ww,w); }

  Generator twisted(Generator s) const { return d_twist[s]; }
  WeylElt twisted(const WeylElt& w) const { return W.translation(w,d_twist); }

  const Twist& twist() const { return d_twist; } // noun "twist"
  void twist(WeylElt& w) const { w=twisted(w); } // verb "twist"

  Twist dual_twist() const; // the twist for the dual twisted Weyl group

  /*!
     \brief Twisted conjugates element |tw| by the generator |s|:
     \f$tw:=s.tw.\delta(s)\f$.
   */
  void twistedConjugate(TwistedInvolution& tw, Generator s) const
  {
    WeylElt& w=tw.contents();
    W.leftMult(w,s);
    W.mult(w,d_twist[s]);
  }
  void twistedConjugate(TwistedInvolution& tw, const WeylWord& ww) const
  {
    for (size_t i=ww.size(); i-->0; )
      twistedConjugate(tw,ww[i]);
  }
  void inverseTwistedConjugate(TwistedInvolution& tw, const WeylWord& ww) const
  {
    for (size_t i=0; i<ww.size(); ++i )
      twistedConjugate(tw,ww[i]);
  }

  TwistedInvolution twistedConjugated(const TwistedInvolution& tw, Generator s)
    const
  {
    TwistedInvolution result=tw; twistedConjugate(result,s); return result;
  }

  void twistedConjugacyClass(TwistedInvolutionList&, const TwistedInvolution&)
    const;

  /*!\brief
    Twisted conjugates element |tw| by |w|: \f$tw:=w.tw.\delta(w^{-1})\f$.
   */
  void twistedConjugate(TwistedInvolution& tw, const WeylElt& w) const;

  bool hasTwistedCommutation(Generator, const TwistedInvolution&) const;

/*!
  \brief Returns the length of |tw| as a twisted involution.
*/
  unsigned long involutionLength(const TwistedInvolution& tw) const;

/*!
  \brief Returns a reduced expression of |tw| as a twisted involution.

  Positive entries correspond to generators of Cayley transforms, while cross
  actions are indicated by making the generator negative by bitwise complement.
*/
  std::vector<signed char>
    involution_expr(TwistedInvolution tw) const; // call by value

  //!\brief Roots that are images of the simple roots under involution of |tw|
  rootdata::RootList simple_images
    (const rootdata::RootSystem& rs, const TwistedInvolution& tw) const;

//!\brief Matrix of involution of |tw| in adjoint coordinates
  latticetypes::LatticeMatrix
    involution_matrix(const rootdata::RootSystem& rs,
		      const TwistedInvolution& tw) const;

}; // |class TwistedWeylGroup|


// A small derived class to allow hash tables of TwistedInvolution values
// (We might have put |Pooltype| and |hashCode| members in WeylElt directly)

struct TI_Entry
  : public weyl::TwistedInvolution
{
  TI_Entry(const weyl::TwistedInvolution& tw): weyl::TwistedInvolution(tw) {}
  TI_Entry(): weyl::TwistedInvolution() {}

  // members required for an Entry parameter to the HashTable template
  typedef std::vector<TI_Entry> Pooltype; // associated storage type
  size_t hashCode(size_t modulus) const;  // hash function
}; // class TI_Entry


}

}

#endif
