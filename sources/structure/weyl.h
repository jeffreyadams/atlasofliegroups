/*
  This is weyl.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017,2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
// Class definitions and function declarations for WeylGroup.

#ifndef WEYL_H  /* guard against multiple inclusions */
#define WEYL_H

#include "../Atlas.h"

#include <cstring>
#include <cassert>
#include <initializer_list>
#include <array>

#include "constants.h"
#include "tags.h"

#include "size.h"   // for stored order of the Weyl group
#include "matrix.h" // for Coxeter matrix


/******** type declarations *************************************************/

namespace atlas {

namespace weyl {


/******** constant declarations *********************************************/

  const unsigned char UndefValue = constants::ucharMax;
  const unsigned long UndefOrder = 0; // is never a valid value


/******** type definitions **************************************************/

  // A mapping between one interpretation of |Generator|s and another
  class Twist // also used under the typedef name |WeylInterface|
  {
    Generator d[constants::RANK_MAX];
  public:
    Twist () // assures definite values are always set
    { std::fill_n(&d[0],constants::RANK_MAX,Generator(~0)); }
    Twist(const ext_gens& orbits);
    Generator& operator[] (Generator i) { return d[i]; }
    const Generator& operator[] (Generator i) const { return d[i]; }

    bool operator!=(const Twist& y) const
    { for (unsigned i=0; i<constants::RANK_MAX; ++i)
	if (d[i]!=y.d[i])
	  return true;
      return false;
    }
    bool operator==(const Twist& y) const { return not operator!=(y); }

  };


/******** function declarations *********************************************/

Twist make_twist(const RootDatum& rd,
		 const WeightInvolution& d);


WeylEltList conjugacy_class(const WeylGroup& W, const WeylElt& w);

/******** main class definitions *******************************************/

class WeylWord : public std::vector<Generator>
{ // often used to accumulate short words, we reserve 32 letters in advance
  using Base = std::vector<Generator>;
 public:
  WeylWord () : Base () { Base::reserve(32); }
  WeylWord (Base&& v) : Base(std::move(v)) {}
  WeylWord (const Base& v) : Base(v) {}
  WeylWord (std::initializer_list<Generator> l) : Base(l.begin(), l.end()) {}
};

// Numbering of minimal length coset representative in some $W_{i-1}\\W_i$.
  using EltPiece = unsigned char;

/*
  Representation of an individual element of a Weyl group.

  The representation is described in detail in the description of the
  class WeylGroup.  An array of RANK_MAX unsigned char, the ith
  representing a shortest length coset representative of a parabolic
  subquotient W_{i-1}\\W_i.

  The class |WeylElt| is also known under the typedef name |TwistedInvolution|
*/
class WeylElt
{
  friend class WeylGroup; // so we can shield implementation from all others

 private:

/*
  A Weyl group element represented as a product of shortest length coset
  representatives for parabolic subquotients $W_{i-1}\\W_i$.

  Entry number $i-1$ is an unsigned char parametrizing the $i$-th coset
  representative $w_i$ for an element of $W_{i-1}\\W_i$.
  Then $w = w_1.w_2...w_n$.
  */
  std::array<EltPiece,constants::RANK_MAX> pieces;

 public:

// constructors and destructors

  // Default constructor; set element to the identity element of W.
  WeylElt() { pieces.fill(0); }

// copy and assignment
// use standard definitions (raw copy of array)

// accessors

  bool operator== (const WeylElt& w) const { return pieces==w.pieces; }
  bool operator!= (const WeylElt& w) const { return pieces!=w.pieces; }

  // lexicograpic ordering by internal parabolic quotient list
  bool operator< (const WeylElt& w) const { return pieces<w.pieces; }

protected: // these are for |WeylGroup| and |TI_Entry|'s eyes only

  // Get an individual |EltPiece|
  EltPiece piece (Generator j) const { return pieces[j]; }

// manipulators
  EltPiece& piece (Generator j)      { return pieces[j]; }

public:
// dummy methods that mark transition of interpretation

  // from WeylElt interpretation to TwistedInvolution
  // nothing: cast to TwistedInvolution is allowed, would do construction if
  // TwistedInvolution were a derived class, and is a no-op in reality

  // from TwistedInvolution interpretation to WeylElt
  // calls of these mark referral to the "base" class of |TwistedInvolution|
  const WeylElt& w() const { return *this; }
  WeylElt& contents() { return *this; }

}; // |class WeylElt|


/*
  The |WeylGroup| class provides supporting context, whose methods allow
  light-weight |WeylElt| objects to be manipulated as Weyl groups elements.
  (The Atlas software uses many such (supporting context,value) class pairs.)

  It uses a variation on Fokko's own and favourite implementation of Coxeter
  groups in terms of transducers. What follows are Fokko's own words:

  [Since this is Fokko's mathematics as well as his code, I have tried to keep
  his voice. Occasional amplifications are in brackets. DV 7/23/06.]

   I have tried to make a careful choice of datatype for the group elements in
   order to strike the right balance between efficiency and generality. This
   has led me to the choice of _fixed size_ arrays of unsigned chars,
   representing the "parabolic subquotient" representation of a given element;
   in other words, in a group of rank n, the first n elements of the array are
   used (but since it is fixed size, we have to allocate it to some constant
   value, namely RANK_MAX.)

  [Meaning of "parabolic subquotient representation": one orders the simple
  generators in an appropriate fashion s_1,...,s_n (henceforth called the
  "internal representation"). Define W_i = <s_1,...,s_i>, so that W_0 = {1}
  subset W_1 subset W_2 subset ... subset W_n = W. Each W_i is a parabolic
  subgroup of W, so its length function is the restriction of that of W. Each
  right coset in W_{i-1}\\W_i has a unique minimal length coset
  representative. Write N_i for the number of cosets (or representatives).
  Then the cardinality of W is N_1.N_2...N_n. More precisely, an element of W
  has a unique factorization as x_1.x_2...x_n, with x_i one of the N_i minimal
  length coset representatives of W_{i-1}\\W_i.

  The key to the implementation is to order the generators in such a way that
  every parabolic subquotient W_{i-1}\\W_i has cardinality fitting into an
  unsigned char: that is, cardinality at most 256. (This is possible for any
  Weyl group of rank at most 128.) A WeylElt is an array of unsigned char of
  size RANK_MAX; the unsigned char in entry i is the number of the factor x_i.

  There is a little more about how the group structure is computed in
  this representation in the description of the class Transducer.  The
  mathematical reference is Fokko du Cloux, "A transducer approach to
  Coxeter groups," J. Symbolic Computation 27 (1999), 311-324.]

   This is not so wasteful as it may sound : of course the
   representation as a single number is more compact, but will overflow
   even on 64-bit machines for some groups of rank <= 16, and much
   sooner of course on 32-bit machines; also it imposes some
   computational overhead for packing and unpacking.

  [Meaning of "representation as a single number:" if the cardinality of W
  fits in an unsigned long, then an element of W may be represented as the
  unsigned long

  x_1 + x_2(N_1) + x_3(N_1N_2) + ... +x_n(N_1...N_n).

  In this formula I have changed the meaning of x_i: now it means the integer
  between 0 and N_i enumerating the coset representative previously called
  x_i. On a 32 bit machine, this representation works for W(E8), but not for
  W(A12) (which has order 13!, which is about 1.5 x 2^{32}).]

   The variable-size structure of an STL vector uses already three pointer
   values for its control structure (the addresses of the start, fill level
   and end of the dynamically allocated storage), and then it requires that
   allocated storage itself. With a special purpose data type this could
   perhaps be simplified to just a pointer (if the size matches capacity
   and is known from the context, as here the semisimple rank of the group)
   but you still have the big time overhead of allocating and deallocating
   memory from the heap, and the obligation to implement your own memory
   management... A lot of work for not "wasting" some byte in a fixed array
   [Paragraph edited for clarity, see revision history for the original; MvL]

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
   group carries a little interface that will translate back and forth from
   the external ordering and the internal one.

  In reality |WeylElt| values are not stored in large quantities in the Atlas
  software, so Fokko's considerations of storage requirements for them are not
  very important. They serve mainly for computational navigation while
  determining the set of twisted involutions in the (twisted) Weyl group, and
  these form a much smaller set for which identification by a single number is
  feasible in all cases. Indeed we even tabulate the full set of twisted
  involutions, which admittedly limits the size of groups that can be handled.

  To speed up some computations we precompute in |d_min_star| for each
  generator the smallest other generator with which it does not commute (or
  the generator itself if there are none).
*/
class WeylGroup
{
  struct Transducer;

  std::vector<Transducer> d_transducer; // list of parabolic quotients
  WeylInterface d_in;  // conversion from external numbering to internal
  WeylInterface d_out; // conversion from internal numbering to external
  std::vector<Generator> d_min_star; // diagram dependent field for efficiency
  std::vector<Generator> upper; // of diagram component

// private member functions
// these interpret Generators and WeylWords internally, so not for public use!

  // first transducer to use when shifting in |s|
  Generator start_gen(Generator s) const { return upper[s]; }

  // transform local generator |g| from transducer |i| to outer numbering
  Generator output_local_gen(Generator i, Generator g) const;

  // auxiliary doing the main work for |inner_mult|
  int transduce(WeylElt& w, Generator start, Generator local_s) const;
  // set |w=ws|, return value is $l(ws)-l(w)\in\{+1,-1}$
  int inner_mult(WeylElt& w, Generator s) const;

  // set |w=wv| where |v| is peice |x| of transducer |i|
  // return value is $l(wv)-l(w)-l(v)\in -2\N$
  int mult_by_piece(WeylElt& w, const WeylElt& v, Generator i) const;

  // set |w=sw|, return value is $l(sw)-l(w)\in\{+1,-1}$
  int inner_left_mult(WeylElt& w, Generator s) const;

  WeylElt inner_gen (Generator i) const; // $s_i$ as Weyl group element
  bool inner_commutes (Generator s, Generator t) const;

  // First generator $<s$ not commuting with |s|, or |s| if none exist; inner
  Generator min_neighbor (Generator s) const { return d_min_star[s]; }

// From this point on each Generator or WeylWord uses external numbering
public:

// constructors and destructors
  WeylGroup(const int_Matrix& cartan); // from Cartan matrix

  WeylGroup(const WeylGroup&) = delete; // Weyl groups should be shared
  WeylGroup(WeylGroup&& W); // but they can be moved

  ~WeylGroup(); // default, but must be declared as it calls |~Transducer|

// accessors

  WeylElt generator (Generator i) const // $s_i$ as Weyl group element
    { return inner_gen(d_in[i]); }

  // set |w=ws|, return length change
  int mult(WeylElt& w, Generator s) const { return inner_mult(w,d_in[s]); }

  // set |w=wv|
  void mult(WeylElt& w, const WeylElt& v) const;

  // set |w=wv|
  void mult(WeylElt&, const WeylWord&) const;

  // set |w=sw| (argument order motivated by modification effect on |w|)
  int left_multiply(WeylElt& w, Generator s) const
  { return inner_left_mult(w,d_in[s]); }
  void left_multiply(WeylElt& w, const WeylWord& ww) const
  {
    for (unsigned int i=ww.size(); i-->0; ) // use letters from right to left
      left_multiply(w,ww[i]);
  }
  void left_multiply(WeylElt& w, const WeylElt& x) const; // |w=xw|

  WeylElt prod(WeylElt w, Generator s) const { mult(w,s); return w; }
  WeylElt prod(Generator s, WeylElt w) const { left_multiply(w,s); return w; }
  WeylElt prod(WeylElt w, const WeylElt& v) const { mult(w,v); return w; }
  WeylElt prod(WeylElt w, const WeylWord& ww) const { mult(w,ww); return w; }
  WeylElt prod(const WeylWord& ww, WeylElt w) const
    { left_multiply(w,ww); return w; }

  /* These additional definitions would be needed if TwistedInvolutions were a
     type distinct from WeylElt (but they are not allowed as it is):

  void left_multiply(TwistedInvolution& w, Generator s) const {
    left_multiplyIn(w.contents(),d_in[s]);
  }
  void left_multiply(TwistedInvolution& w, const WeylWord& ww) const {
    left_multiply(w.contents(),ww);
  }
  void left_multiply(TwistedInvolution& w, const WeylElt& x) const {
    left_multiply(w.contents(),x);
  }
  */


  unsigned int length(const WeylElt&) const;

  // return (only) $l(w*s)-l(w)\in\{ -1, 1 \}$
  int length_change(WeylElt w, Generator s) const
  {
    return inner_mult(w,d_in[s]); // discard copied then modified |w|
  }

  // return $l(s*w)-l(w)$
  int length_change(Generator s,const WeylElt& w) const
  {
    return has_descent(s,w) ? -1 : 1;
  }

  // Conjugate |w| by the generator |s|: |w=sws|.
  void conjugate(WeylElt& w, Generator s) const
    { left_multiply(w,s); mult(w,s); }

  WeylElt inverse(const WeylElt& w) const;
  void invert(WeylElt& w) const { w=inverse(w); }

  /*
   Find a left descent for |w|, or |UndefGenerator| if |w==e|

   In fact the descent whose internal number is lowest is returned, but
   converted to external numbering.
 */
  Generator leftDescent(const WeylElt& w) const;

  WeylElt longest() const;

  unsigned int max_length() const;

  // the order of the Weyl group
  size::Size order() const;

  Generator rank() const // |d_transducer.size()| makes compiler unhappy here
  { return d_min_star.size(); } // and also |upper.size()|

  bool commutes (Generator s, Generator t) const
  { return s==t or inner_commutes(d_in[s],d_in[t]); }

  WeylWord word(const WeylElt& w) const;
  WeylElt element(const WeylWord& ww) const { return prod(WeylElt(),ww); }

  Generator Chevalley_dual(Generator s) const; // conjugation by |longest()|
  Twist Chevalley_twist () const; // combine values for all |s|

  // give representation of |w| as integral number.
  arithmetic::big_int to_big_int(const WeylElt& w) const;

  // inverse operation of |to_big_int|
  WeylElt to_WeylElt(arithmetic::big_int) const;

  bool has_descent(Generator, const WeylElt&) const; // on the left
  bool has_descent(const WeylElt&, Generator) const; // on the right

  // apply automorphism of $(W,S)$ given by |f| in terms of outer numbering
  WeylElt translation(const WeylElt& w, const WeylInterface& f) const;
  // same as |translation| but for $w^{-1}$ instead of $w$; this is faster!
  WeylElt reverse_translation(const WeylElt& w, const WeylInterface& f) const;

  void translate(WeylElt& w, const WeylInterface& i) const
    { w=translation(w,i); }

  // reflection action of Weyl group on a root
  void act(const RootSystem& rd, const WeylElt& w, RootNbr& alpha) const;
  // standard reflection action of Weyl group on weights for a root datum
  template<typename C>
    void act(const RootDatum& rd, const WeylElt& w, matrix::Vector<C>& v)
    const;
  // standard reflection action of Weyl group on coweights for a root datum
  template<typename C>
    void co_act(const RootDatum& rd, matrix::Vector<C>& v, const WeylElt& w)
    const;
  // standard reflection action of Weyl group using a root datum
  void act(const RootDatum& rd, const WeylElt& w, RatWeight& v) const;
  // standard reflection left action of Weyl group using a root datum
  void act(const RootDatum& rd, const WeylElt& w, LatticeMatrix& M) const;

  // same using only lists of simple (co)roots avoiding construction root datum
  template<typename C>
    void act(const PreRootDatum& prd, const WeylElt& w, matrix::Vector<C>& v)
    const;
  void act(const PreRootDatum& prd, const WeylElt& w, RatWeight& v) const;
  void act(const PreRootDatum& prd, const WeylElt& w, LatticeMatrix& M) const;
  // Nondestructive version of |act| method
  Weight
    image_by(const RootDatum& rd, const WeylElt& w, Weight v) const
    { act(rd,w,v); return v; }

  void inverse_act(const RootDatum& rd, const WeylElt& w, Weight& v) const;

  // Nondestructive version of |inverse_act| method
  Weight
    image_by_inverse(const RootDatum& rd, const WeylElt& w, Weight v) const
    { inverse_act(rd,w,v); return v; }

  Coweight
    image_by(const RootDatum& rd, Coweight v, const WeylElt& w) const
    { co_act(rd,v,w); return v; }

// manipulators: nothing can modify the |WeylGroup| itself after construction

}; // |class WeylGroup|


/*
  We have split off in the following class functionality that involves an
  involutive diagram automorphism |delta|. While WeylElt values are assumed to
  live just in the Weyl group W, this class implements operations that
  are more naturally interprests them in $W semidirect Z/2Z$ (the second
  factor acting on the first via |delta|), namely in the non-identity coset
  for $W$.

  This class is not derived from |WeylGroup|, although many methods are just
  handed over to an underlying |WeylGroup|. We assume the latter to be
  constructed beforehand, and store a reference to it. For one thing, this
  makes it possible for two |TwistedWeylGroup|s with different twists to share
  the same |WeylGroup|; this is useful because the Weyl groups of a root datum
  and its dual can be identified, but not their twists.

  Having no privilege with respect to |WeylGroup|, we are totally ignorant of
  its internal numbering.
*/
class TwistedWeylGroup
{
  const WeylGroup& W;  // non-owned reference
  const Twist d_twist; // cannot be reference, if dual is to be constructible

  void operator=(const TwistedWeylGroup&); // forbid assignment
 protected:
  TwistedWeylGroup(const TwistedWeylGroup& g) // only derived classes may copy
   : W(g.W), d_twist(g.d_twist) {} // share |W|, copy |d_twist|
public:
  TwistedWeylGroup(const WeylGroup&, const Twist&);

  // construct the "dual" twisted Weyl group: differs by a dual twist
  TwistedWeylGroup(const TwistedWeylGroup&, tags::DualTag);

  const WeylGroup& weylGroup() const { return W; }
  Generator rank() const { return W.rank(); }

  int mult(WeylElt& w, Generator s) const { return W.mult(w,s); }
  void mult(WeylElt& w, const WeylElt& v) const { W.mult(w,v); }
  void mult(WeylElt& w, const WeylWord& ww) const { W.mult(w,ww); }

  int left_multiply(WeylElt& w, Generator s) const
    { return W.left_multiply(w,s); }
  void left_multiply(WeylElt& w, const WeylWord& ww) const
    { W.left_multiply(w,ww); }
  WeylWord word(const WeylElt& w) const { return W.word(w); }

  WeylElt prod(const WeylElt& w, Generator s) const { return W.prod(w,s); }
  WeylElt prod(Generator s, const WeylElt& w) const { return W.prod(s,w); }
  WeylElt prod(const WeylElt& w, const WeylElt& v) const { return W.prod(w,v); }
  WeylElt prod(const WeylElt& w, const WeylWord& ww) const
    { return W.prod(w,ww); }
  WeylElt prod(const WeylWord& ww,const WeylElt& w) const
    { return W.prod(ww,w); }

  bool has_descent(Generator s, const WeylElt& w) const
    { return W.has_descent(s,w); }
  bool has_descent(const WeylElt& w, Generator s) const
    { return W.has_descent(w,s); }

  Generator twisted(Generator s) const { return d_twist[s]; }
  WeylElt twisted(const WeylElt& w) const { return W.translation(w,d_twist); }

  Generator dual_twisted(Generator s) const
  { return W.Chevalley_dual(d_twist[s]); }
  WeylElt dual_twisted(const WeylElt& w) const;

  const Twist& twist() const { return d_twist; } // noun "twist"
  void twist(WeylElt& w) const { w=twisted(w); } // verb "twist"

  ext_gens twist_orbits () const;
  Twist dual_twist() const; // the twist for the dual twisted Weyl group

  /*
     Twisted conjugate element |tw| by the generator |s|:
     $tw:=s.tw.\delta(s)$. Returns length change, in $\{-2,0,2\}$. For
     consistency, these functions should have |s| resp. |ww| as first argument
   */
  int twistedConjugate(TwistedInvolution& tw, Generator s) const
  {
    WeylElt& w=tw.contents();
    int d = W.left_multiply(w,s);
    return d+W.mult(w,d_twist[s]);
  }
  int twistedConjugate(const WeylWord& ww,TwistedInvolution& tw) const
  {
    int d=0;
    for (unsigned int i=ww.size(); i-->0; )
      d+=twistedConjugate(tw,ww[i]);
    return d;
  }
  int inverseTwistedConjugate(TwistedInvolution& tw, const WeylWord& ww) const
  {
    int d=0;
    for (unsigned int i=0; i<ww.size(); ++i )
      d+=twistedConjugate(tw,ww[i]);
    return d;
  }

  TwistedInvolution twistedConjugated(TwistedInvolution tw, Generator s) const
    { twistedConjugate(tw,s); return tw; }

  void twistedConjugacyClass(TwistedInvolutionList&, const TwistedInvolution&)
    const;

  // Twisted conjugate element |tw| by |w|: $tw:=w.tw.\delta(w^{-1})$.
  void twistedConjugate(TwistedInvolution& tw, const WeylElt& w) const;

  bool hasTwistedCommutation(Generator, const TwistedInvolution&) const;

  // The length of |tw| as a twisted involution.
  unsigned int involutionLength(const TwistedInvolution& tw) const;

/*
  Return a reduced expression of |tw| as a twisted involution.
  The second form assures the lexicographically first possible form is used
*/
   InvolutionWord involution_expr(TwistedInvolution tw) const; // call by value
   InvolutionWord canonical_involution_expr(TwistedInvolution tw) const; // idem

   // extended involution expression for a twist-fixed twisted involution
   InvolutionWord extended_involution_expr(TwistedInvolution tw) const; // idem

  // Roots that are images of the simple roots under the true involution of |tw|
  RootNbrList simple_images
    (const RootSystem& rs, const TwistedInvolution& tw) const;

// Matrix of the true involution of |tw|, in adjoint coordinates
  WeightInvolution
    involution_matrix(const RootSystem& rs, const TwistedInvolution& tw) const;

}; // |class TwistedWeylGroup|


// A small derived class to allow hash tables of TwistedInvolution values
// (We might have put |Pooltype| and |hashCode| members in WeylElt directly)

struct TI_Entry
  : public TwistedInvolution
{
  TI_Entry(const TwistedInvolution& tw): TwistedInvolution(tw) {}
  TI_Entry(): TwistedInvolution() {} // identity TwistedInvolution

  // members required for an Entry parameter to the HashTable template
  typedef std::vector<TI_Entry> Pooltype; // associated storage type
  size_t hashCode(size_t modulus) const;  // hash function
}; // class TI_Entry


} // |namespace weyl|

} // |namespace atlas|

#endif
