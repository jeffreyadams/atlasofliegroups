/*!
\file
\brief Implementation of WeylGroup.

  I have decided to represent elements as fixed-size arrays of
  unsigned characters. This forces expressing things in the standard
  ordering of the generators, and hence to have a small i/o interface
  for resetting the numbering to and from the numbering used by the
  outside world.

  It has seemed to me that this is the best compromise between size of the
  dataype, generality and efficiency.

*/
/*
  This is weyl.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "weyl.h"

#include <algorithm>
#include <set>
#include <stack>

#include "dynkin.h"
#include "setutils.h"
#include "rootdata.h"

/*****************************************************************************

  Implementation of Weyl groups. I have decided to represent elements as
  fixed-size arrays of unsigned characters. This forces to express things
  in the standard ordering of the generators, and hence to have a small
  i/o interface for resetting the numbering to and from the outside world's.

  It has seemed to me that this is the best compromise between size of the
  dataype, generality and efficiency.

******************************************************************************/

namespace atlas {

namespace {

  using namespace latticetypes;
  using namespace setutils;
  using namespace weyl;

  void fillCoxMatrix(LatticeMatrix&, const LatticeMatrix&, const Permutation&);

  EltPiece dihedralMin(const Transducer&, EltPiece, Generator, Generator);

  EltPiece dihedralShift(const Transducer&, EltPiece, Generator, Generator,
			 unsigned long);

  /* This class exists only so that its constructor (for a static object) is
     run exactly once, before any Weyl group needs the identity twist. It
     would not be necessary if Twist had been a class rather than a typedef.
  */
  class IdentityTwist{
  private:
    Twist d_twist;
  public:
    IdentityTwist() {
      for (size_t j = 0; j < constants::RANK_MAX; ++j)
	d_twist[j] = j;
    }
    const Twist& twist() {return d_twist;}
  };

}

/*****************************************************************************

        Chapter I -- The WeylGroup and WeylElt classes

  The Weyl group class is a variation on my favourite implementation in terms
  of transducers.

  I have tried to make a careful choice of datatype for the group elements in
  order to strike the right balance between efficiency and generality. This has
  led me to the choice of _fixed size_ arrays of unsigned chars, representing
  the "parabolic subquotient" representation of a given element; in other
  words, in a group of rank n, the first n elements of the array are used
  (but since it is fixed size, we have to allocate it to some constant value,
  namely RANK_MAX.) This is not so wasteful as it may sound : of course the
  representation as a single number is more compact, but will overflow even
  on 64-bit machines for some groups of rank <= 16, and much sooner of course
  on 32-bit machines; also it imposes some computational overhead for packing
  and unpacking. Any variable-size structure like the STL vector uses already
  three unsigned longs for its control structure (the address of the data, the
  size and the capacity), and then it still has to allocate. This could
  perhaps be simplified to just a pointer (after all the size of the
  allocation is known to the group) but you still have the big overhead of
  allocating and deallocating memory from the heap, and remembering to delete
  the pointers when they go out of scope, or else use autopointers ...

  And if worst comes to worst and one really is worried about a factor 2
  waste for type E8 (which might be significant if one deals with huge
  containers of group elements), then one can still recompile with MAX_RANK=8,
  which will then give a datatype that in 64 bits is the same size as an
  unsigned long.

  Notice that the unsigned char type miraculously suffices for all subquotients
  of all groups up to rank 128 (indeed, the biggest subquotient for B128 is
  of order 256), _provided_ the generators are enumerated in an appropriate
  order. This forces us to do quite a bit of type recognition, which is
  relegated to the dynkin namespace. Because of this reordering, the group
  carries a little interface that will translate back and forth from the
  external ordering and the internal one.

******************************************************************************/

namespace weyl {

/*!
  \brief Build the Weyl group corresponding to the Cartan matrix c,
  and incorporate the possibly given twist (take identity if twist==NULL).

  NOTE : |c| and |twist| use some (consistent) labelling of simple roots,
  but we will determine an internal renumbering making the subquotients small
*/
WeylGroup::WeylGroup(const latticetypes::LatticeMatrix& c, const Twist* twist)
  : d_rank(c.numColumns()) // all other fields are computed later
  , d_order() // this initialises this size::Size value to 1
  , d_maxlength(0)
  , d_longest()
  , d_coxeterMatrix()
  , d_transducer(d_rank)
  , d_twist() // not copied from twist even if non-NULL: renumbering is needed
  , d_in()    // being arrays, |d_in| and |d_out| cannot be initialised
  , d_out()
  , d_min_star(d_rank)

{
  /* analyse the Coxeter matrix */

  setutils::Permutation a;
  dynkin::DynkinDiagram d(c); // make diagram from Cartan matrix
  normalize(a,d);  // find renumbering |a| putting labels in canonical order

  /* now put appropriate permutations into |d_in| and |d_out|, so that
     internal number |j| gives external number |d_out[j]|, and of course
     |d_in[d_out[j]]==j|. By convention of |normalize|, this means |d_out==a|
   */
  for (size_t j = 0; j < d_rank; ++j) {
    d_out[j] = a[j];
    d_in[d_out[j]] = j;
  }


  // the default twist is the identity twist
  if (twist == NULL)
    copy(d_twist,identityTwist()); // no need for renumbering; id=id always
  else
    for (size_t j = 0; j < d_rank; ++j)
      d_twist[j] = d_in[(*twist)[d_out[j]]]; // |j| internal, |twist| external

  // deduce the Coxeter matrix from Cartan matrix |c| and renumbering |a|
  fillCoxMatrix(d_coxeterMatrix,c,a);

  // now construct the transducer
  for (size_t j = 0; j < d_rank; ++j)
    d_transducer[j] = Transducer(d_coxeterMatrix,j);

  // the longest element has the maximal valid value in each of its pieces
  for (size_t j = 0; j < d_rank; ++j)
    d_longest[j] = d_transducer[j].size()-1;

  // its length is obtained by summing the (maximal) lengths of its pieces
  for (size_t j = 0; j < d_rank; ++j)
    d_maxlength += d_transducer[j].maxlength();

  // and the Weyl group size is the product of the numbers of transducer states
  for (size_t j = 0; j < d_rank; ++j) {
    d_order *= d_transducer[j].size();
  }

  // precompute for each |j| the first non-commuting or equal |i<=j|
  for (size_t j = 0; j < d_rank; ++j)
    for (size_t i=0; i<=j; ++i)
      if (d_coxeterMatrix(i,j)!=2)
      {
	d_min_star[j]=i; break;
      }

}

  /* the "dual" Weyl group: the only difference with W is that the twist is
      multiplied by conjugation with the longest element.
  */
WeylGroup::WeylGroup(const WeylGroup& W, tags::DualTag)
  : d_rank(W.d_rank)
  , d_order(W.d_order)
  , d_maxlength(W.d_maxlength)
  , d_longest(W.d_longest)
  , d_coxeterMatrix(W.d_coxeterMatrix)
  , d_transducer(W.d_transducer)
  , d_twist() // cannot initialise here
  , d_in()    // cannot initialise here
  , d_out()   // cannot initialise here
  , d_min_star(W.d_min_star)

{
  memcpy(d_in,W.d_in,sizeof(d_in));
  memcpy(d_out,W.d_out,sizeof(d_out));
  memset(d_twist,0,sizeof(d_twist));

  for (size_t s = 0; s < d_rank; ++s) {
    WeylElt w = d_longest;
    prodIn(w,W.d_twist[s]); // no conversion, |d_twist| is already internal
    prod(w,d_longest); // now |w| represents $w_0 twist_s w_0$ as WeylElt

    // |w| should be some generator |t|; find which, and store it
    Generator t = 0;
    for (; t < d_rank; ++t)
      if (w[t]!=0)
	break;
    d_twist[s] = t;
  }
}

void WeylGroup::swap(WeylGroup& other)

{

  std::swap(d_rank,other.d_rank);
  std::swap(d_order,other.d_order);
  std::swap(d_maxlength,other.d_maxlength);
  std::swap(d_longest,other.d_longest);
  d_coxeterMatrix.swap(other.d_coxeterMatrix);
  d_transducer.swap(other.d_transducer);

  // STL swap doesn't work for arrays!

  Generator tmp[constants::RANK_MAX];

  memcpy(tmp,d_in,constants::RANK_MAX);
  memcpy(d_in,other.d_in,constants::RANK_MAX);
  memcpy(other.d_in,tmp,constants::RANK_MAX);

  memcpy(tmp,d_out,constants::RANK_MAX);
  memcpy(d_out,other.d_out,constants::RANK_MAX);
  memcpy(other.d_out,tmp,constants::RANK_MAX);

  memcpy(tmp,d_twist,constants::RANK_MAX);
  memcpy(d_twist,other.d_twist,constants::RANK_MAX);
  memcpy(other.d_twist,tmp,constants::RANK_MAX);

  d_min_star.swap(other.d_min_star);
}

/******** accessors **********************************************************/

WeylElt WeylGroup::generator (Generator i) const
{
  WeylElt s; // initialise to identity
  s[i]=1;    // then shift to piece #1 in $W_{i-1}\\W_i$, which is $s_i$
  return s;
}

int WeylGroup::prodIn(WeylElt& w, Generator s) const

/*!
  \brief Multiply w on the right by the (internally labelled)
  generator s: w *= s.

  This just means exercising the transducer tables as they were designed to.

  Returns +1 if the length moves up, -1 if the length goes down.

  Amazingly, I could simplify Fokko's original to the code below from. I left
  the test for |UndefGenerator| first, because transduce is more frequent than
  shift. Therefore the following even simpler code could be less efficient:

  Generator t=s;
  for (Generator j=d_rank; j-->0; t=d_transducer[j].out(w[j],t))
  {
    EltPiece wj=w[j]; w[j]=d_transducer[j].shift(wj,s);
    if (wj!=w[j]) return w[j]>wj ? 1 : -1;
  }

  MvL
*/

{
  unsigned long j = d_rank-1; // current transducer

  // in the next loop |j| cannot pass |0| since transducer 0 only has shifts
  for (Generator t; (t=d_transducer[j].out(w[j],s))!=UndefGenerator; s=t)
    --j;

  // now transductions are exhausted and one nontrivial shift remains
  EltPiece wj=w[j]; w[j]=d_transducer[j].shift(wj,s); // assert: |w[j]!=wj|

  return w[j]>wj ? 1 : -1; // no need to use d_length, numeric '>' suffices
}

void WeylGroup::prodIn(WeylElt& w, const WeylWord& v) const

/*!
  \brief Multiply w on the right by v (written in internal generators): w *= v.

  Precondition: v is written in the _internal_ generators.

  Algorithm: do the elementary multiplication by the generators, running
  through v left-to-right.
*/

{
  for (size_t j = 0; j < v.size(); ++j)
    prodIn(w,v[j]);
}

void WeylGroup::leftProdIn(WeylElt& w, Generator s) const

/*!
  \brief Transforms w into s*w.

  Algorithm: note that our transducers are geared towards _right_
  multiplication by a generator. But we note that passing from $w$ to $sw$
  only affects the pieces $x_j$ in $w$ for $t <= j <= s$, where
  |t=min_neighbor(s)| is the first generator that does not commute with |s|
  (remarkably, if $v$ is the product of those pieces, $sv$ does have non-zero
  components only for that set of indices; hard to believe at first but easy
  to prove).
*/

{
  WeylElt sw=generator(s);

  // now compute $sv$ as above
  for (size_t j = min_neighbor(s); j <= s; ++j)
    prodIn(sw,wordPiece(w,j));

  // and copy its relevant pieces into $w$
  for (size_t j = min_neighbor(s); j <= s; ++j)
    w[j] = sw[j];
}

void WeylGroup::prod(WeylElt& w, const WeylElt& v) const

/*!
  \brief Multiply w on the right by v, and put the product in w: w*=v.

  Algorithm: increment w by the various pieces of v, whose reduced
  expressions are available from the transducer.
*/

{
  for (size_t j = 0; j < d_rank; ++j)
    prodIn(w,wordPiece(v,j));
}

void WeylGroup::prod(WeylElt& w, const WeylWord& v) const

/*!
  \brief Multiply w on the right by v, and put the product in w: w*=v.

  Algorithm: do the elementary multiplication by the generators, running
  through v left-to-right.
*/

{
  for (size_t j = 0; j < v.size(); ++j)
    prod(w,v[j]);
}

/*! \brief set |w=xw| */
void WeylGroup::leftMult(WeylElt& w, const WeylElt& x) const
{
  WeylElt xx=x; prod(xx,w); w=xx;
}

WeylElt WeylGroup::inverse(const WeylElt& w) const

/*!
  \brief return inverse of |w|

  Algorithm: read backwards the reduced expression gotten from the
  pieces of w.
*/

{
  WeylElt wi;

  for (size_t j = d_rank; j-->0 ;) {
    const WeylWord& x_ww = wordPiece(w,j);
    for (size_t i = x_ww.size(); i-->0;) {
      prodIn(wi,x_ww[i]);
    }
  }

  return wi;
}

void WeylGroup::twistedConjugate
  (TwistedInvolution& tw, const WeylElt& w) const
{
  WeylElt x=w; prod(x,tw.w());

  // now multiply $x$ by $\delta(w^{-1})$
  for (size_t j = d_rank; j-->0 ;) {
    const WeylWord& wj = wordPiece(w,j);
    for (size_t i = wj.size(); i-->0;)
      prodIn(x,d_twist[wj[i]]);
  }
  tw.contents()=x;
}

void WeylGroup::conjugacyClass(WeylEltList& c, const WeylElt& w)
  const

/*!
  \brief Puts in c the conjugacy class of w.

  Algorithm: straightforward enumeration of the connected component of |w| in
  the graph defined by the operation conjugate or twistedConjugate, using a
  |std::set| structure to record previously encountered elements and a
  |std::stack| to store elements whose neighbors have not yet been generated.

*/

{
  std::set<WeylElt> found;
  std::stack<WeylElt> toDo;

  found.insert(w);
  toDo.push(w);

  while (not toDo.empty()) {

    WeylElt v = toDo.top(); toDo.pop();

    for (Generator s = 0; s < rank(); ++s) {
      conjugate(v,s);
      if (found.insert(v).second) // then it was new, put it on to-do list
	toDo.push(v);
    }
  }

  // now convert set |found| to vector
  c.clear(); c.reserve(found.size());
  std::copy(found.begin(),found.end(),back_inserter(c));
}

void WeylGroup::twistedConjugacyClass
  (TwistedInvolutionList& c, const TwistedInvolution& tw)
  const
/*!
  \brief Puts in c the twistes conjugacy class of w.

  Algorithm: straightforward enumeration of the connected component of |w| in
  the graph defined by the operation conjugate or twistedConjugate, using a
  |std::set| structure to record previously encountered elements and a
  |std::stack| to store elements whose neighbors have not yet been generated.

*/

{
  std::set<TwistedInvolution> found;
  std::stack<TwistedInvolution> toDo;

  found.insert(tw);
  toDo.push(tw);

  while (not toDo.empty()) {

    TwistedInvolution v = toDo.top();
    toDo.pop();

    for (Generator s = 0; s < rank(); ++s) {
	twistedConjugate(v,s);
      if (found.insert(v).second)
	toDo.push(v);
    }
  }

  // now convert set |found| to vector
  c.clear(); c.reserve(found.size());
  std::copy(found.begin(),found.end(),back_inserter(c));

}

bool WeylGroup::hasDescent(Generator s, const WeylElt& w) const

/*!
  \brief Tells whether sw < w.

  Method: we multiply from $s$ to $sw$, at least on the relevant pieces (those
  from |min_neighbor(s)| to |s| inclusive); if any decent occurs we return
  |true|, otherwise |false|. Despite the double loop below, this question is
  resolved in very few operations on the average.
*/

{
  s=d_in[s]; // inner numbering is used below

  WeylElt x = generator(s); // becomes (part of) $sw$, unless early |return|

  for (size_t j = min_neighbor(s); j <= s; ++j)
  {
    const WeylWord& piece=wordPiece(w,j);
    for (size_t i=0; i<piece.size(); ++i)
      if (prodIn(x,piece[i])<0) // multiply and see if a descent occurs
	return true;
  }

  return false; // since only ascents occur, we have $l(sw)>l(w)$
}

Generator WeylGroup::leftDescent(const WeylElt& w) const

/*!

  \brief Returns a left descent generator for |w|, or |UndefGenerator| if
  there is no such (i.e., if $w = e$). In fact this is the index |i| of the
  first piece of |w| that is not 0 (identity), converted to external
  numbering, since the canonical (minimal for ShortLex) expression for |w|
  starts with this letter.
*/

{
  for (Generator i=0; i<d_rank; ++i)
    if (w[i]>0) return d_out[i];

  // if we come here, |w==e|
  return UndefGenerator;
}


/*!
  \brief Tells whether |w| twisted-commutes with |s|: $s.w.\delta(s)=w$

  Precondition: |w| is a twisted involution: $w^{-1}=\delta(w)$. Therefore
  twisted commutation is equivalent to $s.w$ being a twisted involution.

  This is in fact the case if and only if $s.w.\delta(s)$ has the same length
  as $w$, by the following reasoning. Suppose first that $s.w$ is reduced,
  then its twisted inverse $w.\delta(s)$ is reduced as well. Then the only
  possible reduction in $s.w.\delta(s)$ is cancellation of the extremal
  generators; whether this reduction applies is equivalent to having twisted
  commutation. If $s.w$ is not reduced, then neither is $w.delta(s)$, and
  $w'=s.w.\delta(s)$ is a twisted involution not longer than $w$. If it is
  strictly shorter then obviously twisted commutation fails. In the remaining
  case that $l(s.w.\delta(s))=l(w)$, let $v=s.w$ so that $w=s.v$ and
  $w'=v.\delta(s)$ are reduced, but $w.\delta(s)=s.v.\delta(s)$ does reduce,
  which can only be by cancelling the extremal generators: $s.v.\delta(s)=v$
  which implies $w'=w$, and one has twisted commutation.
*/
bool WeylGroup::hasTwistedCommutation(Generator s, const TwistedInvolution& tw)
  const
{
  WeylElt x = tw.w();
  int m = prodIn(x,d_twist[d_in[s]]); // now |x| is $w.delta(s)$

  return (m>0)==hasDescent(s,x); // lengths match iff members are equivalent
}

void WeylGroup::involutionOut
  (weyl::WeylWord& ww, const weyl::TwistedInvolution& tw) const

/*!
  \brief Puts in |ww| a reduced expression of |tw| as a twisted involution.

  Precondition: tw is a twisted involution: $tw^{-1}=\delta(tw)$.

  The argument given under |hasTwistedCommutation| shows that for every
  generator |s| exactly one of $s.tw$ and $s.tw.\delta(s)$ is a twisted
  involution distinct from $tw$, and that if $l(s.tw)<l(tw)$ (for the usual
  length function on the Weyl group) then the length of this new twisted
  involution is less than that of $tw$ (by 1 or 2, respectively). A "reduced
  expression as a twisted involution" for $tw$ is obtained by iterating this
  to bring the length down to $0$, and writing down the generators $s$ used in
  reverse order. Working back from the identity, it can be determined for each
  letter to which if the two types of transformation it corresponds. It is not
  immediately obvious that all such reduced expressions have the same length.

  The code below choses the first possible generator (for the internal
  numbering) at each step, so the reduced expression found can be described as
  backwards-lexicographically first for the internal ordering of the generators
*/

{
  TwistedInvolution x = tw;
  ww.clear();

  for (Generator s = leftDescent(x.w()); s != UndefGenerator;
                 s = leftDescent(x.w())) {
    ww.push_back(s);
    if (hasTwistedCommutation(s,x))
      leftMult(x,s);
    else
      twistedConjugate(x,s);
  }

  // reverse ww
  std::reverse(ww.begin(),ww.end());
}

unsigned long WeylGroup::involutionLength
  (const weyl::TwistedInvolution& tw) const

/*!
  \brief Returns the length of tw as a twisted involution.

  Precondition: tw is a twisted involution;

  Algorithm: this is a simplified version of |involutionOut| that records only
  the length
*/

{
  TwistedInvolution x = tw;
  unsigned long length = 0;

  for (Generator s = leftDescent(x.w()); s != UndefGenerator;
                 s = leftDescent(x.w()),++length) {
    if (hasTwistedCommutation(s,x))
      leftMult(x,s);
    else
      twistedConjugate(x,s);
  }

  return length;
}

unsigned long WeylGroup::length(const WeylElt& w) const

/*!
  \brief Returns the length of w.
*/

{
  unsigned long p = 0;

  for (size_t j = 0; j < d_rank; ++j) {
    p += d_transducer[j].length(w[j]);
  }

  return p;
}

void WeylGroup::out(WeylWord& ww, const WeylElt& w) const

/*!
  \brief Puts in ww an _external_ reduced expression for w.
*/

{
  ww.resize(length(w));
  size_t p = 0;

  for (size_t j = 0; j < d_rank; ++j) {
    const WeylWord& xw = wordPiece(w,j);
    for (size_t i = 0; i < xw.size(); ++i) {
      ww[p] = d_out[xw[i]];
      ++p;
    }
  }

  return;
}

void WeylGroup::outerTwist(Twist& t) const

/*!
  \brief Puts the twist of the outer generators in t.

  This reconstructs the twist with which the Weyl group was constructed
*/

{
  for (size_t s = 0; s < d_rank; ++s)
    t[s] = d_out[d_twist[d_in[s]]];

  return;
}

unsigned long WeylGroup::toUlong(const WeylElt& w) const

/*!
  \brief Return the packed form of w. This will work only if the order of
  the group fits into an |unsigned long|, and should therefore not be used.

  This is the mixed-radix interpretation of the sequence of pieces, where
  piece 0 is the least significant one: if $N_i$ is |d_transducer[i].size()|
  and $a_i\in\{0,\ldots,N_i-1\}$ is the value of piece |i|, the value is
  $a_0+N_0*(a_1+N_1*(a_2+... N_{n-2}*(a_{n-1})...))$
*/

{
  unsigned long u = 0;

  for (size_t j; j-->0; )
    u=u*d_transducer[j].size()+w[j];

  return u;
}

WeylElt WeylGroup::toWeylElt(unsigned long u) const

/*!
  \brief Returns the WeylElt whose packed form is u

  Its pieces for the mixed-radix representation of |u|, which as usual can be
  found starting at the least significant end by repeated Euclidean divisions
*/

{
  WeylElt w;

  for (size_t j = 0; j < d_rank; ++j) {
    w[j] = u%d_transducer[j].size();
    u /= d_transducer[j].size();
  }

  return w;
}

void WeylGroup::translate(WeylElt& w, const WeylInterface& I) const

/*!
\brief Applies to w the homomorphism defined by the generator permutation
  in I.

  Algorithm: we used the standard reduced decomposition of w.
*/

{
  WeylElt wt;

  for (size_t j = 0; j < rank(); ++j) {
    const WeylWord& xw = wordPiece(w,j);
    for (size_t i = 0; i < xw.size(); ++i)
      prodIn(wt,I[xw[i]]);
  }

  w = wt;
}

void WeylGroup::twist(WeylElt& w) const

/*!
  \brief Twists the element w by the involutive diagram automorphism |delta|

  Algorithm: we used the standard reduced decomposition of w.
*/

{
  translate(w,d_twist);
}

/*!
  \brief let |w| act on |v| according to reflection action in root datum |rd|
  Note that rightmost factors act first, as in a product of matrices
*/
void WeylGroup::act(const rootdata::RootDatum& rd,
		    const WeylElt& w,
		    latticetypes::LatticeElt& v) const
{
  for (size_t j = d_rank; j-->0; ) {
    const WeylWord& xw = wordPiece(w,j);
    for (size_t i = xw.size(); i-->0; )
      rd.simpleReflection(v,d_out[xw[i]]);
  }
}

/*!
  \brief nondestructive version of |act| method
*/
latticetypes::LatticeElt
WeylGroup::imageBy(const rootdata::RootDatum& rd,
		   const WeylElt& w,
		   const latticetypes::LatticeElt& v) const
{
  latticetypes::LatticeElt result=v; act(rd,w,result); return result;
}

/* One constructor for WeylElt was not defined in header file:
   construct the element from a Weyl word |ww| in a given Weyl group |W|
*/

WeylElt::WeylElt(const WeylWord& ww, const WeylGroup& W)
{
  memset(d_data,0,sizeof(d_data)); // set to identity
  W.prod(*this,ww); // multiply |ww| into current element
}


} // namespace weyl


/*****************************************************************************

        Chapter II -- The Transducer class

  We implement here the construction of the Transducer tables (all accessors
  are defined in the class definition, and there are are no manipulators).
  This is described in section 4 of Fokko's 1999 paper "Transducer approach.."
  Actually, that paper does almost everything by induction, in particular it
  makes a (somewhat vague) reference to using previously constructed
  transducer tables while bootstrapping the current one; this is not what is
  done here, which proceeds strictly independently of other Transducer tables.
  I think the mention of dihedral groups below replaces the inductive part. MvL

******************************************************************************/

namespace weyl {

Transducer::Transducer(const latticetypes::LatticeMatrix& c, size_t r)
  : d_shift(1), d_out(1) // start with tables of size 1
  , d_length(1,0), d_piece(1,WeylWord()) // with empty word, length 0.

/*!
\brief Constructs subquotient \#r for the Coxeter matrix c.

  This uses the Coxeter matrix only up to index r. In fact we can behave as
  if generator |r| is the final one, since we ignore any higher ones for now.

  Precondition : c is a _normalized_ Coxeter matrix (meaning that all
  the parabolic subquotients W_{r-1}\\W_r are small enough to fit in
  an unsigned char); and r is < rank(c);

  Algorithm :

  The algorithm is a version of my favorite bootstrapping procedure for the
  construction of Weyl groups and parabolic quotients. We start with a
  partially defined automaton containing only one element, and for which all
  shifts by the final generator $r$ are not yet defined (generators $i<r$ give
  transduction of $i$). At each point in time, all shifts for all elements in
  the automaton that do _not_ take the length up are defined; and we maintain
  a queue of elements that may have as yet undefined shifts.

  We start up with one element in the automaton, with just one undefined
  shift, the one by r. Then run through the elements $x$ of the automaton in
  order of generation (which will also be in ShortLex order), and for each as
  yet undefined shift of $x$ by $s$ :

    - add a new element $xs$ (this increases the size of the automaton, during
      the loop!). At this point it is sure that this element is really new,
      and that its normal form is obtained by adding $s$ to that of $x$

    - then find all other elements $x'$ already in the automaton and $t$ such
      that $xs==x't$ (so $x't$ gives a non-canonical but reduced expression
      for $xs$). Do this by trying generators $t\neq s$: if $xst$ goes down
      (has the same length as $x$) then $x'==xst$ gives such a case. The trick
      for this is to look at the orbit of x under the dihedral group <s,t>. In
      the full group, this has necessarily cardinality 2m, with m = m(s,t) the
      coefficient in the Coxeter matrix, and $l(xst)==l(x)$ iff $xs$ is the
      unique elt. of maximal length in the orbit, hence to have this $x$ must
      goes down $m-1$ times when applying successively $t$, $s$, $t$, ... In
      the parabolic quotient, the orbit of the dihedral group (which is not
      reduced to a point) can either have cardinality $2m$ or cardinality $m$,
      and in the latter case it is a string with $m-1$ steps between the
      bottom and the top, with a stationary step at either extreme (to see
      this, note that on one hand each step up in the full group gives a step
      in the quotient that is either up or stationary, while on the other hand
      a stationary step in the quotient causes then next step to be the
      reverse of the previous one). So we have one of the following cases: (1)
      $x$ goes down $m-1$ times; then the image of the orbit has $2m$ elements
      and $xs==x't$ for $x'=x.(ts)^(m-1)$. (2) $x$ goes down $m-2$ times to
      some element $a$ and is then stationary (if $s'\in\{s,t\}$ is the next
      to apply, then $a.s'=g.a$ for some generator $g\in W_{r-1}$); then the
      orbit has $m$ elements, and if $v$ is the alternating word in $\{s,t\}$
      of length $m-2$ not starting with $s'$, so that $a.v=x$, one has
      $v.st=s'vs$ whence $x.st=a.v.st=a.s'vs=g.a.vs=g.xs$ so that $xs$ has a
      transduction for $t$ that outputs the generator $g$. (3) either $x$ goes
      down less than $m-2$ times, or $m-2$ times followed by an upward step;
      then $xst$ goes up.
*/

{
  // first first row of transition and of transduction table

  ShiftRow& firstShift=d_shift[0];
  OutRow& firstOut=d_out[0];

  // all shifts lower than r should be defined
  for (size_t j = 0; j < r; ++j) {
    firstShift[j] = 0; // shift to (current) state 0, i.e., no transition
    firstOut[j] = j;   // but transduction of the same generator
  }
  // note that firstShift[r]==UndefValue==UndefEltPiece and
  // firstOut[r]==UndefValue==UndefGenerator remain from default construction

  // in this loop, the table grows! the loop stops when x overtakes the
  // table size because no more new elements are created.

  for (EltPiece x = 0; x < d_shift.size(); ++x) {
    for (Generator s = 0; s <= r; ++s)
      /* since RANK_MAX<128, |UndefEltPiece| is never a valid Piece number, so
         its presencein a slot in |d_shift| assures that this slot is
         unchanged from its intialisation value
      */
      if (d_shift[x][s] == UndefEltPiece) {

	EltPiece xs = d_shift.size(); // index of state that will be added

	d_shift.push_back(ShiftRow()); // push a row of UndefValue both onto
	d_out.push_back(OutRow());     // |d_shift| and onto |d_out|

	d_shift[x][s] = xs;
	d_shift[xs][s] = x;
	// new element has length length(x)+1
	d_length.push_back(d_length[x]+1);

	WeylWord xs_ww = d_piece[x]; // get normal form for |x|
	xs_ww.push_back(s);          // and append a letter |s|
	d_piece.push_back(xs_ww);    // which gives normal form for |xs|


	// now define the shifts or transductions that do not take xs up

	for (Generator t = 0; t <= r; ++t) {

	  if (t == s)
	    continue;

	  LatticeCoeff m = c(s,t); // coefficient of the Coxeter matrix

	  EltPiece y  = dihedralMin(*this,xs,s,t);
	  unsigned long d = d_length[xs] - d_length[y];

	  if (d < static_cast<unsigned long>(m-1))
	    continue;  // case (3):  $t$ takes $xs$ up, do nothing

	  Generator st[] = {s,t};

	  if (d == static_cast<unsigned long>(m)) {
	    // case (1): there is no transduction
	    // xs.t is computed by shifting up from y the other way around
	    y = dihedralShift(*this,y,st[m%2],st[(m+1)%2],m-1);
	    d_shift[xs][t] = y;
	    d_shift[y][t] = xs;
	    continue;  // next t
	  }

	  // now d == m-1;

	  Generator u = st[(m+1)%2];

	  if (d_shift[y][u] == y) {
	    // case (2): $xs$ is fixed by $t$, outputs the same $g$ as $y$
	    d_shift[xs][t] = xs;
	    d_out[xs][t] = d_out[y][u];
	  }

	  // else case (3) : $xs$ moves up, do nothing
	}
      }
  }
}

}

/*****************************************************************************

        Chapter III -- Functions declared in weyl.h

******************************************************************************/

namespace weyl {

void copy(Twist& dest, const Twist& source)

/*!
  \brief Copies source onto dest.
*/

{
  memcpy(dest,source,constants::RANK_MAX);
  return;
}

const Twist& identityTwist()

/*!
  \brief Returns a constant reference to a trivial twist.
*/

{
  static IdentityTwist t;

  return t.twist();
}

void makeReflections(WeylEltList& refl, const WeylGroup& W)

/*!
  \brief Puts in |refl| the list of all reflections (conjugates of generators)
         in W.

  NOTE: the ordering of the reflections is the ordering induced by our
  operator<, which is not very significative mathematically, but has the
  advantage that the STL search tools may be used.
*/

{
  WeylEltList simple(W.rank());

  // put in simple the set of simple reflections
  for (size_t j = 0; j < simple.size(); ++j)
    simple[j][j] = 1;

  std::set<WeylElt> found;

  for (size_t j = 0; j < simple.size(); ++j)
    // make new class if the reflection is not already found
    if (found.insert(simple[j]).second) {
      WeylEltList c;
      W.conjugacyClass(c,simple[j]);
      for (size_t i = 0; i < c.size(); ++i)
	found.insert(c[i]);
    }

  refl.clear();
  copy(found.begin(),found.end(),back_inserter(refl));
}

}

/*****************************************************************************

        Chapter IV -- Auxiliary functions for this module

******************************************************************************/

namespace {

EltPiece dihedralMin(const Transducer& qa, EltPiece x, Generator s,
		       Generator t)

/*!
  \brief Returns the minimal element in the orbit of x under s and t.

  Precondition : s is in the descent set of x;
*/

{
  Generator u = s;
  Generator v = t;

  EltPiece y = x;

  for (;;) {
    // this is ok even if the shift is still undefined
    if (qa.shift(y,u) >= y)
      return y;
    else
      y = qa.shift(y,u);
    std::swap(u,v);
  }
}

EltPiece dihedralShift(const Transducer& qa, EltPiece x, Generator s,
		       Generator t, unsigned long d)

/*!
  \brief Returns the result of applying s and t alternately to x, for a
  total of d times.
*/

{
  Generator u = s;
  Generator v = t;

  EltPiece y = x;

  for (unsigned long j = 0; j < d; ++j) {
    y = qa.shift(y,u);
    std::swap(u,v);
  }

  return y;
}

void fillCoxMatrix(LatticeMatrix& cox, const LatticeMatrix& cart,
		   const Permutation& a)

/*!
  \brief Fills in the Coxeter matrix cox.

  Precondition: cart is a Cartan matrix; a holds a normalizing permutation
  for cart, such as constructed by normalize(a,d) where d is the Dynkin diagram
  of cart (declared in dynkin.h);

  Postcondition : cox holds the normalized Coxeter matrix corresponding to
  cox and a;
*/

{
  cox.resize(cart.numColumns(),cart.numColumns(),2);

  // fill in the appropriate values

  for (size_t j = 0; j < cox.numColumns(); ++j)
    cox(j,j) = 1;

  for (size_t j = 0; j < cart.numColumns(); ++j)
    for (size_t i = j+1; i < cart.numRows(); ++i)
      if (cart(i,j)!=0) {
	if (cart(i,j)*cart(j,i) == 1) {
	  cox(i,j) = 3;
	  cox(j,i) = 3;
	}
	else if (cart(i,j)*cart(j,i) == 2) {
	  cox(i,j) = 4;
	  cox(j,i) = 4;
	}
	else if (cart(i,j)*cart(j,i) == 3) {
	  cox(i,j) = 6;
	  cox(j,i) = 6;
	}
      }

  // permute

  cox.permute(a);
}

} // namespace

} // namespace atlas
