/*!
\file
\brief Implementation of WeylGroup.

  I have decided to represent elements as fixed-size arrays of
  unsigned characters. This forces expressing things in the standard
  ordering of the generators, and hence to have a small I/O interface
  for resetting the numbering to and from the numbering used by the
  outside world.

  It has seemed to me that this is the best compromise between size of the
  dataype, generality and efficiency.
  [Fokko du Cloux]
*/
/*
  This is weyl.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007--2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "weyl.h"

#include <algorithm>
#include <set>
#include <stack>

#include "ratvec.h"	// to act upon |RatWeight|s
#include "dynkin.h"	// to analyze Cartan matrices
#include "permutations.h"// to hold the result from dynkin
#include "prerootdata.h"// for defining action using only simple (co)roots
#include "rootdata.h"	// also needed for defining action, and deducing twist
#include "lietype.h"  // for |ext_gen|
#include "sl_list.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include <iterator>
#endif

/*****************************************************************************

  Implementation of Weyl groups. I have decided to represent elements as
  fixed-size arrays of unsigned characters. This forces to express things
  in the standard ordering of the generators, and hence to have a small
  i/o interface for resetting the numbering to and from the outside world's.

  It has seemed to me that this is the best compromise between size of the
  dataype, generality and efficiency.

******************************************************************************/

namespace atlas {

  namespace weyl { // constants needed only in this file

  const WeylElt::EltPiece UndefEltPiece = UndefValue;
  const Generator UndefGenerator = UndefValue;

}

//			     auxiliary functions
namespace {


  void fillCoxMatrix(int_Matrix&,
		     const int_Matrix&,
		     const Permutation&);

  WeylElt::EltPiece dihedralMin(const weyl::Transducer&,
				      WeylElt::EltPiece,
				      weyl::Generator,
				      weyl::Generator);

  WeylElt::EltPiece dihedralShift(const weyl::Transducer&,
					WeylElt::EltPiece,
					weyl::Generator,
					weyl::Generator,
					unsigned long);

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
  and incorporate the possibly given twist (take identity if |twist==NULL|).

  NOTE : |c| and |twist| use some (consistent) labelling of simple roots,
  but we will determine an internal renumbering making the subquotients small
*/
WeylGroup::WeylGroup(const int_Matrix& c)
  : d_rank(c.numColumns()) // all other fields are computed later
  , d_order() // this initialises this size::Size value to 1
  , d_maxlength(0)
  , d_longest()
  , d_coxeterMatrix()
  , Chevalley()
  , d_transducer(d_rank)
  , d_in()
  , d_out()
  , d_min_star(d_rank)

{
  /* analyse the Cartan matrix */

  DynkinDiagram d(c); // make diagram from Cartan matrix
  // find renumbering |a| putting labels in canonical order
  Permutation a= dynkin::normalize(d);

  /* now put appropriate permutations into |d_in| and |d_out|, so that
     internal number |j| gives external number |d_out[j]|, and of course
     |d_in[d_out[j]]==j|. By convention of |normalize|, this means |d_out==a|
   */
  for (size_t j = 0; j < d_rank; ++j) {
    d_out[j] = a[j];
    d_in[d_out[j]] = j;
  }

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
  for (size_t j = 0; j < d_rank; ++j)
    d_order *= d_transducer[j].size();

  // precompute for each |j| the first non-commuting or equal generator |i<=j|
  for (size_t j = 0; j < d_rank; ++j)
    for (size_t i=0; i<=j; ++i)
      if (d_coxeterMatrix(i,j)!=2)
      {
	d_min_star[j]=i; break;
      }

  // now our Weyl group is operational for computations

  // precompute the Chevalley involution
  // since we dont have a root datum at hand, we conjugate by |longest()|
  for (Generator s=0; s<rank(); ++s)
  {
    WeylElt w = longest();
    mult(w,s);
    mult(w,longest());

    // |w| should be some generator |t|; find which, and store it
    for (Generator t=0; t<rank(); ++t)
      if (w==generator(t))
      {	Chevalley[s] = t; break; }
  }

} // |WeylGroup::WeylGroup|


/******** accessors **********************************************************/

WeylElt WeylGroup::genIn (Generator i) const
{
  WeylElt s; // initialise to identity
  s[i]=1;    // then shift to piece #1 in $W_{i-1}\\W_i$, which is $s_i$
  return s;
}

/*!
  \brief Multiply w on the right by internally numbered generator s: w *= s.

  Returns +1 if the length moves up, -1 if the length goes down.

  This just means exercising the transducer tables as they were designed to.
*/

/*
  Amazingly, I could simplify Fokko's original to the code below.
  Fokko's original code was (more or less):

  Generator t = s;
  for (unsigned long j = d_rank; j-->0;) {
    if (d_transducer[j].out(w[j],t) == UndefGenerator) {
      // no transduction; shift and terminate
      unsigned long l = d_transducer[j].length(w[j]);
      w[j] = d_transducer[j].shift(w[j],t);
      if (d_transducer[j].length(w[j]) > l)
	return 1;
      else
	return -1;
    }
    else {
      // no shift; transduce and continue
      t = d_transducer[j].out(w[j],t);
    }
  }

  The main simplifications are to avoid calling |d_transducer[j].out| twice in
  succession, and to avoid calling |d_transducer[j].length| at all (since
  comparison of the |EltPiece| values suffices to decide sense of change).

  I left the test for |UndefGenerator| coming first, as transduce is more
  frequent than shift; hence the following even simpler code could be less
  efficient:

  Generator t=s;
  for (Generator j=d_rank; j-->0; t=d_transducer[j].out(w[j],t))
  {
    EltPiece wj=w[j]; w[j]=d_transducer[j].shift(wj,s);
    if (wj!=w[j]) return w[j]>wj ? 1 : -1;
  }

  MvL
*/
int WeylGroup::multIn(WeylElt& w, Generator s) const
{
  unsigned long j = d_rank-1; // current transducer

  // in the next loop |j| cannot pass |0| since transducer 0 only has shifts
  for (Generator t; (t=d_transducer[j].out(w[j],s))!=UndefGenerator; s=t)
  {
    assert(j!=0);
    --j;
  }

  // now transductions are exhausted and one nontrivial shift remains
  WeylElt::EltPiece wj=w[j];
  w[j]=d_transducer[j].shift(wj,s); // assert: |w[j]!=wj|

  return w[j]>wj ? 1 : -1; // no need to use d_length, numeric '>' suffices
}

/*!\brief Multiply w on the right by v (in internal numbering): w *= v.

  Returns nonpositive even value $l(wv)-l(w)-l(v)$

  Precondition: v is written in the _internal_ generators.

  Algorithm: do the elementary multiplication by the generators, running
  through v left-to-right.
*/
int WeylGroup::multIn(WeylElt& w, const WeylWord& v) const

{
  int result=0;
  for (size_t j = 0; j < v.size(); ++j)
    if (multIn(w,v[j])<0) result-=2;

  return result;
}

/*!
  \brief Transforms w into s*w, with s in internal numbering;
  returns $l(sw)-l(w)\in\{+1,-1}$

  Algorithm: note that our transducers are geared towards _right_
  multiplication by a generator. But we note that passing from $w$ to $sw$
  only affects the pieces $x_j$ in $w$ for $t <= j <= s$, where
  |t=min_neighbor(s)| is the first generator that does not commute with |s|
  (remarkably, if $v$ is the product of those pieces, $sv$ does have non-zero
  components only for that set of indices; hard to believe at first but easy
  to prove).
*/
int WeylGroup::leftMultIn(WeylElt& w, Generator s) const
{
  WeylElt sw=genIn(s);
  int l=1; //

  // now compute $sv$ as above, keeping track of any length drop (at most 1)
  for (size_t j = min_neighbor(s); j <= s; ++j)
    l+=multIn(sw,wordPiece(w,j));

  // and copy its relevant pieces into $w$
  for (size_t j = min_neighbor(s); j <= s; ++j)
    w[j] = sw[j];

  return l;
}

/*!
  \brief Multiply w on the right by v, and put the product in w: w*=v.

  Algorithm: increment w by the various pieces of v, whose reduced
  expressions are available from the transducer.
*/
void WeylGroup::mult(WeylElt& w, const WeylElt& v) const
{
  for (size_t j = 0; j < d_rank; ++j)
    multIn(w,wordPiece(v,j));
}

/*!
  \brief Multiply w on the right by v, and put the product in w: w*=v.

  Algorithm: do the elementary multiplication by the generators, running
  through v left-to-right.
*/
void WeylGroup::mult(WeylElt& w, const WeylWord& v) const
{
  for (size_t j = 0; j < v.size(); ++j)
    mult(w,v[j]);
}

/*! \brief set |w=xw| */
void WeylGroup::leftMult(WeylElt& w, const WeylElt& x) const
{
  WeylElt xx=x; mult(xx,w); w=xx;
}

/*!
  \brief return inverse of |w|

  Algorithm: read backwards the reduced expression gotten from the
  pieces of w.
*/
WeylElt WeylGroup::inverse(const WeylElt& w) const
{
  WeylElt wi;

  for (size_t j = d_rank; j-->0 ;) {
    const WeylWord& x_ww = wordPiece(w,j);
    for (size_t i = x_ww.size(); i-->0;) {
      multIn(wi,x_ww[i]);
    }
  }

  return wi;
}

/*!
  \brief Puts into |c| the conjugacy class of |w|.

  Algorithm: straightforward enumeration of the connected component of |w| in
  the graph defined by the operation conjugatee, using a |std::set| structure
  to record previously encountered elements and a |std::stack| to store
  elements whose neighbors have not yet been generated.

*/
void WeylGroup::conjugacyClass(WeylEltList& c, const WeylElt& w) const
{
  std::set<WeylElt> found;
  std::stack<WeylElt,containers::mirrored_simple_list<WeylElt> > toDo;

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

/*!
  \brief Tells whether sw < w.

  Method: we multiply from $s$ to $sw$, at least by the word pieces of |w| at
  the relevant pieces: those from |min_neighbor(s)| to |s| inclusive. If any
  descent occurs we return |true|, otherwise |false|. Despite the double loop
  below, this question is resolved in relatively few operations on the average.
*/
bool WeylGroup::hasDescent(Generator s, const WeylElt& w) const
{
  s=d_in[s]; // inner numbering is used below

  WeylElt x = genIn(s); // element operated upon, starts out as |s|

  for (size_t j = min_neighbor(s); j <= s; ++j)
  {
    const WeylWord& piece=wordPiece(w,j);
    for (size_t i=0; i<piece.size(); ++i)
      if (multIn(x,piece[i])<0) // multiply and see if a descent occurs
	return true;
  }

  return false; // since only ascents occur, we have $l(sw)>l(w)$
}

// same question, but on the right
bool WeylGroup::hasDescent(const WeylElt& w,Generator s) const
{
  s=d_in[s]; // inner numbering is used below
  unsigned int j = d_rank-1; // current transducer

  // in the next loop |j| cannot pass |0| since transducer 0 only has shifts
  for (Generator t; (t=d_transducer[j].out(w[j],s))!=UndefGenerator; s=t)
    --j;

  WeylElt::EltPiece wj=w[j];
  return d_transducer[j].shift(wj,s)<wj;
}


/*!

  \brief Returns a left descent generator for |w|, or |UndefGenerator| if
  there is no such (i.e., if $w = e$). In fact this is the index |i| of the
  first piece of |w| that is not 0 (identity), converted to external
  numbering, since the canonical (minimal for ShortLex) expression for |w|
  starts with this letter. Returning the first generator in external numbering
  would take a bit more time, as we would have to repeatedly use |hasDescent|.
*/
Generator WeylGroup::leftDescent(const WeylElt& w) const
{
  for (Generator i=0; i<d_rank; ++i)
    if (w[i]>0) return d_out[i];

  // if we come here, |w==e|
  return UndefGenerator;
}


/*!
  \brief Returns the length of w.

  This is relatively efficient (compared to |involutionLength|)
*/
unsigned int WeylGroup::length(const WeylElt& w) const
{
  unsigned long l = 0;

  for (size_t i = 0; i < d_rank; ++i)
    l += d_transducer[i].length(w[i]);

  return l;
}

Generator WeylGroup::letter(const WeylElt& w, unsigned int k) const
{
  for (size_t i = 0; i<d_rank; ++i)
    if (k>=d_transducer[i].length(w[i]))
      k -= d_transducer[i].length(w[i]);
    else
      return d_out[wordPiece(w,i)[k]];
  assert(false);
  return Generator(~0);
}

WeylWord WeylGroup::word(const WeylElt& w) const
{
  WeylWord result; result.reserve(length(w));
  for (size_t j = 0; j < d_rank; ++j)
  {
    const WeylWord& xw = wordPiece(w,j);
    for (size_t i = 0; i < xw.size(); ++i)
      result.push_back(d_out[xw[i]]);
  }

  return result;
}

/*!
  \brief Returns the list of all reflections (conjugates of generators).

  NOTE: the ordering of the reflections is the ordering induced by our
  operator<, which is not very significative mathematically, but has the
  advantage that the STL search tools may be used.
*/
WeylEltList WeylGroup::reflections() const
{
  WeylEltList simple; simple.reserve(rank());

  // put in simple the set of simple reflections, along internal numbering
  for (size_t j = 0; j < rank(); ++j)
    simple.push_back(genIn(j));

  std::set<WeylElt> found;

  for (size_t j = 0; j < simple.size(); ++j)
    if (found.insert(simple[j]).second) // then generator itself still absent
    { // generate conjugacy class of generator and insert its elements
      WeylEltList c; conjugacyClass(c,simple[j]);
      found.insert(c.begin(),c.end()); // add new class to result set
    }

  return WeylEltList(found.begin(),found.end()); // convert set to vector
}

/*!
  \brief Return the packed form of w. This will work only if the order of
  the group fits into an |unsigned long|, and should therefore not be used.

  This is the mixed-radix interpretation of the sequence of pieces, where
  piece 0 is the least significant one: if $N_i$ is |d_transducer[i].size()|
  and $a_i\in\{0,\ldots,N_i-1\}$ is the value of piece |i|, the value is
  $a_0+N_0*(a_1+N_1*(a_2+... N_{n-2}*(a_{n-1})...))$
*/
unsigned long WeylGroup::toUlong(const WeylElt& w) const
{
  unsigned long u = 0;

  for (size_t j=d_rank; j-->0; )
    u=u*d_transducer[j].size()+w[j];

  return u;
}

/*!
  \brief Returns the WeylElt whose packed form is u

  Its pieces for the mixed-radix representation of |u|, which as usual can be
  found starting at the least significant end by repeated Euclidean divisions
*/
WeylElt WeylGroup::toWeylElt(unsigned long u) const
{
  WeylElt w;

  for (size_t j = 0; j < d_rank; ++j) {
    w[j] = u%d_transducer[j].size();
    u /= d_transducer[j].size();
  }

  return w;
}

/*!
\brief Applies to |w| the generator permutation in |I|, which should be an
  automorphism of the Dynkin diagram, expressed in terms of outer numbering.

  Algorithm: we use the standard reduced decomposition of |w|, and rebuild the
  return value by repeated right-multiplication. We can do the multiplication
  in the same Weyl group as the decomposition because |I| is supposed to be an
  automorphism; if it weren't we would need a reference to a second Weyl group.
*/
WeylElt WeylGroup::translation(const WeylElt& w, const WeylInterface& f) const
{
  WeylElt result;

  for (size_t j = 0; j < rank(); ++j) {
    const WeylWord& xw = wordPiece(w,j);
    for (size_t i = 0; i < xw.size(); ++i)
      mult(result,f[d_out[xw[i]]]);
  }

  return result;
}

/*
  Let |w| act on |alpha| according to reflection action in root datum |rd|
  Note that rightmost factors act first, as in a product of matrices
*/
void WeylGroup::act(const RootDatum& rd, const WeylElt& w, RootNbr& alpha) const
{
  for (size_t i = d_rank; i-->0; )
    alpha = rd.permuted_root(wordPiece(w,i),alpha);
}

// Let |w| act on |v| according to reflection action in root datum |rd|
template<typename C>
  void WeylGroup::act
    (const RootDatum& rd, const WeylElt& w,  matrix::Vector<C>& v) const
{
  for (size_t i = d_rank; i-->0; )
  {
    const WeylWord& xw = wordPiece(w,i);
    for (size_t j = xw.size(); j-->0; )
      rd.simple_reflect(d_out[xw[j]],v);
  }
}

void WeylGroup::act(const RootDatum& rd, const WeylElt& w, RatWeight& v) const
{ act(rd,w,v.numerator()); }

void WeylGroup::act(const RootDatum& rd, const WeylElt& w, LatticeMatrix& M)
  const
{
  for (size_t i = d_rank; i-->0; )
  {
    const WeylWord& xw = wordPiece(w,i);
    for (size_t j = xw.size(); j-->0; )
      rd.simple_reflect(d_out[xw[j]],M);
  }
}


template<typename C>
  void WeylGroup::act
  (const PreRootDatum& prd, const WeylElt& w, matrix::Vector<C>& v) const
{
  for (size_t i = d_rank; i-->0; )
  {
    const WeylWord& xw = wordPiece(w,i);
    for (size_t j = xw.size(); j-->0; )
      prd.simple_reflect(d_out[xw[j]],v);
  }
}

void WeylGroup::act(const PreRootDatum& prd, const WeylElt& w, RatWeight& v)
  const
{ act(prd,w,v.numerator()); }

void WeylGroup::act(const PreRootDatum& prd, const WeylElt& w, LatticeMatrix& M)
  const
{
  for (size_t i = d_rank; i-->0; )
  {
    const WeylWord& xw = wordPiece(w,i);
    for (size_t j = xw.size(); j-->0; )
      prd.simple_reflect(d_out[xw[j]],M);
  }
}

/*!
  \brief
  Same as |act(rd,inverse(w),v)|, but avoiding computation of |inverse(w)|.
  Here the leftmost factors act first.
*/
void WeylGroup::inverse_act(const RootDatum& rd, const WeylElt& w, Weight& v)
  const
{
  for (size_t i=0; i<d_rank; ++i )
  {
    const WeylWord& xw = wordPiece(w,i);
    for (size_t j=0; j<xw.size(); ++j )
      rd.simple_reflect(d_out[xw[j]],v);
  }
}

/* One constructor for WeylElt was not defined in header file:
   construct the element from a Weyl word |ww| in a given Weyl group |W|
*/
WeylElt::WeylElt(const WeylWord& ww, const WeylGroup& W)
{
  std::memset(d_data,0,sizeof(d_data)); // set to identity
  W.mult(*this,ww); // multiply |ww| into current element
}


/*
               TwistedWeylGroup implementation
*/

TwistedWeylGroup::TwistedWeylGroup
  (const WeylGroup& d_W, const Twist& twist)
  : W(d_W)
  , d_twist(twist)
{
}

WeylElt TwistedWeylGroup::dual_twisted(const WeylElt& w) const
{
  Twist dual_twist;
  for (Generator s=0; s<rank(); ++s)
    dual_twist[s] = W.Chevalley_dual(d_twist[s]);
  return W.translation(w,dual_twist);
}

ext_gens TwistedWeylGroup::twist_orbits ()  const
{
  unsigned int size=0;
  for (weyl::Generator s=0; s<rank(); ++s)
    if (twisted(s)>=s)
      ++size;

  ext_gens result; result.reserve(size);

  for (weyl::Generator s=0; s<rank(); ++s)
    if (twisted(s)==s)
      result.push_back(ext_gen(s));
    else if (twisted(s)>s)
      result.push_back(ext_gen(W.commutes(s,twisted(s)),s,twisted(s)));

  return result;
} // |twist_orbits|


Twist TwistedWeylGroup::dual_twist() const
{
  Twist twist; // "dimensioned" but not initialised
  for (Generator s=0; s<W.rank(); ++s)
    twist[s] = W.Chevalley_dual(twisted(s));
  return twist;
}

/* the "dual" Weyl group: the only difference with W is that the twist is
   multiplied by conjugation with the longest element.
*/
TwistedWeylGroup::TwistedWeylGroup(const TwistedWeylGroup& tW, tags::DualTag)
  : W(tW.W)
  , d_twist(tW.dual_twist())
{}

void TwistedWeylGroup::twistedConjugate // $tw = w.tw.twist(w)^{-1}$
  (TwistedInvolution& tw, const WeylElt& w) const
{
  WeylElt x=w; W.mult(x,tw.w());

  // now multiply $x$ by $\delta(w^{-1})$
  for (size_t i = W.length(w); i-->0 ;)
    mult(x,twisted(W.letter(w,i)));

  tw.contents()=x;
}

/*!
  \brief Puts in c the twistes conjugacy class of w.

  Algorithm: straightforward enumeration of the connected component of |w| in
  the graph defined by the operation twistedConjugate, using a
  |std::set| structure to record previously encountered elements and a
  |std::stack| to store elements whose neighbors have not yet been generated.

*/
void TwistedWeylGroup::twistedConjugacyClass
  (TwistedInvolutionList& c, const TwistedInvolution& tw)
  const
{
  std::set<TwistedInvolution> found;
  std::stack<TwistedInvolution,
	     containers::mirrored_simple_list<TwistedInvolution> > toDo;

  found.insert(tw);
  toDo.push(tw);

  while (not toDo.empty()) {

    TwistedInvolution v = toDo.top();
    toDo.pop();

    for (Generator s=0; s<W.rank(); ++s)
    {
	twistedConjugate(v,s);
      if (found.insert(v).second)
	toDo.push(v);
    }
  }

  // now convert set |found| to vector
  c.clear(); c.reserve(found.size());
  std::copy(found.begin(),found.end(),back_inserter(c));

}


/*!
  \brief Tells whether |w| twisted-commutes with |s|: $s.w.\delta(s)=w$

  Precondition: |w| is a twisted involution: $w^{-1}=\delta(w)$. Therefore
  twisted commutation is equivalent to $s.w$ being a twisted involution.

  This is in fact the case if and only if $s.w.\delta(s)$ has the same
  length as $w$, by the following reasoning. Suppose first that $s.w$ is
  reduced, then its twisted inverse $w.\delta(s)$ is reduced as well. Then
  the only possible reduction in $s.w.\delta(s)$ is cancellation of the
  extremal generators; whether this reduction applies is equivalent to having
  twisted commutation. If $s.w$ is not reduced, then neither is
  $w.\delta(s)$, and $w'=s.w.\delta(s)$ is a twisted involution not
  longer than $w$. If it is strictly shorter then obviously twisted
  commutation fails. In the remaining case that $l(s.w.\delta(s))=l(w)$,
  let $v=s.w$ so that $w=s.v$ and $w'=v.\delta(s)$ are reduced, but
  $w.\delta(s)=s.v.\delta(s)$ does reduce, which can only be by cancelling
  the extremal generators: $s.v.\delta(s)=v$ which implies $w'=w$, and one
  has twisted commutation.
*/
bool TwistedWeylGroup::hasTwistedCommutation
  (Generator s, const TwistedInvolution& tw) const
{
  WeylElt x = tw.w();
  int m = mult(x,d_twist[s]); // now |x| is $w.\delta(s)$

  return (m>0)==W.hasDescent(s,x); // lengths match iff members are equivalent
}


/*!
  Precondition: tw is a twisted involution: $tw^{-1}=\delta(tw)$.

  The argument given under |hasTwistedCommutation| shows that for every
  generator |s| exactly one of $s.tw$ and $s.tw.\delta(s)$ is a twisted
  involution distinct from $tw$, and that if $l(s.tw)<l(tw)$ (for the usual
  length function on the Weyl group) then the length of this new twisted
  involution is less than that of $tw$ (by 1 or 2, respectively). A "reduced
  expression as a twisted involution" for $tw$ is obtained by iterating this
  to bring the length down to $0$. Working back from the identity (so reading
  our expression from right to left), it can be determined for each letter to
  which if the two types of transformation it corresponds; nevertheless, we
  encode which case prevails in the sign of the generator recorded: we
  bitwise-complement for the case of conjugation.

  Although not immediately obvious, all such reduced expressions do have the
  same length.

  The code below chooses the first possible generator (for the internal
  numbering, as returned by |leftDescent|) at each step, so the reduced
  expression found is lexicographically first for the internal renumbering.
*/
InvolutionWord TwistedWeylGroup::involution_expr(TwistedInvolution tw) const
{
  InvolutionWord result; result.reserve(involutionLength(tw));

  for (Generator s = W.leftDescent(tw.w()); s != UndefGenerator;
                 s = W.leftDescent(tw.w()))
    if (hasTwistedCommutation(s,tw))
    {
      result.push_back(s);
      leftMult(tw,s);
    }
    else
    {
      result.push_back(~s);
      twistedConjugate(tw,s);
    }

  return result;
}

// This one trades some efficiency for assurance of external least lex repr
InvolutionWord TwistedWeylGroup::canonical_involution_expr(TwistedInvolution tw)
  const
{
  InvolutionWord result; result.reserve(involutionLength(tw));

  TwistedInvolution delta; // distinguished
  while (tw!=delta)
  {
    Generator s=0;
    while (not hasDescent(s,tw))
      ++s;
    // now |s| is the first twisted right equivalently left descent

    if (hasTwistedCommutation(s,tw))
    {
      result.push_back(s);
      leftMult(tw,s);
    }
    else
    {
      result.push_back(~s);
      twistedConjugate(tw,s);
    }

  } // while(tw!=delta)
  return result;
} // |canonical_involution_expr|

// expression as series of extended-ascents from the trivial involution
// precondition: |TwistedWeylGroup::twisted(tw)==tw|
InvolutionWord TwistedWeylGroup::extended_involution_expr(TwistedInvolution tw)
  const
{
  assert(twisted(tw)==tw);
  ext_gens orbit = twist_orbits();

  InvolutionWord result; result.reserve(involutionLength(tw));

  TwistedInvolution delta; // distinguished
  while (tw!=delta)
  {
    Generator s=0;
    while (not hasDescent(orbit[s].s0,tw))
    {
      ++s;
      assert(s<orbit.size());
    }
    // now |s| is the first right descent

    switch(orbit[s].type)
    {
    case ext_gen::one:
      if (hasTwistedCommutation(orbit[s].s0,tw))
      {
	leftMult(tw,orbit[s].s0); // real
	result.push_back(s);
      }
      else
      {
	twistedConjugate(tw,orbit[s].s0); // complex descent
	result.push_back(~s);
      }
      break;
    case ext_gen::two:
      if (hasTwistedCommutation(orbit[s].s0,tw))
      {
	result.push_back(s); // a double real descent
	leftMult(tw,orbit[s].s0);
	leftMult(tw,orbit[s].s1);
      }
      else
      {
	twistedConjugate(tw,orbit[s].s0);
	if (hasDescent(orbit[s].s1,tw)) // is |s1| a left descent?
	{
	  twistedConjugate(tw,orbit[s].s1);
	  result.push_back(~s); // two-complex descent
	}
	else
	  result.push_back(s); // semi-real descent
      }
      break;
    case ext_gen::three:
      if (hasTwistedCommutation(orbit[s].s0,tw))
      {
	result.push_back(s); // a real-complex descent
	leftMult(tw,orbit[s].s0);
	twistedConjugate(tw,orbit[s].s1);
      }
      else
      {
	twistedConjugate(tw,orbit[s].s0);
	if (hasTwistedCommutation(orbit[s].s1,tw))
	{
	  leftMult(tw,orbit[s].s1);
	  result.push_back(s); // complex-real descent
	}
	else
	{
	  twistedConjugate(tw,orbit[s].s1);
	  twistedConjugate(tw,orbit[s].s0);
	  result.push_back(s); // triple complex descent
	}
      }
      break;
    }
    assert(twisted(tw)==tw); // stay within the twist-fixed subset
  } // while(tw!=delta)
  return result;
}


/*!
  \brief Returns the length of tw as a twisted involution.

  Precondition: tw is a twisted involution;

  Algorithm: this is a simplified version of |involutionOut| that records only
  the length. This statistic plays a predominant role in the kgb and block
  structures; avoid calling this in sorting routines, since it is inefficient
  in such circumstances; instead do with the stored length information there.
*/
unsigned long TwistedWeylGroup::involutionLength
  (const TwistedInvolution& tw) const
{
  TwistedInvolution x = tw;
  unsigned long length = 0;

  for (Generator s = W.leftDescent(x.w()); s != UndefGenerator;
                 s = W.leftDescent(x.w()),++length) {
    if (hasTwistedCommutation(s,x))
      leftMult(x,s);
    else
      twistedConjugate(x,s);
  }

  return length;
}

RootNbrList TwistedWeylGroup::simple_images
(const RootSystem& rs, const TwistedInvolution& tw) const
{
  assert(rank()==rs.rank()); // compatibility of Weyl groups required
  RootNbrList result(rank());
  for (size_t i=0; i<rank(); ++i)
    result[i]=rs.simpleRootNbr(twisted(i));
  WeylWord ww=word(tw.w());
  for (size_t i=ww.size(); i-->0;)
    rs.simple_root_permutation(ww[i]).renumber(result);

  return result;
}

WeightInvolution TwistedWeylGroup::involution_matrix
  (const RootSystem& rs, const TwistedInvolution& tw) const
{
  RootNbrList simple_image = simple_images(rs,tw);
  WeightList b(rank());
  for (size_t i=0; i<rank(); ++i)
    b[i] = rs.root_expr(simple_image[i]);

  return WeightInvolution(b,b.size());
}

} // |namespace weyl|


/*****************************************************************************

        Chapter II -- The Transducer class

  We implement here the construction of the Transducer tables (all accessors
  are defined in the class definition, and there are are no manipulators).
  This is described in section 4 of Fokko's 1999 paper "Transducer approach.."
  Actually, that paper does almost everything by induction, in particular it
  makes a (somewhat vague) reference to using previously constructed
  transducer tables while bootstrapping the current one; this is not what is
  done here, which proceeds strictly independently of other Transducer tables.
  The mention of dihedral groups below replaces the inductive part. [MvL]

******************************************************************************/

namespace weyl {

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
      for this is to look at the orbit of x under the dihedral group $\<s,t>$.
      In the full group, this has necessarily cardinality 2m, with m = m(s,t)
      the coefficient in the Coxeter matrix, and $l(xst)==l(x)$ iff $xs$ is
      the unique elt. of maximal length in the orbit, hence to have this $x$
      must goes down $m-1$ times when applying successively $t$, $s$, $t$, ...
      In the parabolic quotient, the orbit of the dihedral group (which is not
      reduced to a point) can either have cardinality $2m$ or cardinality $m$,
      and in the latter case it is a string with $m-1$ steps between the
      bottom and the top, with a stationary step at either extreme (to see
      this, note that on one hand each step up in the full group gives a step
      in the quotient that is either up or stationary, while on the other hand
      a stationary step in the quotient causes then next step to be the
      reverse of the previous one). So we have one of the following three
      cases: (1) $x$ goes down $m-1$ times; then the image of the orbit has
      $2m$ elements and $xs==x't$ for $x'=x.(ts)^(m-1)$. (2) $x$ goes down
      $m-2$ times to some element $a$ and is then stationary (if
      $s'\in\{s,t\}$ is the next to apply, then $a.s'=g.a$ for some generator
      $g\in W_{r-1}$). In case (2) the orbit has $m$ elements, and if $v$ is
      the alternating word in $\{s,t\}$ of length $m-2$ not starting with
      $s'$, so that $a.v=x$, one has $v.st=s'vs$ whence
      $x.st=a.v.st=a.s'vs=g.a.vs=g.xs$ so that $xs$ has a transduction for $t$
      that outputs the generator $g$. (3) either $x$ goes down less than $m-2$
      times, or $m-2$ times followed by an upward step; then $xst$ goes up.
*/

Transducer::Transducer(const int_Matrix& c, size_t r)
  : d_shift(1), d_out(1) // start with tables of size 1
  , d_length(1,0), d_piece(1,WeylWord()) // with empty word, length 0.

{
  // first row of transition and of transduction table

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

  for (WeylElt::EltPiece x = 0; x < d_shift.size(); ++x)
    for (Generator s = 0; s <= r; ++s)
      /* since RANK_MAX<128, |UndefEltPiece| is never a valid Piece number, so
         its presencein a slot in |d_shift| assures that this slot is
         unchanged from its intialisation value
      */
      if (d_shift[x][s] == UndefEltPiece)
      {

	WeylElt::EltPiece xs =
	  d_shift.size(); // index of state that will be added

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

	  unsigned long m = c(s,t); // coef from Coxeter matrix

	  WeylElt::EltPiece y  = dihedralMin(*this,xs,s,t);
	  unsigned long d = d_length[xs] - d_length[y];

	  if (d < m-1)
	    continue;  // case (3):  $t$ takes $xs$ up, do nothing

	  Generator st[] = {s,t};

	  if (d == m) {
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
      } // |if (..==UndefEltPiece)|, |for|, |for|
} // |Transducer::Transducer|



/*			    Other, small, classes			*/

// extract |Twist| information from a list of "extended generators"
Twist::Twist(const ext_gens& orbits)
{
  std::fill_n(&d[0],constants::RANK_MAX,Generator(~0));
  for (auto it=orbits.begin(); it!=orbits.end(); ++it)
    if (it->length()==1)
      d[it->s0]=it->s0;
    else
    {
      d[it->s0]=it->s1;
      d[it->s1]=it->s0;
    }
}

size_t TI_Entry::hashCode(size_t modulus) const
{
  unsigned int hash=0;
  for (size_t i=constants::RANK_MAX; i-->0; )
    hash = 13*(hash+(*this)[i]);
  return hash & (modulus-1);
}

} // |namespace weyl|

/*****************************************************************************

        Chapter III -- Functions declared in weyl.h

******************************************************************************/

namespace weyl {

/*!\brief Returns the twist defined by |d| relative to |rd|.

  Precondition: |d| is an involution of the root datum |rd|. If not an
  involution of the based datum, an appropriate Weyl group action is applied
 */
Twist make_twist(const RootDatum& rd, const WeightInvolution& d)
{
  RootNbrList simple_image(rd.semisimpleRank());

  for (size_t i = 0; i<simple_image.size(); ++i)
    simple_image[i] = rd.root_index(d*rd.simpleRoot(i));

  rootdata::wrt_distinguished(rd,simple_image); // and forget the Weyl element

  Twist result;

  for (size_t i = 0; i<simple_image.size(); ++i)
    result[i] = rd.simpleRootIndex(simple_image[i]);

  return result;
}


} // |namespace weyl|

/*****************************************************************************

        Chapter IV -- Auxiliary functions for this module

******************************************************************************/

namespace {

/*!
  \brief Returns the minimal element in the orbit of x under s and t.

  Precondition : s is in the descent set of x;
*/
WeylElt::EltPiece dihedralMin(const weyl::Transducer& qa,
				    WeylElt::EltPiece x,
				    weyl::Generator s,
				    weyl::Generator t)
{
  weyl::Generator u = s;
  weyl::Generator v = t;

  WeylElt::EltPiece y = x;

  for (;;) {
    // this is ok even if the shift is still undefined
    if (qa.shift(y,u) >= y)
      return y;
    else
      y = qa.shift(y,u);
    std::swap(u,v);
  }
}


/*!
  \brief Returns the result of applying s and t alternately to x, for a
  total of d times.
*/
WeylElt::EltPiece dihedralShift(const weyl::Transducer& qa,
				      WeylElt::EltPiece x,
				      weyl::Generator s,
				      weyl::Generator t,
				      unsigned long d)
{
  weyl::Generator u = s;
  weyl::Generator v = t;

  WeylElt::EltPiece y = x;

  for (unsigned long j = 0; j < d; ++j) {
    y = qa.shift(y,u);
    std::swap(u,v);
  }

  return y;
}


/*!
  \brief Fills in the Coxeter matrix cox.

  Precondition: cart is a Cartan matrix; a holds a normalizing permutation
  for cart, such as constructed by normalize(a,d) where d is the Dynkin diagram
  of cart (declared in dynkin.h);

  Postcondition : cox holds the normalized Coxeter matrix corresponding to
  cox and a;
*/
void fillCoxMatrix(int_Matrix& cox,
		   const int_Matrix& cart,
		   const Permutation& a)
{
  assert (cart.numRows()==cart.numColumns());
  int_Matrix(cart.numRows(),cart.numRows()) //create matrix
    .swap(cox); // and make |cox| refer to it

  static const int translate[] // from product of cart entries -> cox entry
    = { 2, 3, 4, 6 }; // N.B.  |0<=cart(i,j)*cart(j,i)<=3| always for |i!=j|
  for (size_t i=0; i<cart.numRows(); ++i)
  {
    cox(i,i) = 1;
    for (size_t j=i+1; j<cart.numRows(); ++j)
      cox(i,j) = cox(j,i) = translate[cart(i,j)*cart(j,i)];
  }

  // permute

  a.inv_conjugate(cox);
}


//				Template instantiation

} // |namespace|

namespace weyl {

template
void WeylGroup::act
  (const RootDatum& rd, const WeylElt& w, matrix::Vector<int>& v) const;


}

} // |namespace atlas|
