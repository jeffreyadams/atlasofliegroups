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

        Chapter I -- The WeylGroup class

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

WeylGroup::WeylGroup(const latticetypes::LatticeMatrix& c, const Twist* twist)
  :d_rank(c.numColumns()),
   d_maxlength(0UL),
   d_transducer(d_rank)

/*!
  Synopsis : constructs the Weyl group corresponding to the Cartan matrix c,
  and to the given twist.

  NOTE: a zero value for the Twist pointer indicates that we have to use the
  default (identity) twist.

  NOTE : we need to reorder the vertices carefully so that the subquotients
  fit in an unsigned char!
*/

{
  using namespace constants;
  using namespace dynkin;
  using namespace setutils;

  // analyze the coxeter matrix

  Permutation a;
  DynkinDiagram d(c);
  normalize(a,d);

  // put appropriate permutations in d_in and d_out
  // it should be so that d_in[j] is the _internal_ generator corresponding
  // to the external j, and similarly d_out[j] is the _external_ generator
  // corresponding to the internal j; of course d_in and d_out are inverses
  // of each other.

  for (size_t j = 0; j < d_rank; ++j) {
    d_out[j] = a[j];
    d_in[d_out[j]] = j;
  }
  
  // the default twist is the identity twist

  if (twist == 0)
    copy(d_twist,identityTwist());
  else for (size_t j = 0; j < d_rank; ++j)
    d_twist[j] = d_in[(*twist)[d_out[j]]];

  // fill in coxeter matrix

  fillCoxMatrix(d_coxeterMatrix,c,a);

  // construct the transducer

  for (size_t j = 0; j < d_rank; ++j)
    d_transducer[j] = Transducer(d_coxeterMatrix,j);

  // fill in the order, maxlength and longest element

  for (size_t j = 0; j < d_rank; ++j)
    d_longest[j] = d_transducer[j].size()-1;

  for (size_t j = 0; j < d_rank; ++j)
    d_maxlength += d_transducer[j].maxlength();

  for (size_t j = 0; j < d_rank; ++j) {
    d_order *= d_transducer[j].size();
  }

}

WeylGroup::WeylGroup(const WeylGroup& W, tags::DualTag)
  :d_rank(W.d_rank),
   d_order(W.d_order),
   d_maxlength(W.d_maxlength),
   d_longest(W.d_longest),
   d_coxeterMatrix(W.d_coxeterMatrix),
   d_transducer(W.d_transducer)

/*!
  Synopsis: the only difference with W is that the twist is multiplied by
  conjugation with the longest element.
*/

{
  memcpy(d_in,W.d_in,constants::RANK_MAX);
  memcpy(d_out,W.d_out,constants::RANK_MAX);
  memset(d_twist,0,constants::RANK_MAX);

  for (size_t s = 0; s < d_rank; ++s) {
    WeylElt w = d_longest;
    prodIn(w,W.d_twist[s]);
    prod(w,d_longest);
    Generator t = 0;
    for (; t < d_rank; ++t)
      if (w[t])
	break;
    d_twist[s] = t;
  }
}

void WeylGroup::swap(WeylGroup& other)

{  
  using namespace constants;

  std::swap(d_rank,other.d_rank);
  std::swap(d_order,other.d_order);
  std::swap(d_maxlength,other.d_maxlength);
  std::swap(d_longest,other.d_longest);
  d_coxeterMatrix.swap(other.d_coxeterMatrix);
  d_transducer.swap(other.d_transducer);

  // STL swap doesn't work for arrays!

  Generator tmp[RANK_MAX];

  memcpy(tmp,d_in,constants::RANK_MAX);
  memcpy(d_in,other.d_in,constants::RANK_MAX);
  memcpy(other.d_in,tmp,constants::RANK_MAX);

  memcpy(tmp,d_out,constants::RANK_MAX);
  memcpy(d_out,other.d_out,constants::RANK_MAX);
  memcpy(other.d_out,tmp,constants::RANK_MAX);

  memcpy(tmp,d_twist,constants::RANK_MAX);
  memcpy(d_twist,other.d_twist,constants::RANK_MAX);
  memcpy(other.d_twist,tmp,constants::RANK_MAX);

  return;
}

/******** accessors **********************************************************/

void WeylGroup::conjugacyClass(WeylEltList& c, const WeylElt& w, bool twisted) 
  const

/*!
  \brief Puts in c the conjugacy class of w.

  Algorithm: the orbit is enumerated in a time proportional to (orbit size)*
  rank. We maintain a stack of elements to be examined (initially w) and a
  set of found elements. Then while the stack is not empty :
    - v = stack.top(); pop the stack;
    - insert the various s.v.s in the set; if they are new, push them on
      the stack;

  NOTE: should throw an error::MemoryOverflow in case of memory overrun! 
  However I'm not sure how set and stack manage overflow conditions (if they 
  do at all!)
*/

{
  std::set<WeylElt> found;
  std::stack<WeylElt> toDo;

  found.insert(w);
  toDo.push(w);

  while (!toDo.empty()) {

    WeylElt v = toDo.top();
    toDo.pop();

    for (Generator s = 0; s < rank(); ++s) {
      WeylElt v1 = v;
      if (twisted)
	twistedConjugate(v1,s);
      else
	conjugate(v1,s);
      if (found.insert(v1).second)
	toDo.push(v1);
    }
  }

  c.clear();
  copy(found.begin(),found.end(),back_inserter(c));

  return;
}

bool WeylGroup::hasDescent(Generator s, const WeylElt& w) const

/*!
  \brief Tells whether sw < w.
*/

{
  WeylElt sw = w;
  leftProd(sw,s);

  return length(sw) < length(w);
}

bool WeylGroup::hasTwistedCommutation(Generator s, const WeylElt& w) const

/*!
  \brief Tells whether w twisted-commutes with s.
*/

{
  WeylElt ws = w;
  int m = prodIn(ws,d_twist[d_in[s]]);

  if (m == 1) // ws > w
    return hasDescent(s,ws);
  else
    return not hasDescent(s,ws);
}

void WeylGroup::invert(WeylElt& w) const

/*!
  \brief Transforms w into its inverse.

  Algorithm: read backwards the reduced expression gotten from the
  pieces of w.
*/

{
  WeylElt wi;

  for (size_t j = d_rank; j ;) {
    --j;
    const WeylWord& x_ww = wordPiece(w,j);
    for (size_t i = x_ww.size(); i;) {
      --i;
      prodIn(wi,x_ww[i]);
    }
  }

  w = wi;
}

unsigned long WeylGroup::involutionLength(const weyl::WeylElt& d_w) const

/*!
  \brief Returns the length of d_w as a twisted involution.

  Precondition: d_w is a twisted involution;

  Algorithm: let s be a generator s.t. s.w < w. Then if w.delta(s) = s.w, the 
  reduced expression is s.red(s.w); otherwise it is s.red(s.w.delta(s)).
*/

{
  using namespace weyl;

  WeylElt w = d_w;
  unsigned long length = 0;

  for (Generator s = leftDescent(w); s != UndefGenerator; s = leftDescent(w)) {
    ++length;
    if (hasTwistedCommutation(s,w))
      leftProd(w,s);
    else
      twistedConjugate(w,s);
  }

  return length;
}

void WeylGroup::involutionOut(weyl::WeylWord& ww, const weyl::WeylElt& d_w) 
  const

/*!
  \brief Puts in ww a reduced expression of d_w as a twisted involution.

  Precondition: d_w is a twisted involution;

  Algorithm: let s be a generator s.t. s.w < w. Then if w.delta(s) = s.w, the 
  reduced expression is s.red(s.w); otherwise it is s.red(s.w.delta(s)).
*/

{
  using namespace weyl;

  WeylElt w = d_w;
  ww.clear();

  for (Generator s = leftDescent(w); s != UndefGenerator; s = leftDescent(w)) {
    ww.push_back(s);
    if (hasTwistedCommutation(s,w))
      leftProd(w,s);
    else
      twistedConjugate(w,s);
  }

  // reverse ww
  std::reverse(ww.begin(),ww.end());

  return;
}

Generator WeylGroup::leftDescent(const WeylElt& w) const

/*!
  \brief Returns a left descent generator for w; UndefGenerator if there
  is no such (i.e., if w = e).
*/

{
  if (w == Identity)
    return UndefGenerator;

  for (size_t j = 0; j < d_rank; ++j) {
    const WeylWord& x = wordPiece(w,j);
    if (x.size())
      return d_out[x[0]];
  }

  return UndefGenerator; // this cannot be reached
}

void WeylGroup::leftProdIn(WeylElt& w, Generator s) const

/*!
  \brief Transforms w into s*w.

  Algorithm: note that our transducers are geared towards _right_ 
  multiplication by a generator. But we note that passing from w to sw
  only affects the pieces x_j in w for t <= j <= s, where t is the first
  generator that does not commute with s (remarkably, if v is the product
  of those pieces, sv does have non-zero components only for that set of
  indices; hard to believe at first but easy to prove.)
*/

{
  // find the generator t

  Generator t = 0;

  for(; t < s; ++t)
    if (d_coxeterMatrix(s,t) != 2)
      break;

  WeylElt sw;
  sw[s] = 1;

  for (size_t j = t; j <= s; ++j)
    prodIn(sw,wordPiece(w,j));

  for (size_t j = t; j <= s; ++j)
    w[j] = sw[j];

  return;
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
*/

{
  for (size_t s = 0; s < d_rank; ++s)
    t[s] = d_out[d_twist[d_in[s]]];

  return;
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

  return;
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

  return;
}

int WeylGroup::prodIn(WeylElt& w, Generator s) const

/*!  
  \brief Multiply w on the right by the (internally labelled)
  generator s: w *= s.

  Algorithm: the information for this is exactly contained in the transducer
  tables!

  Returns +1 if the length moves up, -1 if the length goes down.

  NOTE: The only way I see this could be written as an operator is if w 
  remembered the group that it came from!
*/

{
  Generator t = s;

  for (unsigned long j = d_rank; j;) {
    --j;
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

  return 0; // this should never happen!
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

  return;
}

unsigned long WeylGroup::toUlong(const WeylElt& w) const

/*!
  \brief Return the packed form of w (see the explanation in toWeylElt).

  NOTE: will work only if the order of the group fits into an unsigned long,
  i.e., if d_order is not equal to UndefOrder = 0;
*/

{
  unsigned long u = 0;
  unsigned long a = 1;

  for (size_t j = 0; j < d_rank; ++j) {
    u += w[j]*a;
    a *= d_transducer[j].size();
  }

  return u;
}

WeylElt WeylGroup::toWeylElt(unsigned long u) const

/*!
  \brief Returns the WeylElt whose packed form is u.

  Algorithm: the packed form of w is defined to be w_0.a_0 + ... + 
  w_{n-1}.a_{n-1}, where w_0, ... , w_{n-1} are the pieces of w, and
  a_0, ... , a_{n-1} are the cumulative subgroup sizes, defined by
  a_0 = 1, a_j = a_{j-1}*d_transducer[j-1].size() (in other words,
  the w_j are the "digits" of u in the variable base defined by the
  a_j.) Hence unpacking is a series of modular operations.

  NOTE: will work only if the order of the group fits into an unsigned long,
  i.e., if d_order is not equal to UndefOrder = 0;
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
  \brief Twists the element w by the outer involution.

  Algorithm: we used the standard reduced decomposition of w.
*/

{
  translate(w,d_twist);
  return;
}

}

/*****************************************************************************

        Chapter II -- The Transducer class

  ... explain here when it is stable ...

******************************************************************************/

namespace weyl {

Transducer::Transducer(const latticetypes::LatticeMatrix& c, size_t r)

/*!
\brief Constructs subquotient \#r for the Coxeter matrix c. 

  This uses the Coxeter matrix only up to index r.

  Precondition : c is a _normalized_ Coxeter matrix (meaning that all
  the parabolic subquotients W_{r-1}\\W_r are small enough to fit in
  an unsigned char); and r is < rank(c);

  Algorithm : 

  The algorithm is a version of my favorite bootstrapping procedure for the 
  construction of Weyl groups and parabolic quotients. We start with
  a partially defined automaton containing only one element, and for which
  all shifts not by smaller generators are undefined. At each point in time,
  all shifts for all elements in the automaton that do _not_ take the length
  up are defined; and we maintain a queue of elements that may have undefined
  shifts.

  We start up with one element in the automaton, with just one undefined
  shift, the one by r. Then run through the elements of the automaton, and
  for each undefined shift of x, say by s :

    - add a new element xs (this increases the size of the automaton, during
      the loop!)
    - for each generator t : we have to find out if xs.t goes up or down.
      The trick for this is to look at the orbit of x under the dihedral
      group <s,t>. In the full group, this has necessarily cardinality
      2m, with m = m(s,t) the coefficient in the Coxeter matrix, and xs.t
      goes down iff xs is the unique elt. of maximal length in the orbit,
      hence x goes down m-1 times when applying t,s,t, ... In the parabolic
      quotient, we can have either cardinality 2m, or cardinality m; in the
      latter case, the orbit is a string of length m, and xs is maximal iff
      x goes down m-2 times; then xs.t = u.xs, where u is the same generator
      as the one output by the minimal element in the orbit for the one
      among s,t which fixes it.
*/

{
  // initialize

  ShiftRow firstShift;
  OutRow firstOut;

  // all shifts lower than r should be defined

  for (size_t j = 0; j < r; ++j) {
    firstShift[j] = 0;
    firstOut[j] = j;
  }

  d_shift.push_back(firstShift);
  d_out.push_back(firstOut);
  d_length.push_back(0);
  d_piece.push_back(WeylWord()); // the empty WeylWord

  // in this loop, the table grows! the loop stops when x overtakes the
  // table size because no more new elements are created.

  for (EltPiece x = 0; x < d_shift.size(); ++x) {
    for (Generator s = 0; s <= r; ++s)
      // we need to check d_out[x][s] only in type B128, and only for the
      // last element of the last transducer ...
      if ((d_shift[x][s] == UndefEltPiece) and 
	  (d_shift[x][s] == UndefGenerator)) {

	EltPiece xs = d_shift.size();

	d_shift.push_back(ShiftRow());
	d_out.push_back(OutRow());

	d_shift[xs][s] = x;
	d_shift[x][s] = xs;
	// new element has length length(x)+1
	d_length.push_back(d_length[x]+1);
	
	WeylWord xs_ww = d_piece[x];
	xs_ww.push_back(s);
	d_piece.push_back(xs_ww);

	// define the shifts or transductions that do not take xs up

	for (Generator t = 0; t <= r; ++t) {

	  if (t == s)
	    continue;

	  LatticeCoeff m = c(s,t); // coefficient of the Coxeter matrix

	  EltPiece y  = dihedralMin(*this,xs,s,t);
	  unsigned long d = d_length[xs] - d_length[y];

	  if (d < static_cast<unsigned long>(m-1))
	    continue;  // t takes xs up

	  Generator st[] = {s,t};

	  if (d == static_cast<unsigned long>(m)) { 
	    // there is no transduction
	    // xs.t is computed by shifting up from y the other way around
	    y = dihedralShift(*this,y,st[m%2],st[(m+1)%2],m-1);
	    d_shift[xs][t] = y;
	    d_shift[y][t] = xs;
	    continue;  // next t
	  }

	  // now d == m-1;

	  Generator u = st[(m+1)%2];

	  if (d_shift[y][u] == y) { 
	    // xs is fixed by t and outputs the same generator as y
	    d_shift[xs][t] = xs;
	    d_out[xs][t] = d_out[y][u];
	  }
	  
	  // else xs moves up, do nothing
	}
      }
  }
}

}

/*****************************************************************************

        Chapter III -- Functions declared in weyl.h

  ... explain here when it is stable ...

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
  \brief Puts in refl the list of all reflections in W.

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
      W.conjugacyClass(c,simple[j],false); // ordinary conjugacy class!
      for (size_t i = 0; i < c.size(); ++i)
	found.insert(c[i]);
    }

  refl.clear();
  copy(found.begin(),found.end(),back_inserter(refl));

  return;
}

}

/*****************************************************************************

        Chapter IV -- Auxiliary functions for this module

  ... explain here when it is stable ...

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
      if (cart(i,j)) {
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

  return;
}

}

}
