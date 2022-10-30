
/*
  This is weyl.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007--2017, 2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  Implementation of WeylGroup.

  I have decided to represent elements as fixed-size arrays of
  unsigned characters. This forces expressing things in the standard
  ordering of the generators, and hence to have a small I/O interface
  for resetting the numbering to and from the numbering used by the
  outside world.

  It has seemed to me that this is the best compromise between size of the
  dataype, generality and efficiency.
  [Fokko du Cloux]
*/

#include "weyl.h"

#include <algorithm>
#include <set>

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

const EltPiece UndefEltPiece = UndefValue;
const Generator UndefGenerator = UndefValue;

//			     auxiliary functions
namespace {


  void fillCoxMatrix(int_Matrix&,
		     const int_Matrix&,
		     const WeylInterface&);

} // anonymous |namespace|

/*****************************************************************************

        Chapter 1 -- The WeylGroup::Transducer class

  The |Transducer| is used to implement right multiplication action of simple
  reflections on a Weyl group modulo (to the left) a parabolic subgroup. In the
  cases where right multiplication by |s| stabilises the parabolic coset, a new
  generator |t| is produced (to the left) to represent the commutation with the
  chosen (minimal length coset representative |w|, so $w*s=t*w|. Here |t| is
  always a generator of the parabolic subgroup to the left.

  In the notation from the description of the class WeylGroup, there will be
  one Transducer object for each parabolic subquotient W_{r-1}\\W_r. List the
  shortest length coset representatives for this subquotient as
  x_0,...,x_{N_r-1}. Recall that the simple roots were ordered to guarantee
  that N_r-1 fits in an unsigned char, so each coset representative can be
  indexed by an unsigned char. We wish to compute the product x_i.s_j for j
  between 1 and r. The key theoretical fact about multiplication is that there
  are two mutually exclusive possibilities:

      x_i.s_j = x_{i'}  (some i' ne i)

  OR

      x_i.s_j = s_k.x_i  (some k < r).

  The first possibility is called _transition_ and the second _transduction_.
  (Confusingly Fokko's 1999 paper interchanges these terms at their definition,
  but their usual meaning and the sequel makes clear that this was an error).

  The Transducer has tables to describe the two cases. the first table
  |d_shift| describes the transistions, namely |d_shift[i][j]==i'| in the
  first case; the cases that are transductions can be distinguished from these
  by the fact that |d_shift[i][j]==i|. In these cases, the value |k| emitted
  by the transduction is stored in |d_out[i][j]|, which otherwise contains the
  value |UndefGenerator|
*/

struct WeylGroup::Transducer
{
  // there will be one such object for each $r\in\{1,2,\ldots,n\}$
  // but $r$ is not explicitly stored in the |Tranducer| object

  // One row of a transducer table for a Weyl group.
  using RowBase = std::array<unsigned char,constants::RANK_MAX>;
  // this used to be a defined class contained in |Transducer|
  // that allowed default initialisation to Undef values; now done explicitly

  struct elem_info
  { RowBase shift; // Right multiplication by $s_j$ transitions to |shift[j]|
    RowBase out; // Right multiplication by $s_j$ may transduce |out[j]|
    unsigned char length; // length of |piece| in quotient Bruhat order
    Generator right; // rightmost factor of chosen word for this piece

    elem_info(unsigned short l)
      : shift(), out(), length(l), right(UndefGenerator)
    { shift.fill(UndefEltPiece); out.fill(UndefGenerator); }
  };

  using PieceIndex = unsigned short; // used to index letters inside a piece
  // as piece length cannot exceed number of states, |unsigned char| would do

  std::vector<elem_info> elt;

// constructors and destructors
  Transducer() {}

  Transducer(const int_Matrix& Coxeter, Generator s);

  ~Transducer() {}

// accessors

  static constexpr int max_piece_length = // for piece word reversal
    std::max(57ul,2*constants::RANK_MAX-1); // E8 needs 57, BCn needs 2n-1

  // Length of minimal coset representative x.
  unsigned int length(EltPiece x) const { return elt[x].length; }


/*
  Maximal length of minimal coset representatives.

  This is the number of positive roots for the Levi subgroup L_r, minus
  the number of positive roots for L_{r-1}.
*/
  unsigned int max_length() const { return elt.back().length; }

/*
  Simple reflection t (strictly preceding s) so that xs = tx, if any

  In case of a transition within the quotient, this returns |UndefGenerator|.
*/
  Generator out(EltPiece x, Generator s) const { return elt[x].out[s]; }

/*
  Right coset x' defined by x' = xs.

  When x' is not equal to s, this is an equality of minimal coset
  representatives. When x'=x, the equation for minimal coset representatives
  is out(x,s).x = x.s.
*/
  EltPiece shift(EltPiece x, Generator s) const { return elt[x].shift[s]; }

  Generator unshift(EltPiece& x) const; // assuming |x>0| to lowering shift

  // Number of cosets W_{r-1}\\W_r.
  unsigned int size() const { return elt.size(); } // must be at most 256


  // Reduced decomposition in W (or W_r) of minimal coset representative x.
  WeylWord word_of_piece(EltPiece x) const;

private:
  EltPiece
  dihedralMin(EltPiece,weyl::Generator,weyl::Generator) const;
  EltPiece
  dihedralShift(EltPiece,weyl::Generator,weyl::Generator,unsigned int)
    const;

  // this class should have no manipulators!
}; // |class WeylGroup::Transducer|


/*****************************************************************************

  We implement here the construction of the Transducer tables (all accessors
  are defined in the class definition, and there are are no manipulators).
  This is described in section 4 of Fokko's 1999 paper "Transducer approach.."
  Actually, that paper does almost everything by induction, in particular it
  makes a (somewhat vague) reference to using previously constructed
  transducer tables while bootstrapping the current one; this is not what is
  done here, which proceeds strictly independently of other Transducer tables.
  The mention of dihedral groups below replaces the inductive part. [MvL]

******************************************************************************/

/*
  Construct subquotient \#r for the Coxeter matrix c.

  This uses the Coxeter matrix only up to index r. In fact we can behave as
  if generator |r| is the final one, since we ignore any higher ones for now.

  Precondition : c is a _normalized_ Coxeter matrix (meaning that all
  the parabolic subquotients W_{r-1}\\W_r are small enough to fit in
  an unsigned char); and r is < rank(c);

  Algorithm :

  The algorithm is a version of my [Fokko] favorite bootstrapping procedure
  for the construction of Weyl groups and parabolic quotients. We start with a
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

WeylGroup::Transducer::Transducer(const int_Matrix& c, Generator r)
  : elt() // with empty word, length 0.

{
  // first row of transition and of transduction table
  elt.emplace_back(0); // length 0, empty |piece|, all fields Undef values

  // all shifts lower than |r| are transductions of unchanged generator
  for (Generator j = 0; j < r; ++j)
  {
    elt[0].shift[j]=0; // shift to self, no transition
    elt[0].out[j]=j;   // transduction of unchanged generator
  }
  // |elt[0].shift[r]=UndefEltPiece;| was set by |elem_info| constructor
  // |elt[0].out[r]  =UndefGenerator;| idem

  // in this loop, the |elt| table grows! the loop stops when x overtakes the
  // table size because no more new elements are created.

  for (EltPiece x = 0; x < elt.size(); ++x)
    for (Generator s = 0; s <= r; ++s)
      /* since RANK_MAX<128, |UndefEltPiece| is never a valid Piece number, so
         its presence in a slot in |d_shift| assures that this slot is
         unchanged from its intialisation value
      */
      if (elt[x].shift[s] == UndefEltPiece)
      {

	const EltPiece xs = elt.size(); // piece that will be added
	elt.emplace_back(elt[x].length+1);

	// only now can we fix (non transient) references
	auto& base = elt[x];
	auto& top  = elt.back();

	top.right=s; // last letter in word that leads to this |top|

	// |shift| and |out| fields of |top| are currently set to Undef values
	base.shift[s] = xs;
	top.shift[s] = x;

	// now define the shifts or transductions that do not take |xs| upwards

	for (Generator t = 0; t <= r; ++t)
	{
	  if (t == s)
	    continue;

	  const EltPiece b  = dihedralMin(xs,s,t);
	  const auto& bot = elt[b];
	  const unsigned int d = top.length - bot.length;
	  const unsigned int m = c(s,t); // positive coef from Coxeter matrix
	  const Generator st[] = {s,t}; // convenience

	  if (d == m)
	  { // case (1): there is no transduction
	    // xs.t is computed by shifting up from y the other way around
	    const EltPiece y = dihedralShift(b,st[m%2],st[1-m%2],m-1);
	    top.shift[t] = y;
	    elt[y].shift[t] = xs;
	  }
	  else if (d == m-1)
	  {
	    Generator u = st[1-m%2];

	    if (bot.shift[u] == b)
	    { // case (2): $xs$ is fixed by $t$, outputs the same $g$ as $y$
	      top.shift[t] = xs;
	      top.out[t]   = bot.out[u];
	    }
	  }
	  else assert(d<m-1); // case (3):  $t$ takes $xs$ up, do nothing
	} // |for(t)|
      } // |if (..==UndefEltPiece)|
    // |for(s)|
  // |for(x)|
} // |Transducer::Transducer|

// temporary implementation; see is |piece| can be replaced with final letter
Generator WeylGroup::Transducer::unshift(EltPiece& x) const
{
  auto s = elt[x].right;
  assert(s!=UndefGenerator); // or simpler |x>0|
  x = shift(x,s);
  return s;
}


WeylWord WeylGroup::Transducer::word_of_piece(EltPiece x) const
{
  unsigned int k=length(x);
  std::vector<Generator> result(k);
  while (x>0)
    result[--k] = unshift(x);
  return WeylWord{std::move(result)};
}

/*****************************************************************************

        Chapter 2 -- The WeylGroup and WeylElt classes

  The Weyl group class is a variation on my favourite implementation in terms
  of transducers.

  I have tried to make a careful choice of datatype for the group elements in
  order to strike the right balance between efficiency and generality. This
  has led me to the choice of _fixed size_ arrays of unsigned chars,
  representing the "parabolic subquotient" representation of a given element;
  in other words, in a group of rank n, the first n elements of the array are
  used (but since it is fixed size, we have to allocate it to some constant
  value, namely RANK_MAX.) This is not so wasteful as it may sound : of course
  the representation as a single number is more compact, but will overflow
  even on 64-bit machines for some groups of rank <= 16, and much sooner of
  course on 32-bit machines; also it imposes some computational overhead for
  packing and unpacking. Any variable-size structure like the STL vector uses
  already three pointers for its control structure (the address of the data,
  the size and the capacity), and then it still has to allocate. This could
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

/*
  Build the Weyl group corresponding to the Cartan matrix |c|,
  and incorporate the possibly given twist (take identity if |twist==nullptr|).

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

  DynkinDiagram d(c); // make diagram from Cartan matrix and classify

/*
  now put appropriate permutations into |d_in| and |d_out|, so that internal
  number |i| gives external number |d_out[i]|, and |d_in[d_out[i]]==i|.
*/
  for (const auto& comp : d.components())
    if (std::strchr("BCD",comp.type)==nullptr)
      for (unsigned i=0; i<comp.position.size(); ++i)
	d_out[comp.offset+i]=comp.position[i];
    else // reverse in types BCD
    {
      unsigned last=comp.offset+comp.position.size()-1;
      for (unsigned i=0; i<comp.position.size(); ++i)
	d_out[last-i]=comp.position[i];
    }

  for (Generator i=0; i<d_rank; ++i)
    d_in[d_out[i]] = i;

  // deduce the Coxeter matrix from Cartan matrix |c| and renumbering |a|
  fillCoxMatrix(d_coxeterMatrix,c,d_out);

  // now construct the transducer
  // each one gets to know the Coxeter matrix and its place in the table
  for (Generator j = 0; j < d_rank; ++j)
    d_transducer[j] = Transducer(d_coxeterMatrix,j);

  // the longest element has the maximal valid value in each of its pieces
  for (Generator j = 0; j < d_rank; ++j)
    d_longest.piece(j) = d_transducer[j].size()-1;

  // its length is obtained by summing the (maximal) lengths of its pieces
  for (Generator j = 0; j < d_rank; ++j)
    d_maxlength += d_transducer[j].max_length();

  // and the Weyl group size is the product of the numbers of transducer states
  for (Generator j = 0; j < d_rank; ++j)
    d_order *= d_transducer[j].size();

  // precompute for each |j| the first non-commuting or equal generator |i<=j|
  for (Generator j = 0; j < d_rank; ++j)
    for (Generator i=0; i<=j; ++i)
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

// these methods use the default, but that requires |Transducer| type complete
WeylGroup::WeylGroup(WeylGroup&& W) = default;
WeylGroup::~WeylGroup() = default;

/******** accessors **********************************************************/

WeylElt WeylGroup::genIn (Generator i) const
{
  WeylElt s;    // initialise all pieces to 0 (identity)
  s.piece(i)=1; // then shift to piece #1 in $W_{i-1}\\W_i$, which is $s_i$
  return s;
}

/*
  Multiply |w| on the right by internally numbered generator |s|: |w *= s|.

  Returns +1 if the length moves up, -1 if the length goes down.

  This just means exercising the transducer tables as they were designed to.
*/

/*
  Amazingly, I could simplify Fokko's original to the code below.
  Fokko's original code was (more or less):

  Generator t = s;
  for (Generator j = d_rank; j-->0;) {
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
int WeylGroup::inner_mult(WeylElt& w, Generator s) const
{
  Generator j = d_rank-1; // current transducer

  // in next loop |j| cannot pass |0| since |d_transducer[0]| only has shifts
  for (Generator t; (t=d_transducer[j].out(w.piece(j),s))!=UndefGenerator; s=t)
  {
    assert(j!=0);
    --j;
  }

  // now transductions are exhausted and one nontrivial shift remains
  EltPiece wj=w.piece(j);
  w.piece(j)=d_transducer[j].shift(wj,s); // assert: |w[j]!=wj|

  return w.piece(j)>wj ? 1 : -1; // no need to use length, numeric '>' is OK
}

/*
  Multiply |w| on the right by piece |i| of |v| (in internal numbering)

  Returns nonpositive even value $l(wv)-l(w)-l(v)$

  Precondition: |v| is written in the _internal_ generators.

  Algorithm: do the elementary multiplication by the generators, running
  through |v| left-to-right.
*/
int WeylGroup::mult_by_piece(WeylElt& w, const WeylElt& v, Generator i) const

{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal
  EltPiece x = v.piece(i);
  const auto& tr = d_transducer[i];
  int result = -static_cast<int>(tr.elt[x].length);
  // final result takes -2 for each shortening |inner_mult|, 0 otherwise

  unsigned int k=0;
  while (x>0)
    stack[k++]=tr.unshift(x);
  while (k-->0)
    result+=inner_mult(w,stack[k]);

  return result;
}

/*
  Transform |w| into |s*w|, with |s| in internal numbering;
  returns legth change $l(sw)-l(w)\in\{+1,-1}$

  Note that our transducers are geared towards _right_ multiplication by
  a generator, so this left multiplication is less straightforward.
  Nonetheless there is considerably better than the obvious algorithm of
  starting with a |WeylElt| for the generator |s|, and right-multiplying it
  successively by all letters of a word for |w|. In fact it turns out that
  comparing $w$ with $sw$, only pieces $x_j$ can change for $t <= j <= s$,
  where |t=min_neighbor(s)| is the first generator that does not commute
  with |s|; due to the numbering of diagrams used there are at most 3 such
  $j$, and when there are 3 (for non linear Dynkin diagrams) they occur
  only for small |s|, where the pieces are relatively short. What we are
  claiming is not just that right-multiplying the letters from the words
  for pieces of |w| before |t| all transduce without change across piece |s|
  and doing so reconstruct the same pieces again (this is easy to see), but
  also that when inserting the letters from the pieces $x_j$ that do need to
  be taken into account, only the pieces for the same range of |j| values
  are affected (meaning that none of the |inner_mult| calls invoked by the
  code below will transduce leftwards of piece |t|). Formulated like this it
  is hard to believe at first, but it is easy to prove; for instance, if one
  decomposes the expression for $w$ given by its pieces as $w_1w_2$ where
  $w_1$ is the part in the parabolic factor that commutes with $s$, then
  $sw=sw_1w_2=w_1sw_2$ can only have reductions inside the part $sw_2$,
  leaving the left part $w_1$ unchanged. These observations justify that
  the two loops below only run over the mentioned interval of values |j|.
*/
int WeylGroup::leftMultIn(WeylElt& w, Generator s) const
{
  WeylElt sw=genIn(s);
  int l=1; //

  // now compute $sv$ as above, keeping track of any length drop (at most 1)
  for (Generator j = min_neighbor(s); j <= s; ++j)
    l+=mult_by_piece(sw,w,j);

  // and copy its relevant pieces into $w$
  for (Generator j = min_neighbor(s); j <= s; ++j)
    w.piece(j) = sw.piece(j);

  return l;
}

WeylWord WeylGroup::word_of_piece(const WeylElt& w, Generator j) const
{ return d_transducer[j].word_of_piece(w.piece(j)); }

/*
  Multiply |w| on the right by |v|, and put the product in |w|: |w*=v|.

  Algorithm: increment |w| by the various pieces of |v|, whose reduced
  expressions are available from the transducer.
*/
void WeylGroup::mult(WeylElt& w, const WeylElt& v) const
{
  for (Generator j = 0; j < d_rank; ++j)
    mult_by_piece(w,v,j);
}

/*
  Multiply |w| on the right by |v|, and put the product in |w|: |w*=v|.

  Algorithm: do the elementary multiplication by the generators, running
  through v left-to-right.
*/
void WeylGroup::mult(WeylElt& w, const WeylWord& v) const
{
  for (unsigned int j = 0; j < v.size(); ++j)
    mult(w,v[j]);
}

// Set |w=xw|
void WeylGroup::leftMult(WeylElt& w, const WeylElt& x) const
{
  WeylElt xx=x; mult(xx,w); w=xx;
}

/*
  Inverse of |w|

  Algorithm: read backwards the reduced expression gotten from the
  pieces of |w|.
*/
WeylElt WeylGroup::inverse(const WeylElt& w) const
{
  WeylElt wi;

  for (Generator j = d_rank; j-->0 ;)
  { EltPiece x = w.piece(j);
    const auto& tr = d_transducer[j];
    while (x>0)
      inner_mult(wi,tr.unshift(x));
  }

  return wi;
}

/*
  Put into |c| the conjugacy class of |w|.

  Algorithm: straightforward enumeration of the connected component of |w| in
  the graph defined by the operation |conjugate|, using a |std::set| structure
  to record previously encountered elements and a |containers::queue| to store
  elements whose neighbors have not yet been generated.

*/
void WeylGroup::conjugacyClass(WeylEltList& c, const WeylElt& w) const
{
  std::set<WeylElt> found { w };
  containers::queue<WeylElt> to_do { w };

  while (not to_do.empty())
  {
    WeylElt v = to_do.front(); to_do.pop();

    for (Generator s = 0; s < rank(); ++s)
    {
      conjugate(v,s);
      if (found.insert(v).second) // then it was new, put it on to-do list
	to_do.push(v);
    }
  }

  // now convert set |found| to vector
  c.reserve(found.size());
  c.assign(found.begin(),found.end());
}

/*
  Whether |sw < w|.

  Method: we multiply from $s$ to $sw$, at least by the word pieces of |w| at
  the relevant pieces: those from |min_neighbor(s)| to |s| inclusive. If any
  descent occurs we return |true|, otherwise |false|. Despite the double loop
  below, this question is resolved in relatively few operations on the average.
*/
bool WeylGroup::hasDescent(Generator s, const WeylElt& w) const
{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal

  s=d_in[s]; // inner numbering is used below
  WeylElt sw = genIn(s); // element operated upon, starts out as |s|

  for (Generator j = min_neighbor(s); j <= s; ++j)
  {
    EltPiece x = w.piece(j);
    const auto& tr = d_transducer[j];

    unsigned int k=0;
    while (x>0)
      stack[k++]=tr.unshift(x);
    while (k-->0)
      if (inner_mult(sw,stack[k])<0) // multiply and see if a descent occurs
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
  for (Generator t; (t=d_transducer[j].out(w.piece(j),s))!=UndefGenerator; s=t)
    --j;

  EltPiece wj=w.piece(j);
  return d_transducer[j].shift(wj,s)<wj;
}


/*
  Return a left descent generator for |w|, or |UndefGenerator| if
  there is no such (i.e., if $w = e$). In fact this is the index |i| of the
  first piece of |w| that is not 0 (identity), converted to external
  numbering, since the canonical (minimal for ShortLex) expression for |w|
  starts with this letter. Returning the first generator in external numbering
  would take a bit more time, as we would have to repeatedly use |hasDescent|.
*/
Generator WeylGroup::leftDescent(const WeylElt& w) const
{
  for (Generator i=0; i<d_rank; ++i)
    if (w.piece(i)>0) return d_out[i];

  // if we come here, |w==e|
  return UndefGenerator;
}


/*
  Return the length of |w|.

  This is relatively efficient (compared to |involutionLength|)
*/
unsigned int WeylGroup::length(const WeylElt& w) const
{
  unsigned int l = 0;

  for (Generator i = 0; i < d_rank; ++i)
    l += d_transducer[i].length(w.piece(i));

  return l;
}

WeylWord WeylGroup::word(const WeylElt& w) const
{
  unsigned int k=length(w);
  std::vector<Generator> result(k);
  for (Generator i = d_rank; i-->0; )
  {
    EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];
    while (x>0)
      result[--k] = d_out[tr.unshift(x)];
  }
  assert(k==0);

  return WeylWord{std::move(result)};
}

/*
  Return the packed form of |w|

  This is the mixed-radix interpretation of the sequence of pieces, where
  piece 0 is the least significant one: if $N_i$ is |d_transducer[i].size()|
  and $a_i\in\{0,\ldots,N_i-1\}$ is the value of piece |i|, the value is
  $a_0+N_0*(a_1+N_1*(a_2+... N_{n-2}*(a_{n-1})...))$
*/
arithmetic::big_int WeylGroup::to_big_int(const WeylElt& w) const
{
  arithmetic::big_int u(0);

  for (Generator j=d_rank; j-->0; )
    (u*=static_cast<int>(d_transducer[j].size()))+=static_cast<int>(w.piece(j));

  return u;
}

/*
  Return the WeylElt whose packed form is |u|

  Its pieces for the mixed-radix representation of |u|, which as usual can be
  found starting at the least significant end by repeated Euclidean divisions
*/
WeylElt WeylGroup::toWeylElt(arithmetic::big_int u) const
{
  WeylElt w;

  for (Generator j = 0; j < d_rank; ++j)
    w.piece(j) = u.shift_modulo(static_cast<int>(d_transducer[j].size()));

  return w;
}

/*
  Apply to |w| the generator permutation in |outer_f|, which should be an
  automorphism of the Dynkin diagram, expressed in terms of outer numbering.

  Algorithm: after transforming |outer| we use the standard reduced
  decomposition of |w|, and rebuild the return value by repeated
  right-multiplication. We can do the multiplication in the same Weyl group as
  the decomposition because |f| is supposed to be an automorphism; if it weren't
  we would need a reference to a second Weyl group.
*/
WeylElt WeylGroup::translation
  (const WeylElt& w, const WeylInterface& outer_f) const
{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal

  WeylInterface inner_f;
  for (Generator i = 0; i<rank(); ++i)
    inner_f[i] = d_in[outer_f[d_out[i]]];

  WeylElt result; // start out with identity
  for (Generator i=0; i<rank(); ++i)
  {
    EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    unsigned int k=0;
    while (x>0)
      stack[k++]=tr.unshift(x);
    while (k-->0)
      inner_mult(result,inner_f[stack[k]]);
  }

  return result;
}

WeylElt WeylGroup::reverse_translation
  (const WeylElt& w, const WeylInterface& outer_f) const
{
  WeylInterface inner_f;
  for (Generator i = 0; i<rank(); ++i)
    inner_f[i] = d_in[outer_f[d_out[i]]];

  WeylElt result; // start out with identity
  for (Generator i=rank(); i-->0; )
  {
    EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      inner_mult(result,inner_f[tr.unshift(x)]);
  }

  return result;
}

/*
  Let |w| act on |alpha| according to reflection action in root datum |rd|
  Note that rightmost factors act first, as in a product of matrices
*/
void
  WeylGroup::act(const RootSystem& rd, const WeylElt& w, RootNbr& alpha) const
{
  for (Generator i = d_rank; i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];
    while (x>0)
      rd.simple_reflect_root(d_out[tr.unshift(x)],alpha);
  }
}

// Let |w| act on weight |v| according to reflection action in root datum |rd|
template<typename C>
  void WeylGroup::act
    (const RootDatum& rd, const WeylElt& w,  matrix::Vector<C>& v) const
{
  for (Generator i = d_rank; i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];
    while (x>0)
      rd.simple_reflect(d_out[tr.unshift(x)],v);
  }
}

// Let |w| act on coweight |v| according to reflection action in root datum |rd|
template<typename C>
  void WeylGroup::co_act
    (const RootDatum& rd,  matrix::Vector<C>& v, const WeylElt& w) const
{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal

  for (Generator i = 0; i<d_rank; ++i)
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    unsigned int k=0;
    while (x>0)
      stack[k++]=tr.unshift(x);
    while (k-->0)
      rd.simple_coreflect(v,d_out[stack[k]]);
  }
}

void WeylGroup::act(const RootDatum& rd, const WeylElt& w, RatWeight& v) const
{ act(rd,w,v.numerator()); }

void WeylGroup::act(const RootDatum& rd, const WeylElt& w, LatticeMatrix& M)
  const
{
  for (Generator i = d_rank; i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      rd.simple_reflect(d_out[tr.unshift(x)],M);
  }
}


template<typename C>
  void WeylGroup::act
    (const PreRootDatum& prd, const WeylElt& w, matrix::Vector<C>& v) const
{
  for (Generator i = d_rank; i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      prd.simple_reflect(d_out[tr.unshift(x)],v);
  }
}

void
  WeylGroup::act(const PreRootDatum& prd, const WeylElt& w, RatWeight& v) const
{ act(prd,w,v.numerator()); }

void WeylGroup::act(const PreRootDatum& prd, const WeylElt& w, LatticeMatrix& M)
  const
{
  for (Generator i = d_rank; i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      prd.simple_reflect(d_out[tr.unshift(x)],M);
  }
}

/*
  Same as |act(rd,inverse(w),v)|, but avoiding computation of |inverse(w)|.
  Here the leftmost factors act first.
*/
void
  WeylGroup::inverse_act(const RootDatum& rd, const WeylElt& w, Weight& v) const
{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal
  for (Generator i=0; i<d_rank; ++i )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    unsigned int k=0;
    while (x>0)
      stack[k++]=tr.unshift(x);
    while (k-->0)
      rd.simple_reflect(d_out[stack[k]],v);
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
{}

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

/*
  the "dual" Weyl group: the only difference with W is that the twist is
  multiplied by conjugation with the longest element.
*/
TwistedWeylGroup::TwistedWeylGroup(const TwistedWeylGroup& tW, tags::DualTag)
  : W(tW.W)
  , d_twist(tW.dual_twist())
{}

void TwistedWeylGroup::twistedConjugate // $tw = w.tw.twist(w)^{-1}$
  (TwistedInvolution& tw, const WeylElt& w) const
{
  WeylElt x=w; W.mult(x,tw.w()); W.mult(x,W.reverse_translation(w,d_twist));

  tw.contents()=x;
}

/*
  Put into |c| the twisted conjugacy class of |w|.

  Algorithm: straightforward enumeration of the connected component of |w| in
  the graph defined by the operation |twistedConjugate|, using a |std::set|
  to record previously encountered elements and a |containers::queue| to store
  elements whose neighbors have not yet been generated.
*/
void TwistedWeylGroup::twistedConjugacyClass
  (TwistedInvolutionList& c, const TwistedInvolution& tw)
  const
{
  std::set<TwistedInvolution> found { tw } ;
  containers::queue<TwistedInvolution> to_do { tw };

  while (not to_do.empty())
  {
    TwistedInvolution v = to_do.front();
    to_do.pop();

    for (Generator s=0; s<W.rank(); ++s)
    {
      twistedConjugate(v,s);
      if (found.insert(v).second)
	to_do.push(v);
    }
  }

  // now convert set |found| to vector
  c.reserve(found.size());
  c.assign(found.begin(),found.end());
}


/*
  Tell whether |w| twisted-commutes with |s|: $s.w.\delta(s)=w$

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


/*
  Precondition: |tw| is a twisted involution: $tw^{-1}=\delta(tw)$.

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


/*
  The length of |tw| as a twisted involution.

  Precondition: |tw| is a twisted involution;

  Algorithm: this is a simplified version of |involutionOut| that records only
  the length. This statistic plays a predominant role in the kgb and block
  structures; avoid calling this in sorting routines, since it is inefficient
  in such circumstances; instead do with the stored length information there.
*/
unsigned int TwistedWeylGroup::involutionLength
  (const TwistedInvolution& tw) const
{
  TwistedInvolution x = tw;
  unsigned int length = 0;

  for (Generator s = W.leftDescent(x.w()); s != UndefGenerator;
                 s = W.leftDescent(x.w()),++length)
    if (hasTwistedCommutation(s,x))
      leftMult(x,s);
    else
      twistedConjugate(x,s);

  return length;
}

RootNbrList
  TwistedWeylGroup::simple_images
     (const RootSystem& rs, const TwistedInvolution& tw) const
{
  assert(rank()==rs.rank()); // compatibility of Weyl groups required
  auto w=tw.w();
  RootNbrList result(rank());
  for (Generator i=0; i<rank(); ++i)
  {
    result[i]=rs.simpleRootNbr(twisted(i));
    W.act(rs,w,result[i]);
  }

  return result;
}

WeightInvolution TwistedWeylGroup::involution_matrix
  (const RootSystem& rs, const TwistedInvolution& tw) const
{
  RootNbrList simple_image = simple_images(rs,tw);
  WeightList b(rank());
  for (Generator i=0; i<rank(); ++i)
    b[i] = rs.root_expr(simple_image[i]);

  return WeightInvolution(b,b.size());
}



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
    hash = 13*(hash+piece(i));
  return hash & (modulus-1);
}


/*****************************************************************************

        Chapter III -- Functions declared in weyl.h

******************************************************************************/

/*
  Return the twist defined by |d| relative to |rd|.

  Precondition: |d| is an involution of the root datum |rd|. If not an
  involution of the based datum, an appropriate (but not recorded) Weyl group
  action is applied to make simple roots map to simple roots.
 */
Twist make_twist(const RootDatum& rd, const WeightInvolution& d)
{
  RootNbrList simple_image(rd.semisimple_rank());

  for (Generator i = 0; i<simple_image.size(); ++i)
    simple_image[i] = rd.root_index(d*rd.simpleRoot(i));

  rootdata::wrt_distinguished(rd,simple_image); // and forget the Weyl element

  Twist result;

  for (Generator i = 0; i<simple_image.size(); ++i)
    result[i] = rd.simpleRootIndex(simple_image[i]);

  return result;
}

/*
  Return the minimal element in the orbit of |x| under |s| and |t|.

  Precondition : |s| is in the descent set of |x|;
*/
EltPiece weyl::WeylGroup::Transducer::dihedralMin
  (EltPiece x, weyl::Generator s, weyl::Generator t) const
{
  weyl::Generator u = s; weyl::Generator v = t; // to be swapped below

  EltPiece y = x;

  for (;;)
  { // this is ok even if the shift is still undefined:
    if (shift(y,u) >= y)
      return y;
    else
      y = shift(y,u);
    std::swap(u,v);
  }
}

/*
  Return the result of applying |s| and |t| alternately to |x|, for a
  total of |d| times.
*/
EltPiece weyl::WeylGroup::Transducer::dihedralShift
(EltPiece x,
 weyl::Generator s,
 weyl::Generator t,
 unsigned int d) const
{
  weyl::Generator u = s;
  weyl::Generator v = t;

  EltPiece y = x;

  for (unsigned int j = 0; j < d; ++j)
  {
    y = shift(y,u);
    std::swap(u,v);
  }

  return y;
}


/*****************************************************************************

        Chapter IV -- Auxiliary functions for this module

******************************************************************************/

namespace {




/*
  Fill in the Coxeter matrix |cox|.

  Precondition: cart is a Cartan matrix; a holds a normalizing permutation
  for cart, such as constructed by normalize(a,d) where d is the Dynkin diagram
  of cart (declared in dynkin.h);

  Postcondition : cox holds the normalized Coxeter matrix corresponding to
  cox and a;
*/
void fillCoxMatrix(int_Matrix& cox,
		   const int_Matrix& cart,
		   const WeylInterface& a)
{
  assert (cart.numRows()==cart.numColumns());
  int_Matrix(cart.numRows(),cart.numRows()) //create matrix
    .swap(cox); // and make |cox| refer to it

  static const int translate[] // from product of cart entries -> cox entry
    = { 2, 3, 4, 6 }; // N.B.  |0<=cart(i,j)*cart(j,i)<=3| always for |i!=j|
  for (Generator i=0; i<cart.numRows(); ++i)
  {
    cox(i,i) = 1;
    unsigned ai=a[i]; // index corresponding to (inner) |i| in (outer) |cart|
    for (Generator j=i+1; j<cart.numRows(); ++j)
      cox(i,j) = cox(j,i) = translate[cart(ai,a[j])*cart(a[j],ai)];
  }

}

} // anonymous |namespace|

//				Template instantiation

template
void WeylGroup::act
  (const RootDatum& rd, const WeylElt& w, matrix::Vector<int>& v) const;

template
void WeylGroup::co_act
  (const RootDatum& rd, matrix::Vector<int>& v, const WeylElt& w) const;


} // |namespace weyl|

} // |namespace atlas|
