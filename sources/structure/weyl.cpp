
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

  struct elem_info
  { unsigned char length; // length of |piece| in quotient Bruhat order
    Generator right; // rightmost factor of chosen word for this piece

    elem_info(unsigned short l): length(l), right(UndefGenerator) {}
  };

  using PieceIndex = unsigned short; // used to index letters inside a piece
  // as piece length cannot exceed number of states, |unsigned char| would do

  // data fields
  Generator offset; // first generator in this Dynkin diagram component
  Generator limit; // |shift| and |out| only valid for |Generator| below |limit|
  std::vector<elem_info> elt;
  matrix::Matrix<Generator> table;

// constructors and destructors
  Transducer() {}

  Transducer(lietype::TypeLetter type, Generator offset, Generator s);

  ~Transducer() {}

// accessors

  static constexpr int max_piece_length = // for piece word reversal
    std::max(57ul,2*constants::RANK_MAX-1); // E8 needs 57, BCn needs 2n-1

  // conversion into numbering for this component only and back again
  Generator local(Generator s) const { return s-offset; }
  Generator unlocal(Generator s) const { return s+offset; }

  // Length of minimal coset representative x.
  unsigned int length(EltPiece x) const { return elt[x].length; }


/*
  Maximal length of minimal coset representatives.

  This is the number of positive roots for the Levi subgroup L_r, minus
  the number of positive roots for L_{r-1}.
*/
  unsigned int max_length() const { return elt.back().length; }

  bool has_shift(EltPiece x, Generator s) const
  { assert(x<elt.size()); assert(s<limit);
    return table(x,s)<elt.size();
  }
/*
  Simple reflection t (strictly preceding s) so that xs = tx, if any

  In case of a transition within the quotient, this returns |UndefGenerator|.
*/
  Generator out(EltPiece x, Generator s) const
  { assert(not has_shift(x,s)); return table(x,s)-elt.size(); }

/*
  Right coset x' defined by x' = xs.

  When x' is not equal to s, this is an equality of minimal coset
  representatives. When x'=x, the equation for minimal coset representatives
  is out(x,s).x = x.s.
*/
  EltPiece shift(EltPiece x, Generator s) const
  { assert(has_shift(x,s)); return table(x,s); }

  Generator unshift(EltPiece& x) const; // assuming |x>0| to lowering shift

  // Number of cosets W_{r-1}\\W_r.
  unsigned int size() const { return elt.size(); } // must be at most 256


  // Reduced decomposition in W (or W_r) of minimal coset representative x.
  WeylWord word_of_piece(EltPiece x) const;

  // this class should have no manipulators!
}; // |class WeylGroup::Transducer|

/* Since the |Transducer| construction assumes effectively a connected Dynkin
   diagram (the |Generator| values used in its table are made relative to the
   current component) and a standard ordering for each type is used, we no
   longer pass a Coxeter matrix to the constructor, but just a (simple) type.
   The following function computes Coxeter matrix entries for any type, in the
   order we use, namely Bourbaki except for reversal in types B,C,D
*/

unsigned int Coxeter_entry(lietype::TypeLetter type, Generator i, Generator j)
{
  if (i>j)
    std::swap(i,j); // ensure |i<=j|, the matrix is symmetric
  if (std::strchr("DE",type)==nullptr) // whether diagram is linear
    switch(j-i) // absolute distance
    {
    default: return 2; // commutation if at least 2 apart
    case 0: return 1; // this case is not really used
    case 1: return
	type=='A' or
	(std::strchr("BC",type)!=nullptr and i>0) or
	(type=='F' and i!=1)
	? 3
	: type=='G' ? 6
	: 4; // cases BC and $(i,j)=(0,1)$ or F and $(i,j)=(1,2)$
    }

  // types DE
  if (i==0)
    return j==2 ? 3 : 2;
  else if (type=='E' and i==1)
    return j==3 ? 3 : 2;
  else
    return j-i==1 ? 3 : 2;
}

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

WeylGroup::Transducer::Transducer
  (lietype::TypeLetter type, Generator offset, Generator r)
    // here |r| is relative to |offset|
    : offset(offset), limit(r+1), elt(), table()
{
  // One row of a transducer table for a Weyl group.
  using RowBase = std::array<unsigned char,constants::RANK_MAX>;
  struct entry
  { RowBase shift; // Right multiplication by $s_j$ transitions to |shift[j]|
    RowBase out; // Right multiplication by $s_j$ may transduce |out[j]|

    entry() { shift.fill(UndefEltPiece); out.fill(UndefGenerator); }
  };
  std::vector<entry> tab; // local that will determine |table| at the end

  // first row of transition and of transduction table
  elt.emplace_back(0); // length 0, empty |piece|, all fields Undef values
  tab.emplace_back();

  // all shifts lower than |r| are transductions of unchanged generator
  for (Generator i = 0; i < r; ++i)
  {
    tab[0].shift[i]=0; // shift to self, no transition
    tab[0].out[i]=i;   // transduction of unchanged generator
  }
  // |tab[0].shift[r]=UndefEltPiece;| was set by |entry| constructor
  // |tab[0].out[r]  =UndefGenerator;| idem

  // In this loop, the |elt| and |tab| tables grow! The loop stops when |x|
  // overtakes the table size because no more new elements are created.

  for (EltPiece x = 0; x < elt.size(); ++x)
    for (Generator s = 0; s <= r; ++s)
      /* since RANK_MAX<128, |UndefEltPiece| is never a valid Piece number, so
         its presence in a slot in |d_shift| assures that this slot is
         unchanged from its intialisation value
      */
      if (tab[x].shift[s] == UndefEltPiece)
      {

	const EltPiece xs = elt.size(); // piece that will be added
	elt.emplace_back(elt[x].length+1);
	tab.emplace_back();

	auto& top  = tab.back(); // or |tab[xs]|

	elt.back().right=s; // last letter in word that leads to the |top|

	// |shift| and |out| fields of |top| are currently set to Undef values
	tab[x].shift[s] = xs;
	top.shift[s] = x;

	// now define the shifts or transductions that do not take |xs| upwards

	for (Generator t = 0; t <= r; ++t)
	{
	  if (t == s)
	    continue;

	  EltPiece b = x;
	  { // set |b| to minimun of dihedral orbit for |s| and |t|
	    EltPiece next;
	    while (true)
	    { // following test is OK even if |next==UndefEltPiece|
	      if ((next=tab[b].shift[t]) < b)
		b = next;
	      else break;
	      // now do the same for the other generator
	      if ((next=tab[b].shift[s]) < b)
		b = next;
	      else break;
	    } // repeat until |break|
	  }
	  const auto& bot = tab[b];
	  const unsigned int d = elt[xs].length - elt[b].length;
	  const unsigned int m = Coxeter_entry(type,s,t);
	  const Generator st[] = {s,t}; // convenience

	  if (d == m)
	  { // case (1): there is no transduction
	    // xs.t is computed by shifting up from y the other way around
	    EltPiece y = b;
	    { // perform |m-1| upward steps, starting with |st[m%2]|
	      const Generator u=st[m%2], v=st[1-m%2];
	      unsigned int d = m-1; // number of steps left; nonzero
	      while(true)
	      {
		assert(tab[y].shift[u]>y and tab[y].shift[u]!=UndefEltPiece);
		y = tab[y].shift[u];
		if (--d==0)
		  break;
		// now do the same for the other generator
		assert(tab[y].shift[v]>y and tab[y].shift[v]!=UndefEltPiece);
		y = tab[y].shift[v];
		if (--d==0)
		  break;
	      }
	    }
	    top.shift[t] = y;
	    tab[y].shift[t] = xs;
	  }
	  else if (d == m-1)
	  {
	    Generator u = st[1-m%2];

	    if (bot.shift[u] == b)
	    { // case (2): $xs$ is fixed by $t$
	      // and upon receiving |t| outputs the same $g$ as |b| for |u|
	      top.shift[t] = xs;
	      top.out[t]   = bot.out[u];
	    }
	    // else case (3) with orbit of size $2m$, do nothing
	  }
	  else assert(d<m-1); // case (3):  $t$ takes $xs$ up, do nothing
	} // |for(t)|
      } // |if (..==UndefEltPiece)|
    // |for(s)|
  // |for(x)|
  assert(tab.size() + r < std::numeric_limits<Generator>::max());

  table=matrix::Matrix<Generator>(tab.size(),r+1);
  for (EltPiece i = 0; i < tab.size(); ++i)
    for (Generator j=0; j<=r; ++j)
      if (tab[i].out[j] == UndefGenerator)
	assert(tab[i].shift[j]!=i),table(i,j)=tab[i].shift[j];
      else
	assert(tab[i].out[j]<r),table(i,j) = tab.size() + tab[i].out[j];
} // |Transducer::Transducer|

Generator WeylGroup::Transducer::unshift(EltPiece& x) const
{
  auto s = elt[x].right;
  assert(s!=UndefGenerator); // or simpler |x>0|
  x = shift(x,s);
  assert (s<limit);
  return s;
}


WeylWord WeylGroup::Transducer::word_of_piece(EltPiece x) const
{
  unsigned int k=length(x);
  std::vector<Generator> result(k);
  while (x>0)
    result[--k] = unlocal(unshift(x));
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

  Actually our code now folds the shift and transduce values, which are
  mututally exclusive, into a single table entry of type |unsigned char|. That
  means that the coset size _plus_ the number of distinct values that can be
  output should not exceed 256 for any transducer. The latter number equals the
  number r of the transducer within its Dynkin diagram component (as it can
  output at most r-1); for the final transducer of E8 we get 240+7=247 which
  fits; in classical types the rank limit gets lowered to 86 (the final
  transducer for C86 has 2*86-1=171 states and can output vales 0..84 for
  171+85=256 possible entries. In practice this still is largely sufficent.

******************************************************************************/

/*
  Build the Weyl group corresponding to the Cartan matrix |c|,
  and incorporate the possibly given twist (take identity if |twist==nullptr|).

  NOTE : |c| and |twist| use some (consistent) labelling of simple roots,
  but we will determine an internal renumbering making the subquotients small
*/
WeylGroup::WeylGroup(const int_Matrix& c)
  : d_transducer(c.numColumns())
  , d_in()
  , d_out()
  , d_min_star(d_transducer.size())
  , upper(d_transducer.size())

{
  const unsigned int rank = d_transducer.size();

  /* analyse the Cartan matrix */

  DynkinDiagram d(c); // make diagram from Cartan matrix and classify

/*
  now put appropriate permutations into |d_in| and |d_out|, so that internal
  number |i| gives external number |d_out[i]|, and |d_in[d_out[i]]==i|.
*/
  for (const auto& comp : d.components())
  {
    Generator last=comp.offset+comp.position.size()-1;
    for (unsigned i=0; i<comp.position.size(); ++i)
      upper[comp.offset+i] = last;
    if (std::strchr("BCD",comp.type)==nullptr)
      for (unsigned i=0; i<comp.position.size(); ++i)
	d_out[comp.offset+i]=comp.position[i];
    else // reverse in types BCD
      for (unsigned i=0; i<comp.position.size(); ++i)
	d_out[last-i]=comp.position[i];
  }

  for (Generator i = 0; i < rank; ++i)
    d_in[d_out[i]] = i;

  // now construct the transducers
  // each one gets to know its local Coxeter matrix and its place in component
  for (const auto& comp : d.components())
    for (unsigned i=0; i<comp.rank(); ++i)
      d_transducer[comp.offset+i] = Transducer(comp.type,comp.offset,i);

  // precompute for each |j| the first non-commuting or equal generator |i<=j|
  for (Generator j = 0; j < rank; ++j)
  { d_min_star[j]=j; // value if loop below does not |break|
    for (Generator i=0; i<j; ++i)
      if (not inner_commutes(i,j))
      {
	d_min_star[j]=i; break;
      }
  }

  // now our Weyl group is operational for computations

} // |WeylGroup::WeylGroup|

// these methods use the default, but that requires |Transducer| type complete
WeylGroup::WeylGroup(WeylGroup&& W) = default;
WeylGroup::~WeylGroup() = default;

/******** accessors **********************************************************/

WeylElt WeylGroup::inner_gen (Generator i) const
{
  WeylElt s;    // initialise all pieces to 0 (identity)
  s.piece(i)=1; // then shift to piece #1 in $W_{i-1}\\W_i$, which is $s_i$
  return s;
}


Generator WeylGroup::output_local_gen(Generator i, Generator g) const
{
  return d_out[d_transducer[i].unlocal(g)];
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
  for (Generator j = rank(); j-->0;) {
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

  This simplification led to the following main loop:
  for (Generator t; (t=d_transducer[i].out(w.piece(i),s))!=UndefGenerator; s=t)
  {
    assert(i!=0); // since |d_transducer[0]| only has shifts
    --i;
  }

  I left the test for |UndefGenerator| coming first, as transduce is more
  frequent than shift; hence the following even simpler code could be less
  efficient:

  Generator t=s;
  for (Generator j=rank(); j-->0; t=d_transducer[j].out(w[j],t))
  {
    EltPiece wj=w[j]; w[j]=d_transducer[j].shift(wj,s);
    if (wj!=w[j]) return w[j]>wj ? 1 : -1;
  }

  Finally the tables for |out| and |shift| were folded into a single table, with
  an added predicate |has_shift| to distinguish whic case applies, and the
  |transduce| code was forced to test this first, leading to the code below.

  MvL
*/

// auxiliary; return value tells whether a descent occurred (in final shift)
inline int WeylGroup::transduce(WeylElt& w, Generator i, Generator s) const
{
  while(true) // loop breaks once we hit a possible shift
  {
    const auto wi = w.piece(i);
    if (d_transducer[i].has_shift(wi,s))
      // no need to use length for return value, numeric '<' is OK:
      return (w.piece(i)=d_transducer[i].shift(wi,s))<wi ? -1 : 1;
    assert(i>0); // since |d_transducer[0]| only has shifts
    s=d_transducer[i--].out(wi,s);
  }
}

int WeylGroup::inner_mult(WeylElt& w, Generator s) const
{
  return transduce(w,start_gen(s),d_transducer[s].local(s));
}

/*
  Multiply |w| on the right by piece |i| of |v| (in internal numbering)

  Returns nonpositive even value $l(wv)-l(w)-l(v)$

  Algorithm: do the elementary multiplication by the generators, running
  through |v| left-to-right.
*/
int WeylGroup::mult_by_piece(WeylElt& w, const WeylElt& v, Generator i) const

{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal
  EltPiece x = v.piece(i);
  const auto& tr = d_transducer[i];
  const Generator start = start_gen(i);
  int result = -static_cast<int>(tr.elt[x].length);
  // final result takes -2 for each shortening |inner_mult|, 0 otherwise

  unsigned int k=0;
  while (x>0)
    stack[k++] = tr.unshift(x);
  while (k-->0)
    result+=transduce(w,start,stack[k]);

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
  WeylElt sw=inner_gen(s);
  int l=1; //

  // now compute $sv$ as above, keeping track of any length drop (at most 1)
  for (Generator i = min_neighbor(s); i <= s; ++i)
    l+=mult_by_piece(sw,w,i);

  // and copy its relevant pieces into $w$
  for (Generator i = min_neighbor(s); i <= s; ++i)
    w.piece(i) = sw.piece(i);

  return l;
}

/*
  Multiply |w| on the right by |v|, and put the product in |w|: |w*=v|.

  Algorithm: multiply |w| by the various pieces of |v|.
*/
void WeylGroup::mult(WeylElt& w, const WeylElt& v) const
{
  for (Generator i = 0; i < rank(); ++i)
    mult_by_piece(w,v,i);
}

/*
  Multiply |w| on the right by |v|, and put the product in |w|: |w*=v|.

  Algorithm: do the elementary multiplication by the generators, running
  through v left-to-right.
*/
void WeylGroup::mult(WeylElt& w, const WeylWord& v) const
{
  for (unsigned int i = 0; i < v.size(); ++i)
    mult(w,v[i]);
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

  for (Generator i = rank(); i-->0 ;)
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];
    while (x>0)
      inner_mult(wi,tr.unlocal(tr.unshift(x))); // maybe better expand
  }

  return wi;
}

WeylElt WeylGroup::longest() const
{
  WeylElt result;
  // the longest element has the maximal valid value in each of its pieces
  for (Generator i = 0; i < rank(); ++i)
    result.piece(i) = d_transducer[i].size()-1;
  return result;
}

unsigned int WeylGroup::max_length() const
{
  unsigned int result=0;
  for (Generator i = 0; i < rank(); ++i)
    result += d_transducer[i].max_length();
  return result;
}

size::Size WeylGroup::order() const
{
  size::Size result; // this initialises the value to 1

  // and the Weyl group size is the product of the numbers of transducer states
  for (Generator i = 0; i < rank(); ++i)
    result *= d_transducer[i].size();

  return result;
}

// whether |s| and |t| are distinct and commute
bool WeylGroup::inner_commutes (Generator s, Generator t) const
{
  if (upper[s]!=upper[t]) // whether distinct diagram components
    return true;
  if (s>t)
    std::swap(s,t); // ensure |s<t|
  const Transducer& tr = d_transducer[t]; // piece 1 of |tr| represents |t|
  s = tr.local(s);
  return not tr.has_shift(1,s); // in which case always |tr.out(1,s)==s|
}

Generator WeylGroup::Chevalley_dual(Generator s) const
{
  const WeylElt w0 = longest();
  // since we dont have a root datum at hand, we conjugate by |longest()|
  WeylElt w = prod(prod(w0,s),w0);

  // |w| should be some generator |t|; find which, and store it
  Generator t;
  for (t = 0; t < rank(); ++t)
    if (w==generator(t))
      return t;
  assert(false); // we should never come here
  return UndefGenerator; // to keep compiler happy
}

Twist WeylGroup::Chevalley_twist () const
{
  Twist result;
  // since we dont have a root datum at hand, we conjugate by |longest()|
  for (Generator s = 0; s < rank(); ++s)
    result[s] = Chevalley_dual(s);
  return result;
}

/*
  Put into |c| the conjugacy class of |w|.

  Algorithm: straightforward enumeration of the connected component of |w| in
  the graph defined by the operation |conjugate|, using a |std::set| structure
  to record previously encountered elements and a |containers::queue| to store
  elements whose neighbors have not yet been generated.

*/
WeylEltList conjugacy_class(const WeylGroup& W,const WeylElt& w)
{
  std::set<WeylElt> found { w };
  containers::queue<WeylElt> to_do { w };

  while (not to_do.empty())
  {
    WeylElt v = to_do.front(); to_do.pop();

    for (Generator s = 0; s < W.rank(); ++s)
    {
      W.conjugate(v,s);
      if (found.insert(v).second) // then it was new, put it on to-do list
	to_do.push(v);
    }
  }

  // now convert set |found| to vector
  return WeylEltList(found.begin(),found.end());
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
  WeylElt sw = inner_gen(s); // element operated upon, starts out as |s|
  const Generator start=start_gen(s);

  for (Generator i = min_neighbor(s); i <= s; ++i)
  {
    EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    unsigned int k=0;
    while (x>0)
      stack[k++]=tr.unshift(x);
    while (k-->0)
      if (transduce(sw,start,stack[k])<0) // multiply and see if a descent occurs
	return true;
  }

  return false; // since only ascents occur, we have $l(sw)>l(w)$
}

// same question, but on the right; this happens to be never used in Atlas
#if 0 // remove this line (and matchinig #endif) if ever an instance is needed
bool WeylGroup::hasDescent(const WeylElt& w,Generator s) const
{
  s=d_in[s]; // inner numbering is used below
  Generator i = start_gen(s); // rightmost relevant transducer
  s = d_transducer[s].local(s); // convert |s| to index within diagram component

  while (not d_transducer[i].has_shift(w.piece(i),s))
  {
    s = d_transducer[i].out(w.piece(i),s);
    assert(i!=0); // since |d_transducer[0]| only has shifts
    --i;
  }

  EltPiece wi=w.piece(i);
  return d_transducer[i].shift(wi,s)<wi;
}
#endif

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
  for (Generator i = 0; i < rank(); ++i)
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

  for (Generator i = 0; i < rank(); ++i)
    l += d_transducer[i].length(w.piece(i));

  return l;
}

WeylWord WeylGroup::word(const WeylElt& w) const
{
  unsigned int k=length(w);
  std::vector<Generator> result(k);
  for (Generator i = rank(); i-->0; )
  {
    EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];
    while (x>0)
      result[--k] = output_local_gen(i,tr.unshift(x));
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

  for (Generator i = rank(); i-->0; )
    (u*=static_cast<int>(d_transducer[i].size()))+=static_cast<int>(w.piece(i));

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

  for (Generator i = 0; i < rank(); ++i)
    w.piece(i) = u.shift_modulo(static_cast<int>(d_transducer[i].size()));

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

  Note that unlike other places where we multiply by the letters obtained from
  repeated |unshift| of a piece of |w|, we must apply |unlocal| to them, and
  then later |local| (inside |inner_mult|), because |inner_f| uses that.
*/
WeylElt WeylGroup::translation
  (const WeylElt& w, const WeylInterface& outer_f) const
{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal

  WeylInterface inner_f;
  for (Generator i = 0; i < rank(); ++i)
    inner_f[i] = d_in[outer_f[d_out[i]]];

  WeylElt result; // start out with identity
  for (Generator i = 0; i < rank(); ++i)
  {
    EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    unsigned int k=0;
    while (x>0)
      stack[k++] = inner_f[tr.unlocal(tr.unshift(x))];
    while (k-->0)
      inner_mult(result,stack[k]);
  }

  return result;
}

// doing the same while reversing the order is easier: no stack is needed
WeylElt WeylGroup::reverse_translation
  (const WeylElt& w, const WeylInterface& outer_f) const
{
  WeylInterface inner_f;
  for (Generator i = 0; i < rank(); ++i)
    inner_f[i] = d_in[outer_f[d_out[i]]];

  WeylElt result; // start out with identity
  for (Generator i = rank(); i-->0; )
  {
    EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      inner_mult(result,inner_f[tr.unlocal(tr.unshift(x))]);
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
  for (Generator i = rank(); i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];
    while (x>0)
      rd.simple_reflect_root(output_local_gen(i,tr.unshift(x)),alpha);
  }
}

// Let |w| act on weight |v| according to reflection action in root datum |rd|
template<typename C>
  void WeylGroup::act
    (const RootDatum& rd, const WeylElt& w,  matrix::Vector<C>& v) const
{
  for (Generator i = rank(); i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];
    while (x>0)
      rd.simple_reflect(output_local_gen(i,tr.unshift(x)),v);
  }
}

// Let |w| act on coweight |v| according to reflection action in root datum |rd|
template<typename C>
  void WeylGroup::co_act
    (const RootDatum& rd,  matrix::Vector<C>& v, const WeylElt& w) const
{
  Generator stack[Transducer::max_piece_length]; // working memory for reversal

  for (Generator i = 0; i < rank(); ++i)
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    unsigned int k=0;
    while (x>0)
      stack[k++] = tr.unshift(x);
    while (k-->0)
      rd.simple_coreflect(v,output_local_gen(i,stack[k]));
  }
}

void WeylGroup::act(const RootDatum& rd, const WeylElt& w, RatWeight& v) const
{ act(rd,w,v.numerator()); }

void WeylGroup::act(const RootDatum& rd, const WeylElt& w, LatticeMatrix& M)
  const
{
  for (Generator i = rank(); i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      rd.simple_reflect(output_local_gen(i,tr.unshift(x)),M);
  }
}


template<typename C>
  void WeylGroup::act
    (const PreRootDatum& prd, const WeylElt& w, matrix::Vector<C>& v) const
{
  for (Generator i = rank(); i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      prd.simple_reflect(output_local_gen(i,tr.unshift(x)),v);
  }
}

void
  WeylGroup::act(const PreRootDatum& prd, const WeylElt& w, RatWeight& v) const
{ act(prd,w,v.numerator()); }

void WeylGroup::act(const PreRootDatum& prd, const WeylElt& w, LatticeMatrix& M)
  const
{
  for (Generator i = rank(); i-->0; )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    while (x>0)
      prd.simple_reflect(output_local_gen(i,tr.unshift(x)),M);
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
  for (Generator i = 0; i < rank(); ++i )
  { EltPiece x = w.piece(i);
    const auto& tr = d_transducer[i];

    unsigned int k=0;
    while (x>0)
      stack[k++]=tr.unshift(x);
    while (k-->0)
      rd.simple_reflect(output_local_gen(i,stack[k]),v);
  }
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
  for (Generator s = 0; s < rank(); ++s)
    dual_twist[s] = W.Chevalley_dual(d_twist[s]);
  return W.translation(w,dual_twist);
}

ext_gens TwistedWeylGroup::twist_orbits ()  const
{
  unsigned int size=0;
  for (weyl::Generator s = 0; s < rank(); ++s)
    if (twisted(s)>=s)
      ++size;

  ext_gens result; result.reserve(size);

  for (weyl::Generator s = 0; s < rank(); ++s)
    if (twisted(s)==s)
      result.push_back(ext_gen(s));
    else if (twisted(s)>s)
      result.push_back(ext_gen(W.commutes(s,twisted(s)),s,twisted(s)));

  return result;
} // |twist_orbits|


Twist TwistedWeylGroup::dual_twist() const
{
  Twist twist; // "dimensioned" but not initialised
  for (Generator s = 0;  s < W.rank(); ++s)
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

    for (Generator s = 0; s < W.rank(); ++s)
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
  for (Generator i = 0; i < rank(); ++i)
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
  for (Generator i = 0; i < rank(); ++i)
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


//			     Template instantiation

template
void WeylGroup::act
  (const RootDatum& rd, const WeylElt& w, matrix::Vector<int>& v) const;

template
void WeylGroup::co_act
  (const RootDatum& rd, matrix::Vector<int>& v, const WeylElt& w) const;


} // |namespace weyl|

} // |namespace atlas|
