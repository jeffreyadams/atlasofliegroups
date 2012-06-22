/*!
\file
\brief Implementation of the class KLContext.

  This module contains code for the computation of the Kazhdan-Lusztig
  polynomials for a given block of representations. We have taken the radical
  approach of not using the Bruhat ordering at all, just ordering by length
  instead, and coping with the ensuing appearance of zero polynomials. It
  is expected that the simplification thus achieved will more than outweigh
  the additional polynomials computed.

  The general scheme is fairly similar to the one in Coxeter: there is a
  "KLSupport" structure, that holds the list of primitive pairs that makes it
  possible to read the d_kl list, plus some additional lists that allow for
  a fast primitivization algorithm, for instance; there are two main lists,
  d_kl (filled in for all primitive pairs), and d_mu (filled in only for
  non-zero mu coefficients.)
*/
/*
  This is kl.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "kl.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <string>
#include <algorithm> // for |std::binary_search|

// #include <cassert>
#include "/usr/include/assert.h"
#include <set>
#include <map>
#include <stack>
#include <stdexcept>

#include "hashtable.h"
#include "kl_error.h"
#include "wgraph.h"	// for the |wGraph| function


// extra defs for windows compilation -spc
#ifdef WIN32
#include <iterator>
#endif

/*
  [Fokko's original description, which referred to a slighly older
  version of the computation. Fokko implemented the change from using
  extremal pairs to the slightly large set of primitive pairs, but did
  not change all of the code comments. I am trying to do that now. DV
  8/26/06.]

  This module contains code for the computation of the Kazhdan-Lusztig
  polynomials for a given block of representations. We have taken the radical
  approach of not using the Bruhat ordering at all, just ordering by length
  instead, and coping with the ensuing appearance of zero polynomials. It
  is expected that the simplification thus achieved will more than outweigh
  the additional polynomials computed.

  The general scheme is fairly similar to the one in Coxeter: there is a
  "KLSupport" structure, that holds the list of extremal pairs that makes it
  possible to read the d_kl list, plus some additional lists that allow for
  a fast extremalization algorithm, for instance; there are two main lists,
  d_kl (filled in for all extremal pairs), and d_mu (filled in only for
  non-zero mu coefficients.)
*/

namespace atlas {

namespace kl {

  /*!
\brief Polynomial 0, which is stored as a vector of size 0.
  */
  const KLPol Zero;

  /*! \brief Polynomial 1.q^0. */
  const KLPol One(0,KLCoeff(1)); // Polynomial(d,1) gives 1.q^d.

// we wrap |KLPol| into a class |KLPolEntry| that can be used in a |HashTable|

/* This associates the type |KLStore| as underlying storage type to |KLPol|,
   and adds the methods |hashCode| (hash function) and |!=| (unequality), for
   use by the |HashTable| template.
 */
class KLPolEntry : public KLPol
{
public:
  // constructors
  KLPolEntry() : KLPol() {} // default constructor builds zero polynomial
  KLPolEntry(const KLPol& p) : KLPol(p) {} // lift polynomial to this class

  // members required for an Entry parameter to the HashTable template
  typedef KLStore Pooltype;		   // associated storage type
  size_t hashCode(size_t modulus) const; // hash function

  // compare polynomial with one from storage
  bool operator!=(Pooltype::const_reference e) const;

}; // |class KLPolEntry|


namespace helper {

  typedef hashtable::HashTable<KLPolEntry,KLIndex> KLHashStore;

  class Helper
  : public KLContext
  {
    // a table for KL polynomial lookup, which is absent in KLContext
    KLHashStore d_hashtable;
    bool d_verbose;

  // Members serving for statistics

    size_t prim_size;            // number of pairs stored in d_prim
    size_t nr_of_prim_nulls;     // number of zeroes suppressed from d_prim

  public:

    // constructors and destructors
    Helper(const KLContext&); // assumes an existing bare-bones |KLContext|

    ~Helper()
    {
      if(d_verbose)
      {
	std::cout << "\nNumber of primitive pairs stored:     "
		  << prim_size << ".\n";
	std::cout << "Number of unrecorded primitive pairs: "
		  << nr_of_prim_nulls << '.' << std::endl;
      }
    }

    //accessors

    size_t firstDirectRecursion(BlockElt y) const;

    std::pair<size_t,size_t> firstImaginaryReal(BlockElt x, BlockElt y) const;

    size_t firstNiceAscent(BlockElt x,BlockElt y) const;

    BlockEltPair inverseCayley(size_t s, BlockElt y) const
    {
      return block().inverseCayley(s,y);
    }

    /* Declare that member function klPol from base class KLContext is also
       considered. This is necessary since it would otherwise be shadowed
       rather than overloaded by a new definition with different signature
       in the Helper class.
    */
    using KLContext::klPol;

    // look up a polynomial in |klv|, for |x| in range |p_begin| to |p_end|
    // this one is used by the Thicket class, which uses |KLRow| internally
    KLPolRef klPol(BlockElt x, BlockElt y,
		   KLRow::const_iterator klv,
		   PrimitiveRow::const_iterator p_begin,
		   PrimitiveRow::const_iterator p_end) const;

    // this variant is for |newRecursionRow|, which uses std::vector<KLPol>
    const KLPol& klPol(BlockElt x, BlockElt y,
		       const std::vector<KLPol>& klv,
		       const PrimitiveRow& pr,
		       unsigned int lwb // no need to look below this index
		       ) const;

    inline bool ascentMu // mu(x,y) when s is ascent for x and descent for y
      (BlockElt x, BlockElt y, size_t s) const;

    inline MuCoeff goodDescentMu
      (BlockElt x, BlockElt y, size_t s) const;

    MuCoeff lengthOneMu(BlockElt x, BlockElt y) const;

    MuCoeff type2Mu(BlockElt x, BlockElt y) const;
    // manipulators
    void rebuild_index() { d_hashtable.reconstruct(); }

    void fill(BlockElt y, bool verbose);

    void fillKLRow(BlockElt y);

    void fillMuRow(BlockElt y);

    void directRecursion(BlockElt y, size_t s);

    void recursionRow(std::vector<KLPol> & klv,
		      const PrimitiveRow& e, BlockElt y, size_t s);

    void muCorrection(std::vector<KLPol>& klv,
		      const PrimitiveRow& e,
		      BlockElt y, size_t s);

    void newRecursion(BlockElt y);

    void newRecursionRow(std::vector<KLPol> & klv,
			 const PrimitiveRow& e, BlockElt y);
    void muNewFormula(std::vector<KLPol>& klv, const PrimitiveRow& e,
		      size_t j, BlockElt y, size_t s);

    void writeRow(const std::vector<KLPol>& klv,
		  const PrimitiveRow& e, BlockElt y);


  }; // class Helper

} // namespace helper
} // namespace kl

/*****************************************************************************

        Chapter I -- Methods of the KLContext and KLPolEntry classes.

 *****************************************************************************/

/* methods of KLContext */


namespace kl {

KLContext::KLContext(const Block_base& b)
  : klsupport::KLSupport(b) // construct unfilled support object from block
  , fill_limit(0)
  , d_prim()
  , d_kl()
  , d_mu()
  , d_store(2)
  , d_zero(0)
  , d_one(1)
{
  d_store[d_zero]=Zero;
  d_store[d_one]=One;

}

/******** copy, assignment and swap ******************************************/


/*!
  \brief Copy constructor. Used only to copy-construct base for Helper class
*/
KLContext::KLContext(const KLContext& other)
  : klsupport::KLSupport(other)
  , fill_limit(other.fill_limit)
  , d_prim(other.d_prim)
  , d_kl(other.d_kl)
  , d_mu(other.d_mu)
  , d_store(other.d_store)
  , d_zero(other.d_zero)
  , d_one(other.d_one)
{}



void KLContext::swap(KLContext& other) // used to extract base out of Helper
{
  klsupport::KLSupport::swap(other); // swap base (support) objects
  std::swap(fill_limit,other.fill_limit);

  d_prim.swap(other.d_prim);

  d_kl.swap(other.d_kl);
  d_mu.swap(other.d_mu);

  d_store.swap(other.d_store);  // this puts the Helper store into base object
  std::swap(d_zero,other.d_zero);
  std::swap(d_one,other.d_one);
}

// If table extension fails, we need to recover the old part of the tables
void KLContext::partial_swap(KLContext& other)
{
  // this is only called both to copy old tables into the helper object |other|,
  // and in case of an overflow to recover those old tables form |other|
  klsupport::KLSupport::swap(other); // swap base (support) objects

  d_prim.swap(other.d_prim);
  d_kl.swap(other.d_kl);
  d_mu.swap(other.d_mu);
  d_store.swap(other.d_store);
  std::swap(d_zero,other.d_zero);
  std::swap(d_one,other.d_one);

  if (other.fill_limit<fill_limit) // we are initially filling helper object
  {
    other.fill_limit=fill_limit; // for fill limit, copy to other but keep ours
    d_store.reserve(other.d_store.size()); // a stupid way to record old size
  }
  else // then we are restoring after a throw; |other| will soon disappear
  {
    d_prim.resize(fill_limit); // truncate, freeing occupied memory
    d_kl.resize(fill_limit);
    d_mu.resize(fill_limit);
    d_store.resize(other.d_store.capacity()); // truncate to old size
  }

}

/******** accessors **********************************************************/


/*!
  \brief Returns the Kazhdan-Lusztig-Vogan polynomial P_{x,y}

  Precondition: x and y are smaller than size();

  Explanation: since d_store holds all polynomials for primitive
  pairs (x,y), this is basically a lookup function. While x is not
  primitive w.r.t. y, it moves x up, using the "easy" induction
  relations. At that point, we look x up in the primitive list for
  y. If it is found, we get a pointer to the result.  If it is not
  found, the polynomial is zero.
*/
KLPolRef KLContext::klPol(BlockElt x, BlockElt y) const
{
  const PrimitiveRow& pr = d_prim[y];
  const KLRow& klr = d_kl[y];

  x=primitivize(x,descentSet(y));
  if (x>y) return Zero; // includes case x==blocks::UndefBlock
  PrimitiveRow::const_iterator xptr =
    std::lower_bound(pr.begin(),pr.end(),x);
  if (xptr == pr.end() or *xptr != x) // not found
    return Zero;
  return d_store[klr[xptr - pr.begin()]];
}


/*!
  \brief Returns mu(x,y).

  Explanation: it is guaranteed that all the x'es such that mu(x,y) != 0
  occur in d_mu[y] (and in fact, that only those occur.) So it is a simple
  matter of looking up x.
*/
MuCoeff KLContext::mu(BlockElt x, BlockElt y) const
{
  const MuRow& mr = d_mu[y];

  std::vector<BlockElt>::const_iterator xloc=
    std::lower_bound(mr.first.begin(),mr.first.end(),x);

  if (xloc==mr.first.end() or *xloc!=x)
    return 0; // x not found in mr

  return mr.second[xloc-mr.first.begin()];
}

/* The following two methods were moved here from the Helper class, since
   they turn out to be useful even when no longer constructing the KLContext
*/

/*!
  \brief Puts in e the list of all x extremal w.r.t. y.

  Explanation: this means that either x = y, or length(x) < length(y), and
  every descent for y is a descent for x.  Or:  asc(x)\cap desc(y)=\emptyset
  Here descent means "in the tau invariant" (possibilities C-, ic, r1, r2).
*/
void
KLContext::makeExtremalRow(PrimitiveRow& e, BlockElt y)
  const
{
  BitMap b(size());
  b.fill(0,lengthLess(length(y))); // start with all elements < y in length
  b.insert(y);                     // and y itself

  // extremalize (filter out those that are not extremal)
  extremalize(b,descentSet(y)); // KLSupport::extremalize does the real work

  // copy from bitset b to list e
  e.reserve(e.size()+b.size()); // ensure tight fit after copy
  std::copy(b.begin(),b.end(),back_inserter(e));

}


/*!
  \brief Puts in e the list of all x primitive w.r.t. y.

  Explanation: this means that either x = y, or length(x) < length(y),
  and every descent for y is either a descent, or an imaginary type II
  ascent for x.
*/
void
KLContext::makePrimitiveRow(PrimitiveRow& e, BlockElt y)
  const
{
  BitMap b(size());
  b.fill(0,lengthLess(length(y))); // start with all elements < y in length
  b.insert(y);                     // and y itself

  // primitivize (filter out those that are not primitive)
  primitivize(b,descentSet(y));
  // copy to list
  e.reserve(e.size()+b.size()); // ensure tight fit after copy
  std::copy(b.begin(),b.end(),back_inserter(e));

}


/******** manipulators *******************************************************/

/*!
  \brief Fills (or extends) the KL- and mu-lists.

  Explanation: this is the main function in this module; all the work is
  deferred to the Helper class.
*/
void KLContext::fill(BlockElt y, bool verbose)
{

  if (y<fill_limit)
    return; // tables present already sufficiently large for |y|

#ifndef VERBOSE
  verbose=false; // if compiled for silence, force this variable
#endif

  if (verbose)
    std::cerr << "computing Kazhdan-Lusztig polynomials ..." << std::endl;

  helper::Helper help(*this); // make helper, with empty base KLContext
  if (fill_limit>0)
  {
    partial_swap(help);   // move previously computed tables into helper
    help.rebuild_index(); // since we just swapped the referred-to |KLStore|
  }

  try
  {
    help.fill(y,verbose); // takes care of filling embedded |KLSupport| too
    swap(help); // swap base object, including the |KLSupport| and |KLStore|

    fill_limit = y+1;

    if (verbose)
      std::cerr << "done" << std::endl;
  }
  catch (std::bad_alloc)
  { // roll back, and transform failed allocation into MemoryOverflow
    std::cerr << "\n memory full, KL computation abondoned." << std::endl;
    partial_swap(help); // restore (only) the previous contents from helper
    throw error::MemoryOverflow();
  }

}

BitMap KLContext::primMap (BlockElt y) const
{
  BitMap b(size()); // block-size bitmap

  // start with all elements < y in length
  b.fill(0,lengthLess(length(y)));
  b.insert(y);   // and y itself

  // primitivize (filter out those that are not primitive)
  primitivize(b,descentSet(y));

  // now b holds a bitmap indicating primitive elements for y

  // our result will be a bitmap of that size
  BitMap result (b.size()); // initiallly all bits are cleared

 // the list of primitive elements with nonzero polynomials at y
  const PrimitiveRow& row=d_prim[y];

  // traverse |b|, and for its elements that occur in, set bits in |result|

  size_t position=0; // position among set bits in b (avoids using b.position)
  size_t j=0; // index into row;
  for (BitMap::iterator it=b.begin(); it(); ++position,++it)
    if (*it==row[j]) // look if |it| points to current element of row
    {
      result.insert(position); ++j; // record position and advance in row
      if (j==row.size()) break;     // stop when row is exhausted
    }

  return result;
}


/* methods of KLPolEntry */


/*!
  \brief calculate a hash value in [0,modulus[, where modulus is a power of 2

  The function is in fact evaluation of the polynomial (with coefficients
  interpreted in Z) at the point 2^21+2^13+2^8+2^5+1=2105633, which can be
  calculated quickly (without multiplications) and which gives a good spread
  (which is not the case if 2105633 is replaced by a small number, because
  the evaluation values will not grow fast enough for low degree
  polynomials!).

*/
inline size_t KLPolEntry::hashCode(size_t modulus) const
{ const KLPol& P=*this;
  if (P.isZero()) return 0;
  polynomials::Degree i=P.degree();
  size_t h=P[i]; // start with leading coefficient
  while (i-->0) h= ((h<<21)+(h<<13)+(h<<8)+(h<<5)+h+P[i]) & (modulus-1);
  return h;
}

bool KLPolEntry::operator!=(KLPolEntry::Pooltype::const_reference e) const
{
  if (degree()!=e.degree()) return true;
  if (isZero()) return false; // since degrees match
  for (polynomials::Degree i=0; i<=degree(); ++i)
    if ((*this)[i]!=e[i]) return true;
  return false; // no difference found
}

} // namespace kl


/*****************************************************************************

        Chapter II -- Methods of the Helper class.

 *****************************************************************************/

namespace kl {
  namespace helper {

    // main constructor

Helper::Helper(const KLContext& kl)
  : KLContext(kl.block())  // re-construct empty KLContext base object
  , d_hashtable(d_store) // hash table identifies helper base object's d_store
  , d_verbose(false)
  , prim_size(0), nr_of_prim_nulls(0)
{
}

/******** accessors **********************************************************/


/*!
  \brief Returns the first descent generator that is not real type II

  Explanation: those are the ones that give a direct recursion formula for the
  K-L basis element. Explicitly, we search for a generator |s| such that
  |descentValue(s,y)| is either |DescentStatus::ComplexDescent| or
  |DescentStatus::RealTypeI|. If no such generator exists, we return |rank()|.
*/
size_t Helper::firstDirectRecursion(BlockElt y) const
{
  const DescentStatus& d = descent(y);

  for (size_t s = 0; s < rank(); ++s) {
    DescentStatus::Value v = d[s];
    if (DescentStatus::isDirectRecursion(v))
      return s;
  }

  return rank();
} // |Helper::firstDirectRecursion|

/*
  Returns the first real nonparity ascent for y that is a complex ascent, or
  imaginary type 2, or compact imaginary for x.

  Explanation: those are the ones that give the best new recursion formula for
  the K-L polynomial

  If no such generator exists, we return |rank()|.
*/
size_t Helper::firstNiceAscent(BlockElt x,BlockElt y) const
{
  const DescentStatus& dx = descent(x);
  const DescentStatus& dy = descent(y);

  for (size_t s = 0; s < rank(); ++s)
    if (dy[s]==DescentStatus::RealNonparity)
      { DescentStatus::Value vx = dx[s];
	if (vx==DescentStatus::ComplexAscent or
	    vx==DescentStatus::ImaginaryCompact or
	    vx==DescentStatus::ImaginaryTypeII)
	  return s;
      }
  return rank(); // failure
} // |Helper::firstNiceAscent|

/*
  Preconditions:
  * all descents for y are of type r2 (firstDirectRecursion has failed)
  * x is extremal for y (none of the descents for y are ascents for x)
  * none of the rn ascents for y is C+, ic or i2 for x (firstNiceAscent failed)

  Returns the first pair (s,t) such that
  1) (s,t) is (rn,r2) for y;
  2) (s,t) is (i1,ic) for x;
  3) (s,t) is (i1,i1/2) for s.x.
  Since the statuses of t for x and s.x differ, s must be ajdacent to t in the
  Dynkin diagram (but this is not tested or used explicitly here)

  Such a pair can be used to compute P_{x,y}+P_{sx,y} using s,
  then P_{sx,y} using t, and so P_{x,y}.

  The test that t is ic for x is omitted, since x and s.x being related by an
  imaginary cross action are in the same fiber, so t is imaginary for x if it
  is so for s.x, and noncompact for x is ruled out by extremality precondition

  If no such pair exists or (rank,*) is returned. Under the given
  preconditions, such failure allows concluding that P_{x,y}=0.
*/
std::pair<size_t,size_t> Helper::firstImaginaryReal(BlockElt x,
						    BlockElt y) const
{
  const DescentStatus& dx = descent(x);
  const DescentStatus& dy = descent(y);

  size_t s,t, r=rank();

  for (s =0; s<r; ++s)
    if (dy[s]==DescentStatus::RealNonparity and
	dx[s]==DescentStatus::ImaginaryTypeI)
    { BlockElt sx = cross(s,x);
      const DescentStatus& dsx = descent(sx);
      for (t = 0; t<r; ++t)
	if (dy[t]==DescentStatus::RealTypeII)
	  if (dsx[t]==DescentStatus::ImaginaryTypeI or
	      dsx[t]==DescentStatus::ImaginaryTypeII)
	    return std::make_pair(s,t);
    }
  return std::make_pair(r,0); // failure

} // |Helper::firstImaginaryReal|

/*!
  \brief Returns the Kazhdan-Lusztig polynomial for x corresponding to
  the given row.

  Precondition: |klv| holds the tail of the set of primitive Kazhdan-Lusztig
  polynomials for |y|, enough to find the required one by elementary lookup;
  |[p_begin,p_end[| is the corresponding range of primitive elements.

  Algorithm: primitivize |x| with respect to the descents in |y|; if a real
  nonparity situation is encountered, return |Zero|; otherwise look up the
  primitive |x| in the range and return the corresponding element from |klv|
*/
KLPolRef Helper::klPol(BlockElt x, BlockElt y,
		       KLRow::const_iterator klv,
		       PrimitiveRow::const_iterator p_begin,
		       PrimitiveRow::const_iterator p_end) const
{
  BlockElt xp = primitivize(x,descentSet(y));

  if (xp>y) return Zero; // includes case |xp==blocks::UndefBlock|
  PrimitiveRow::const_iterator xptr =
    std::lower_bound(p_begin,p_end,xp);
  if (xptr == p_end or *xptr != xp)
    return Zero; // missing may mean equal length, or computed but found Zero
  return d_store[klv[xptr-p_begin]];
}

const KLPol& Helper::klPol(BlockElt x, BlockElt y,
			   const std::vector<KLPol>& klv,
			   const PrimitiveRow& pr,
			   unsigned int lwb // no need to look below this index
			   ) const
{
  BlockElt xp = primitivize(x,descentSet(y));
  if (xp>y or (length(xp)==length(y) and xp!=y))
    return Zero; // includes case |xp==blocks::UndefBlock|

  PrimitiveRow::const_iterator xptr =
    std::lower_bound(pr.begin()+lwb,pr.end(),xp);
  assert(xptr != pr.end() and *xptr == xp); // here all primitives are stored
  return klv[xptr-pr.begin()];
}



/*!
  \brief Computes whether |mu(x,y)==1| in a good ascent situation

  Preconditions: $l(y)>0$; $l(x)=l(y)-1$; $s$ is an ascent for $x$ w.r.t. $y$

  Explanation: this is the situation where |mu(x,y)| is directly expressible.
  The point is that |x| will ascend to an element of the same length as |y|,
  so that the corresponding K-L polynomial is zero unless |x| ascends to |y|.
*/
inline bool
Helper::ascentMu(BlockElt x, BlockElt y, size_t s) const
{
  switch (descentValue(s,x))
  {
  case DescentStatus::ComplexAscent:   return cross(s,x) == y;
  case DescentStatus::RealNonparity:   return false;
  case DescentStatus::ImaginaryTypeI:  return cayley(s,x).first == y;
  case DescentStatus::ImaginaryTypeII: return inverseCayley(s,y).first == x;
  default: // this cannot happen
    assert(false);
    return false; // keep compiler happy
  }
}

/*!
  \brief Computes $\mu(x,y)$ in the special case that $l(y)-l(x) = 1$.

  Preconditions: $l(y) > 0$; $l(x) = l(y)-1$; the $\mu$-rows for all $y$
  of smaller lengths have already been filled in.

  Explanation: these are the $\mu$-values that can be nonzero without being in
  an extremal situation. The value can be obtained in all cases as the
  constant term of |klPol(x,y)|, but if x is not extremal it is easier to
  check, using ascentMu, whether y is in some $s$-up-set of x
*/
MuCoeff Helper::lengthOneMu(BlockElt x, BlockElt y) const
{
  unsigned int s = ascent_descent(x,y);
  if (s<rank())
    return ascentMu(x,y,s) ? MuCoeff(1) : MuCoeff(0);

  // For extremal x, just lookup the KL polynomial
  KLPolRef p=klPol(x,y);
  return p.isZero() ? MuCoeff(0) : p[0];
}


/*!
  \brief Gets |mu(x,y)| by a good descent recursion.

  Precondition: $l(y)>0$; $l(x)=l(y)-1$; all previous mu-rows have been filled
  in; |s| is a good descent for |y|;

  Explanation: a good descent for |y| is a descent that is neither real type
  II nor imaginary compact (so it is either a complex or real type I); these
  are the cases where |mu(x,y)| is directly expressed by recursion.
*/
inline MuCoeff
Helper::goodDescentMu(BlockElt x, BlockElt y, size_t s) const
{
  BlockElt y1;

  if (descentValue(s,y) == DescentStatus::ComplexDescent)
    y1 = cross(s,y);
  else // s is real type I for y, chooise one of double-valued inverse Cayley
    y1 = inverseCayley(s,y).first;

  switch (descentValue(s,x))
  {
  case DescentStatus::ImaginaryCompact:  return MuCoeff(0);
  case DescentStatus::ComplexDescent:    return mu(cross(s,x),y1);
  case DescentStatus::RealTypeI:
    {
      BlockEltPair x1 = inverseCayley(s,x);
      return mu(x1.first,y1)+mu(x1.second,y1);
    }
  case DescentStatus::RealTypeII:  return mu(inverseCayley(s,x).first,y1);

  default: // this cannot happen
    assert(false);
    return MuCoeff(0); // keep compiler happy
  }
}


/******** manipulators *******************************************************/


/*!
  \brief Dispatches the work of filling the KL- and mu-lists.

  The current implementation blindly assumes that the tables are empty
  initially. While this is true in the way it is now used, it should really be
  rewritten to take into account the possibility of initial |fill_limit>0|
*/
void Helper::fill(BlockElt last_y, bool verbose)
{
  // make sure the support (base of base) is filled
  klsupport::KLSupport::fill();

  // resize the outer lists to the block size
  d_prim.resize(size());
  d_kl.resize(size());
  d_mu.resize(size());

  // fill the lists
  size_t minLength = length(0);
  size_t maxLength = length(last_y<size() ? last_y : size()-1);

  // do the minimal length cases; they come first in the enumeration
  for (BlockElt y = 0; y < lengthLess(minLength+1); ++y) {
    d_prim[y].push_back(y); // singleton list for this row
    ++prim_size;
    // the K-L polynomial is 1
    d_kl[y].push_back(d_one);
    // there are no mu-coefficients
  }

  //set timers for KL computation
  d_verbose=verbose; // inform our destructor of |verbose| setting
  std::time_t time0;
  std::time(&time0);
  std::time_t time;

  std::ifstream statm("/proc/self/statm"); // open file for memory status

  // do the other cases
  for (size_t l=minLength+1; l<=maxLength; ++l)
  {
    BlockElt y_limit = l<maxLength ? lengthLess(l+1) : last_y+1;
    for (BlockElt y=lengthLess(l); y<y_limit; ++y)
    {
      if (verbose)
	std::cerr << y << "\r";

      try
      {
	fillKLRow(y);
      }
      catch (kl_error::KLError& e)
      {
	std::ostringstream os;
	os << "negative coefficient in P_{" << e.x << ',' << e.y
	   << "} at line " << e.line << '.';
	throw std::runtime_error(os.str()); // so that realex may catch it
      }
      fillMuRow(y);
    }

    // now length |l| is completed
    if (verbose)
    {
      size_t p_capacity // currently used memory for polynomials storage
	=d_hashtable.capacity()*sizeof(KLIndex)
	+ d_store.capacity()*sizeof(KLPol);
      for (size_t i=0; i<d_store.size(); ++i)
	p_capacity+= (d_store[i].degree()+1)*sizeof(KLCoeff);

      std::time(&time);
      double deltaTime = difftime(time, time0);
      std::cerr << "t="    << std::setw(5) << deltaTime
		<< "s. l=" << std::setw(3) << l // completed length
		<< ", y="  << std::setw(6)
		<< lengthLess(l+1)-1 // last y value done
		<< ", polys:"  << std::setw(11) << d_store.size()
		<< ", pmem:" << std::setw(11) << p_capacity
		<< ", mat:"  << std::setw(11) << prim_size
		<<  std::endl;

      unsigned int size, resident,share,text,lib,data;
      statm.seekg(0,std::ios_base::beg);
      statm >> size >> resident >> share >> text >> lib >> data;
      std::cerr << "Current data size " << data*4 << "kB (" << data*4096
		<< "), resident " << resident*4 << "kB, total "
		<< size*4 << "kB.\n";
    }

  } // for (l=min_length+1; l<=max_Length; ++l)

  if (verbose)
  {
    std::time(&time);
    double deltaTime = difftime(time, time0);
    std::cerr << std::endl;
    std::cerr << "Total elapsed time = " << deltaTime << "s." << std::endl;
    std::cerr << d_store.size() << " polynomials, "
	      << prim_size << " matrix entries."<< std::endl;

    std::cerr << std::endl;
  }

}


/*!
  \brief Fills in the row for y in the KL-table.

  Precondition: all lower rows have been filled; y is of length > 0;
  R-packets are consecutively numbered;

  Explanation: this function actually fills out the "row" of y (the set of
  all P_{x,y} for x<y, which is actually more like a column)
*/
void Helper::fillKLRow(BlockElt y)
{
  if (d_kl[y].size()>0) // then row has already been filled
    return;
  size_t s = firstDirectRecursion(y);
  if (s<rank())  // a direct recursion was found, use it for y, for all x
    directRecursion(y,s);
  else
    newRecursion(y); // other cases, with formulas depending on position of x
}

/*!
  \brief Fills in the row for y using a direct recursion.

  Precondition: s is either a complex, or a real type I descent generator for
  y.

  The real work is done by the recursionRow function, that will be used also
  in the real type II descents.
*/
void Helper::directRecursion(BlockElt y, size_t s)
{
  std::vector<KLPol> klv;
  PrimitiveRow e;

  // put result of recursion formula in klv
  makeExtremalRow(e,y);
  recursionRow(klv,e,y,s);

  // write result
  writeRow(klv,e,y);

}

/*!
  \brief Puts into klv the right-hand side of the recursion formula for y
  corresponding to the descent s.

  Precondition: s is either a complex, or a real type I descent for y.

  Explanation: the shape of the formula is:

    P_{x,y} = (c_s.c_{y1})-part - correction term

  where y1 = cross(s,y) when s is complex for y, one of the two elements in
  inverseCayley(s,y) when s is real. The (c_s.c_{y1})-part depends on the
  status of x w.r.t. s (we look only at extremal x, so we know it is a
  descent). The correction term, coming from $\sum_z mu(z,y1)c_z$, is handled
  by |muCorrection|; the form of the summation depends only on |y1| (which it
  recomputes), but involves polynomials $P_{x,z}$ that depend on $x$ as well.
*/
void Helper::recursionRow(std::vector<KLPol>& klv,
			  const PrimitiveRow& e,
			  BlockElt y,
			  size_t s)
{
  klv.resize(e.size());

  BlockElt sy =
    descentValue(s,y) == DescentStatus::ComplexDescent ? cross(s,y)
    : inverseCayley(s,y).first;  // s is real type I for y here, ignore .second

  size_t i; // keep outside for error reporting
  try {

    // last K-L polynomial is 1
    klv.back() = One;

    // the following loop could be run in either direction: no dependency.
    // however it is natural to take x descending from y (exclusive)
    for (i=klv.size()-1; i-->0; )
    {
      BlockElt x = e[i];
      switch (descentValue(s,x))
      {
      case DescentStatus::ImaginaryCompact:
	{ // (q+1)P_{x,sy}
	  klv[i] = klPol(x,sy);
	  klv[i].safeAdd(klv[i],1);
	}
	break;
      case DescentStatus::ComplexDescent:
	{ // P_{sx,sy}+q.P_{x,sy}
	  BlockElt sx = cross(s,x);
	  klv[i] = klPol(sx,sy);
	  klv[i].safeAdd(klPol(x,sy),1);
	}
	break;
      case DescentStatus::RealTypeI:
	{ // P_{sx.first,sy}+P_{sx.second,sy}+(q-1)P_{x,sy}
	  BlockEltPair sx = inverseCayley(s,x);
	  klv[i] = klPol(sx.first,sy);
	  klv[i].safeAdd(klPol(sx.second,sy));
	  klv[i].safeAdd(klPol(x,sy),1);
	  klv[i].safeSubtract(klPol(x,sy));
	}
	break;
      case DescentStatus::RealTypeII:
	{ // P_{x_1,y_1}+qP_{x,sy}-P_{s.x,sy}
	  BlockElt sx = inverseCayley(s,x).first;
	  klv[i] = klPol(sx,sy);
	  klv[i].safeAdd(klPol(x,sy),1);
	  klv[i].safeSubtract(klPol(cross(s,x),sy));
	}
	break;
      default: // this cannot happen
	assert(false);
	break;
      }
    }
  }
  catch (error::NumericUnderflow& err){
    throw kl_error::KLError(e[i],y,__LINE__,
			    static_cast<const KLContext&>(*this));
  }

  // do mu-correction
  muCorrection(klv,e,y,s);

} // |Helper::recursionRow|

/*!
  \brief Subtracts from all polynomials in |klv| the correcting terms in the
  K-L recursion.

  Precondtion: |klv| already contains, for all $x$ that are primitive w.r.t.
  |y| in increasing order, the terms in $P_{x,y}$ corresponding to
  $c_s.c_{y1}$, whery |y1| is $s.y$ if |s| is a complex descent, and |y1| is
  an inverse Cayley transform of |y| if |s| is real type I or II.
  The mu-table and KL-table have been filled in for elements of length < l(y).

  Explanation: the recursion formula is of the form:
  $$
    lhs = c_s.c_{y1} - \sum_{z} mu(z,y1)c_z
  $$
  where |z| runs over the elements $< y1$ such that |s| is a descent for |z|.
  Here $lhs$ stands for $c_y$ when |s| is a complex descent or real type I for
  |y|, and for $c_{y}+c_{s.y}$ when |s| is real type II; however it plays no
  part in this function that only subtracts $\mu$-terms.

  We construct a loop over |z| first, before traversing |klv| (the test for
  $z<y1$ is absent, but $\mu(z,y1)\neq0$ implies $z\leq y1$, and $z=y1$ will
  be rejected by the descent condition). The chosen loop order allows fetching
  $\mu(z,y1)$ only once, and terminating the scan of |klv| once its values |x|
  become too large to produce a non-zero $P_{x,z}$.
*/
void Helper::muCorrection(std::vector<KLPol>& klv,
			  const PrimitiveRow& e,
			  BlockElt y, size_t s)
{
  BlockElt sy =
    descentValue(s,y) == DescentStatus::ComplexDescent ? cross(s,y)
    : inverseCayley(s,y).first;  // s is real type I for y here, ignore .second

  const MuRow& mrow = d_mu[sy];
  size_t l_y = length(y);

  for (size_t i = 0; i<mrow.first.size(); ++i)
  {

    BlockElt z = mrow.first[i];
    size_t l_z = length(z);

    DescentStatus::Value v = descentValue(s,z);
    if (not DescentStatus::isDescent(v))
      continue;

    MuCoeff mu = mrow.second[i]; // mu!=MuCoeff(0)

    polynomials::Degree d = (l_y-l_z)/2; // power of q used in the loops below

    size_t j; // define outside for error reporting
    if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
      try {
	for (j = 0; j < e.size(); ++j)
	{
	  BlockElt x = e[j];
	  if (length(x) > l_z) break; // once reached, no more terms for |z|

	  KLPolRef pol = klPol(x,z);
	  klv[j].safeSubtract(pol,d); // subtract q^d.P_{x,z} from klv[j]
	} // for (j)
      }
      catch (error::NumericUnderflow& err){
	throw kl_error::KLError(e[j],y,__LINE__,
				static_cast<const KLContext&>(*this));
      }
    else // mu!=MuCoeff(1)
      try {
	for (j = 0; j < e.size(); ++j)
	{
	  BlockElt x = e[j];
	  if (length(x) > l_z) break; // once reached, no more terms for |z|

	  KLPolRef pol = klPol(x,z);
	  klv[j].safeSubtract(pol,d,mu); // subtract q^d.mu.P_{x,z} from klv[j]
	} // for {j)
      }
      catch (error::NumericUnderflow& err){
	throw kl_error::KLError(e[j],y,__LINE__,
				static_cast<const KLContext&>(*this));
      }

  } // for (i)

} // |Helper::muCorrection|


/*!
  \brief Stores into |klv[j]| the $\mu$-sum appearing a new K-L recursion.

  Precondition: |pr| is the primitive row for |y|, $s$ is real nonparity for
  $y$ and either C+ or imaginary for $x=pr[j]$ (those are the cases for which
  the formula is used; the status w.r.t. $x$ is not actually used by the
  code), and for all $k>j$ one already has stored $P_{pr[k],y}$ in |klv[k]|.

  The mu-table and KL-table have been filled in for elements of length < l(y),
  so that for $z<y$ we can call |klPol(x,z)|.

  Explanation: the various recursion formulas involve a sum:
  $$
    \sum_{x<z<y} mu(z,y) q^{(l(y)-l(z)+1)/2}P_{x,z}
  $$
  where in addition to the condition given, |s| must be a descent for |z|.
  Since mu(z,y) cannot be nonzero unless z is primitive (indeed unless it is
  either extremal or primitive and of length l(y)-1, since in all other cases
  the recursion formulas show that $P_{z,y}$ cannot attain the maximal
  authorised degree $(l(y)-l(z)-1)/2$), so we can loop over elements of |pr|.

  We construct a loop over |z|.  The test for
  $z<y$ is absent, but $\mu(z,y)\neq0$ implies $z\leq y$. The chosen
  loop order allows fetching
  $\mu(z,y)$ only once, and terminating the scan of |klv| once its values |x|
  become too large to produce a non-zero $P_{x,z}$.

  We can't use d_mu(y), which hasn't yet been written, so mu(z,y) is extracted
  manually from the appropriate klv[k].
*/
void Helper::muNewFormula(std::vector<KLPol>& klv, const PrimitiveRow& pr,
			  size_t j, BlockElt y, size_t s)

{
  klv[j]=Zero;

  size_t l_y = length(y);
  // should iterate over z in the extremal row pr for y, maybe
  // decreasing z, only odd length differences, stopping above length(x)

  BlockElt x = pr[j];

  try
  {
    for (size_t k = pr.size()- 2; length(pr[k]) > length(x) ; --k)
    {
      BlockElt z = pr[k];
      size_t l_z = length(z);

      unsigned int d2=l_y - l_z +1; // twice the exponent of $q$ in the formula
      if (d2 % 2 !=0) continue; // which must be even

      DescentStatus::Value v = descentValue(s,z);
      if (not DescentStatus::isDescent(v))  continue;

      unsigned int d = d2/2; // power of q used in the
      KLPolRef mupol = klv[k]; // this fetches P_{z,y}
      if (mupol.degree() != d-1) continue;

      // now we have a contribution with nonzero $\mu$
      MuCoeff mu = mupol[d-1];
      KLPolRef pol = klPol(x,z);

      if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
	klv[j].safeAdd(pol,d); // add q^d.P_{x,z} to klv[j]
      else // mu!=MuCoeff(1)
	klv[j].safeAdd(pol,d,mu); // add q^d.mu.P_{x,z} to klv[j]

    } // for (k)
  }
  catch (error::NumericOverflow& e){
    throw kl_error::KLError(x,y,__LINE__,
			    static_cast<const KLContext&>(*this));
  }
} //muNewFormula

/*!
  \brief Fills in the row for y in the absence of complex descents for y.

  Precondition: every simple for y is a complex ascent, or imaginary, or real

  The real work is done by the newRecursionRow function.
*/
void Helper::newRecursion(BlockElt y)
{
  std::vector<KLPol> klv;
  PrimitiveRow pr;

  // list in pr the x primitive for y
  // (any ascents for x that are descents for y must be imaginary type II)
  makePrimitiveRow(pr,y);
  // put result of recursion formula in klv
  newRecursionRow(klv,pr,y);

  // write result
  writeRow(klv,pr,y);

} // |Helper::newRecursion|


/*!
  \brief Puts in klv[i] the right-hand side of a recursion formula for
  P_{e[i],y} corresponding to some complex ascent s for x=e[i] that is real
  (necessarily nonparity) for y.

  Precondition: every simple root is for y either a complex ascent or
  imaginary or real (no complex descents for y). (split 1)

  Explanation: the shape of the formula is roughly

    P_{x,y} = P_{x',y} +  correction term

  where x' is an ascent of x. The correction term, coming from $\sum_z
  mu(z,y1)c_z$, is handled by |newMuCorrection|.
*/

// called only when y corresponds to theta-stable q with l split.
void Helper::newRecursionRow(std::vector<KLPol>& klv,
			     const PrimitiveRow& pr,
			     BlockElt y)
{
  klv.resize(pr.size());
  // last K-L polynomial is 1
  klv.back() = One;
  size_t j = klv.size()-1; ; // declare outside try block for error reporting
  try {
    while (j-->0) {
      BlockElt x = pr[j];
      unsigned int s= ascent_descent(x,y);
      if (s<rank()) // a primitive element that is not extremal; easy case
      { // equation (1.9b) in recursion.pdf
	assert(descentValue(s,x)==DescentStatus::ImaginaryTypeII);
	BlockEltPair p = cayley(s,x);
	klv[j] = klPol(p.first,y,klv,pr,j+1);
	klv[j].safeAdd(klPol(p.second,y,klv,pr,j+1));
	continue; // done with |x|, go on to the next
      }

      // now x is extremal for y. By (split 1) and Lemma 3.1 of recursion.pdf
      // this implies that if x<y in the Bruhat order, there is at least one s
      // real for y that is a true ascent (not rn) of x and therefore rn for y
      // we first hope that at least one of them is not i1 for x

      // we first seek a real nonparity ascent for y that is C+,ic or i2 for x
      s = firstNiceAscent(x,y);
      if (s < rank()) // there is such an ascent s
      {
	// start setting klv[j] to the expression (3.4) in recursion.pdf
	muNewFormula(klv,pr,j,y,s);

	if (descentValue(s,x)==DescentStatus::ComplexAscent)
	{ // use equations (3.3a)=(3.4)
	  BlockElt sx = cross(s,x);
	  klv[j].safeSubtract(klPol(sx,y,klv,pr,j+1),1);
	  // subtract qP_{sx,y} from mu terms
	} // ComplexAscent case

	else if (descentValue(s,x)==DescentStatus::ImaginaryTypeII)
	{ // use equations (3.3a)=(3.4)
	  BlockEltPair p = cayley(s,x);
	  KLPol pol = klPol(p.first,y,klv,pr,j+1);
	  pol.safeAdd(klPol(p.second,y,klv,pr,j+1));
	  klv[j].safeAdd(pol);
	  klv[j].safeSubtract(pol,1); //now we've added (1-q)(P_{sx,y}
	  klv[j].safeDivide(2);   //this may throw
	} // ImaginaryTypeII case

	else if (descentValue(s,x)==DescentStatus::ImaginaryCompact)
	{ // here s is a emph{descent} for x, which causes an extra unknown
	  // leading (if nonzero) term to appear in addition to (3.4), giving
	  // rise to equation (3.7). Yet we can determine the quotient by q+1.
	  klv[j].safeQuotient(length(y)-length(x));
	} // ImaginaryCompact case

	// no 'else'; we've handled all possible NiceAscents
      } // NiceAscent case

      else
      {
	/* just setting klv[j]=Zero; won't do here, even in C2. We need to use
	   idea on p. 8 of recursion.pdf. This means: find s and t, both real
	   for y and imaginary for x, moreover repectively nonparity and
	   parity (r2) for y, repectively i1 and compact for x, while t is
	   noncompact for s.x (the imaginary cross of x), which implies t is
	   adjacent to s. Then we can compute P_{s.x,y} using t (an easy
	   recursion, (1.9) but for t, expresses it as sum of one or two
	   already computed polynomials), and for the sum P_{sx,y}+P_{x,y} we
	   have a formula (3.6) of the kind used for NiceAscent, and it
	   suffices to subtract P_{s.x,y} from it.

	   If no such s,t exist then we may conclude x is not Bruhat below y,
	   so P_{x,y}=0.
	*/
	std::pair<size_t,size_t> st = firstImaginaryReal(x,y);
	s = st.first;
	unsigned int t = st.second;

	if (s < rank())
	{
	  muNewFormula(klv,pr,j,y,s);

	  //subtract (q-1)P_{xprime,y} from terms of expression (3.4)
	  BlockElt xprime = cayley(s,x).first;
	  const KLPol& P_xprime_y =  klPol(xprime,y,klv,pr,j+1);
	  klv[j].safeAdd(P_xprime_y);
	  klv[j].safeSubtract(P_xprime_y,1);

	  //now klv[j] holds P_{x,y}+P_{s.x,y}
	  //compute P_{s.x,y} using t

	  BlockEltPair sx_up_t = cayley(t,cross(s,x));

	  klv[j].safeSubtract(klPol(sx_up_t.first,y,klv,pr,j+1));
	  if (sx_up_t.second != blocks::UndefBlock)
	    klv[j].safeSubtract(klPol(sx_up_t.second,y,klv,pr,j+1));
	}
	else // |firstImaginaryReal| found nothing
	  klv[j]=Zero;
      } // end of no NiceAscent case
    } // while (j-->0)
  }
  catch (error::NumericUnderflow& err) // repackage error, reporting x,y
  {
    throw kl_error::KLError(pr[j],y,__LINE__,
			    static_cast<const KLContext&>(*this));
  }
} // |Helper::newRecursionRow|


/*!
  \brief Writes down row y in d_kl and d_prim.

  Precondition: The parallel pair (klv,er) records the polynomials for at
  least all the extremal values $x$ for $y$, and at most for all primitive
  values $x$ for $y$. So when $x=er[i]$ then $P_{x,y}=klv[i]$. In practice
  |er| will contain either all extremal elements (when called from
  |directRecursion|) or all primitive elements (for |newRecursion|).

  This function writes out these data to |d_prim[y]| and |d_kl[y]|,
  transformed as follows: (1) rather than storing polynomials from |klv| (or
  others computed here), these are looked up in |d_hashtable| and the index is
  stored; (2) when a polynomial turns out to be 0, nothing is recorded either
  in |d_prim[y]| or in |d_kl[y]|; (3) for primitive elements not present in
  |er|, the polynomial is computed here on-the-fly (using an imaginary type II
  ascent that exists in this case) and then stored along with those from |er|.

  Case (3) will not apply if |er| already contains all primitive elements, and
  for that case this function could be considerably simplified, but it works
  well as is, so we didn't write a simplified version.
 */
void Helper::writeRow(const std::vector<KLPol>& klv,
		      const PrimitiveRow& er, BlockElt y)
{
  PrimitiveRow pr;
  makePrimitiveRow(pr,y);

  KLRow klr(pr.size()); // nonzero primitive entries (indexes into d_store)
  PrimitiveRow nzpr(pr.size()); // columns of nonzero prim. entries
  KLRow::iterator new_pol = klr.end();
  PrimitiveRow::iterator new_nzpr = nzpr.end();
  const PrimitiveRow::iterator nzpr_end = nzpr.end();

  // set stops in pr at elements that occur in er
  std::vector<size_t> stop(er.size()+1);
  stop[0] = 0;

  for (size_t j = 0; j < er.size(); ++j) {
    size_t epos = std::lower_bound(pr.begin(),pr.end(),er[j]) - pr.begin() + 1;
    stop[j+1] = epos;
  }

  // write polynomials, beginning with largest x (which is y).
  for (size_t j = er.size(); j-->0;)
  {
    // insert extremal polynomial for er[j]
    if (not klv[j].isZero()) {
      *--new_nzpr = er[j];
      *--new_pol  = d_hashtable.match(klv[j]); // look up klv[j], store index
    }
    // do the others (primitive but not extremal ones), down to |stop[j]|
    for (size_t i = stop[j+1]-1; i--> stop[j];)
    {
      unsigned int s = ascent_descent(pr[i],y);
      assert (s<rank()); // not extremal
      BlockEltPair xs = cayley(s,pr[i]);
      assert(xs.second != blocks::UndefBlock); // must be imaginary type II
      KLPol pol = klPol(xs.first,y,new_pol,new_nzpr,nzpr_end);
      pol.safeAdd(klPol(xs.second,y,new_pol,new_nzpr,nzpr_end));

      if (not pol.isZero()) // it could be primitive non-extremal yet zero
      {
	*--new_nzpr = pr[i];
	*--new_pol  = d_hashtable.match(pol);
      }
    }
  }

  // commit
  d_prim[y].reserve(klr.end() - new_pol);
  d_kl[y].reserve(klr.end() - new_pol);

  copy(new_nzpr,nzpr.end(),back_inserter(d_prim[y]));
  copy(new_pol,klr.end(),back_inserter(d_kl[y]));

  prim_size        += nzpr.end()-new_nzpr;
  nr_of_prim_nulls += new_nzpr  -nzpr.begin(); // measure unused space

} // |Helper::writeRow|

/*!
  \brief Fills in the row for y in the mu-table.

  Precondition: the row for y in the KL-table has been filled; length(y) > 0;

  Explanation: for the elements of length < length(y) - 1, mu(x,y) can
  be non-zero only if x is extremal w.r.t. y; so we run through d_kl[y],
  and look at the cases where the polynomial is of degree (1/2)(l(y)-l(x)-1)
  (the biggest possible). For the elements of colength 1, in the classical
  case we always had mu(x,y)=1. Here that's not true anymore, zero might
  be a possibility, and also larger values I believe; but at any rate the
  mu-values can be computed in terms of other mu-values.

  NOTE: we are not using the hasse-list here, although it is probably a
  good idea to compute that; that will reduce the mu-computation.
*/
void Helper::fillMuRow(BlockElt y)
{
  const PrimitiveRow& e = d_prim[y]; // list of nonzero polynomials

  size_t ly = length(y);
  if (ly==0) // we are in fact never called for |y| values of length 0
    return;  // but this is prudent, since next loop would fail for |ly==0|

  PrimitiveRow::const_iterator start= e.begin();
  // traverse lengths of opposite parity, up to ly-3
  for (size_t lx=(ly-1)%2,d = (ly-1)/2; d>0; --d,lx+=2) {// d=(ly-1-lx)/2

    PrimitiveRow::const_iterator stop =
      std::lower_bound(start,e.end(),lengthLess(lx+1));
    for (start= std::lower_bound(start,stop,lengthLess(lx));
	 start<stop; ++start) {
      BlockElt x = *start;
      KLIndex klp = d_kl[y][start-e.begin()];

      KLPolRef p=d_store[klp];
      if (p.degree()== d) { // then we have found a mu-coefficient for x
	d_mu[y].first.push_back(x);
	d_mu[y].second.push_back(p[d]);
      }
    }
  }

  // do cases of length ly-1
  BlockElt x_begin = lengthLess(ly-1);
  BlockElt x_end = lengthLess(ly);

  for (BlockElt x = x_begin; x < x_end; ++x) {
    MuCoeff mu = lengthOneMu(x,y);
    if (mu!=MuCoeff(0)) {
      d_mu[y].first.push_back(x);
      d_mu[y].second.push_back(mu);
    }
  }

} // |Helper::fillMuRow|

} // namespace helper
} // namespace kl

/*****************************************************************************

        Chapter V -- Functions declared in kl.h

 *****************************************************************************/

namespace kl {


/*!
  \brief Puts in wg the W-graph for this block.

  Explanation: the W-graph is a graph with one vertex for each element of the
  block; the corresponding descent set is the tau-invariant, i.e. the set of
  generators s that are either complex descents, real type I or II, or
  imaginary compact. Let x < y in the block such that mu(x,y) != 0, and
  descent(x) != descent(y). Then there is an edge from x to y unless
  descent(x) is contained in descent(y), and an edge from y to x unless
  descent(y) is contained in descent(x). Note that the latter containment
  always holds when the length difference is > 1, so that in that case there
  will only be an edge from x to y (the edge must be there because we already
  assumed that the descent sets were not equal.) In both cases, the
  coefficient corresponding to the edge is mu(x,y).

  NOTE: if I'm not mistaken, the edgelists come already out sorted.
*/
void wGraph(wgraph::WGraph& wg, const KLContext& klc)
{
  wg.reset();
  wg.resize(klc.size());

  // fill in descent sets
  for (BlockElt y = 0; y < klc.size(); ++y)
    wg.descent(y) = klc.descentSet(y);

  // fill in edges and coefficients
  for (BlockElt y = 0; y < klc.size(); ++y) {
    const RankFlags& d_y = wg.descent(y);
    const MuRow& mrow = klc.muRow(y);
    for (size_t j = 0; j < mrow.first.size(); ++j) {
      BlockElt x = mrow.first[j];
      const RankFlags& d_x = wg.descent(x);
      if (d_x == d_y)
	continue;
      MuCoeff mu = mrow.second[j];
      if (klc.length(y) - klc.length(x) > 1) { // add edge from x to y
	wg.edgeList(x).push_back(y);
	wg.coeffList(x).push_back(mu);
	continue;
      }
      // if we get here, the length difference is 1
      if (not d_y.contains(d_x)) { // then add edge from x to y
	wg.edgeList(x).push_back(y);
	wg.coeffList(x).push_back(mu);
      }
      if (not d_x.contains(d_y)) { // then add edge from y to x
	wg.edgeList(y).push_back(x);
	wg.coeffList(y).push_back(mu);
      }
    }
  }

} // |wGraph|

} // namespace kl
} // namespace atlas
