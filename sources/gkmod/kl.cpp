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
  Copyright 2012 David Vogan, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "kl.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <string>
#include <algorithm> // for |std::lower_bound|

#include <sys/time.h>
#include <sys/resource.h> // for getrusage in verbose

#include <cassert>
#include <set>  // for |down_set|
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


/*****************************************************************************

        Chapter I -- Public methods of the KLPolEntry and KLContext classes.

 *****************************************************************************/

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
  while (i-->0) h= (h<<21)+(h<<13)+(h<<8)+(h<<5)+h+P[i];
  return h & (modulus-1);
}

bool KLPolEntry::operator!=(KLPolEntry::Pooltype::const_reference e) const
{
  if (degree()!=e.degree()) return true;
  if (isZero()) return false; // since degrees match
  for (polynomials::Degree i=0; i<=degree(); ++i)
    if ((*this)[i]!=e[i]) return true;
  return false; // no difference found
}

/* methods of KLContext */


KLContext::KLContext(const Block_base& b)
  : klsupport::KLSupport(b) // construct unfilled support object from block
  , fill_limit(0)
  , d_prim()
  , d_kl()
  , d_mu()
  , d_store(2)
{
  // make sure the support (base class) is filled
  klsupport::KLSupport::fill();

  d_store[d_zero]=Zero; // ensure these polynomials are present
  d_store[d_one]=One;   // at expected indices, even if maybe absent in |d_kl|
}

/******** copy, assignment and swap ******************************************/


/******** accessors **********************************************************/


/*!
  \brief Returns the Kazhdan-Lusztig-Vogan polynomial $P_{x,y}$

  Precondition: row $y$ is completely computed, stored in |d_kl[y]|.

  Since |d_kl| holds all polynomials for primitive pairs $(x,y)$, this is
  basically a lookup function. While |x| is not primitive for |y|, move |x| up
  (done in |primitivize|). If this has made |x>y| (in particular if it has
  made |x==UndefBlock|, which might even be its initial value) return a zero
  polynomial. Otherwise look up $x$ in the primitive list for $y$; if found,
  use offset of result to find polynomial in |d_store|, if not found, the
  polynomial is zero. Always returns a value from |d_store|, maybe |d_zero|.
  */
KLPolRef KLContext::klPol(BlockElt x, BlockElt y) const
{
  const KLRow& klr = d_kl[y];
  x=primitivize(x,descentSet(y));
  if (x>=y) return d_store[x==y ? d_one : d_zero];

  const PrimitiveRow& pr = d_prim[y];
  PrimitiveRow::const_iterator xptr =
    std::lower_bound(pr.begin(),pr.end(),x);
  return d_store[xptr == pr.end() or *xptr != x ? d_zero
		 : klr[xptr - pr.begin()]];
}

// The same, but just return the index into |d_store| that gives $P_{x,y}$
KLIndex KLContext::KL_pol_index(BlockElt x, BlockElt y) const
{
  const KLRow& klr = d_kl[y];
  x=primitivize(x,descentSet(y));
  if (x>=y) return x==y ? d_one : d_zero;

  const PrimitiveRow& pr = d_prim[y];
  PrimitiveRow::const_iterator xptr =
    std::lower_bound(pr.begin(),pr.end(),x);
  return xptr == pr.end() or *xptr != x ? d_zero : klr[xptr - pr.begin()];
}

// an auxiliary needed to be passed in a call of |std::lower_bound|
bool mu_entry_compare(const std::pair<BlockElt,MuCoeff>& x,
		      const std::pair<BlockElt,MuCoeff>& y)
{ return x.first<y.first; }

/*!
  \brief Returns mu(x,y).

  Explanation: it is guaranteed that all the x'es such that mu(x,y) != 0
  occur in d_mu[y] (and in fact, that only those occur.) So it is a simple
  matter of looking up x. We can say 0 without lookup in some easy cases.
*/
MuCoeff KLContext::mu(BlockElt x, BlockElt y) const
{
  unsigned int lx=length(x),ly=length(y);
  if (ly<=lx or (ly-lx)%2==0)
    return MuCoeff(0);
  const MuRow& mr = d_mu[y];
  std::pair<BlockElt,MuCoeff> x0(x,0); // second component is unused
  MuRow::const_iterator xloc=
    std::lower_bound(mr.begin(),mr.end(),x0,&mu_entry_compare);

  if (xloc==mr.end() or xloc->first!=x)
    return MuCoeff(0); // x not found in mr

  return mr[xloc-mr.begin()].second;
}

/*
  Returns the list of all x extremal w.r.t. y.

  Explanation: this means that length(x) < length(y), and every descent
  for y is either a descent for x.  Or:  asc(x)\cap desc(y)=\emptyset
  Here descent means "in the tau invariant" (possibilities C-, ic, r1, r2).
*/
PrimitiveRow KLContext::extremalRow(BlockElt y)
  const
{
  BitMap b(size());
  b.fill(0,lengthLess(length(y)));  // start with all elements < y in length
  filter_extremal(b,descentSet(y)); // filter out those that are not extremal

  return PrimitiveRow(b.begin(),b.end()); // convert to list
}


/*
  Returns the list of all x primitive w.r.t. y.

  Explanation: this means that length(x) < length(y), and every descent
  for y is either a descent, or an imaginary type II ascent for x.
*/
PrimitiveRow KLContext::primitiveRow(BlockElt y) const
{
  BitMap b(size());
  b.fill(0,lengthLess(length(y)));   // start with all elements < y in length
  filter_primitive(b,descentSet(y)); // filter out those that are not primitive

  return PrimitiveRow(b.begin(),b.end());
}


/******** manipulators *******************************************************/

/*!
  \brief Fills (or extends) the KL- and mu-lists.
*/
void KLContext::fill(BlockElt y, bool verbose)
{
  if (y<fill_limit)
    return; // tables present already sufficiently large for |y|

#ifndef VERBOSE
  verbose=false; // if compiled for silence, force this variable
#endif

  try
  {
    d_prim.resize(y+1);
    d_kl.resize(y+1);
    d_mu.resize(y+1);
    if (verbose)
    {
      std::cerr << "computing Kazhdan-Lusztig polynomials ..." << std::endl;
      verbose_fill(y);
      std::cerr << "done" << std::endl;
    }
    else
      silent_fill(y);

    fill_limit = y+1; // commit extension of tables
  }
  catch (std::bad_alloc)
  { // roll back, and transform failed allocation into MemoryOverflow
    std::cerr << "\n memory full, KL computation abondoned." << std::endl;
    d_prim.resize(fill_limit);
    d_kl.resize(fill_limit);
    d_mu.resize(fill_limit); // truncate to previous contents
    throw error::MemoryOverflow();
  }

}

BitMap KLContext::primMap (BlockElt y) const
{
  BitMap b(size()); // block-size bitmap

  // start with all elements < y in length
  b.fill(0,lengthLess(length(y)));
  b.insert(y);   // and y itself

  filter_primitive(b,descentSet(y)); // filter out those that are not primitive

  // now b holds a bitmap indicating primitive elements for y

  // our result will be a bitmap of that capacity
  BitMap result (b.size()); // initiallly all bits are cleared

 // the list of primitive elements with nonzero polynomials at y
  const PrimitiveRow& row=d_prim[y];

  // traverse |b|, and for elements that occur in |row|, set bits in |result|

  size_t position=0; // position among set bits in b (avoids using b.position)
  size_t j=0; // index into row;
  for (BitMap::iterator it=b.begin(); j<row.size() and it(); ++position,++it)
    if (*it==row[j]) // look if |*it| occurs in |row| (indexes nonzero element)
    {
      result.insert(position); ++j; // record its position and advance in row
    }

  return result;
}


/*****************************************************************************

        Chapter II -- Private methods used during construction

 *****************************************************************************/

/*!
  \brief Returns the first descent generator that is not real type II

  Explanation: these are the ones that give a direct recursion formula for the
  K-L basis element. Explicitly, we search for a generator |s| such that
  |descentValue(s,y)| is either |DescentStatus::ComplexDescent| or
  |DescentStatus::RealTypeI|. If no such generator exists, we return |rank()|.
*/
weyl::Generator KLContext::firstDirectRecursion(BlockElt y) const
{
  const DescentStatus& d = descent(y);
  weyl::Generator s;
  for (s=0; s<rank(); ++s)
    if (DescentStatus::isDirectRecursion(d[s]))
      break;

  return s;

} // |KLContext::firstDirectRecursion|

/*
  Returns the first real nonparity ascent for y that is a complex ascent, or
  imaginary type 2, or compact imaginary for x.

  Explanation: those are the ones that give a nice new recursion formula for
  the K-L polynomial

  If no such generator exists, we return |rank()|.
*/
weyl::Generator KLContext::first_nice_and_real(BlockElt x,BlockElt y) const
{
  const DescentStatus& dx = descent(x);
  const DescentStatus& dy = descent(y);

  weyl::Generator s;
  for (s=0; s<rank(); ++s)
    if (dy[s]==DescentStatus::RealNonparity)
      { DescentStatus::Value vx = dx[s];
	if (vx==DescentStatus::ComplexAscent or
	    vx==DescentStatus::ImaginaryTypeII or
	    vx==DescentStatus::ImaginaryCompact)
	  break; // and return the current |s|
      }
  return s;

} // |KLContext::first_nice_and_real|

/*
  Preconditions:
  * all descents for y are of type r2 (firstDirectRecursion has failed)
  * x is extremal for y (none of the descents for y are ascents for x)
  * none of the rn ascents for y is C+, ic or i2 for x (so
    |first_nice_and_real| failed to find anything)

  Returns the first pair (s,t) such that
  1) (s,t) is (rn,r2) for y;
  2) (s,t) is (i1,ic) for x;
  3) (s,t) is (i1,i1/2) for s.x.
  Since the statuses of t for x and s.x differ, s must be ajdacent to t in the
  Dynkin diagram (but this is not tested or used explicitly here)

  Such a pair can be used to compute P_{x,y}+P_{s.x,y} using s,
  then P_{s.x,y} using t, and so P_{x,y}.

  The test that t is ic for x is omitted, since x and s.x being related by an
  imaginary cross action are in the same fiber, so t is imaginary for x if it
  is so for s.x, and noncompact for x is ruled out by extremality precondition

  If no such pair exists, (rank,*) is returned. Under the given preconditions,
  such failure allows concluding that P_{x,y}=0.

  The function may also return (s,rank), indicating that a good s was found,
  but that s.x was undefined (not in the partial block) so that no t could
  even be searched for. In this case one always has P_{s.x,y}=0, so the above
  method can still be used, and indeed simplifies by not needing t.
*/
std::pair<weyl::Generator,weyl::Generator>
   KLContext::first_endgame_pair(BlockElt x, BlockElt y) const
{
  const DescentStatus& dx = descent(x);
  const DescentStatus& dy = descent(y);

  weyl::Generator s,t, r=rank();

  for (s=0; s<r; ++s)
    if (dy[s]==DescentStatus::RealNonparity and
	dx[s]==DescentStatus::ImaginaryTypeI)
    { BlockElt sx = cross(s,x);
      if (sx==UndefBlock) // cross image might be outside partial block
	return std::make_pair(s,r); // cannot and need not search for |t|
      const DescentStatus& dsx = descent(sx);
      for (t = 0; t<r; ++t)
	if (dy[t]==DescentStatus::RealTypeII)
	  if (dsx[t]==DescentStatus::ImaginaryTypeI or
	      dsx[t]==DescentStatus::ImaginaryTypeII)
	    return std::make_pair(s,t);
    }
  return std::make_pair(r,0); // failure

} // |KLContext::first_endgame_pair|

// A convenience method that is "derived" from (the non-ancestor) |Block|
inline BlockEltPair KLContext::inverseCayley(size_t s, BlockElt y) const
{ return block().inverseCayley(s,y); }

// compute the down-set of $y$, the non-extremal $x$ with $\mu(x,y)\neq0$
std::set<BlockElt> KLContext::down_set(BlockElt y) const
{
  std::set<BlockElt> result;

  for (RankFlags::iterator it=descentSet(y).begin(); it(); ++it)
    switch (descentValue(*it,y))
    {
    case DescentStatus::ComplexDescent: result.insert(cross(*it,y));
      break;
    case DescentStatus::RealTypeI:
      {
	BlockEltPair sy = inverseCayley(*it,y);
	result.insert(sy.first); result.insert(sy.second);
      }
      break;
    case DescentStatus::RealTypeII:
      result.insert(inverseCayley(*it,y).first);
      break;
    default: // |case DescentStatus::ImaginaryCompact| nothing
      break;
    }
  return result;

} // |KLContext::down_set|

/*!
  \brief Returns the Kazhdan-Lusztig polynomial for x corresponding to
  the given row.

  Precondition: |klv| holds the tail of the set of primitive Kazhdan-Lusztig
  polynomials for |y|, enough to find the required one by elementary lookup;
  |[p_begin,p_end[| is the corresponding range of primitive elements.

  Algorithm: primitivize |x| with respect to the descents in |y|; if a real
  nonparity situation is encountered, return |Zero|; otherwise look up the
  primitive |x| in the range and return the corresponding element from |klv|.

  Like the basic |klPol|, this will return zero when |x==UndefBlock|.
*/
KLPolRef KLContext::klPol(BlockElt x, BlockElt y,
			  KLRow::const_iterator klv,
			  PrimitiveRow::const_iterator p_begin,
			  PrimitiveRow::const_iterator p_end) const
{
  x = primitivize(x,descentSet(y));

  if (x>=y) return d_store[x==y ? d_one : d_zero];
  PrimitiveRow::const_iterator xptr =
    std::lower_bound(p_begin,p_end,x);
  return d_store[xptr == p_end or *xptr != x ? d_zero : klv[xptr-p_begin]];
}


// private manipulators

/*!
  \brief Fills in the row for y in the KL-table.

  Precondition: all lower rows have been filled

  Row of $y$ is the set of all $P_{x,y}$ for $x<y$; actually more like a column
*/
size_t KLContext::fillKLRow(BlockElt y, KLHash& hash)
{
  size_t sparseness=0; // number of entries saved by suppressing zero polys
  if (d_kl[y].size()>0)
    return 0; // row has already been filled
  weyl::Generator s = firstDirectRecursion(y);
  if (s<rank())  // a direct recursion was found, use it for |y|, for all |x|
  {
    std::vector<KLPol> klv;
    PrimitiveRow e = extremalRow(y); // we compute for |x| extremal only

    recursionRow(klv,e,y,s); // compute all polynomials for these |x|
    // write result
    sparseness += writeRow(klv,e,y,hash);
  }
  else // we must use an approach that distinguishes on |x| values
  {
    KLRow klv;
    PrimitiveRow pr = primitiveRow(y); // here we do all |x| primitive for |y|
    // (any ascents for x that are descents for y must be imaginary type II)
    newRecursionRow(klv,pr,y,hash); // put result of recursion formula in klv
    sparseness +=
      remove_zeros(klv,pr,y); // write |d_prim[y]|, |d_kl[y]|, suppressing zeros
  }
  return sparseness;
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
void KLContext::recursionRow(std::vector<KLPol>& klv,
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

    // the following loop could be run in either direction: no dependency.
    // however it is natural to take |x| descending from |y| (exclusive)
    for (i=e.size(); i-->0; )
    {
      BlockElt x = e[i]; // extremal for $y$, so $s$ is descent for $x$
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
	  KLPolRef Pxsy = klPol(x,sy);
	  klv[i].safeAdd(Pxsy,1);
	  klv[i].safeSubtract(Pxsy);
	}
	break;
      case DescentStatus::RealTypeII:
	{ // P_{sx,sy}+qP_{x,sy}-P_{s.x,sy}
	  BlockElt sx = inverseCayley(s,x).first;
	  klv[i] = klPol(sx,sy);
	  klv[i].safeAdd(klPol(x,sy),1);
	  klv[i].safeSubtract(klPol(cross(s,x),sy));
	}
	break;
      default: assert(false); // this cannot happen
      }
    }
  }
  catch (error::NumericUnderflow& err){
    throw kl_error::KLError(e[i],y,__LINE__,
			    static_cast<const KLContext&>(*this));
  }

  muCorrection(klv,e,y,s); // subtract mu-correction from all of |klv|

} // |KLContext::recursionRow|

/*!
  \brief Subtracts from all polynomials in |klv| the correcting terms in the
  K-L recursion.

  Precondtion: |klv| already contains, for all $x$ that are primitive w.r.t.
  |y| in increasing order, the terms in $P_{x,y}$ corresponding to
  $c_s.c_{y'}$, whery |y'| is $s.y$ if |s| is a complex descent, and |y'| is
  an inverse Cayley transform of |y| if |s| is real type I.
  The mu-table and KL-table have been filled in for elements of length < l(y).

  Explanation: the recursion formula is of the form:
  $$
    lhs = c_s.c_{y'} - \sum_{z} mu(z,y')c_z
  $$
  where |z| runs over the elements $< y'$ such that |s| is a descent for |z|.
  Here $lhs$ stands for $c_y$ when |s| is a complex descent or real type I for
  |y|, and for $c_{y}+c_{s.y}$ when |s| is real type II; however it plays no
  part in this function that only subtracts $\mu$-terms.

  The element $y'$ is called |sy| in the code below.

  We construct a loop over |z| first, before traversing |klv| (the test for
  $z<sy$ is absent, but $\mu(z,sy)\neq0$ implies $z<sy$ (strict, as mu(sy,sy)
  is 0; in any case no coefficient for |sy| is stored in |d_mu[sy]|, and
  moreover $z=sy$ would be rejected by the descent condition). The choix have
  the out loop over $z$ and the inner loop over $x$ (i.e., over |klv|) allows
  fetching $\mu(z,sy)$ only once, and terminating each scan of |klv| once its
  values |x| become too large to produce a non-zero $P_{x,z}$. (In fact we
  stop once $l(x)=l(z)$, and separately consider the case $x=z$.) Either
  direction of the loop on $z$ would work, but taking it decreasing is more
  natural; we keep track of the index |zi| at which $x=z$ occurs, if it does.

  Elements of length at least $l(sy)=l(y)-1$ on the list |e| are always
  rejected, so the tail of |e| never reached.
 */
void KLContext::muCorrection(std::vector<KLPol>& klv,
			     const PrimitiveRow& e,
			     BlockElt y, size_t s)
{
  BlockElt sy =
    descentValue(s,y) == DescentStatus::ComplexDescent ? cross(s,y)
    : inverseCayley(s,y).first;  // s is real type I for y here, ignore .second

  const MuRow& mrow = d_mu[sy];
  size_t l_y = length(y);

  size_t zi=e.size(); // should satisfy |e[zi]==z| whenever such |zi| exists

  size_t j; // define outside for error reporting
  try {
    for (size_t i = mrow.size(); i-->0; ) // loop over |z| decreasing from |sy|
    {
      BlockElt z = mrow[i].first;
      DescentStatus::Value v = descentValue(s,z);
      if (not DescentStatus::isDescent(v))
	continue;

      size_t l_z = length(z);
      while (zi>0 and e[zi-1]>=z)
	--zi; // ensure |e[k]>=z| if and only if |k>=iz|

      MuCoeff mu = mrow[i].second; // mu!=MuCoeff(0)

      polynomials::Degree d = (l_y-l_z)/2; // power of q used in the loops below

      assert( zi==e.size() or e[zi]>z or
	      (klv[zi].degree()==d and klv[zi][d]==mu) );

      if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
	for (j = 0; j < e.size(); ++j)
	{
	  BlockElt x = e[j];
	  if (length(x) >= l_z) break; // once reached, $x=z$ is only case left

	  KLPolRef pol = klPol(x,z);
	  klv[j].safeSubtract(pol,d); // subtract q^d.P_{x,z} from klv[j]
	} // for (j)
      else // mu!=MuCoeff(1)
	for (j = 0; j < e.size(); ++j)
	{
	  BlockElt x = e[j];
	  if (length(x) >= l_z) break; // once reached, $x=z$ is only case left

	  KLPolRef pol = klPol(x,z);
	  klv[j].safeSubtract(pol,d,mu); // subtract q^d.mu.P_{x,z} from klv[j]
	} // for {j)

      if (zi<e.size() and e[zi]==z) // handle final term |x==z|
      {
	assert( klv[zi].degree()==d and klv[zi][d]==mu );
	klv[zi].safeSubtract(KLPol(d,mu)); // subtract off the term $mu.q^d$
      }

    } // for (i)
  }
  catch (error::NumericUnderflow& err){
    throw kl_error::KLError(e[j],y,__LINE__,
			    static_cast<const KLContext&>(*this));
  }

} // |KLContext::muCorrection|

/*!
  \brief Writes down row y in d_kl and d_prim.

  Precondition: The parallel pair (klv,er) records the polynomials for at
  least all the extremal values $x$ for $y$, and at most for all primitive
  values $x$ for $y$. So when $x=er[i]$ then $P_{x,y}=klv[i]$. In practice
  |er| will contain either all extremal elements (when called from
  |recursionRow|) or all primitive elements (for |newRecursionRow|).

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
size_t KLContext::writeRow(const std::vector<KLPol>& klv,
			   const PrimitiveRow& er, BlockElt y,
			   KLHash& hash)
{
  PrimitiveRow pr = primitiveRow(y);

  PrimitiveRow nzpr(pr.size()); // columns of the nonzero primimitive entries
  KLRow KL(pr.size()); // nonzero primitive entries (indexes into d_store)

  PrimitiveRow::iterator nzpr_p=nzpr.end();
  KLRow::iterator KL_p=KL.end();

  BitMap mu_elements(pr.size()); // flags |i| with |mu(pr[i],y)>0|

  unsigned int ly = length(y);

  size_t j= er.size()-1; // points to current extremal element

  for (size_t i = pr.size(); i-->0; )
    if (j<er.size() and pr[i]==er[j]) // must test for underflow |j|
    { // extremal element; use stored polynomial
      const KLPol& Pxy=klv[j--];
      if (not Pxy.isZero())
      {
	*--nzpr_p = pr[i];
        *--KL_p = hash.match(Pxy);
	unsigned int lx=length(pr[i]);
	if ((ly-lx)%2>0 and Pxy.degree()==(ly-lx)/2) // in fact |(ly-lx-1)/2|
	  mu_elements.insert(i);
      }
    }
    else // insert a polynomial for primitive non-extremal pr[i] if nonzero
    {
      unsigned int s = ascent_descent(pr[i],y);
      assert(descentValue(s,pr[i])==DescentStatus::ImaginaryTypeII);
      BlockEltPair xs = cayley(s,pr[i]);
      KLPol Pxy = klPol(xs.first,y,KL_p,nzpr_p,nzpr.end()); // look up using KL
      Pxy.safeAdd(klPol(xs.second,y,KL_p,nzpr_p,nzpr.end()));
      if (not Pxy.isZero())
      {
	*--nzpr_p = pr[i];
        *--KL_p = hash.match(Pxy);
      } // no need to check for |mu| here; |down_set| covers possible cases
    }

  if (ly==0)
    return 0; // nothing left to do at minimal length

  std::set<BlockElt> downs = down_set(y);

  // commit
  d_prim[y] = PrimitiveRow(nzpr_p,nzpr.end()); // copy shifting, from |nzpr_p|
  d_kl[y] = KLRow(KL_p,KL.end());

  d_mu[y].reserve(mu_elements.size()+downs.size());

  for (BitMap::iterator it=mu_elements.begin(); it(); ++it)
  {
    BlockElt x=pr[*it];
    KLPolRef Pxy = klPol(x,y,KL_p,nzpr_p,nzpr.end());
    assert(not Pxy.isZero());
    d_mu[y].push_back(std::make_pair(x,Pxy[Pxy.degree()]));
  }
  for (std::set<BlockElt>::iterator it=downs.begin(); it!=downs.end(); ++it)
  {
    d_mu[y].push_back(std::make_pair(*it,MuCoeff(1)));
    for (size_t i=d_mu[y].size()-1; i>0 and d_mu[y][i-1].first>*it; --i)
      std::swap(d_mu[y][i-1],d_mu[y][i]); // insertion-sort
  }

  return nzpr_p  -nzpr.begin(); // measure unused space

} // |KLContext::writeRow|

// this method is called instead of |writeRow| in cases involving new recursion
size_t KLContext::remove_zeros(const KLRow& klv,
			       const PrimitiveRow& pr, BlockElt y)
{
  PrimitiveRow nzpr(pr.size()); // columns of the nonzero primimitive entries
  KLRow KL(pr.size()); // nonzero primitive entries (indexes into d_store)

  PrimitiveRow::iterator nzpr_p=nzpr.end();
  KLRow::iterator KL_p=KL.end();

  for (size_t i = pr.size(); i-->0; )
    if (not isZero(klv[i]))
    {
      *--nzpr_p = pr[i];
      *--KL_p = klv[i];
    }


  // commit
  d_prim[y] = PrimitiveRow(nzpr_p,nzpr.end()); // copy shifting, from |nzpr_p|
  d_kl[y] = KLRow(KL_p,KL.end());

  return nzpr_p  -nzpr.begin(); // measure unused space

} // |KLContext::remove_zeros|

/*
  Puts in klv[i] the polynomial P_{e[i],y} for every primtitve x=pr[i],
  computed by a recursion formula for those |y| admitting no direct recursion.

  Precondition: every simple root is for y either a complex ascent or
  imaginary or real (no complex descents for y). (split 1)

  In fact real type 1 descents for |y| don't occur, but this is not used.

  From the precondition we get: for each extremal |x| for |y|, there either
  exists a true ascent |s| that is real for |y|, necessarily nonparity because
  |x| is extremal (split 3), or we are assured that $P_{x,y}=0$.

  Here there is a recursion formula of a somewhat opposite nature than in the
  case of direct recursion. The terms involving $P_{x',y}$ where $x'$ are in
  the up-set of |x| appear in what is most naturally the left hand side of the
  equation, while the sum involving |mu| values appears on the right (3.2). As
  a consequence, the |mu| terms will be computed first, and then modifications
  involving such $P_{x',y}$ and subtraction are applied. However if |s| is
  type 'i1' for |x| the left hand side has (apart from $P_{x,y}$) another term
  $P_{s.x,y}$ for the imaginary cross image $s.x$, and so $P_{x,y}$ cannot be
  directly obtained in this manner; we shall avoid this case if we can.

  An additional case where |s| is 'ic' for |x| and 'rn' for |y| can be handled
  similarly, and we do so if the occasion presents itself. Here there are no
  terms from the up-set of |x| so the right hand side is just the mu terms,
  but the left hand side is $(q+1)P_{x,y}$ truncated to terms of degree at
  most $(l(y)-l(x)-1)/2$, from which we can recover $P_{x,y}$ by |safeDivide|.

  All in all the following cases are handled easily: |x| is (primitive but)
  not extremal, |x| has some |s| which is 'rn' for |y| and one of 'C+', 'i2'
  or 'ic' for |x| (|first_nice_and_real|). If no such |s| exists we can almost
  conclude $P_{x,y}=0$, but need to handle an exceptional "endgame" situation
  in which the formula for an 'i1' ascent can be exploited in spite of the
  presence of $P_{s.x,y}$, because that term can be computed on the fly.

  The sum involving mu, produced by |muNewFormula|, has terms involving
  $P_{x,u}\mu(u,y}$, so when doing a downward loop over |x| it pays to keep
  track of the previous |u| with nonzero $\mu(u,y)$.

  This code gets executed for |y| that are of minimal length, in which case
  it only contributes $P_{y,y}=1$; the |while| loop will be executed 0 times.
*/
void KLContext::newRecursionRow
( KLRow& klv,
  const PrimitiveRow& pr, // primitive elements of length less than |length(y)|
  BlockElt y,
  KLHash& hash)
{
  klv.resize(pr.size());
  KLRow::iterator kl_p=klv.end();
  PrimitiveRow::const_iterator pr_p=pr.end();

  unsigned int l_y = length(y);

  MuRow mu_y; mu_y.reserve(lengthLess(length(y))); // a very gross estimate

  // start off |mu_y| with ones for |down_set(y)|, not otherwise computed
  std::set<BlockElt> downs = down_set(y);
  for (std::set<BlockElt>::iterator it=downs.begin(); it!=downs.end(); ++it)
    mu_y.push_back(std::make_pair(*it,MuCoeff(1)));

  size_t j = klv.size(); // declare outside try block for error reporting
  try {
    while (--kl_p,--pr_p, j-->0) // |*kl_p=klv[j]|, |*pr_p=pr[j]|
    {
      BlockElt x = pr[j];

      unsigned int s= ascent_descent(x,y);
      if (s<rank()) // a primitive element that is not extremal; easy case
      { // equation (1.9b) in recursion.pdf
	assert(descentValue(s,x)==DescentStatus::ImaginaryTypeII);
	BlockEltPair p = cayley(s,x);
	KLPol pol = klPol(p.first,y,kl_p,pr_p,pr.end());
	pol.safeAdd(klPol(p.second,y,kl_p,pr_p,pr.end()));
	klv[j] = hash.match(pol);
	continue; // done with |x|, go on to the next
      }

      unsigned int l_x = length(x);

      // now x is extremal for y. By (split 1) and Lemma 3.1 of recursion.pdf
      // this implies that if x<y in the Bruhat order, there is at least one s
      // real for y that is a true ascent (not rn) of x and therefore rn for y
      // we first hope that at least one of them is not i1 for x

      // we first seek a real nonparity ascent for y that is C+,i2 or ic for x
      s = first_nice_and_real(x,y);
      if (s < rank()) // there is such an ascent s
      {
	// start setting |pol| to the expression (3.4) in recursion.pdf
	KLPol pol = muNewFormula(x,y,s,mu_y);

	switch (descentValue(s,x))
	{
	case DescentStatus::ComplexAscent:
	{ // use equations (3.3a)=(3.4)
	  BlockElt sx = cross(s,x);
	  pol.safeSubtract(klPol(sx,y,kl_p,pr_p,pr.end()),1);
	  // subtract qP_{sx,y} from mu terms
	} // ComplexAscent case
	break;

	case DescentStatus::ImaginaryTypeII:
	{ // use equations (3.3a)=(3.5)
	  BlockEltPair p = cayley(s,x);
	  KLPol sum = klPol(p.first,y,kl_p,pr_p,pr.end());
	  sum.safeAdd(klPol(p.second,y,kl_p,pr_p,pr.end()));
	  pol.safeAdd(sum);
	  pol.safeSubtract(sum,1); //now we've added (1-q)(P_{x',y}+P_{x'',y})
	  pol.safeDivide(2);   //this may throw
	} // ImaginaryTypeII case
	break;

	case DescentStatus::ImaginaryCompact:
	  /* here s is a emph{descent} for x, which causes an extra unknown
	     leading (if nonzero) term to appear in addition to (3.4), giving
	     rise to equation (3.7). Yet we can determine the quotient by q+1.
	  */
	  pol.safeQuotient(length(y)-length(x));
	  break;

	default: assert(false); //we've handled all possible NiceAscents
	}
	klv[j] = hash.match(pol);
	if ((l_y-l_x)%2!=0 and pol.degree()==(l_y-l_x)/2)
	  mu_y.push_back(std::make_pair(x,pol[pol.degree()]));

      } // end of |first_nice_and_real| case

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
	std::pair<size_t,size_t> st = first_endgame_pair(x,y);
	if ((s=st.first) < rank())
	{
	  KLPol pol = muNewFormula(x,y,s,mu_y);

	  //subtract (q-1)P_{xprime,y} from terms of expression (3.4)
	  BlockElt xprime = cayley(s,x).first;
	  const KLPol& P_xprime_y =  klPol(xprime,y,kl_p,pr_p,pr.end());
	  pol.safeAdd(P_xprime_y);
	  pol.safeSubtract(P_xprime_y,1);

	  //now klv[j] holds P_{x,y}+P_{s.x,y}

	  unsigned int t=st.second;

	  if (t<rank()) // nothing to subtract if $s.x$ not in partial block
	  {
	    //compute P_{s.x,y} using t
	    BlockEltPair sx_up_t = cayley(t,cross(s,x));

	    // any |UndefBlock| component of |sx_up_t| will contribute $0$
	    pol.safeSubtract(klPol(sx_up_t.first,y,kl_p,pr_p,pr.end()));
	    pol.safeSubtract(klPol(sx_up_t.second,y,kl_p,pr_p,pr.end()));

	  }

	  klv[j] = hash.match(pol);
	  if ((l_y-l_x)%2!=0 and pol.degree()==(l_y-l_x)/2)
	    mu_y.push_back(std::make_pair(x,pol[pol.degree()]));
	}
	else // |first_endgame_pair| found nothing
	  klv[j]=d_zero;
      } // end of no NiceAscent case
    } // while (j-->0)
  }
  catch (error::NumericUnderflow& err) // repackage error, reporting x,y
  {
    throw kl_error::KLError(pr[j],y,__LINE__,
			    static_cast<const KLContext&>(*this));
  }

  d_mu[y]=MuRow(mu_y.rbegin(),mu_y.rend()); // reverse to a tight copy

  for (unsigned int k=mu_y.size()-downs.size(); k<mu_y.size(); ++k)
  { // successively insertion-sort the down-set x's into previous list
    for (size_t i=k; i>0 and d_mu[y][i-1].first>d_mu[y][i].first; --i)
      std::swap(d_mu[y][i-1],d_mu[y][i]);
  }

} // |KLContext::newRecursionRow|

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

  We construct a loop over |z|. The test for $z<y$ is absent, but implied by
  $\mu(z,y)\neq0$; the $\mu(\cdot,y)$ information is passed in the |mu_y|
  argument to this method. The chosen loop order allows fetching $\mu(z,y)$
  only once, and terminating the scan of |klv| once its values |x| become too
  large to produce a non-zero $P_{x,z}$.
*/
KLPol KLContext::muNewFormula
  (BlockElt x, BlockElt y, size_t s, const MuRow& mu_y)
{
  KLPol pol=Zero;

  unsigned int lx=length(x), ly = length(y);

  try
  {
    for (MuRow::const_iterator it=mu_y.begin(); it!=mu_y.end(); ++it)
    {
      BlockElt z = it->first; // a block element with $\mu(z,y)\neq0$
      unsigned int lz = length(z);
      if (lz<=lx)
	break; // length |z| decreases, and |z==x| must be excluded, so stop
      if (not DescentStatus::isDescent(descentValue(s,z))) continue;

      // now we have a true contribution with nonzero $\mu$
      unsigned int d = (ly - lz +1)/2; // power of $q$ used in the formula
      MuCoeff mu = it->second;
      KLPolRef Pxz = klPol(x,z); // which is known because $z<y$

      if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
	pol.safeAdd(Pxz,d); // add $q^d.P_{x,z}$ to |pol|
      else // mu!=MuCoeff(1)
	pol.safeAdd(Pxz,d,mu); // add $q^d.\mu(z,y).P_{x,z}$ to |pol|

    } // for (k)
  }
  catch (error::NumericOverflow& e){
    throw kl_error::KLError(x,y,__LINE__,
			    static_cast<const KLContext&>(*this));
  }

  return pol;
} // |KLContext::muNewFormula|


void KLContext::silent_fill(BlockElt last_y)
{
  try
  {
    KLHash hash(d_store); // (re-)construct a hastable for polynomial storage
    // fill the lists
    for (BlockElt y=fill_limit; y<=last_y; ++y)
      fillKLRow(y,hash);
    // after all rows are done the hash table is freed, only the store remains
  }
  catch (kl_error::KLError& e)
  {
    std::ostringstream os;
    os << "negative coefficient in P_{" << e.x << ',' << e.y
       << "} at line " << e.line << '.';
    throw std::runtime_error(os.str()); // so that atlas may catch it
  }
}

/*
  New routine that does verbose filling of existing |KLContext| object
*/
void KLContext::verbose_fill(BlockElt last_y)
{
  try
  {
    KLHash hash(d_store);

    size_t minLength = length(fill_limit); // length of first new |y|
    size_t maxLength = length(last_y<size() ? last_y : size()-1);

    //set timers for KL computation
    std::time_t time0;
    std::time(&time0);
    std::time_t time;

    struct rusage usage; //holds Resource USAGE report
    size_t storesize = 0; // previous size of d_store
    size_t polsize = 0; // running total of sum of (polynomial degrees+1)

    size_t nr_of_prim_nulls = 0, prim_size = 0;

    for (size_t l=minLength; l<=maxLength; ++l) // by length for progress report
    {
      BlockElt y_start = l==minLength ? fill_limit : lengthLess(l);
      BlockElt y_limit = l<maxLength ? lengthLess(l+1) : last_y+1;
      for (BlockElt y=y_start; y<y_limit; ++y)
      {
	std::cerr << y << "\r";

	nr_of_prim_nulls += fillKLRow(y,hash);
	prim_size += d_prim[y].size();
      }

      // now length |l| is completed
      size_t p_capacity // currently used memory for polynomials storage
	= hash.capacity()*sizeof(KLIndex) + d_store.capacity()*sizeof(KLPol);
      for (size_t i=storesize; i<d_store.size(); ++i)
	polsize+= (d_store[i].degree()+1)*sizeof(KLCoeff);
      storesize = d_store.size(); // avoid recounting polynomials!
      p_capacity += polsize;

      std::cerr // << "t="    << std::setw(5) << deltaTime << "s.
	<< "l=" << std::setw(3) << l // completed length
	<< ", y="  << std::setw(6)
	<< lengthLess(l+1)-1 // last y value done
	<< ", polys:"  << std::setw(11) << d_store.size()
	<< ", mat:"  << std::setw(11) << prim_size
	<<  std::endl;
      unsigned cputime, resident; //memory usage in megabytes
      if(getrusage(RUSAGE_SELF, &usage) != 0)
	std::cerr << "getrusage failed" << std::endl;
      resident = usage.ru_maxrss/1024; //largest so far??
#ifdef __APPLE__
      resident = resident/1024;
#endif
      cputime = usage.ru_utime.tv_sec;
      std::cerr << "CPU time = " << std::setw(5) << cputime
		<< " secs, Max res size="
		<< std::setw(5) << resident << "MB, pmem="
		<< std::setw(6) << p_capacity/1048576 << "MB, matmem="
		<< std::setw(6) << prim_size*sizeof(KLIndex)/1048576
		<< "MB \n";

    } // for (l=min_length+1; l<=max_Length; ++l)

    std::time(&time);
    double deltaTime = difftime(time, time0);
    std::cerr << std::endl;
    std::cerr << "Total elapsed time = " << deltaTime << "s." << std::endl;
    std::cerr << d_store.size() << " polynomials, "
	      << prim_size << " matrix entries."<< std::endl;

    std::cerr << "Number of unrecorded primitive pairs: "
	      << nr_of_prim_nulls << '.' << std::endl;
    std::cerr << std::endl;

  }
  catch (kl_error::KLError& e)
  {
    std::ostringstream os;
    os << "negative coefficient in P_{" << e.x << ',' << e.y
       << "} at line " << e.line << '.';
    throw std::runtime_error(os.str()); // so that atlas may catch it
  }

}



/*****************************************************************************

        Chapter V -- Functions declared in kl.h

 *****************************************************************************/


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
    for (size_t j = 0; j < mrow.size(); ++j) {
      BlockElt x = mrow[j].first;
      const RankFlags& d_x = wg.descent(x);
      if (d_x == d_y)
	continue;
      MuCoeff mu = mrow[j].second;
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

} // |namespace kl|
} // |namespace atlas|
