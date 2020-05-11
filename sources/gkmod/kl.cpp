/*
  This is kl.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright 2012 David Vogan, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  Implementation of the class KL_table.

  This module contains code for the computation of the Kazhdan-Lusztig
  polynomials for a given block of representations. We have taken the radical
  approach of not using the Bruhat ordering at all, just ordering by length
  instead, and coping with the ensuing appearance of zero polynomials. It
  is expected that the simplification thus achieved will more than outweigh
  the additional polynomials computed.

  The general scheme is fairly similar to the one in Coxeter: there is a
  "KLSupport" structure, that holds the list of primitive pairs that makes it
  possible to read the |d_KL| list, plus some additional lists that allow for
  a fast primitivization algorithm, for instance; there are two main lists,
  |d_KL| (filled in for all primitive pairs), and |d_mu| (filled in only for
  non-zero mu coefficients.)
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
  possible to read the |d_KL| list, plus some additional lists that allow for
  a fast extremalization algorithm, for instance; there are two main lists,
  |d_KL| (filled in for all extremal pairs), and |d_mu| (filled in only for
  non-zero mu coefficients.)
*/

namespace atlas {

namespace kl {

// Polynomial 0, which is stored as a vector of size 0.
  const KLPol Zero;

// Polynomial $1.q^0$.
  const KLPol One(0,KLCoeff(1)); // since |Polynomial(d,1)| gives |1.q^d|.

/*****************************************************************************

        Chapter I -- Public methods of the KLPolEntry and KL_table classes.

 *****************************************************************************/

/* methods of KLPolEntry */


/*
  Calculate a hash value in [0,modulus[, where modulus is a power of 2

  The function is in fact evaluation of the polynomial (with coefficients
  interpreted in $\Z$) at the point $2^{21}+2^{13}+2^8+2^5+1=2105633$, which can
  be calculated quickly (without multiplications) and which gives a good spread
  (which is not the case if 2105633 is replaced by a small number, because the
  evaluation values will not grow fast enough for low degree polynomials!).
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

/* methods of KL_table */


KL_table::KL_table(const Block_base& b)
  : klsupport::KLSupport(b) // construct unfilled support object from block
  , d_holes(b.size()) // start with ambition to fill everything
  , d_KL(b.size()) // create empty slots for whole block; doesn't cost much
  , d_mu(b.size())
  , d_store(2)
{
  d_holes.fill();
  d_store[d_zero]=Zero; // ensure these polynomials are present
  d_store[d_one]=One;   // at expected indices, even if maybe absent in |d_KL|
}

/******** copy, assignment and swap ******************************************/


/******** accessors **********************************************************/


/*
  Return the Kazhdan-Lusztig-Vogan polynomial $P_{x,y}$

  Here column $y$ must have been completely computed, and stored in |d_KL[y]|.

  Since |d_KL| holds all polynomials for primitive pairs $(x,y)$, this is
  basically a lookup function. While |x| is not primitive for |y|, move |x| up
  (done in |primitivize|). If this has made |x>y| (in particular if it has
  made |x==UndefBlock|, which might even be its initial value) return a zero
  polynomial. Otherwise look up $x$ in the primitive list for $y$; if found,
  use offset of result to find polynomial in |d_store|, if not found, the
  polynomial is zero. Always returns a value from |d_store|, maybe |d_zero|.
  */
KLPolRef KL_table::KL_pol(BlockElt x, BlockElt y) const
{
  x=primitivize(x,descentSet(y));
  if (x>=y) return d_store[x==y ? d_one : d_zero];

  KL_pair target(x,d_zero); // provide dummy second component for search
  const auto& kl_col = d_KL[y];
  auto xptr = std::lower_bound(kl_col.cbegin(),kl_col.cend(),target);
  return d_store[xptr == kl_col.cend() or xptr->x != x ? d_zero : xptr->P];
}

// The same, but just return the index into |d_store| that gives $P_{x,y}$
KLIndex KL_table::KL_pol_index(BlockElt x, BlockElt y) const
{
  x=primitivize(x,descentSet(y));
  if (x>=y) return x==y ? d_one : d_zero;

  KL_pair target(x,d_zero); // provide dummy second component for search
  const auto& kl_col = d_KL[y];
  auto xptr = std::lower_bound(kl_col.cbegin(),kl_col.cend(),target);
  return xptr == kl_col.cend() or xptr->x != x ? d_zero : xptr->P;
}

/*
  Return $\mu(x,y)$. Since |d_mu[y]| list the $x$'es such that $\mu(x,y)\neq 0$
  in increasing otder, this is a simple matter of looking up $x$. We can say 0
  without lookup in some easy cases.
*/
MuCoeff KL_table::mu(BlockElt x, BlockElt y) const
{
  unsigned int lx=length(x),ly=length(y);
  if (ly<=lx or (ly-lx)%2==0)
    return MuCoeff(0);
  const Mu_column& mc = d_mu[y];
  auto xloc= std::lower_bound(mc.cbegin(),mc.cend(),Mu_pair{x,MuCoeff(0)});

  if (xloc==mc.cend() or xloc->x!=x)
    return MuCoeff(0); // x not found in mr

  return mc[xloc-mc.cbegin()].coef;
}

/*
  Return the list of all |x| extremal w.r.t. |y|.

  Explanation: this means that |length(x) < length(y)|, and every descent for |y|
  is also a descent for |x|, in other words $asc(x)\cap desc(y)=\emptyset$.
  Here descent means "in the $\tau$ invariant" (possibilities C-, ic, r1, r2).
*/
PrimitiveColumn KL_table::extremal_column(BlockElt y)
  const
{
  BitMap b(size());
  b.fill(0,lengthLess(length(y)));  // start with all elements < y in length
  filter_extremal(b,descentSet(y)); // filter out those that are not extremal

  return PrimitiveColumn(b.begin(),b.end()); // convert to vector
}


/*
  Return the list of all |x| primitive w.r.t. |y|.

  Explanation: this means that |length(x) < length(y)|, and every descent
  for |y| is either a descent, or an imaginary type II ascent for |x|.
*/
PrimitiveColumn KL_table::primitive_column(BlockElt y) const
{
  BitMap b(size());
  b.fill(0,lengthLess(length(y)));   // start with all elements < y in length
  filter_primitive(b,descentSet(y)); // filter out those that are not primitive

  return PrimitiveColumn(b.begin(),b.end());
}


/******** manipulators *******************************************************/


KLHash KL_table::pol_hash () { return KLHash(d_store); }

// Fill (or extend) the KL- and mu-lists.
void KL_table::fill(BlockElt y, bool verbose)
{
  if (y<first_hole())
    return; // tables present already sufficiently large for |y|

#ifndef VERBOSE
  verbose=false; // if compiled for silence, force this variable
#endif

  try
  {
    if (verbose)
    {
      std::cerr << "computing Kazhdan-Lusztig polynomials ..." << std::endl;
      verbose_fill(y);
      std::cerr << "done" << std::endl;
    }
    else
      silent_fill(y);
  }
  catch (std::bad_alloc)
  { // roll back, and transform failed allocation into MemoryOverflow
    std::cerr << "\n memory full, KL computation abondoned.\n";
    for (auto it = d_holes.begin(); it() and *it<=y; ++it)
    { // remove any partially written columns
      d_KL[*it].clear();
      d_mu[*it].clear();
    }
    throw error::MemoryOverflow();
  }

}

BitMap KL_table::primMap (BlockElt y) const
{
  BitMap b(size()); // block-size bitmap

  // start with all elements < y in length
  b.fill(0,lengthLess(length(y)));
  b.insert(y);   // and y itself

  filter_primitive(b,descentSet(y)); // filter out those that are not primitive

  // now b holds a bitmap indicating primitive elements for y

  // our result will be a bitmap of that capacity
  BitMap result (b.size()); // initiallly all bits are cleared

 // the list of pairs for primitive elements with nonzero polynomials at y
  const auto& col = d_KL[y];

  // traverse |b|, and for elements that occur in |col|, set bits in |result|

  size_t position=0; // position among set bits in b (avoids using b.position)
  size_t j=0; // index into |col|;
  for (BitMap::iterator it=b.begin(); j<col.size() and it(); ++position,++it)
    if (*it==col[j].x) // look if |*it| occurs in |col| (indexes nonzero element)
    {
      result.insert(position); ++j; // record its position and advance in |col|
    }

  return result;
}


/*****************************************************************************

        Chapter II -- Private methods used during construction

 *****************************************************************************/

/*
  Return the first descent generator that is not real type II

  Explanation: these are the ones that give a direct recursion formula for the
  K-L basis element. Explicitly, we search for a generator |s| such that
  |descentValue(s,y)| is either |DescentStatus::ComplexDescent| or
  |DescentStatus::RealTypeI|. If no such generator exists, we return |rank()|.
*/
weyl::Generator KL_table::firstDirectRecursion(BlockElt y) const
{
  const DescentStatus& d = descent(y);
  weyl::Generator s;
  for (s=0; s<rank(); ++s)
    if (DescentStatus::isDirectRecursion(d[s]))
      break;

  return s;

} // |KL_table::firstDirectRecursion|

/*
  Return the first real nonparity ascent for y that is a complex ascent, or
  imaginary type 2, or compact imaginary for x.

  Explanation: those are the ones that give a nice new recursion formula for
  the K-L polynomial

  If no such generator exists, we return |rank()|.
*/
weyl::Generator KL_table::first_nice_and_real(BlockElt x,BlockElt y) const
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

} // |KL_table::first_nice_and_real|

/*
  Preconditions:
  * all descents for y are of type r2 (firstDirectRecursion has failed)
  * x is extremal for y (none of the descents for y are ascents for x)
  * none of the rn ascents for y is C+, ic or i2 for x (so
    |first_nice_and_real| failed to find anything)

  Return the first pair (s,t) such that
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
   KL_table::first_endgame_pair(BlockElt x, BlockElt y) const
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

} // |KL_table::first_endgame_pair|

// A convenience method that is "derived" from (the non-ancestor) |Block|
inline BlockEltPair KL_table::inverse_Cayley(weyl::Generator s, BlockElt y) const
{ return block().inverseCayley(s,y); }


/*
  Return the Kazhdan-Lusztig polynomial for x corresponding to the given column.

  Precondition: |klv| holds the tail of the set of primitive Kazhdan-Lusztig
  polynomials for |y|, enough to find the required one by elementary lookup;
  |[p_begin,p_end[| is the corresponding range of primitive elements.

  Algorithm: primitivize |x| with respect to the descents in |y|; if a real
  nonparity situation is encountered, return |Zero|; otherwise look up the
  primitive |x| in the range and return the corresponding |KLIndex| value.

  Like the basic |KL_pol|, this will return zero when |x==UndefBlock|.
*/
KLPolRef KL_table::KL_pol(BlockElt x, BlockElt y,
			  KL_column::const_iterator p_begin,
			  KL_column::const_iterator p_end) const
{
  x = primitivize(x,descentSet(y));

  if (x>=y) return d_store[x==y ? d_one : d_zero];

  KL_pair target(x,d_zero);
  auto xptr = std::lower_bound(p_begin,p_end,target);
  return d_store[xptr == p_end or xptr->x != x ? d_zero : xptr->P];
}


// private manipulators

/*
  Fill the column for |y| in the KL-table.

  Precondition: all lower rows have been filled

  Column of $y$ is the set of all $P_{x,y}$ for $x<y$
*/
size_t KL_table::fill_KL_column(BlockElt y, KLHash& hash)
{
  size_t sparseness=0; // number of entries saved by suppressing zero polys
  if (d_KL[y].size()>0)
    return 0; // column has already been filled
  weyl::Generator s = firstDirectRecursion(y);
  if (s<rank())  // a direct recursion was found, use it for |y|, for all |x|
  {
    std::vector<KLPol> klv;
    PrimitiveColumn e = extremal_column(y); // we compute for |x| extremal only

    recursion_column(klv,e,y,s); // compute all polynomials for these |x|
    // write result
    sparseness += write_column(klv,e,y,hash);
  }
  else // we must use an approach that distinguishes on |x| values
  {
    KL_column klv;
    auto pc = primitive_column(y); // here we do all |x| primitive for |y|

    klv.reserve(pc.size());
    for (BlockElt x : pc)
      klv.emplace_back(x,d_zero);

    // (any ascents for x that are descents for y must be imaginary type II)
    new_recursion_column(klv,y,hash); // put new-recursion result into |klv|
    sparseness +=
      remove_zeros(klv,y); // write |d_KL[y]|, suppressing zeros
  }
  return sparseness;
}

/*
  Put into |klv[i]| the the right-hand sides of the recursion formulae for the
  elements |x=e[i]| with |y|, corresponding to the descent |s| for |y|. Here |e|
  contains the block elements extremal for |y| (so in particular their length is
  less than that of |y| and |s| is a descent for all of them), and |s| is either
  a complex, or a real type I descent for |y|. The formula takes the form

    P_{x,y} = (c_s.c_{y'})-part - correction term

  where y' = cross(s,y) when s is complex for y, one of the two elements in
  inverseCayley(s,y) when s is real. The (c_s.c_{y'})-part depends on what kind
  of descent |s| is for |x|. The correction term comes from $\sum_z mu(z,y1)c_z$
  and is handled by |muCorrection|; the form of the summation depends only on
  |y1| (which is recomputed here), but the summation itself involves polynomials
  $P_{x,z}$ that depend on $x$ as well.
*/
void KL_table::recursion_column(std::vector<KLPol>& klv,
				const PrimitiveColumn& e, // extremals for y
				BlockElt y,
				weyl::Generator s)
{
  klv.resize(e.size());

  BlockElt sy =
    descentValue(s,y) == DescentStatus::ComplexDescent ? cross(s,y)
    : inverse_Cayley(s,y).first;  // s is real type I for y here, ignore .second

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
	{ // $(q+1)P_{x,sy}$
	  klv[i] = KL_pol(x,sy);
	  klv[i].safeAdd(klv[i],1); // mulitply by $1+q$
	}
	break;
      case DescentStatus::ComplexDescent:
	{ // $P_{sx,sy}+q.P_{x,sy}$
	  BlockElt sx = cross(s,x);
	  klv[i] = KL_pol(sx,sy);
	  klv[i].safeAdd(KL_pol(x,sy),1);
	}
	break;
      case DescentStatus::RealTypeI:
	{ // $P_{sx.first,sy}+P_{sx.second,sy}+(q-1)P_{x,sy}$
	  BlockEltPair sx = inverse_Cayley(s,x);
	  klv[i] = KL_pol(sx.first,sy);
	  klv[i].safeAdd(KL_pol(sx.second,sy));
	  KLPolRef Pxsy = KL_pol(x,sy);
	  klv[i].safeAdd(Pxsy,1);
	  klv[i].safeSubtract(Pxsy); // subtraction must be last
	}
	break;
      case DescentStatus::RealTypeII:
	{ // $P_{sx,sy}+qP_{x,sy}-P_{s.x,sy}$
	  BlockElt sx = inverse_Cayley(s,x).first;
	  klv[i] = KL_pol(sx,sy);
	  klv[i].safeAdd(KL_pol(x,sy),1);
	  klv[i].safeSubtract(KL_pol(cross(s,x),sy)); // subtraction must be last
	}
	break;
      default: assert(false); // this cannot happen
      }
      // for now |klv[i].degree()| might be one notch too high, which will be
      // corrected in |mu_correction|; also |assert| there are based on this one
      assert(klv[i].isZero() or 2*klv[i].degree()<length(y)-length(x) or
	     (2*klv[i].degree()==length(y)-length(x) and
	      klv[i][klv[i].degree()]==mu(x,sy)
	     ));
    } // |for (i=e.size()-->0)|

  }
  catch (error::NumericUnderflow& err){
    throw kl_error::KLError(e[i],y,__LINE__,
			    static_cast<const KL_table&>(*this));
  }

  mu_correction(klv,e,y,s); // subtract mu-correction from all of |klv|

} // |KL_table::recursion_column|

/*
  Subtract from all polynomials in |klv| the correcting terms in the
  K-L recursion.

  Precondtion: |klv| already contains, for all $x$ that are extremal for |y|,
  which are listed in |e| in increasing order, the terms in $P_{x,y}$
  corresponding to $c_s.c_{y'}$, whery |y'| is $s.y$ if |s| is a complex
  descent, and |y'| is an inverse Cayley transform of |y| if |s| is real type I.
  The mu-table and KL-table have been filled in for elements of length < l(y).

  The recursion formula is of the form:
  $$
    lhs = c_s.c_{y'} - \sum_{z} mu(z,y')c_z
  $$
  where |z| runs over the elements $< y'$ such that |s| is a descent for |z|.
  Here $lhs$ stands for $c_y$ when |s| is a complex descent or real type I for
  |y|, and for $c_{y}+c_{s.y}$ when |s| is real type II; however it plays no
  part in this function that only subtracts $\mu$-terms.


  The element $y'$ is called |sy| in the code below.

  We construct a loop over |z| first, before traversing |klv| (the test for
  $z<sy$ is absent, but $\mu(z,sy)\neq0$ implies $z<sy$ (strict, as mu(sy,sy) is
  0; in any case no coefficient for |sy| is stored in |d_mu[sy]|, and moreover
  $z=sy$ would be rejected by the descent condition). The choix have the out
  loop over $z$ and the inner loop over $x$ (i.e., over |klv|) allows fetching
  $\mu(z,sy)$ only once, and terminating each scan of |klv| once its values |x|
  become too large to produce a non-zero $P_{x,z}$. (In fact we stop once
  $l(x)=l(z)$, and separately consider the case $x=z$.) Either direction of the
  loop on $z$ would work, but taking it decreasing is more natural; we keep
  track of the index |zi| at which $z$ occurs in |e|, if it does.

  Elements of length at least $l(sy)=l(y)-1$ on the list |e| are always
  rejected, so the tail of |e| never reached.
 */
void KL_table::mu_correction(std::vector<KLPol>& klv,
			     const PrimitiveColumn& e,
			     BlockElt y, weyl::Generator s)
{
  BlockElt sy =
    descentValue(s,y) == DescentStatus::ComplexDescent ? cross(s,y)
    : inverse_Cayley(s,y).first;  // s is real type I for y here, ignore .second

  const Mu_column& mcol = d_mu[sy];
  size_t ly = length(y);

  size_t inx_z=e.size(); // should satisfy |e[inx_z]==z| whenever that exists

  size_t j; // define outside for error reporting
  try {
    for (auto it = mcol.rbegin(); it!=mcol.rend(); ++it) // makes |z| decreasing
      if (DescentStatus::isDescent(descentValue(s,it->x))) // |s| descent for |z|
      {
	BlockElt z = it->x;
	MuCoeff mu = it->coef; // $\mu(z,sy)$, which is nonzero

	size_t lz = length(z);
	polynomials::Degree d = (ly-lz)/2; // power of |q| used below

	if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
	  for (j=0; j<e.size(); ++j)
	  {
	    BlockElt x = e[j];
	    if (length(x) >= lz) break; // once reached, $x=z$ is only case left
	    KLPolRef pol = KL_pol(x,z);
	    klv[j].safeSubtract(pol,d); // subtract q^d.P_{x,z} from klv[j]
	  } // |for (j)|
	else // (rare) case that |mu>1|
	  for (j=0; j<e.size(); ++j)
	  {
	    BlockElt x = e[j];
	    if (length(x) >= lz) break; // once reached, $x=z$ is only case left
	    KLPolRef pol = KL_pol(x,z);
	    klv[j].safeSubtract(pol,d,mu); // subtract q^d.mu.P_{x,z} from klv[j]
	  } // |for (j)|

	while (inx_z>0 and e[inx_z-1]>=z)
	  --inx_z; // ensure |e[k]>=z| if and only if |k>=inx_z|

	if (inx_z<e.size() and e[inx_z]==z) // handle final term |x==z|
	{ // none of the larger |z| should have altered leading coefficient
	  assert( klv[inx_z].degree()==d and klv[inx_z][d]==mu );
	  klv[inx_z].safeSubtract(KLPol(d,mu)); // subtract off the term $mu.q^d$
	}

      } // |for (it->reverse(mcol))|
  }
  catch (error::NumericUnderflow& err){
    throw kl_error::KLError(e[j],y,__LINE__,
			    static_cast<const KL_table&>(*this));
  }

} // |KL_table::mu_correction|

/*
  Write down column |y| in |d_KL|.

  Precondition: The parallel pair (klv,er) records the polynomials for at
  least all the extremal values $x$ for $y$, and at most for all primitive
  values $x$ for $y$. So when $x=er[i]$ then $P_{x,y}=klv[i]$. In practice
  |er| will contain either all extremal elements (when called from
  |recursion_column|) or all primitive elements (for |new_recursion_column|).

  This function writes out these data to |d_KL[y]|, transformed as follows:
  (1) rather than storing polynomials from |klv| (or others computed here),
  these are looked up in |d_hashtable| and the index is stored; (2) when a
  polynomial turns out to be 0, nothing is recorded in |d_KL[y]|; (3) for
  primitive elements not present in |er|, the polynomial is computed here
  on-the-fly (using an imaginary type II ascent that exists in this case) and
  then stored along with those from |er|.

  Case (3) will not apply if |er| already contains all primitive elements, and
  for that case this function could be considerably simplified, but it works
  well as is, so we didn't write a simplified version.
 */
size_t KL_table::write_column(const std::vector<KLPol>& klv,
			      const PrimitiveColumn& ec, BlockElt y,
			      KLHash& hash)
{
  auto pc = primitive_column(y); // the elements for which we must write an entry

  KL_column nz_KL(pc.size()); // columns of (primitive entry,nonzero KL poly)

  auto nz_KL_p = nz_KL.end();

  Mu_list mu_pairs; // those |x| with |mu(x,y)>0|

  unsigned int ly = length(y);

  size_t j= ec.size()-1; // points to current extremal element

  for (size_t i = pc.size(); i-->0; )
    if (j<ec.size() and pc[i]==ec[j]) // must test for underflow |j|
    { // extremal element; use stored polynomial
      const KLPol& Pxy=klv[j--];
      if (not Pxy.isZero())
      {
	*--nz_KL_p = KL_pair(pc[i],hash.match(Pxy));
	unsigned int lx=length(pc[i]);
	if (ly==lx+2*Pxy.degree()+1) // in particular parities |lx|, |ly| differ
	  mu_pairs.emplace_front(pc[i],MuCoeff(Pxy[Pxy.degree()]));
      }
    }
    else // insert a polynomial for primitive non-extremal pc[i] if nonzero
    {
      unsigned int s = ascent_descent(pc[i],y);
      assert(descentValue(s,pc[i])==DescentStatus::ImaginaryTypeII);
      BlockEltPair xs = cayley(s,pc[i]);
      KLPol Pxy = KL_pol(xs.first ,y,nz_KL_p,nz_KL.end()); // look up using KL
      Pxy.safeAdd(KL_pol(xs.second,y,nz_KL_p,nz_KL.end()));
      if (not Pxy.isZero())
	*--nz_KL_p = KL_pair(pc[i],hash.match(Pxy));
      // no need to check for |mu| here: |down_set| has the only possible cases
    }

  if (ly==0)
    return 0; // nothing left to do at minimal length

  Mu_list downs;
  for (BlockElt x : down_set(block(),y))
    downs.emplace_back(x,MuCoeff(1));
  // These $x$s are non-extremal for $y$, yet have $\mu(x,y)=1\neq0$

  // some elements in |mu_pairs| may have same length as |downs|: merge is needed
  mu_pairs.merge(std::move(downs)); // need not call |unique|: sets are disjoint

  // commit
  d_KL[y].assign(nz_KL_p,nz_KL.end()); // copy, shifting, from |nz_KL_p| onward
  d_mu[y].assign(mu_pairs.wcbegin(),mu_pairs.wcend()); // convert to vector

  return nz_KL_p - nz_KL.begin(); // measure unused space

} // |KL_table::write_column|

// this method is called instead of |write_column| when new recursion applies
size_t KL_table::remove_zeros(const KL_column& klv, BlockElt y)
{
  KL_column nz_KL; // columns of (primitive entry,nonzero KL poly)

  for (const auto& pair : klv)
    if (not isZero(pair.P))
      nz_KL.push_back(pair);

  // commit
  d_KL[y] = std::move(nz_KL);

  return klv.size()-nz_KL.size(); // measure unused space

} // |KL_table::remove_zeros|

/*
  Puts in klv[i] the polynomial P_{e[i],y} for every primtitve x=pc[i],
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
void KL_table::new_recursion_column
( KL_column& klv, // with entries primitive elements of length |< length(y)|
  BlockElt y,
  KLHash& hash)
{
  auto kl_p = klv.end();

  unsigned int l_y = length(y);

  Mu_list mu_pairs; // those |x| with |mu(x,y)>0|

  // start off |mu_pairs| with ones for |down_set(y)|, not otherwise computed
  for (BlockElt x : down_set(block(),y))
    mu_pairs.emplace_back(x,MuCoeff(1)); // initial part |mu_pairs| is increasing

  // remainder of |mu_pairs| will be decreasing by |x|; it does not matter that
  // all |mu_pairs| is not decreasing by |x|, but it must be decreasing by length
  const auto downs_end = mu_pairs.end(); // record separation for final sorting

  size_t j = klv.size(); // declare outside try block for error reporting
  try {
    while (--kl_p, j-->0) // |*kl_p=klv[j]|
    {
      BlockElt x = klv[j].x;

      unsigned int s= ascent_descent(x,y);
      if (s<rank()) // a primitive element that is not extremal; easy case
      { // equation (1.9) in recursion.pdf
	assert(descentValue(s,x)==DescentStatus::ImaginaryTypeII);
	BlockEltPair p = cayley(s,x);
	KLPol pol = KL_pol(p.first,y,kl_p,klv.end());
	pol.safeAdd(KL_pol(p.second,y,kl_p,klv.end()));
	klv[j].P = hash.match(pol);
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
	KLPol pol = mu_new_formula(x,y,s,mu_pairs);

	switch (descentValue(s,x))
	{
	case DescentStatus::ComplexAscent:
	{ // use equations (3.3a)=(3.4)
	  BlockElt sx = cross(s,x);
	  pol.safeSubtract(KL_pol(sx,y,kl_p,klv.end()),1);
	  // subtract qP_{sx,y} from mu terms
	} // ComplexAscent case
	break;

	case DescentStatus::ImaginaryTypeII:
	{ // use equations (3.3a)=(3.5)
	  BlockEltPair p = cayley(s,x);
	  KLPol sum = KL_pol(p.first,y,kl_p,klv.end());
	  sum.safeAdd(KL_pol(p.second,y,kl_p,klv.end()));
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
	  pol.safe_quotient_by_1_plus_q(length(y)-length(x));
	  break;

	default: assert(false); //we've handled all possible NiceAscents
	}
	klv[j].P = hash.match(pol);
	if (l_y==l_x+2*pol.degree()+1)
	  mu_pairs.emplace_back(x,pol[pol.degree()]);

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
	  KLPol pol = mu_new_formula(x,y,s,mu_pairs);

	  //subtract (q-1)P_{xprime,y} from terms of expression (3.4)
	  BlockElt xprime = cayley(s,x).first;
	  const KLPol& P_xprime_y =  KL_pol(xprime,y,kl_p,klv.end());
	  pol.safeAdd(P_xprime_y);
	  pol.safeSubtract(P_xprime_y,1);

	  //now klv[j] holds P_{x,y}+P_{s.x,y}

	  unsigned int t=st.second;

	  if (t<rank()) // nothing to subtract if $s.x$ not in partial block
	  {
	    //compute P_{s.x,y} using t
	    BlockEltPair sx_up_t = cayley(t,cross(s,x));

	    // any |UndefBlock| component of |sx_up_t| will contribute $0$
	    pol.safeSubtract(KL_pol(sx_up_t.first,y,kl_p,klv.end()));
	    pol.safeSubtract(KL_pol(sx_up_t.second,y,kl_p,klv.end()));

	  }

	  klv[j].P = hash.match(pol);
	  if (l_y==l_x+2*pol.degree()+1)
	    mu_pairs.emplace_back(x,pol[pol.degree()]);
	} // |if (endgame_pair(x,y)) |
	else // |first_endgame_pair| found nothing
	  klv[j].P=d_zero;
      } // end of no NiceAscent case
    } // while (j-->0)
  }
  catch (error::NumericUnderflow& err) // repackage error, reporting x,y
  {
    throw kl_error::KLError(klv[j].x,y,__LINE__,
			    static_cast<const KL_table&>(*this));
  }

  {
    Mu_list downs; // set apart initial part which is increasing
    downs.splice(downs.begin(),mu_pairs,mu_pairs.begin(),downs_end);
    mu_pairs.reverse(); // remainder was decreasing, so make it increasing
    mu_pairs.merge(std::move(downs)); // fusion with initial increasing part
  }
  d_mu[y].assign(mu_pairs.wcbegin(),mu_pairs.wcend());

} // |KL_table::new_recursion_column|

/*
  Store into |klv[j]| the $\mu$-sum appearing in a new K-L recursion.

  Here |pc| is the primitive column for |y|, $s$ is real nonparity for $y$ and
  either C+ or imaginary for $x=pc[j]$ (those are the cases for which the
  formula is used; the status w.r.t. $x$ is not actually used by the code), and
  for all $k>j$ one already has stored $P_{pc[k],y}$ in |klv[k]|.

  The mu-table and KL-table have been filled in for elements of length < l(y),
  so that for $z<y$ we can call |KL_pol(x,z)|.

  The various recursion formulas involve a sum:
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
KLPol KL_table::mu_new_formula
  (BlockElt x, BlockElt y, weyl::Generator s, const Mu_list& mu_y)
{
  KLPol pol=Zero;

  unsigned int lx=length(x), ly = length(y);

  try
  {
    for (const auto& pair : mu_y) // a block element with $\mu(z,y)\neq0$
    {
      BlockElt z = pair.x;
      unsigned int lz = length(z);
      if (lz<=lx)
	break; // length |z| decreases, and |z==x| must be excluded, so stop
      if (not DescentStatus::isDescent(descentValue(s,z))) continue;

      // now we have a true contribution with nonzero $\mu$
      unsigned int d = (ly - lz +1)/2; // power of $q$ used in the formula
      MuCoeff mu = pair.coef;
      KLPolRef Pxz = KL_pol(x,z); // which is known because $z<y$

      if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
	pol.safeAdd(Pxz,d); // add $q^d.P_{x,z}$ to |pol|
      else // mu!=MuCoeff(1)
	pol.safeAdd(Pxz,d,mu); // add $q^d.\mu(z,y).P_{x,z}$ to |pol|

    } // for (k)
  }
  catch (error::NumericOverflow& e){
    throw kl_error::KLError(x,y,__LINE__,
			    static_cast<const KL_table&>(*this));
  }

  return pol;
} // |KL_table::muNewFormula|


void KL_table::silent_fill(BlockElt last_y)
{
  try
  {
    KLHash hash(d_store); // (re-)construct a hastable for polynomial storage
    // fill the lists
    for (auto it = d_holes.begin(); it() and *it<=last_y; ++it)
    {
      fill_KL_column(*it,hash);
      d_holes.remove(*it);
    }
    // after all columns are done the hash table is freed, only the store remains
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
  New routine that does verbose filling of existing |KL_table| object
*/
void KL_table::verbose_fill(BlockElt last_y)
{
  try
  {
    KLHash hash(d_store);

    size_t minLength = length(first_hole()); // length of first new |y|
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
      BlockElt y_start = l==minLength ? first_hole() : lengthLess(l);
      BlockElt y_limit = l<maxLength ? lengthLess(l+1) : last_y+1;
      for (BlockElt y=y_start; y<y_limit; ++y)
      {
	std::cerr << y << "\r";

	nr_of_prim_nulls += fill_KL_column(y,hash);
	prim_size += d_KL[y].size();
	d_holes.remove(y);
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

void KL_table::swallow (KL_table&& sub, const BlockEltList& embed, KLHash& hash)
{
#ifndef NDEBUG
  check_sub(sub,embed);
#endif
  std::vector<KLIndex> poly_trans(sub.d_store.size());;
  for (KLIndex i=0; i<sub.d_store.size(); ++i)
    poly_trans[i]=hash.match(sub.d_store[i]); // should also extend |d_store|

  for (BlockElt z=0; z<sub.block().size(); ++z)
    if (d_holes.isMember(embed[z]))
      {
	for (auto& entry : sub.d_KL[z])
	{
	  entry.x=embed[entry.x];      // renumber block elements
	  entry.P=poly_trans[entry.P]; // renumber polynomial indices
	}
	d_KL[embed[z]] = std::move(sub.d_KL[z]);

	for (auto& entry : sub.d_mu[z])
	  entry.x = embed[entry.x]; // renumber block elements (coef unchanged)
	d_mu[embed[z]] = std::move(sub.d_mu[z]);

	d_holes.remove(embed[z]);
      }
}


/*****************************************************************************

        Chapter V -- Functions declared in kl.h

 *****************************************************************************/


/*
  Return the W-graph for this block.

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

  The edge lists are constructed in already sorted order: for a given element,
  the outgoing edges to smaller elements are first constructed when |y| equals
  that element (and they come in increasing order because |mcol| has its
  first components (|x|) increasing), then to larger elements when the given
  element occurs as |x| for another as |y|; the |y| are always increasing.

*/
wgraph::WGraph wGraph(const KL_table& kl_tab)
{
  wgraph::WGraph wg(kl_tab.rank(),kl_tab.size());

  // fill in descent sets, edges and coefficients
  for (BlockElt y = 0; y < kl_tab.size(); ++y)
  {
    const RankFlags& d_y = kl_tab.descentSet(y);
    wg.descent_sets[y] = d_y;
    const Mu_column& mcol = kl_tab.mu_column(y);
    for (size_t j = 0; j < mcol.size(); ++j)
    {
      BlockElt x = mcol[j].x;
      assert(x<y); // this is a property of |mu_column|
      const RankFlags& d_x = kl_tab.descentSet(x);
      if (d_x == d_y)
	continue;
      MuCoeff mu = mcol[j].coef;
      if (kl_tab.length(y) - kl_tab.length(x) > 1)
      { // nonzero $\mu$, unequal descents, $l(x)+1<l(y)$: edge from $x$ to $y$
	wg.oriented_graph.edgeList(x).push_back(y);
	wg.coefficients[x].push_back(mu);
	continue;
      }
      // now length difference is 1: edges except to a larger descent set
      if (not d_y.contains(d_x)) // then add edge from $x$ to $y$
      {
	wg.oriented_graph.edgeList(x).push_back(y);
	wg.coefficients[x].push_back(mu);
      }
      if (not d_x.contains(d_y)) // then add edge from $y$ to $x$
      {
	wg.oriented_graph.edgeList(y).push_back(x);
	wg.coefficients[y].push_back(mu);
      }
    }
  }
  return wg;
} // |wGraph|

} // |namespace kl|
} // |namespace atlas|
