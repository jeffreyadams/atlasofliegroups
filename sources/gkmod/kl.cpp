/*
  This is kl.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2012 David Vogan
  Copyright (C) 2005-2020, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  Implementation of the class |KL_table|.

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

  Since |d_KL| holds all polynomials for primitive pairs $(x,y)$, this is just a
  lookup function. Find the index |inx| of the primitivisation of |x| in the
  column |kl_col==d_kl[y]| (done in using quick lookup in |prim_index|). If this
  has made |inx| out of bounds (in particular if |inx==UndefBlock|) return a
  zero polynomial. Otherwise fetch from |d_KL| and |d_store|.
  */
KLPolRef KL_table::KL_pol(BlockElt x, BlockElt y) const
{
  if (x==UndefBlock) // partial blocks can cause this in many ways
    return d_store[d_zero];
  const auto& kl_col = d_KL[y];
  unsigned int inx = prim_index(x,descent_set(y));

  if (inx>=kl_col.size()) // l(x)>=l(y), includes case x==-1: no primitivization
    return d_store[inx==self_index(y) ? d_one : d_zero];
  return d_store[kl_col[inx]];
}

// The same, but just return the index into |d_store| that gives $P_{x,y}$
KLIndex KL_table::KL_pol_index(BlockElt x, BlockElt y) const
{
  if (x==UndefBlock) // partial blocks can cause this in many ways
    return d_zero;
  const auto& kl_col = d_KL[y];
  unsigned int inx = prim_index(x,descent_set(y));

  if (inx>=kl_col.size()) // l(x)>=l(y), includes case |inx==(unsigned)-1|
    return inx==self_index(y) ? d_one : d_zero;
  return kl_col[inx];
}

/*
  Return $\mu(x,y)$.

  This function is not used internally, so we are sure all tables are computed.
  We prefer fast look-up in |d_KL| (via |klPol|) over binary search in |d_mu|.
*/
MuCoeff KL_table::mu(BlockElt x, BlockElt y) const
{
  KLPolRef p = KL_pol(x,y);
  return p.isZero() or 2*p.degree()+1+length(x)<length(y) ? MuCoeff(0)
						: MuCoeff(p[p.degree()]);
}


/*
  Return the list of all |x| primitive w.r.t. |y|.

  Explanation: this means that |length(x) < length(y)|, and every descent
  for |y| is either a descent, or an imaginary type II ascent for |x|.
*/
BitMap KL_table::primitives (BlockElt y) const
{
  const BlockElt limit = length_floor(y);
  const RankFlags desc_y = descent_set(y);
  BitMap result(limit);
  BlockElt x=limit;
  while (prim_back_up(x,desc_y))
    result.insert(x);
  return result;
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

BitMap KL_table::prim_map (BlockElt y) const
{
  // the vector of polynomial indices at primitive elements |x|, all with |y|
  const auto& col = d_KL[y];
  BitMap result(col.size());
  for (unsigned int i=0;  i<col.size(); ++i)
    result.set_to(i,col[i]!=d_zero);

  return result;
}


/*****************************************************************************

        Chapter II -- Private methods used during construction

 *****************************************************************************/

/*
  Return the first generator that is a strict descent other than real type II

  These are the ones that give a direct recursion formula for the K-L basis
  element. Concretely, |first_direct_recursion(y)| returns the first generator
  |s| such that |descentValue(s,y)| is either |DescentStatus::ComplexDescent| or
  |DescentStatus::RealTypeI|. If no such generator exists, it returns |rank()|.
*/
weyl::Generator KL_table::first_direct_recursion(BlockElt y) const
{
  const DescentStatus& d = descent(y);
  weyl::Generator s;
  for (s=0; s<rank(); ++s)
    if (DescentStatus::isDirectRecursion(d[s])) // complex descent or real type1
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
  When this routing is called, it has been observed that:
  * any true descent for y is of type r2 (|first_direct_recursion| has failed)
  * |x| is extremal for |y| (none of the weak descents for y are ascents for x)
  * none of the rn ascents for y is C+, ic or i2 for x (|first_nice_and_real|
    failed to find anything)

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


// private manipulators

/*
  Fill the column for |y| in the KL-table.

  Precondition: all lower rows have been filled

  Column of $y$ is the set of all $P_{x,y}$ for $x<y$
*/
void KL_table::fill_KL_column(BlockElt y, KLHash& hash)
{
  if (d_KL[y].size()>0)
    return; // column has already been filled

  prepare_prim_index(descent_set(y)); // so looking up |KL_pol(x,y)| will be OK

  weyl::Generator s = first_direct_recursion(y);
  if (s<rank())  // a direct recursion was found, use it for |y|, for all |x|
  {
    std::vector<KLPol> klv = // compute all polynomials for extremal |x| for |y|
      recursion_column(y,s);
    // write result
    complete_primitives(klv,y,hash); // add primitives; store in |d_KL|
  }
  else // we must use an approach that distinguishes on |x| values
    d_KL[y] = new_recursion_column(y,hash); // compute and install column
} // |KL_table::fill_KL_column|

/*
  Put into |klv[x]| the the right-hand sides of the recursion formulae for the
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
std::vector<KLPol> KL_table::recursion_column (BlockElt y,weyl::Generator s)
{
  std::vector<KLPol> klv(block().size()+1); // storage, indexed by |BlockElt|
  klv[y]=One; // ensure diagonal entry of it taken to be $P_{y,y}=1$

  const RankFlags desc_y = descent_set(y);
  const BlockElt sy =
    descent_value(s,y) == DescentStatus::ComplexDescent ? cross(s,y)
    : inverse_Cayley(s,y).first;  // s is real type I for y here, ignore .second

  // make increasing list of all extremal elements shorter than |y|
  BlockEltList extremals; extremals.reserve(col_size(y)); // more than enough
  for (BlockElt x=0; x<length_floor(y); ++x)
    if (is_extremal(x,desc_y))
      extremals.push_back(x);

  BlockElt x = // declaration outside is needed for error reporting
    length_floor(y);
  try {

    // while more natural to do descending |x| doing forward avoids |std:reverse|
    for (auto it=extremals.cbegin(); it!=extremals.cend(); ++it)
    { // now |x| is extremal for $y$, so $s$ is descent for $x$
      KLPol& Pxy = klv[x=*it];
      switch (descent_value(s,x))
      {
      case DescentStatus::ImaginaryCompact:
	{ // $(q+1)P_{x,sy}$
	  Pxy = KL_pol(x,sy);
	  Pxy.safeAdd(Pxy,1); // mulitply by $1+q$
	}
	break;
      case DescentStatus::ComplexDescent:
	{ // $P_{sx,sy}+q.P_{x,sy}$
	  BlockElt sx = cross(s,x);
	  Pxy = KL_pol(sx,sy);
	  Pxy.safeAdd(KL_pol(x,sy),1);
	}
	break;
      case DescentStatus::RealTypeI:
	{ // $P_{sx.first,sy}+P_{sx.second,sy}+(q-1)P_{x,sy}$
	  BlockEltPair sx = inverse_Cayley(s,x);
	  Pxy = KL_pol(sx.first,sy);
	  Pxy.safeAdd(KL_pol(sx.second,sy));
	  KLPolRef Pxsy = KL_pol(x,sy);
	  Pxy.safeAdd(Pxsy,1);
	  Pxy.safeSubtract(Pxsy); // subtraction must be last
	}
	break;
      case DescentStatus::RealTypeII:
	{ // $P_{sx,sy}+qP_{x,sy}-P_{s.x,sy}$
	  BlockElt sx = inverse_Cayley(s,x).first;
	  Pxy = KL_pol(sx,sy);
	  Pxy.safeAdd(KL_pol(x,sy),1);
	  Pxy.safeSubtract(KL_pol(cross(s,x),sy)); // subtraction must be last
	}
	break;
      default: assert(false); // this cannot happen
      }
      // for now |Pxy.degree()| might be one notch too high, which will be
      // corrected in |mu_correction|; also |assert| there are based on this one
      assert(Pxy.isZero() or 2*Pxy.degree()<length(y)-length(x) or
	     (2*Pxy.degree()==length(y)-length(x) and
	      Pxy[Pxy.degree()]==mu(x,sy)
	     ));
    } // |for (i=e.size()-->0)|

  }
  catch (error::NumericUnderflow& err){
    throw kl_error::KLError(x,y,__LINE__,
			    static_cast<const KL_table&>(*this));
  }

  // now subtract mu-corrections from all of |klv|
  mu_correction(klv,extremals,desc_y,sy,s);
  return klv;

} // |KL_table::recursion_column|

/*
  Subtract from all polynomials in |klv| the correcting terms in the
  K-L recursion.

  When we call |mu_correction|, the polynomial |klv[x]| already contains, for
  all $x$ that are extremal for |y| (the members of |e|), the terms in $P_{x,y}$
  corresponding to $c_s.c_{y'}$, where |y'| is an |s| descent of |y| as before.
  The tables |d_KL| and |d_mu| have been filled in for elements of length < l(y).

  The recursion formula is of the form:
  $$
    lhs = c_s.c_{y'} - \sum_{z} mu(z,y')c_z
  $$
  where $y'$ is the |s|-descent of |y| passed as argument |sy|, with the sum
  over |z| runing over the elements $< y'$ such that |s| is a descent for |z|.
  (Here $lhs$ stands for $c_y$ when |s| is a complex descent or real type I for
  |y|, and for $c_{y}+c_{s.y}$ when |s| is real type II; however it plays no
  part in this function that only subtracts $\mu$-terms.)

  We construct a loopfirst over those |z| for which $\mu(z,y')$ is nonzero
  (which implies $z<y'$) and for which |s| is a descent, before traversing |e|
  for the values of |x| for which |klv[x]| needs correction. This allows
  fetching $\mu(z,sy)$ only once, and terminating each inner loop once |x|
  becomes too large to produce a non-zero $P_{x,z}$. (In fact we stop once
  $l(x)=l(z)$, and separately consider the possibility $x=z$ with $P_{x,z}=1$.)
  Either direction of the loop on $z$ would work, but taking it decreasing is
  more natural.
 */
void KL_table::mu_correction(std::vector<KLPol>& klv,
			     const BlockEltList& extremals,
			     RankFlags desc_y, BlockElt sy, weyl::Generator s)
{
  const Mu_column& mcol = d_mu[sy];
  size_t ly = length(sy)+1; // the length of |y|, otherwise |y| is not used here

  BlockElt xx=UndefBlock; // define outside for error reporting
  try {
    for (auto it = mcol.rbegin(); it!=mcol.rend(); ++it) // makes |z| decreasing
      if (DescentStatus::isDescent(descent_value(s,it->x))) // descent for |z|?
      {
	BlockElt z = it->x;
	MuCoeff mu = it->coef; // $\mu(z,sy)$, which is nonzero

	size_t lz = length(z);
	polynomials::Degree d = (ly-lz)/2; // power of |q| used below

	if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
	  for (BlockElt x : extremals)
	  {
	    if (length(x)>=lz)
	      break; // now only case |x==z| possibly needs consideration
	    KLPolRef pol = KL_pol(x,z);
	    klv[xx=x].safeSubtract(pol,d); // subtract $q^d.P_{x,z}$ from klv[x]
	  }
	else // (rare) case that |mu>1|
	  for (BlockElt x : extremals)
	  {
	    if (length(x)>=lz)
	      break; // now only case |x==z| possibly needs consideration
	    KLPolRef pol = KL_pol(x,z);
	    klv[xx=x].safeSubtract(pol,d,mu); // subtract $q^d.mu.P_{x,z}$
	  }

	if (is_extremal(z,desc_y)) // then handle final term |x==z|
	{ // none of the larger |z| should have altered the leading coefficient
	  assert( klv[z].degree()==d and klv[z][d]==mu );
	  klv[xx=z].safeSubtract(KLPol(d,mu)); // subtract off the term $mu.q^d$
	}

      } // |for (it->reverse(mcol))| |if(isDescent(descentValue(s,it->x))|
  }
  catch (error::NumericUnderflow& err){
    BlockElt y = block().unique_ascent(s,sy); // reconstruct |y| uniquely
    throw kl_error::KLError(xx,y,__LINE__,
			    static_cast<const KL_table&>(*this));
  }

} // |KL_table::mu_correction|

/* A method that takes a row |klv| of completed KL polynomials, computed by
   |recursion_column| at |y| and extremal elements |x| listed in |ext|, and
   transfers them to the main storage structures. Its tasks are

   - generate the list of all primitve elements for |y|, which contains |ext|
   - for each primitive element |x|, if it is extremal just look up $P_{x,y}$
     as |klv[x]|, if |x| is primitive but not extremal, compute that polynomial
     (as sum of two $P_{x',y}$ in the same column); hash and store the result
   - record nonzero $P_{x,y}$ as $(x,P)$ and similarly any non-zero $\mu(x,y)$

   For the latter point there are two categories of |x|: the extremal ones
   (which can conveniently be handled in the loop over |x|), and those found
   by a (complex or real) descent from |y| itself (they have $\mu(x,y)=1$).
   The latter are of length one less than |y| (but there can be extremal |x|
   of that length as well with nonzero mu), and are primitive only in the real
   type 2 case; we must treat them outside the loop over primitive elements.
 */
void KL_table::complete_primitives(std::vector<KLPol>& klv,
				     BlockElt y, KLHash& hash)
{
  KL_column& KL = d_KL[y]; // the column that we must write to
  KL.resize(col_size(y)); // create slots for all pertinent elements |x|

  Mu_list mu_pairs; // those |x| with |mu(x,y)>0|
  const unsigned int ly = length(y);
  const RankFlags desc_y = descent_set(y);

  auto KL_it = KL.rbegin(); // prepare for writing |KL| backwards
  // traverse primitives for |y| of length |y| less than |ly| backwards
  for(BlockElt x=length_floor(y); prim_back_up(x,desc_y);++KL_it)
    if (is_extremal(x,desc_y))
    { // extremal element for |y|; use polynomial from vector passed to us
      const KLPol& Pxy = klv[x];
      *KL_it = hash.match(Pxy);
      unsigned int lx = length(x);
      if (ly==lx+2*Pxy.degree()+1) // in particular parities |lx|, |ly| differ
	mu_pairs.emplace_front(x,MuCoeff(Pxy[Pxy.degree()]));
    }
    else // must insert a polynomial for primitive non-extramal |x|
    {
      unsigned int s = ascent_descent(x,y);
      assert(descent_value(s,x)==DescentStatus::ImaginaryTypeII);
      BlockEltPair xs = cayley(s,x);
      KLPol Pxy = KL_pol(xs.first,y); // look up P_{x',y} in current row, above
      Pxy.safeAdd(KL_pol(xs.second,y)); // current point, and P_{x'',y} as well
      *KL_it = hash.match(Pxy); // add poly at primitive non-extremal x
    }
  assert(KL_it==KL.rend());

  Mu_list downs;
  auto ds = down_set(block(),y);
  for (auto it=ds.begin(); not ds.at_end(it); ++it)
    downs.emplace_back(*it,MuCoeff(1));
  // These $x$s are non-extremal for $y$, yet have $\mu(x,y)=1\neq0$

  // some elements in |mu_pairs| may have same length as |downs|: merge is needed
  mu_pairs.merge(std::move(downs)); // need not call |unique|: sets are disjoint

  // commit
  d_mu[y].assign(mu_pairs.wcbegin(),mu_pairs.wcend()); // convert to vector
} // |KL_table::complete_primitives|

/*
  Compute polynomials $P_{x,y}$ for all $x$ of length less than and primitive
  for |y|, look them up and return a vector of their indices in |d_store|.

  These KL polynomials are computed by a recursion formula designed for those
  elements |y| for which the direct recursion does not apply.

  When we come here, every simple root |s| is for |y| either a complex ascent or
  imaginary or real (so there are no complex descents for |y|). label:(split 1)

  In fact real type 1 descents for |y| don't occur, but this is not used.

  From that condition we get: for each extremal |x| for |y|, there either exists
  a true ascent |s| that is real for |y|, necessarily nonparity because |x| is
  extremal (split 3), or we are assured that $P_{x,y}=0$.

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

  The sum involving mu, produced by |mu_new_formula|, has terms involving
  $P_{x,u}\mu(u,y}$, so when doing a downward loop over |x| it pays to keep
  track of the previous |u| with nonzero $\mu(u,y)$.

  This code gets executed for |y| that are of minimal length, in which case
  it only contributes $P_{y,y}=1$; the |while| loop will be executed 0 times.
*/
std::vector<KLIndex> KL_table::new_recursion_column(BlockElt y, KLHash& hash)
{
  const unsigned int l_y = length(y);
  const auto desc_y = descent_set(y);

  std::vector<KLIndex> cur_col(nr_of_primitives(desc_y)+1,d_zero);
  cur_col[self_index(y)]=d_one; // everything above will remain |d_zero|
  auto KL_y = [this,&cur_col,desc_y] (BlockElt x) -> KLPol
    { return d_store[cur_col[prim_index(x,desc_y)]]; };

  Mu_list mu_pairs; // those |x| with |mu(x,y)>0|
  auto ds = down_set(block(),y);
  // start off |mu_pairs| with ones for |down_set(y)|, not otherwise computed
  for (auto it=ds.begin(); not ds.at_end(it); ++it)
    mu_pairs.emplace_back(*it,MuCoeff(1)); // initial part |mu_pairs|, increasing

  // remainder of |mu_pairs| will be decreasing by |x|; it does not matter that
  // all |mu_pairs| is not decreasing by |x|, but it must be decreasing by length
  const auto downs_end = mu_pairs.end(); // record separation for final sorting

  BlockElt x = length_less(l_y);
  try {
    auto col_it = // |*col_it| will be entry for $P_{x,y}$
      cur_col.begin()+col_size(y); // the upper bound for written part
    while (prim_back_up(x,desc_y)) // reverse loop through primitive elements
    {
      --col_it; // synchronise with |x|
      unsigned int s= ascent_descent(x,y);
      if (s<rank()) // a primitive element that is not extremal; easy case
      { // equation (1.9) in recursion.pdf
	assert(descent_value(s,x)==DescentStatus::ImaginaryTypeII);
	BlockEltPair p = cayley(s,x);
	KLPol Pxy = KL_y(p.first);
	Pxy.safeAdd(KL_y(p.second));
	*col_it = hash.match(Pxy);
	continue; // done with |x|, go on to the next
      }

      unsigned int l_x = length(x);

      /* now |x| is extremal for |y|. By (split 1) and Lemma 3.1 of recursion.pdf
         this implies that if $x<y$ in the Bruhat order, there is at least one
         |s| real for |y| that is a strict ascent (not rn) for |x| and therefore
         rn for |y|; we first hope that at least one of them is not i1 for |x|
      */
      // first seek a real nonparity ascent for |y| that is C+,i2 or ic for |x|
      s = first_nice_and_real(x,y);
      if (s < rank()) // there is such an ascent s
      {
	// start setting |pol| to the expression (3.4) in recursion.pdf
	KLPol pol = mu_new_formula(x,y,s,mu_pairs);

	switch (descent_value(s,x))
	{
	case DescentStatus::ComplexAscent: // use equations (3.3a)=(3.4)
	  pol.safeSubtract(KL_y(cross(s,x)),1); // subtract qP_{sx,y}
	break;

	case DescentStatus::ImaginaryTypeII:
	{ // use equations (3.3a)=(3.5)
	  BlockEltPair p = cayley(s,x);
	  KLPol sum = KL_y(p.first);
	  sum.safeAdd(KL_y(p.second));
	  pol.safeAdd(sum);
	  pol.safeSubtract(sum,1); //now we've added (1-q)(P_{x',y}+P_{x'',y})
	  pol.safeDivide(2);   //this could throw, but should not
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
	*col_it = hash.match(pol); // record definitive value $P_{x,y}$
	if (l_y==l_x+2*pol.degree()+1)
	  mu_pairs.emplace_back(x,pol[pol.degree()]);

      } // end of |first_nice_and_real| case

      else // there is no Weyl group generator "nice for |x| and real for |y|"
      {
      /*
        The need for the new recusion and the absence of "nice and real"
        generators almost implies $P_{x,y}=0$, but not quite; already in the
        case of C2 there are exceptions. To find them we need to use the idea
        described on p. 8 of recursion.pdf: find $s$ and $t$, both real for $y$
        and imaginary for $x$, moreover being repectively nonparity and parity
        (r2) for $y$ while being repectively i1 and compact for x, while
        moreover $t$ is noncompact for $s.x$ (the imaginary cross image of $x$),
        which can only happen when |t| is adjacent in the Dynkin diagram to $s$.
        If such $(s,t)$ exist, then we can compute $P_{s.x,y}$ using $t$ (since
        an easy recursion, (1.9) but for $t$, expresses it as sum of one or two
        already computed polynomials), while for the sum $P_{sx,y}+P_{x,y}$ we
        have a formula (3.6) of the kind used for NiceAscent; it then suffices
        to compute that formula and subtract $P_{s.x,y}$ from it.

	Finally if no such $(s,t)$ exist, then we have exhausted all
	possibilities where $x$ is below $y$ in the Bruhat order, so we may
	validly conclude that $P_{x,y}=0$.
      */
	auto st = first_endgame_pair(x,y);
	if ((s=st.first) < rank())
	{
	  KLPol pol = mu_new_formula(x,y,s,mu_pairs);

	  //subtract (q-1)P_{xprime,y} from terms of expression (3.4)
	  const auto& P_xprime_y = KL_y(cayley(s,x).first);
	  pol.safeAdd(P_xprime_y);
	  pol.safeSubtract(P_xprime_y,1);

	  //now |pol| holds P_{x,y}+P_{s.x,y}

	  weyl::Generator t = st.second;

	  if (t<rank()) // nothing to subtract if $s.x$ not in partial block
	  {
	    //compute P_{s.x,y} using t
	    BlockEltPair sx_up_t = cayley(t,cross(s,x));

	    // any |UndefBlock| component of |sx_up_t| will contribute $0$
	    pol.safeSubtract(KL_y(sx_up_t.first));
	    pol.safeSubtract(KL_y(sx_up_t.second));
	  }

	  *col_it = hash.match(pol); // record definitive value $P_{x,y}$
	  if (l_y==l_x+2*pol.degree()+1)
	    mu_pairs.emplace_back(x,pol[pol.degree()]);
	} // |if (endgame_pair(x,y)) |
	else // |first_endgame_pair| found nothing
	  assert(*col_it==d_zero); // just check unchanged since initialised
      } // end of no NiceAscent case
    } // while (j-->0)
  }
  catch (error::NumericUnderflow& err) // repackage error, reporting x,y
  {
    throw kl_error::KLError(x,y,__LINE__,
			    static_cast<const KL_table&>(*this));
  }

  {
    Mu_list downs; // set apart initial part which is increasing
    downs.splice(downs.begin(),mu_pairs,mu_pairs.begin(),downs_end);
    mu_pairs.reverse(); // remainder was decreasing, so make it increasing
    mu_pairs.merge(std::move(downs)); // fusion with initial increasing part
  }
  d_mu[y].assign(mu_pairs.wcbegin(),mu_pairs.wcend());

  // slice off the relevant initial part of |cur_col|, for elements in |prims|
  return std::vector<KLIndex>(cur_col.begin(), cur_col.begin()+col_size(y));
} // |KL_table::new_recursion_column|

/*
  Compute the $\mu$-sum appearing in a new K-L recursion.

  Here $s$ is real nonparity for $y$ and either C+ or imaginary for $x$ (those
  are the cases for which the formula is used; the status with respect to $x$ is
  not actually used by the code). The list of $\mu$ values for |y| for elements
  down to $x$ is given as |mu_y|. For $z<y$ we can safely call |KL_pol(x,z)|.

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
      if (not DescentStatus::isDescent(descent_value(s,z))) continue;

      // now we have a true contribution with nonzero $\mu$
      unsigned int d = (ly - lz +1)/2; // power of $q$ used in the formula
      MuCoeff mu = pair.coef;
      KLPolRef Pxz = KL_pol(x,z); // we can look this up because $z<y$

      if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible
	pol.safeAdd(Pxz,d); // add $q^d.P_{x,z}$ to |pol|
      else // mu!=MuCoeff(1)
	pol.safeAdd(Pxz,d,mu); // add $q^d.\mu(z,y).P_{x,z}$ to |pol|

    } // |for (pair : mu_y)|
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
    KLHash hash(d_store,4);

    size_t minLength = length(first_hole()); // length of first new |y|
    size_t maxLength = length(last_y<size() ? last_y : size()-1);

    //set timers for KL computation
    std::time_t time0;
    std::time(&time0);
    std::time_t time;

    struct rusage usage; //holds Resource USAGE report
    size_t storesize = 0; // previous size of d_store
    size_t polsize = 0; // running total of sum of (polynomial degrees+1)

    size_t kl_size = 0;

    for (size_t l=minLength; l<=maxLength; ++l) // by length for progress report
    {
      BlockElt y_start = l==minLength ? first_hole() : length_less(l);
      BlockElt y_limit = l<maxLength ? length_less(l+1) : last_y+1;
      for (BlockElt y=y_start; y<y_limit; ++y)
      {
	std::cerr << y << "\r";

	fill_KL_column(y,hash);
	kl_size += d_KL[y].size();
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
	<< y_limit-1 // last y value done
	<< ", polys:"  << std::setw(11) << d_store.size()
	<< ", mat:"  << std::setw(11) << kl_size
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
		<< std::setw(6) << kl_size*sizeof(KLIndex)/1048576
		<< "MB \n";

    } // for (l=min_length+1; l<=max_Length; ++l)

    std::time(&time);
    double deltaTime = difftime(time, time0);
    std::cerr << std::endl;
    std::cerr << "Total elapsed time = " << deltaTime << "s." << std::endl;
    std::cerr << d_store.size() << " polynomials, "
	      << kl_size << " matrix entries."<< std::endl;

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
    if (not sub.d_holes.isMember(z) and d_holes.isMember(embed[z]))
    { // then transfer |sub.d_KL[z]| and |sub.d_mu[z]| to new block
      RankFlags desc = sub.descent_set(z);
      prepare_prim_index(desc); // first make sure |KLSuport| is ready for |z|
      auto sub_prims = sub.primitives(z);
      auto prims = primitives(embed[z]);
      // we need to convert these |BitMap|s to vectors
      BlockEltList sub_pc(sub_prims.begin(),sub_prims.end());
      BlockEltList pc(prims.begin(),prims.end());
      assert(sub.d_KL[z].size()==sub_pc.size());
      assert(desc == descent_set(embed[z]));
      d_KL[embed[z]].resize(pc.size(),d_zero); // default to |d_zero|
      for (unsigned int i=0; i<sub_pc.size(); ++i)
      {
	unsigned int new_i = prim_index(embed[sub_pc[i]],desc);
	assert(sub.prim_index(sub_pc[i],desc)==i); // |sub_pc[i]| is primitive
	assert(prim_index(pc[new_i],desc)==new_i); // |pc[new_i]| is primitive
	d_KL[embed[z]][new_i] = poly_trans[sub.d_KL[z][i]];
      }

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
    const RankFlags& d_y = kl_tab.descent_set(y);
    wg.descent_sets[y] = d_y;
    const Mu_column& mcol = kl_tab.mu_column(y);
    for (size_t j = 0; j < mcol.size(); ++j)
    {
      BlockElt x = mcol[j].x;
      assert(x<y); // this is a property of |mu_column|
      const RankFlags& d_x = kl_tab.descent_set(x);
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
