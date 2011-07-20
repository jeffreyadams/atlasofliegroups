/*!
\file
\brief Implementation of the functions in namespace kltest.

These functions are designed to test a mathematical assertion used in
the implementation of the KL algorithm.
*/
/*
  This is kltest.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "kltest.h"

#include <iostream>
#include <cassert>

#include "interactive.h"
#include "ioutils.h"
#include "kgb.h"
#include "kl.h"
#include "prettyprint.h"
#include "permutations.h"
#include "weyl.h"

namespace atlas {

namespace {

  void involutionList(TwistedInvolutionList&, const KGB&);

  std::ostream& printBasePts(std::ostream&, const TwistedInvolutionList&,
			     const KGBEltList&, const KGB&);

// the following class was copied integrally from kgb.cpp, helper namespace
class InvolutionCompare
{
private:
  const TwistedWeylGroup* d_W;
public:
  explicit InvolutionCompare(const TwistedWeylGroup& W):d_W(&W) {}

  // one should have a < b iff
  // (a) involutionLength(a) < involutionLength(b) or
  // (b) involutionLengths are equal and length(a) < length (b) or
  // (c) both lengths are equal and a < b
  bool operator()
   (const TwistedInvolution& a, const TwistedInvolution& b) const
  {
    const WeylGroup& W=d_W->weylGroup();
    if      (d_W->involutionLength(a) != d_W->involutionLength(b))
      return d_W->involutionLength(a) <  d_W->involutionLength(b) ;
    else if (W.length(a.w()) != W.length(b.w()))
      return W.length(a.w()) <  W.length(b.w());
    else
      return a < b;
  }
}; // |class InvolutionCompare|

}

/*****************************************************************************

        Chapter I -- Functions declared in kltest.h

******************************************************************************/

namespace kltest {

/*!
  \brief Checks whether the conjectural basepoint in each R-packet is
  independent of the choice of reduced expression.

  Explanation: let w be a twisted involution, written as an involution-reduced
  expression w = s_1 ... s_p (here s_j means either conjugation or commuting
  multiplication). Then the conjecture is that if x is large in the
  fundamental fiber, the element s_1 ... s_p.x, where the action is
  cross-action for a conjugation, cayley transform for a commuting
  multiplication, is independent of the choice of the reduced expression.

  Algorithm: traverse involutions w in increasing order, so descents of w are
  visited before w. Denoting for v<w by x_v the base point for v, we may check
  at w that for all descents sw of w, the KGB element s.x_{sw} is the same; it
  will then define the element x_w.
*/
bool checkBasePoint(const KGB& kgb)
{
#ifdef VERBOSE
  std::cerr << "entering checkBasePoint ..." << std::endl;
#endif

  const TwistedWeylGroup& W = kgb.twistedWeylGroup();
  InvolutionCompare comp(W);
  TwistedInvolutionList wl;
  involutionList(wl,kgb);

  KGBEltList basepts;

  for (size_t x0 = 0; x0 < kgb.size() and kgb.length(x0) == 0; ++x0)
  { { size_t s; // check if x0 is large
      for (s=0; s<kgb.rank(); ++s) // check simple roots
      {
	gradings::Status::Value v = kgb.status(s,x0);
	if (v == gradings::Status::ImaginaryCompact) // none should be compact
	  break;
      }
      if (s<kgb.rank()) // previous loop was broken out of
	continue;
    }

    // if we get here, x0 is large
    // check that basepoint is well-defined
    basepts.assign(wl.size(),UndefKGB);
    basepts[0] = x0;
    for (size_t w_pos = 1; w_pos < wl.size(); ++w_pos)
    {
#ifdef VERBOSE
      std::cerr << w_pos << "\r";
#endif
      const TwistedInvolution& tw = wl[w_pos];

      for (size_t s=0; s<kgb.rank(); ++s)
	if (W.weylGroup().hasDescent(s,tw.w()))  // try all descents
	{
	  KGBElt sx_sw; // will hold candidate basepoint at w
	  if (W.hasTwistedCommutation(s,tw)) // descent is inverse Cayley
	  {
	    TwistedInvolution sw = W.prod(s,tw);
	    size_t sw_pos = // locate index of |sw| using binary search
              std::lower_bound(wl.begin(),wl.end(),sw,comp) -  wl.begin();
	    sx_sw = kgb.cayley(s,basepts[sw_pos]);
	    if (sx_sw == UndefKGB) // Cayley undefined; should not happen
	      return false;
	  }
	  else
	  {
	    TwistedInvolution sw = W.twistedConjugated(tw,s);
	    size_t sw_pos = // locate index of |sw| using binary search
              std::lower_bound(wl.begin(),wl.end(),sw,comp) -  wl.begin();
	    sx_sw = kgb.cross(s,basepts[sw_pos]);
	  }

	  if (basepts[w_pos] == UndefKGB) // x_w is new
	  {
	    basepts[w_pos] = sx_sw;
	    // check if basepoint is large (all simple-imaginaries noncompact)
	    // this is automatic according to David
	    // ... not done yet, need to think a little ... [Fokko]
	  }
	  else if (sx_sw != basepts[w_pos]) // inconsistent base point
	    return false;
	}
    } // |for (w_pos)|
#ifdef VERBOSE
    std::cerr << std::endl;
#endif
  }

#ifdef VERBOSE
  std::cerr << "done" << std::endl;
#endif
  if (basepts.size()==wl.size())
  {
    ioutils::OutputFile file;
    printBasePts(file,wl,basepts,kgb);
  }
  else if (basepts.size()==0)
    std::cout << "No basepoints found" << std::endl;
  return true;
}

void dualityPermutation(permutations::Permutation& a, const kl::KLContext& klc)

{}

bool dualityVerify(const kl::KLContext& klc, const kl::KLContext& dual_klc)
{
  std::vector<BlockElt> dual =
    blocks::dual_map(klc.block(),dual_klc.block());

  permutations::Permutation inv_dual
    (permutations::Permutation(dual.begin(),dual.end()),-1);

  BlockElt block_size=klc.size();

  for (BlockElt x=0; x<block_size; ++x)
  {
#ifdef VERBOSE
    std::cerr << x <<'\r';
#endif
    BlockElt limit=dual_klc.lengthLess(dual_klc.length(dual[x]));
    for (BlockElt dx=0; dx<limit; ++dx)
    {
      assert(dx<inv_dual.size());
      kl::KLPol sum=klc.klPol(x,inv_dual[dx]);

      for (size_t l=klc.length(inv_dual[dx]); l-->klc.length(x)+1; )
      {
	kl::KLPol s=kl::Zero;
	for (BlockElt y=klc.lengthLess(l); y<klc.lengthLess(l+1); ++y)
	  s+= klc.klPol(x,y) * dual_klc.klPol(dx,dual[y]);
	sum.subtract_from(s);
      }

      sum.subtract_from(dual_klc.klPol(dx,dual[x])); // contribution |y==x|

      if (not sum.isZero())
      {
	std::cerr << "at (" << x << ',' << dx << "), non-zero entry :";
	prettyprint::printPol(std::cerr,sum,"q"); std::cerr << std::endl;
	return false;
      }
    }
  }
  return true;
}

} // namespace kltest

/*****************************************************************************

        Chapter I -- Functions local to this module

******************************************************************************/

namespace {

std::ostream&
printBasePts(std::ostream& strm, const TwistedInvolutionList& wl,
	     const KGBEltList& bp, const KGB& kgb)

/*
  Synopsis: outputs the list of basepoints to strm.

  Precondition: wl contains the list of twisted involutions, bp the list of
  basepoints (in the same order.)
*/

{
  using namespace prettyprint;

  for (size_t j = 0; j < wl.size(); ++j) {
    strm << "(";
    printWeylElt(strm,wl[j].w(),kgb.weylGroup());
    strm << ",";
    printInvolution
      (strm,wl[j],kgb.twistedWeylGroup()); // added by jda: red. form of invol.
    strm << "," << bp[j] << ")";
    strm << std::endl;
  }

  return strm;
}

/*
  Synopsis: put in wl the list of twisted involutions that appear in kgb.

  The result is equivalent to the internal d_tau in KGB.
*/
void involutionList(TwistedInvolutionList& wl, const KGB& kgb)
{
  wl.push_back(kgb.involution(0));

  for (size_t x = 0; x < kgb.size(); ++x)
    if (kgb.involution(x) != wl.back())
      wl.push_back(kgb.involution(x));
}


}

}
