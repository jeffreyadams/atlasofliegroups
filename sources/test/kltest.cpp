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

  See file main.cpp for full copyright notice
*/

#include "kltest.h"

#include <iostream>

#include "interactive.h"
#include "ioutils.h"
#include "kgb.h"
#include "kl.h"
#include "prettyprint.h"
#include "setutils.h"
#include "weyl.h"

namespace atlas {

namespace {

  void involutionList(weyl::TwistedInvolutionList&, const kgb::KGB&);

  std::ostream& printBasePts(std::ostream&, const weyl::TwistedInvolutionList&,
			     const kgb::KGBEltList&, const kgb::KGB&);

  // the following class is copied integrally from kgb.cpp, helper namespace
class InvolutionCompare {
private:
  const weyl::WeylGroup* d_W;
public:
  explicit InvolutionCompare(const weyl::WeylGroup& W):d_W(&W) {}

  // one should have a < b iff
  // (a) involutionLength(a) < involutionLength(b) or
  // (b) involutionLengths are equal and length(a) < length (b) or
  // (c) both lengths are equal and a < b
  bool operator()
   (const weyl::TwistedInvolution& a, const weyl::TwistedInvolution& b) const
  {
    if      (d_W->involutionLength(a) != d_W->involutionLength(b))
      return d_W->involutionLength(a) < d_W->involutionLength(b) ;
    else if (d_W->length(a.w()) != d_W->length(b.w()))
      return d_W->length(a.w()) < d_W->length(b.w());
    else
      return a < b;
  }
};

}

/*****************************************************************************

        Chapter I -- Functions declared in kltest.h

  ... explain here when it is stable ...

******************************************************************************/

namespace kltest {

bool checkBasePoint(const kgb::KGB& kgb)

/*!
  \brief Checks whether the conjectural basepoint in each R-packet is
  independent of the choice of reduced expression.

  Explanation: let w be a twisted involution, written as an involution-reduced
  expression w = s_1 ... s_p (here s_j means either conjugation or commuting
  multiplication). Then the conjecture is that if x is fundamental, the element
  s_1 ... s_p.x, where the action is cross-action for a conjugation, cayley
  transform for a commuting multiplication, is independent of the choice of
  the reduced expression.

  Algorithm: denote x_w the conjectural element defined for w. Then the
  conjecture amounts to saying that for any descent s for w, s.x_{sw} is
  independent of the choice of s.
*/

{
  using namespace gradings;
  using namespace ioutils;
  using namespace kgb;
  using namespace weyl;

#ifdef VERBOSE
  std::cerr << "entering checkBasePoint ..." << std::endl;
#endif

  const WeylGroup& W = kgb.weylGroup();
  InvolutionCompare comp(W);
  TwistedInvolutionList wl;
  involutionList(wl,kgb);

  KGBEltList basepts;

  for (size_t x0 = 0; x0 < kgb.size() and kgb.length(x0) == 0; ++x0) {
    // check if x0 is large
    for (size_t s = 0; s < kgb.rank(); ++s) {
      Status::Value v = kgb.status(s,x0);
      if (v == Status::ImaginaryCompact)
	goto nextx0;
    }
    // if we get here, x0 is large
    // check that basepoint is well-defined
    basepts.assign(wl.size(),UndefKGB);
    basepts[0] = x0;
    for (size_t w_pos = 1; w_pos < wl.size(); ++w_pos) {
#ifdef VERBOSE
    std::cerr << w_pos << "\r";
#endif
      const TwistedInvolution& tw = wl[w_pos];
      WeylWord w_red;
      W.involutionOut(w_red,tw);
      for (size_t s = 0; s < kgb.rank(); ++s)
	if (W.hasDescent(s,tw.w())) {
	  TwistedInvolution sw = tw;
	  KGBElt sx_sw;
	  if (W.hasTwistedCommutation(s,tw)) {
	    W.leftMult(sw,s);
	    size_t sw_pos = std::lower_bound(wl.begin(),wl.end(),sw,comp) -
	      wl.begin();
	    sx_sw = kgb.cayley(s,basepts[sw_pos]);
	    if (sx_sw == UndefKGB)
	      return false;
	  } else {
	    W.twistedConjugate(sw,s);
	    size_t sw_pos = std::lower_bound(wl.begin(),wl.end(),sw,comp) -
	      wl.begin();
	    sx_sw = kgb.cross(s,basepts[sw_pos]);
	  }
	  if (basepts[w_pos] == UndefKGB) { // x_w is new
	    basepts[w_pos] = sx_sw;
	    // check if basepoint is large
	    // this is automatic according to David
	    // ... not done yet, need to think a little ...
	  } else {
	    if (sx_sw != basepts[w_pos])
	      return false;
	  }
	}
    }
#ifdef VERBOSE
    std::cerr << std::endl;
#endif
 nextx0:
    continue;
  }

#ifdef VERBOSE
  std::cerr << "done" << std::endl;
#endif

  OutputFile file;
  printBasePts(file,wl,basepts,kgb);

  return true;
}

void dualityPermutation(setutils::Permutation& a, const kl::KLContext& klc)

{}

}

/*****************************************************************************

        Chapter I -- Functions local to this module

  ... explain here when it is stable ...

******************************************************************************/

namespace {

std::ostream&
printBasePts(std::ostream& strm, const weyl::TwistedInvolutionList& wl,
	     const kgb::KGBEltList& bp, const kgb::KGB& kgb)

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
    printInvolution(strm,wl[j],kgb.weylGroup());  // added by jda: reduced form of the involution
    strm << "," << bp[j] << ")";
    strm << std::endl;
  }

  return strm;
}

void involutionList(weyl::TwistedInvolutionList& wl, const kgb::KGB& kgb)

/*
  Synopsis: put in wl the list of twisted involutions that appear in kgb.

  The result is equivalent to the internal d_tau in KGB.
*/

{
  using namespace weyl;

  wl.push_back(kgb.involution(0));

  for (size_t x = 0; x < kgb.size(); ++x)
    if (kgb.involution(x) != wl.back())
      wl.push_back(kgb.involution(x));

  return;
}


}

}
