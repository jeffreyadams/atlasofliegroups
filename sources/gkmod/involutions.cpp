/*!
\file
\brief Implementation of the class InvolutionSet.
*/
/*
  This is involutions.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "involutions.h"

#include <map>
#include <set>
#include <stack>

#include "complexredgp.h"
#include "latticetypes.h"
#include "rootdata.h"


/*
   In this file, the unadorned term involution always means twisted involution
*/

namespace atlas {

  // namespace {
  namespace involutions {
  namespace helper {

  using namespace involutions;

      /*!
\brief Derived class of InvolutionSet, to carry out the construction
of InvolutionSet.
      */
  class Helper:public InvolutionSet {

  private:

    /*! \brief Entry |d_toDualWeyl[s]| gives generator whose inner
representation in $W$ coincides with that of |s| in the dual Weyl group
    */
    weyl::WeylInterface d_toDualWeyl;

  public:

    // constructors and destructors
    explicit Helper(complexredgp::ComplexReductiveGroup&);

    virtual ~Helper() {}

    // manipulators
    void fill(const complexredgp::ComplexReductiveGroup&);

    void fillCartan(const complexredgp::ComplexReductiveGroup&);

    void fillInvolutions(const weyl::TwistedWeylGroup&);

    void fillDualInvolutions(const weyl::TwistedWeylGroup&);

    void weylCorrelation(const complexredgp::ComplexReductiveGroup&);
 };

  }

  }


/*****************************************************************************

        Chapter I -- The InvolutionSet class

******************************************************************************/

namespace involutions {

  using namespace atlas::involutions::helper;

InvolutionSet::InvolutionSet()

{}


/*!
\brief Constructor; computes the set of twisted involutions for G.

  NOTE: the set of involutions is actually constructed twice: once by
  fillCartan(), and then here. This is unfortunate, but negligible in the
  grand scheme of things, and having G remember its involution set runs
  against the principle of keeping the ComplexReductiveGroup class fairly
  small. Also this should be easy enough to change in the future if desired.
*/
InvolutionSet::InvolutionSet(complexredgp::ComplexReductiveGroup& G)
{
  Helper help(G);
  swap(help);
}

/******** copy, assignment and swap *****************************************/
void InvolutionSet::swap(InvolutionSet& other)

{
  std::swap(d_size,other.d_size);
  std::swap(d_rank,other.d_rank);
  d_action.swap(other.d_action);
  d_cartan.swap(other.d_cartan);
  d_involution.swap(other.d_involution);
  d_dualInvolution.swap(other.d_dualInvolution);

  return;
}

/******** assignment *******************************************************/

/*!
\brief Index of w in d_involution.

  Precondition: w is a twisted involution for W;

  Algorithm: find a reduced expression of w as an involution, and follow the
  action pointers.
*/
size_t InvolutionSet::involutionNbr(const weyl::TwistedInvolution& w,
				    const weyl::TwistedWeylGroup& W) const
{
  std::vector<signed char> ww=W.involution_expr(w);
  size_t x = 0; // initial, distinguished, involution

  for (size_t i = ww.size(); i-->0; )
    x = action(ww[i]>=0 ? ww[i] : ~ww[i],x);

  return x;
}

} // |namespace involutions|

/*****************************************************************************

        Chapter II -- The Helper class

******************************************************************************/

// namespace {
  namespace involutions {
  namespace helper {

Helper::Helper(complexredgp::ComplexReductiveGroup& G)
  : InvolutionSet()
  , d_toDualWeyl()
{
  d_size = G.numInvolutions();
  d_rank = G.semisimpleRank();

  weylCorrelation(G); // fills |d_toDualWeyl|
  fill(G); // completes filling the base object
}

/******** manipulators *******************************************************/


/*!
\brief Fills the tables.

  Precondition: size and rank have been set.
*/
void Helper::fill(const complexredgp::ComplexReductiveGroup& G)
{

  const weyl::TwistedWeylGroup& W = G.twistedWeylGroup();

  // fill the action and involution tables
  fillInvolutions(W);

  // fill in the dual involution table
  fillDualInvolutions(W);

  // fill in the cartan table
  fillCartan(G);

  return;
}


/*!
\brief Fills the Cartan table.

  Precondition: the action and involution tables have been filled;

  Explanation: it is important that we use the same numbering of Cartan
  subgroups as in G. The algorithm is to use the representative of Cartan
  \#j returned by G.twistedInvolution(j), locate that in the involution set,
  and then number its cross-orbit with j's.
*/
void Helper::fillCartan(const complexredgp::ComplexReductiveGroup& G)
{
  const weyl::TwistedWeylGroup& W = G.twistedWeylGroup();

  d_cartan.resize(d_size);

  std::set<size_t> found;
  std::stack<size_t> toDo;

  for (size_t cn = 0; cn < G.numCartanClasses(); ++cn) {

    const weyl::TwistedInvolution& w = G.twistedInvolution(cn);
    size_t x0 = involutionNbr(w,W);

    // find cross-orbit of x0
    found.clear();

    found.insert(x0);
    toDo.push(x0);

    while (!toDo.empty()) {

      size_t x = toDo.top();
      toDo.pop();

      for (weyl::Generator s = 0; s < d_rank; ++s) {
	const weyl::TwistedInvolution& w = involution(x);
	if (W.hasTwistedCommutation(s,w))
	  continue;
	size_t sx = action(s,x);
	if (found.insert(sx).second)
	  toDo.push(sx);
      }
    }

    // write result
    std::set<size_t>::iterator found_end = found.end();
    for (std::set<size_t>::iterator i = found.begin(); i != found_end; ++i)
      d_cartan[*i] = cn;

  }


}

/*!
\brief Fills the action and involution tables.

  Precondition: size and rank have been set; W is the relevant Weyl group.

  Explanation: action(s,w) is s.w if s and w twisted-commute, s.w.twist(s)
  otherwise.

  Algorithm: we fill the table in order of increasing involution length, by
  the naive algorithm of looking at all the slots which have not yet been
  filled, computing the result and looking it up in a set of elements for the
  next length.
*/
void Helper::fillInvolutions(const weyl::TwistedWeylGroup& W)
{
  d_action.resize(d_rank);
  for (size_t s = 0; s < d_action.size(); ++s)
    d_action[s].resize(d_size,UndefInvolution);

  d_involution.resize(d_size);

  std::map<weyl::TwistedInvolution,size_t> found;
  size_t nextLength = 1;
  size_t firstNew = 1;

  for (size_t x = 0; x < d_size; ++x) {

    if (x == nextLength) { // update
      nextLength += found.size();
      found.clear();
    }

    for (size_t s = 0; s < d_rank; ++s) {

      if (action(s,x) != UndefInvolution)
	continue;

      const weyl::TwistedInvolution& w = involution(x);
      weyl::TwistedInvolution sw = w;

      if (W.hasTwistedCommutation(s,w)) { // action is product
	W.leftMult(sw,s);
      } else { // action is twisted commutation
	W.twistedConjugate(sw,s);
      }

      if (found.insert(std::make_pair(sw,firstNew)).second) {
	// found a new element
	d_involution[firstNew] = sw;
	d_action[s][x] = firstNew;
	d_action[s][firstNew] = x;
	++firstNew;
      } else { // sw is already known
	size_t sx = found.find(sw)->second;
	d_action[s][x] = sx;
	d_action[s][sx] = x;
      }
    }
  }
}

/*!
\brief Fills the dual involution table.

  Precondition: the action and involution tables have been filled; d_toDualWeyl
  is set;
*/
void Helper::fillDualInvolutions(const weyl::TwistedWeylGroup& tW)
{
  const weyl::WeylGroup& W=tW.weylGroup();
  d_dualInvolution.resize(d_size);

  for (size_t i = 0; i < d_size; ++i)
  {
    weyl::WeylElt w = involution(i).w();
    tW.twist(w);
    weyl::WeylElt v = W.longest();
    W.mult(v,w);
    W.invert(v);
    d_dualInvolution[i] =
      weyl::TwistedInvolution(W.translation(v,d_toDualWeyl));
  }
}

/*!
\brief Fills in d_toDualWeyl.

  Explanation: |d_toDualWeyl[s]| is the outer representation of the generator
  that in the Weyl group has the same inner representation as |s| has in the
  dual Weyl group. Thus we can move elements from the Weyl group to the dual
  Weyl group by translating them via |d_toDualWeyl| and then interpreting the
  resulting inner representation in the dual Weyl group (translation operates
  on inner representations but uses a table defined in terms of outer
  representations).
*/
void Helper::weylCorrelation(const complexredgp::ComplexReductiveGroup& G)
{
  using namespace latticetypes;
  using namespace rootdata;
  using namespace weyl;

  const WeylGroup& W = G.weylGroup();

  // make dual Weyl group
  LatticeMatrix c;
  cartanMatrix(c,G.rootDatum());
  c.transpose();   // transposing the Cartan matrix may change Weyl group
  WeylGroup dW(c); // twist is irrelevant here

  // fill in d_toDualWeyl
  for (size_t s = 0; s < W.rank(); ++s) {
    weyl::WeylElt w=dW.generator(s); // converts |s| to inner numbering |dW|
    weyl::WeylWord ww=W.word(w); // interpret |w| in |dW|; gives singleton
    d_toDualWeyl[s] = ww[0];
  }
}

} // namespace helper

} // namespace involutions

} // namespace atlas
