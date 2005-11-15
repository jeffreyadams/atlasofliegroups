/*
  This is kl.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include "kl.h"

#include <cassert>
#include <map>
#include <stack>

#include "bitmap.h"
#include "blocks.h"

/*
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

namespace {

  using namespace kl;

  size_t firstAscent(const descents::DescentStatus&,
		     const descents::DescentStatus&, size_t);

  size_t goodAscent(const descents::DescentStatus&,
		    const descents::DescentStatus&, size_t);

  struct MuCompare {
    bool operator() (const MuData& lhs, const MuData& rhs) const {
      return lhs.first < rhs.first;
    }
  };

  class Thicket;

  class Helper:public KLContext {

    friend class Thicket;

  private:

  typedef std::set<KLPol>::iterator KLPtr;

  public:

    // constructors and destructors
    Helper() {}

    Helper(const KLContext&);

    virtual ~Helper() {}

    //accessors
    MuCoeff ascentMu(size_t, size_t, size_t) const;

    blocks::BlockEltPair cayley(size_t s, size_t y) const {
      return d_support->block().cayley(s,y);
    }

    size_t cross(size_t s, size_t y) const {
      return d_support->block().cross(s,y);
    }

    const descents::DescentStatus& descent(size_t y) const {
      return d_support->block().descent(y);
    }

    descents::DescentStatus::Value descentValue(size_t s, size_t y) const {
      return d_support->descentValue(s,y);
    }

    size_t dualOrbit(size_t y) const {
      return d_support->block().y(y);
    }

    size_t firstDirectRecursion(size_t) const;

    MuCoeff goodDescentMu(size_t, size_t, size_t) const;

    blocks::BlockEltPair inverseCayley(size_t s, size_t y) const {
      return d_support->block().inverseCayley(s,y);
    }

    MuCoeff lengthOneMu(size_t, size_t) const;

    size_t orbit(size_t y) const {
      return d_support->block().x(y);
    }

    MuCoeff recursiveMu(size_t, size_t) const;

    MuCoeff type2Mu(size_t, size_t) const;

    // manipulators
    void completePacket(size_t);

    void directRecursion(size_t, size_t);

    virtual void fill();

    void fillKLRow(size_t);

    void fillMuRow(size_t);

    void fillThickets(size_t);

    void muCorrection(std::vector<KLPol>&, size_t, size_t);

    void recursionRow(std::vector<KLPol>&, size_t, size_t);
  };

  class Thicket {

  public:

    struct Edge;
    typedef std::vector<Edge> EdgeList;

  private:

    std::vector<size_t> d_vertices;
    std::vector<EdgeList> d_edges;
    std::vector<size_t> d_xlist;
    Helper* d_helper;
     
  public:

    // constructors and destructors
    Thicket() {}

    Thicket(Helper&, size_t);

    ~Thicket() {}

    // accessors
    size_t cross(size_t s, size_t y) const {
      return d_helper->cross(s,y);
    }

    const descents::DescentStatus& descent(size_t y) const {
      return d_helper->descent(y);
    }

    descents::DescentStatus::Value descentValue(size_t s, size_t y) const {
      return d_helper->descentValue(s,y);
    }

    const EdgeList& edgeList(size_t j) const {
      return d_edges[j];
    }

    const klsupport::ExtremalRow& extremalRow(size_t j) const {
      return d_helper->extremalRow(d_vertices[j]);
    }

    size_t nonExtremal(size_t) const;

    size_t rank() const {
      return d_helper->rank();
    }

    size_t size() const {
      return d_vertices.size();
    }

    const std::vector<size_t>& vertices() const {
      return d_vertices;
    }

    // manipulators
    void fill();

    void fillXList();

    KLRow& klRow(size_t y) {
      return d_helper->d_kl[y];
    }

    std::set<KLPol>& store() {
      return d_helper->d_store;
    }

  };

  struct Thicket::Edge {
    size_t y;
    size_t s;
    std::vector<KLPol> recursion;
    // constructors and destructors
    Edge() {}
    Edge(size_t y, size_t s, const std::vector<KLPol>& r)
      :y(y), s(s), recursion(r) {}
  };

  class ThicketIterator {

  private:
    typedef Thicket::EdgeList EdgeList;

    std::stack<size_t> d_stack;
    bitmap::BitMap d_done;
    EdgeList::const_iterator d_current;
    EdgeList::const_iterator d_currentEnd;
    size_t d_pos;

    const Thicket* d_thicket;

  public:

    // constructors and destructors
    ThicketIterator(const Thicket&, size_t);

    // accessors
    bool operator() () const {
      return not d_stack.empty();
    }

    size_t operator* () const {
      return d_pos;
    }

    const Thicket::Edge& edge() const {
      return *d_current;
    }

    // manipulators
    ThicketIterator& operator++();
  };

}

/*****************************************************************************

        Chapter I -- The KLContext class.

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {

KLContext::KLContext(klsupport::KLSupport& kls)
  :d_support(&kls)

{
  d_store.insert(One);
}

/******** copy, assignment and swap ******************************************/
void KLContext::swap(KLContext& other)

{
  d_state.swap(other.d_state);

  std::swap(d_support,other.d_support);

  d_kl.swap(other.d_kl);
  d_mu.swap(other.d_mu);

  d_store.swap(other.d_store);

  return;
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/
void KLContext::fill()

/*
  Synopsis: fills the kl- and mu-lists.

  Explanation: this is the main function in this module; all the work is
  deferred to the Helper class.
*/

{
  if (d_state.test(KLFilled))
    return;

  Helper help(*this);

  help.fill();
  swap(help);

  d_state.set(KLFilled);

  return;
}

const KLPol& KLContext::klPol(size_t x, size_t y)

/*
  Synopsis: returns the Kazhdan-Lusztig-Vogan polynomial P_{x,y}

  Precondition: x and y are smaller than size();

  Explanation: this is basically a lookup function, but it still has to
  do some work. Basically, while x is not extremal w.r.t. y, it moves x up,
  using the "easy" induction relations. At that point, we look x (or rather,
  the finite set of extremal x'es that we have reduced to) up in the extremal 
  list for y. Those that are not found have a zero polynomial. From the others,
  we get the result.

  NOTE: when the result is an actual sum, we add it to d_store, so that the
  answer can still be given as a reference.

  NOTE: this is a lazy recursive implementation. It could suffer from
  efficiency problems, and has also the drawback of adding not only the
  final result, but also a number of intermediate polynomials to the store.
*/

{
  using namespace blocks;
  using namespace descents;
  using namespace klsupport;

  if (length(x) >= length(y))
    return x == y ? One : Zero;

  const Block& block = d_support->block();

  const DescentStatus& d_y = block.descent(y);
  const DescentStatus& d_x = block.descent(x);

  size_t s = goodAscent(d_x,d_y,rank());

  if (s != rank()) // found a good ascent
    switch (block.descentValue(s,x)) {
    case DescentStatus::ComplexAscent: {
      size_t x1 = block.cross(s,x);
      return klPol(x1,y);
    }
      break;
    case DescentStatus::RealNonparity: {
      return Zero;
    }
      break;
    case DescentStatus::ImaginaryTypeI: {
      size_t x1 = block.cayley(s,x).first;
      return klPol(x1,y);
    }
      break;
    default: // this cannot happen
      assert(false);
      break;
    }

  s = firstAscent(d_x,d_y,rank());

  if (s != rank()) { // have type II reduction
    BlockEltPair x1 = block.cayley(s,x);
    KLPol p = klPol(x1.first,y);
    p.safeAdd(klPol(x1.second,y));
    d_store.insert(p);
    return *d_store.find(p);
  }

  // if we reach this point, x is extremal w.r.t. y

  const ExtremalRow& e = extremalRow(y);
  if (not binary_search(e.begin(),e.end(),x)) // x is not found
    return Zero;
  // get position of x in e
  size_t j = lower_bound(e.begin(),e.end(),x) - e.begin();
  KLPtr p = d_kl[y][j];
  if (isZero(p))
    return Zero;
  else
    return *p;
}

MuCoeff KLContext::mu(size_t x, size_t y) const

/*
  Synopsis: returns mu(x,y).

  Explanation: it is guaranteed that all the x'es such that mu(x,y) != 0
  occur in d_mu[y] (and in fact, that only those occur.) So it is a simple
  matter of looking up x.
*/

{
  const MuRow& mr = d_mu[y];

  if (not std::binary_search(mr.begin(),mr.end(),std::make_pair(x,0),
			     MuCompare()))
    return 0;

  return std::lower_bound(mr.begin(),mr.end(),std::make_pair(x,0), 
			  MuCompare())->second;
}

}

/*****************************************************************************

        Chapter II -- The Helper class.

  ... explain here when it is stable ...

 *****************************************************************************/

namespace {

Helper::Helper(const KLContext& kl)
  :KLContext(kl)

{}

/******** accessors **********************************************************/
MuCoeff Helper::ascentMu(size_t x, size_t y, size_t s) const

/*
  Synopsis: computes mu(x,y) in a good ascent situation.

  Preconditions: l(y) > 0; l(x) = l(y)-1; s is an ascent for x w.r.t. y;

  Explanation: this is the situation where mu(x,y) is directly expressible.
  The point is that x will ascend to an element of the same length as y,
  so that the corresponding k-l polynomial is zero unless x ascends to y.
*/

{    
  using namespace blocks;
  using namespace descents;

  switch (descentValue(s,x)) {
  case DescentStatus::ComplexAscent: {
    size_t x1 = cross(s,x);
    return x1 == y ? 1 : 0;
  }
    break;
  case DescentStatus::RealNonparity: {
    return 0;
  }
    break;
  case DescentStatus::ImaginaryTypeI: {
    size_t x1 = cayley(s,x).first;
    return x1 == y ? 1 : 0;
  }
    break;
  case DescentStatus::ImaginaryTypeII: {
    BlockEltPair x1 = cayley(s,x);
    // mu(x,y) = mu(x1.first,y) + mu(x1.second,y)
    return (x1.first == y or x1.second == y) ? 1 : 0;
  }
    break;
  default: // this cannot happen
    assert(false);
    break;
  }

  return UndefMuCoeff;
}

size_t Helper::firstDirectRecursion(size_t y) const

/*
  Synopsis: returns the first descent generator that is not real type II

  Explanation: those are the ones that give a direct recursion formula for
  the k-l basis element.
*/

{
  using namespace descents;

  const descents::DescentStatus& d = descent(y);

  for (size_t s = 0; s < rank(); ++s) {
    DescentStatus::Value v = d[s];
    if (DescentStatus::isDirectRecursion(v))
      return s;
  }

  return rank();
}

MuCoeff Helper::goodDescentMu(size_t x, size_t y, size_t s) const

/*
  Synopsis: gets mu(x,y) by a good descent recursion.

  Precondition: length(y) > 0; length(x) = length(y)-1; all previous
  mu-rows have been filled in; s is a good descent for y;

  Explanation: a good descent for y is a descent that is not real type II,
  and also not imaginary compact; these are the cases where mu(x,y) is 
  directly expressed by recursion.
*/

{    
  using namespace blocks;
  using namespace descents;

  size_t y1;

  if (descentValue(s,y) == DescentStatus::ComplexDescent) {
    y1 = cross(s,y);
  } else { // s is real type I for y
    y1 = inverseCayley(s,y).first;
  }

  switch (descentValue(s,x)) {
  case DescentStatus::ImaginaryCompact: {
    return 0;
  }
  case DescentStatus::ComplexDescent: {
    size_t x1 = cross(s,x);
    return mu(x1,y1);
  }
  case DescentStatus::RealTypeI: {
    BlockEltPair x1 = inverseCayley(s,x);
    return mu(x1.first,y1)+mu(x1.second,y1);
  }
    break;
  case DescentStatus::RealTypeII: {
    size_t x1 = inverseCayley(s,x).first;
    return mu(x1,y1);
  }
    break;
  default: // this cannot happen
    assert(false);
    break;
  }

  return UndefMuCoeff;
}

MuCoeff Helper::lengthOneMu(size_t x, size_t y) const

/*
  Synopsis: computes mu(x,y) for l(y)-l(x) = 1.

  Preconditions: l(y) > 0; l(x) = l(y)-1; the mu-rows for y of smaller
  lengths have already been filled in.

  Explanation: these are the mu-values that can come up in possibly 
  non-extremal situations. They can be computed recursively by
  formulas using only other mu-values of the same kind.
*/

{
  // look if x is extremal w.r.t. y
  size_t s = firstAscent(descent(x),descent(y),rank());

  if (s != rank())
    return ascentMu(x,y,s);

  // if we get here, x is extremal w.r.t. y
  return recursiveMu(x,y);
}

MuCoeff Helper::recursiveMu(size_t x, size_t y) const

/*
  Synopsis: gets mu(x,y) by direct recursion.

  Precondition: length(y) > 0; length(x) = length(y)-1; all previous
  mu-rows have been filled in.

  Explanation: the mu-coefficients with length difference 1 are not
  guaranteed to be in the k-l table, so they have to be computed directly.
  This is the "hard case" in lengthOneMu().
*/

{  
  using namespace blocks;
  using namespace descents;

  size_t s = firstDirectRecursion(y);

  if (s != rank()) // the answer is another mu
    return goodDescentMu(x,y,s);

  return type2Mu(x,y);
}

MuCoeff Helper::type2Mu(size_t x, size_t y) const

/*
  Synopsis: gets mu(x,y) by type II recursion.

  Precondition: length(y) > 0; length(x) = length(y)-1; all previous
  mu-rows have been filled in; x is extremal w.r.t. y, and all descents
  of y are real type II;

  Explanation: in this situation, we have recursion formulas that will
  yield mu(x,y)+mu(x,s.y). It is known that proceeding in this way we must
  end up with an y' for which x is not extremal.

  Algorithm: we keep a stack of recursion formulas of the above kind, walking
  through the orbit of y under the simple descents, until we reach an y
  that (a) has a non-type-II descent, or (b) for which x is non-extremal.
*/

{
  using namespace descents;

  std::map<size_t,size_t> rel;
  std::stack<size_t> toDo;

  rel.insert(std::make_pair(y,rank()));
  toDo.push(y);

  size_t s;
  size_t y1 = y;
  MuCoeff mu = 0;

  while (not toDo.empty()) {
    y1 = toDo.top();
    toDo.pop();
    if (y1 < y) { // mu(x,y1) can be gotten recursively
      mu = Helper::mu(x,y1);
      goto unwind;
    }
    s = firstAscent(descent(x),descent(y1),rank());
    if (s != rank()) { // done
      mu = ascentMu(x,y1,s);
      goto unwind;
    }
    for (s = 0; s < rank(); ++s) {
      DescentStatus::Value v = descentValue(s,y1);
      if (not DescentStatus::isDescent(v))
	continue;
      if (DescentStatus::isDirectRecursion(v)) {
	mu = goodDescentMu(x,y1,s);
	goto unwind;
      }
      // at this point s is either real type II or imaginary compact
      if (v == DescentStatus::ImaginaryCompact)
	continue; 
      // at this point, there is no hope to resolve the situation at y1
      {
	size_t y2 = cross(s,y1);
	if (rel.insert(std::make_pair(y2,s)).second) // new element
	  toDo.push(y2);
      }
    }
  }

 unwind:
  while (y1 != y) {
    std::map<size_t,size_t>::iterator i = rel.find(y1);
    s = i->second;
    y1 = cross(s,y1);
    // goodDescent() is mu + mu(x,y2)
    mu = goodDescentMu(x,y1,s) - mu;
  }

  return mu;
}

/******** manipulators *******************************************************/
void Helper::completePacket(size_t y)

/*
  Synopsis: finishes the filling of the R-packet starting at y.

  Precondition: all the rows in the packet capable of a direct recursion
  are filled; at least one row is of this form;

  Explanation: what we do here, is completing the rows that can be completed
  form the already filled ones, using real type II recursions. Whatever is left
  will constitute "thickets", and will be dealt with by another function.

  Algorithm: we traverse the set of y1 in the packet that can be obtained from
  one of the filled ones through real type II cross-actions.
*/

{
  using namespace descents;
  using namespace klsupport;

  size_t o = orbit(y);
  std::stack<size_t> filled;
  std::set<size_t> empty;

  for (size_t y1 = y; orbit(y1) == o; ++y1) {
    if (d_kl[y1].size() != 0)
      filled.push(y1);
    else
      empty.insert(y1);
  }

  while (not filled.empty()) {

    size_t y1 = filled.top();
    filled.pop();

    for (size_t s = 0; s < rank(); ++s) {
      if (descentValue(s,y1) != DescentStatus::RealTypeII)
	continue;
      size_t y2 = cross(s,y1);
      std::set<size_t>::iterator i = empty.find(y2);
      if (i != empty.end()) { // found a new row
	// fill row y2
	std::vector<KLPol> klr;
	recursionRow(klr,y2,s);
	const ExtremalRow& e = d_support->extremalRow(y2);
	d_kl[y2].resize(e.size(),d_store.end());
	// klr[j] is P_{x,y2}+P_{x,y1}, for x = e[j]
	for (size_t j = 0; j < klr.size(); ++j) {
	  size_t x = e[j];
	  KLPol pol = klr[j];
	  pol.safeSubtract(klPol(x,y1));
	  if (not pol.isZero()) {
	    d_store.insert(pol);
	    d_kl[y2][j] = d_store.find(pol);
	  }
	}
	// update empty and filled
	empty.erase(i);
	filled.push(y2);
      }
    }
  }

  if (not empty.empty()) // we are not done
    fillThickets(y);

  return;
}

void Helper::directRecursion(size_t y, size_t s)

/*
  Synopsis: fills in the row for y using a direct recursion.

  Precondition: s is either a complex, or a real type I descent generator for 
  y.

  The real work is done by the recursionRow function, that will be used also
  in the real type II descents.
*/

{  
  std::vector<KLPol> klr;

  // put result of recursion formula in klr
  recursionRow(klr,y,s);

  // write result
  for (size_t j = 0; j < klr.size(); ++j)
    if (not klr[j].isZero()) {
      d_store.insert(klr[j]);
      d_kl[y][j] = d_store.find(klr[j]);
    }

  return;
}

void Helper::fill()

/*
  Synopsis: dispatches the work of filling the kl- and mu-lists.
*/

{
  // make sure the support is filled
  d_support->fill();

  // resize the lists
  d_kl.resize(d_support->size());
  d_mu.resize(d_support->size());

  // fill the lists
  size_t y = 0;

  // do the length zero cases; they come first in the enumeration
  for (; length(y) == 0; ++y) {
    // the k-l polynomial is 1
    KLPtr klp = d_store.find(One);
    d_kl[y].push_back(klp);
    // there are no mu-coefficients
  }

  // do the other cases
  for (; y < d_kl.size(); ++y) {
    fillKLRow(y);
    fillMuRow(y);
  }

  return;
}

void Helper::fillKLRow(size_t y)

/*
  Synopsis: fills in the row for y in the kl-table.

  Precondition: all lower rows have been filled; y is of length > 0;
  R-packets are consecutively numbered;

  Explanation: this function actually fills out the whole R-packet to
  which y belongs (unless it returns immediately). This is done in two
  stages: first, we fill the rows for which there is a direct recursion
  relation. The second stage, if any, is forwarded to the completePacket
  function.
*/

{
  using namespace klsupport;

  if (d_kl[y].size()) // row has already been filled
    return;

  size_t o = orbit(y);
  bool done = true;

  // fill in the direct recursions in the R-packet
  for (size_t y1 = y; orbit(y1) == o; ++y1) {
    size_t s = firstDirectRecursion(y1);
    if (s != rank()) { // direct recursion
      directRecursion(y1,s);
    } else
      done = false;
  }

  // deal with the others
  if (not done)
    completePacket(y);

  return;
}

void Helper::fillMuRow(size_t y)

/*
  Synopsis: fills in the row for y in the mu-table.

  Precondition: the row for y in the kl-table has been filled; length(y) > 0;

  Explanation: for the elements of length < lenght(y) - 1, mu(x,y) can
  be non-zero only if x is extremal w.r.t. y; so we run through d_kl[y],
  and loook at the cases where the polynomial is of degree (1/2)(l(y)-l(x)-1)
  (the biggest possible). For the elements of colength 1, in the classical
  case we always had mu(x,y)=1. Here that's not true anymore, zero might
  be a possibility, and also larger values I believe; but at any rate the
  mu-values can be computed in terms of other mu-values.

  NOTE: we are not using the hasse-list here, although it is probably a
  good idea to compute that; that will reduce the mu-computation.
*/

{
  using namespace klsupport;

  const ExtremalRow& e = extremalRow(y);

  size_t ly = length(y);

  // do cases of length < ly-1
  for (size_t j = 0; j < e.size(); ++j) {
    size_t x = e[j];
    size_t lx = length(x);
    if (lx >= ly-1)
      break;
    // check that lx and ly have opposite parity
    if ((lx&1ul) == (ly&1ul))
      continue;
    size_t d = (ly-lx-1)/2;
    KLPtr klp = d_kl[y][j];
    if (isZero(klp))
      continue;
    if (klp->degree() < d)
      continue;
    // if we get here, we found a mu-coefficient for x
    d_mu[y].push_back(std::make_pair(x,(*klp)[d]));
  }

  // do cases of length ly-1
  size_t x_begin = d_support->lengthLess(ly-1);
  size_t x_end = d_support->lengthLess(ly);

  for (size_t x = x_begin; x < x_end; ++x) {
    MuCoeff mu = lengthOneMu(x,y);
    if (mu)
      d_mu[y].push_back(std::make_pair(x,mu));
  }

  return;
}

void Helper::fillThickets(size_t y)

/*
  Synopsis: finishes the filling of the R-packet starting at y.

  Precondition: all the rows in the packet that are not part of a "thicket"
  have already been filled.

  Explanation: a thicket is an element that has only real type II descents,
  and that moreover is such that all elements in its conected component for
  the relation of being connected by a sequence of real type II cross-actions,
  is of the same form (for instance, it happens quite often that the principal
  series representations are a thicket.) Then we must use the "structural fact"
  in Lemma 6.2 of David's Park City notes: for each x lower than y, there is
  an element y' of the thicket for which x has an ascent (or if there is no 
  such y', the k-ll polynomial is zero.)

  Algorithm: (a) for each y' in the R-packet that has not already been filled,
  we determine the thicket of y' via a traversal algorithm, as usual (b) we
  fill in all the P_{x,y} in that thicket by a downwards recursion. In order
  to do that, we put the union of all the x'es in the extremal lists in a 
  common ordered list. Then we move down the list; for each x, by assumption
  there is an y' in the thicket for which it is not extremal, and starting
  from there we can fill in P_{x,y} in all the rows where it needs to be
  filled.
  
  NOTE: thickets are always small, at the very worst a few hundred elements,
  so we don't have to worry about efficiency here.

  NOTE: perhaps it might be expected that there is usually (always?) a "large"
  element in the thicket, that would provide ascents for all x'es. So it might
  be worthwhile to sort the elements in the thicket in order of number of
  descents.
*/

{
  using namespace descents;
  using namespace klsupport;

  size_t o = orbit(y);

  for (size_t y1 = y; orbit(y1) == o; ++y1)
    if (d_kl[y1].size() == 0) {
      Thicket thicket(*this,y1);
      thicket.fill();
    }

  return;
}

void Helper::muCorrection(std::vector<KLPol>& klp, size_t y, size_t s)

/*
  Synopsis: subtracts from klp the correcting terms in the k-l recursion.

  Precondtion: klp contains the terms corresponding to c_s.c_y, for the x that
  are extremal w.r.t. y; the mu-table and kl-table has been filled in for 
  elements of length <= y.

  Explanation: the recursion formula is of the form:
  
    lhs = c_s.c_{y1} - \sum_{z} mu(z,y1)c_z

  where z runs over the elements < y s.t. s is a descent for z, y1 is s.y,
  and lhs is c_y when s is a complex descent or real type I for y, and 
  c_{y}+c_{s.y} when s is real type II.
*/

{
  using namespace descents;
  using namespace klsupport;
  using namespace polynomials;

  size_t y1;

  if (descentValue(s,y) == DescentStatus::ComplexDescent) {
    y1 = cross(s,y);
  } else { // s is real for y
    y1 = inverseCayley(s,y).first;
  }

  const MuRow& mrow = d_mu[y1];
  const ExtremalRow& e = extremalRow(y);
  size_t l_y = length(y);

  for (size_t i = 0; i < mrow.size(); ++i) {

    size_t z = mrow[i].first;
    size_t l_z = length(z);

    DescentStatus::Value v = descentValue(s,z);
    if (not DescentStatus::isDescent(v))
      continue;

    MuCoeff mu = mrow[i].second;

    for (size_t j = 0; j < e.size(); ++j) {
      size_t x = e[j];
      if (length(x) > l_z)
	break;
      // subtract x^d.mu.P_{x,z} from klp[j], where d = 1/2(l(y)-l(z))
      const KLPol& pol = klPol(x,z);
      Degree d = (l_y-l_z)/2;
      klp[j].safeSubtract(pol,d,mu);
    }

  }

  return;
}

void Helper::recursionRow(std::vector<KLPol>& klr, size_t y, size_t s)

/*
  Synopsis: puts in klr the right-hand side of the recursion formula for y
  corresponding to the descent s.

  Precondition: s is either a complex, or a real type I or type II descent 
  generator for y.

  Explanation: the shape of the formula is: 

    P_{x,y} = (c_s.c_{y1})-part - correction term

  where y1 = cross(s,y) when s is complex for y, one of the two elements in
  inverseCayley(s,y) when s is real. The (c_s.c_{y1})-part depends on the 
  status of x w.r.t. s (we look only at extremal x, so we know it is a 
  descent); the correction term, coming from \sum_z \mu(z,y1)c_z, depends only 
  on y1.
*/

{  
  using namespace blocks;
  using namespace descents;
  using namespace klsupport;

  size_t y1;

  if (descentValue(s,y) == DescentStatus::ComplexDescent) {
    y1 = cross(s,y);
  } else { // s is real type I or type II for y
    y1 = inverseCayley(s,y).first;
  }

  const ExtremalRow& e = d_support->extremalRow(y);
  d_kl[y].resize(e.size(),d_store.end());

  klr.resize(e.size());

  for (size_t j = 0; j < klr.size()-1; ++j) {
    size_t x = e[j];
    switch (descentValue(s,x)) {
    case DescentStatus::ImaginaryCompact: {
      // (q+1)P_{x,y1}
      klr[j] = klPol(x,y1);
      klr[j].safeAdd(klr[j],1);
    }
      break;
    case DescentStatus::ComplexDescent: {
      size_t x1 = cross(s,x);
      // P_{x1,y1}+q.P_{x,y1}
      klr[j] = klPol(x1,y1);
      klr[j].safeAdd(klPol(x,y1),1);
    }
      break;
    case DescentStatus::RealTypeI: {
      BlockEltPair x1 = inverseCayley(s,x);
      // P_{x1.first,y1}+P_{x1.second,y1}+(q-1)P_{x,y1}
      klr[j] = klPol(x1.first,y1);
      klr[j].safeAdd(klPol(x1.second,y1));
      klr[j].safeAdd(klPol(x,y1),1);
      klr[j].safeSubtract(klPol(x,y1));
    }
      break;
    case DescentStatus::RealTypeII: {
      size_t x1 = inverseCayley(s,x).first;
      // P_{x_1,y_1}+qP_{x,y1}-P_{s.x,y1}
      klr[j] = klPol(x1,y1);
      klr[j].safeAdd(klPol(x,y1),1);
      klr[j].safeSubtract(klPol(cross(s,x),y1));
    }
      break;
    default: // this cannot happen
      assert(false);
      break;
    }
  }

  // last k-l polynomial is 1
  klr.back() = One;

  // do mu-correction
  muCorrection(klr,y,s);

  return;
}

}

/*****************************************************************************

        Chapter III -- The Thicket class

  ... explain here when it is stable ...

 *****************************************************************************/

namespace {

Thicket::Thicket(Helper& h, size_t y)
  :d_helper(&h)

/*
  Synopsis: constructs a thicket from y.

  Explanation: the thicket of y is the connected component of y in the graph
  whose edges are the real type II cross-actions. It is therefore contained
  in the R-packet of y.
*/

{
  using namespace descents;

  typedef std::map<size_t,EdgeList> M;
  typedef M::iterator MI;

  std::stack<size_t> toDo;
  M thicket;

  toDo.push(y);
  thicket.insert(std::make_pair(y,EdgeList()));

  while (not toDo.empty()) {
    size_t y1 = toDo.top();
    toDo.pop();
    EdgeList& e1 = thicket.find(y1)->second;
    for (size_t s = 0; s < rank(); ++s) {
      if (descentValue(s,y1) != DescentStatus::RealTypeII)
	continue;
      size_t y2 = cross(s,y1);
      std::pair<MI,bool> ins = thicket.insert(std::make_pair(y2,EdgeList()));
      if (ins.second) { // y2 was new
	EdgeList& e2 = ins.first->second;
	e1.push_back(Edge(y2,s,std::vector<KLPol>()));
	e2.push_back(Edge(y1,s,std::vector<KLPol>()));
	d_helper->recursionRow(e1.back().recursion,y2,s);
	d_helper->recursionRow(e2.back().recursion,y1,s);
	toDo.push(y2);
      }
    }
  }

  // write out result
  d_vertices.resize(thicket.size());
  d_edges.resize(thicket.size());
  size_t j = 0;

  for (MI i = thicket.begin(); i != thicket.end(); ++i) {
    d_vertices[j] = i->first;
    d_edges[j].swap(i->second);
    ++j;
  }

  // revert to relative addresses in EdgeLists
  for (size_t j = 0; j < d_edges.size(); ++j) {
    EdgeList& el = d_edges[j];
    for (size_t i = 0; i < el.size(); ++i)
      // replace el[i].y by its position in d_vertices
      el[i].y = std::lower_bound(d_vertices.begin(),d_vertices.end(),el[i].y)
	- d_vertices.begin();
  }
}

/******** accessors **********************************************************/
size_t Thicket::nonExtremal(size_t x) const

/*
  Synopsis: returns the position in d_vertices of the first y that is not 
  extremal w.r.t. x, size() if there is no such y.
*/

{
  for (size_t j = 0; j < d_vertices.size(); ++j) {
    size_t y = d_vertices[j];
    size_t s = firstAscent(descent(x),descent(y),rank());
    if (s != rank()) // y was found
      return j;
  }

  return size();
}

/******** manipulators *******************************************************/

void Thicket::fill()

/*
  Synopsis: fills in the k-l polynomials for the elements of the thicket.
*/

{
  using namespace klsupport;

  fillXList();

  // initialize rows, fill in last element
  for (size_t j = 0; j < size(); ++j) {
    size_t y = d_vertices[j];
    const ExtremalRow& er = extremalRow(j);
    klRow(y).resize(er.size(),store().end());
    klRow(y).back() = store().find(One);
  }

  // fill in rows
  for (size_t j = d_xlist.size(); j;) {
    --j;
    size_t x = d_xlist[j];
    // find a vertex s.t. x is not extremal
    size_t iy = nonExtremal(x);
    if (iy == size()) // do nothing; polynomials are already zero
      continue;
    for (ThicketIterator i(*this,iy); i(); ++i) {
      size_t pos = *i;
      size_t y = d_vertices[pos];
      const ExtremalRow& er = extremalRow(pos);
      if (not std::binary_search(er.begin(),er.end(),x))
	continue;
      // if we get here, we need to compute the polynomial
      size_t xpos = std::lower_bound(er.begin(),er.end(),x) - er.begin();
      const Edge& e = i.edge();
      KLPol pol = e.recursion[xpos];
      size_t sy = cross(e.s,y);
      pol.safeSubtract(d_helper->klPol(x,sy));
      if (not pol.isZero()) {
	store().insert(pol);
	klRow(y)[xpos] = store().find(pol);
      }
    }
  }

  return;
}

void Thicket::fillXList()

/*
  Synopsis: puts in d_xlist the union of the extremal lists for the elements
  of the thicket, top element excluded.
*/

{
  using namespace klsupport;

  std::set<size_t> xs;

  for (size_t j = 0; j < size(); ++j) {
    const ExtremalRow& e = extremalRow(j);
    for (size_t i = 0; i < e.size()-1; ++i)
      xs.insert(e[i]);
  }

  d_xlist.reserve(xs.size());
  copy(xs.begin(),xs.end(),std::back_inserter(d_xlist));

  return;
}

}

/*****************************************************************************

        Chapter IV -- The ThicketIterator class

  ... explain here when it is stable ...

 *****************************************************************************/

namespace {

ThicketIterator::ThicketIterator(const Thicket& th, size_t j)
  : d_thicket(&th)

{
  d_stack.push(j);
  d_done.resize(th.size());
  d_done.insert(j);
  d_current = th.edgeList(j).begin();
  d_currentEnd = th.edgeList(j).end();
  d_pos = j;
}

/******** manipulators *******************************************************/
ThicketIterator& ThicketIterator::operator++ ()

/*
  Synopsis: pre-increment operator.

  Incrementing the iterator means going to the next new element in the
  current edgelist; if there is no such, popping the stack and finding new
  elements in the new edgelist; if the stack is empty, return.
*/

{
  for (; d_current != d_currentEnd; ++d_current) {
    size_t j = d_current->y;
    if (not d_done.isMember(j)) { // j is new
      d_pos = j;
      d_stack.push(j);
      d_done.insert(j);
      return *this;
    }
  }

  // if we get here, we have exhausted the current row

  while (not d_stack.empty()) {
    size_t j = d_stack.top();
    d_stack.pop();
    d_current = d_thicket->edgeList(j).begin();
    d_currentEnd = d_thicket->edgeList(j).end();
    for (; d_current != d_currentEnd; ++d_current) {
      size_t i = d_current->y;
      if (not d_done.isMember(i)) { // i is new
	d_pos = i;
	d_stack.push(i);
	return *this;
      }
    }
  }

  // if we get here, we have exhausted the thicket
  d_pos = d_thicket->size();

  return *this;
}

}

/*****************************************************************************

        Chapter V -- Functions local to this module

  ... explain here when it is stable ...

 *****************************************************************************/

namespace {

size_t firstAscent(const descents::DescentStatus& d1,
		   const descents::DescentStatus& d2, size_t rank)

/*
  Synopsis: returns the first s that is an ascent for d1, and a descent for d2;
  rank if there is none such.

  Precondition: rank is the number of valid fields in d1 and d2.
*/

{

  using namespace descents;

  for (size_t s = 0; s < rank; ++s) {
    if ((DescentStatus::isDescent(d2[s])) and 
	(not DescentStatus::isDescent(d1[s])))
      return s;
  }

  return rank;
}

size_t goodAscent(const descents::DescentStatus& d1,
		  const descents::DescentStatus& d2, size_t rank)

/*
  Synopsis: returns the first s that is a non-ImaginaryTypeII ascent for
  d1, and a descent for d2; rank if there is none such.

  Precondition: rank is the number of valid fields in d1 and d2.
*/

{
  using namespace descents;

  for (size_t s = 0; s < rank; ++s) {
    if ((not DescentStatus::isDescent(d2[s])) or 
	(DescentStatus::isDescent(d1[s])))
      continue;
    if (d1[s] != DescentStatus::ImaginaryTypeII)
      return s;
  }

  return rank;
}

}

}
