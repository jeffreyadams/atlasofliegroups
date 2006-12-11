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

  See file main.cpp for full copyright notice
*/

#include "kl.h"

#ifdef VERBOSE
#include <iostream>
#include <ctime>
#endif

#include <cassert>
#include <map>
#include <stack>
#include <stdexcept>

#include "bitmap.h"
#include "blocks.h"
#include "error.h"
#include "kl_error.h"

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

  std::ostream& operator<<
    (std::ostream& out,  const kl::KLPol& p)
  { return prettyprint::printPol(out,p,"X"); }

  namespace kl {

    void alert(KLIndex i) // illegal index into hash table
    {
#ifdef VERBOSE
      std::cerr << "Illegal KL pointer: " << i <<".\n";
#endif
      assert(false); // bomb out here
    }

    /* methods of KLPolEntry */

KLPolEntry::KLPolEntry(KLPolRef p) // extract polynomial from PolRef
      : KLPol(p.freeze()) {}

  /*!
    \brief calculate a hash value in [0,modulus[, where modulus is a power of 2

    The function is in fact evaluation of the polynomial (with coefficients
    interpreted in Z) at the point 2^21+2^13+2^8+2^5+1=2105633, which can be
    calculated quickly (without multiplications) and which gives a good spread
    (which is not the case if 8481 is replaced by a small number, because the
    evaluation values will not groz fast enough for low degree polynomials!).

  */
size_t KLPolEntry::hashCode(KLIndex modulus) const
  { const polynomials::Polynomial<KLCoeff>& P=*this;
    if (P.isZero()) return 0;
    size_t i=P.degree();
    KLIndex h=P[i];
    while (i-->0) h= ((h<<21)+(h<<13)+(h<<8)+(h<<5)+h+P[i]) & (modulus-1);
    return h;
  }

bool KLPolEntry::operator!=(KLPolRef e) const
 {
   if (degree()!=e.degree()) return true;
   if (isZero()) return false; // since degrees match
   for (polynomials::Degree i=0; i<=degree(); ++i)
     if ((*this)[i]!=e[i]) return true;
   return false; // no difference found
 }

    /* methods of KLPool */

void KLPool::push_back(const KLPol& p)
  {
    index.push_back(pool.size());
    pool.push_back(KLCoeff::raw_val(p.degree()+1));
    if (not p.isZero()) // allows writing i<=p.degree() below
      for (size_t i=0; i<=p.degree(); ++i) pool.push_back(p[i]);
  }

KLPolRef KLPool::operator[] (KLIndex i) const // get polynomial reference
  {
    size_t ii=index[i]; // index into pool
    size_t l=static_cast<unsigned char>(pool[ii++]);

    return KLPolRef(&pool[ii],0,l);
  }




  namespace helper {

  void pause() {} // possible breakpoint for debugger

  using namespace kl;

  size_t firstAscent(const descents::DescentStatus&,
		     const descents::DescentStatus&, size_t);

  size_t goodAscent(const descents::DescentStatus&,
		    const descents::DescentStatus&, size_t);

  struct MuCompare {
    bool operator() (const MuData& lhs, const MuData& rhs) const {
      return lhs.first < rhs.first;
    }
  }; // struct MuCompare

  class Thicket;

  class Helper:public KLContext {

    friend class Thicket;

  private:

  public:

    // constructors and destructors
    Helper() {}

    Helper(const KLContext&);

    virtual ~Helper() {}

    //accessors
    MuCoeff ascentMu(size_t x, size_t y, size_t s) const;

    /*!
 \brief Cayley transform of block element y through simple root s.
    */
    blocks::BlockEltPair cayley(size_t s, size_t y) const {
      return d_support->block().cayley(s,y);
    }

    /*!
 \brief Cross action of simple root s on block element y.
    */
    size_t cross(size_t s, size_t y) const {
      return d_support->block().cross(s,y);
    }

    /*!
 \brief Returns vector of RANK_MAX unsigned char; entry s gives
 descent status of simple root s for block element y.
    */
    const descents::DescentStatus& descent(size_t y) const {
      return d_support->block().descent(y);
    }

    /*!
\brief Unsigned char whose value gives the descent status of
simple root s for block element y.
    */
    descents::DescentStatus::Value descentValue(size_t s, size_t y) const {
      return d_support->descentValue(s,y);
    }

    /*!
\brief Second coordinate (corresponding to K^vee orbit on G^vee/B^vee)
of pair of integers specifying block element y.
    */
    size_t dualOrbit(size_t y) const {
      return d_support->block().y(y);
    }

    size_t firstDirectRecursion(size_t y) const;

    MuCoeff goodDescentMu(size_t x, size_t y, size_t s) const;

    blocks::BlockEltPair inverseCayley(size_t s, size_t y) const {
      return d_support->block().inverseCayley(s,y);
    }

    /* Declare that member function klPol from base class KLContext is also
       considered. This is necessary since it would otherwise be shadowed
       rather than overloaded by a new definition with different signature
       in the Helper class.
    */
    using KLContext::klPol;

    KLPolRef klPol(size_t x, size_t y, KLRow::const_iterator klv,
		   klsupport::PrimitiveRow::const_iterator p_begin,
		   klsupport::PrimitiveRow::const_iterator p_end) const;

    MuCoeff lengthOneMu(size_t x, size_t y) const;

    void makeExtremalRow(klsupport::PrimitiveRow& e, size_t y) const;

    void makePrimitiveRow(klsupport::PrimitiveRow& e, size_t y) const;

    /*!
\brief First coordinate (corresponding to K orbit on G/B) of pair of
integers specifying block element y.
    */
    size_t orbit(size_t y) const {
      return d_support->block().x(y);
    }

    MuCoeff recursiveMu(size_t x, size_t y) const;

    MuCoeff type2Mu(size_t x, size_t y) const;

    // manipulators
    void completePacket(size_t y);

    void directRecursion(size_t y, size_t s);

    virtual void fill();

    void fillKLRow(size_t y);

    void fillMuRow(size_t y);

    void fillThickets(size_t y);

    void muCorrection(std::vector<KLPol>& klv,
		      const klsupport::PrimitiveRow& e,
		      size_t y, size_t s);

    void recursionRow(std::vector<KLPol> & klv,
		      const klsupport::PrimitiveRow& e, size_t y, size_t s);

    void writeRow(const std::vector<KLPol>& klv,
		  const klsupport::PrimitiveRow& e, size_t s);
  }; // class Helper

    /*!
\brief Collection of block elements y_j of the same length, differing
by type II real cross actions, and having no other descents.

A pair (y_j, s x y_j) (with s type II real) is an Edge of the
thicket. Such an s is a descent for both y_j and s x y_j.  The pair is
the image of the Cayley transform (type II imaginary) of a single
element y_0, of length one less.  The first class of KL recursion
relations used in a Thicket is this:

P_{x,y_j} + P_{x,s x y_j} = [formula involving various P_{? ,
y_0}].

The right sides here are known by induction on y, and recorded in the
recursion data member of the Edge.  We can therefore compute P_{x,y}
as soon as we know P_{x,y'} for a single element y' in the
thicket. For that (given x) we find essentially three possibilities:

1) x is equal to y', so P_{x,y'} = 1.

2) There is an s that is a good ascent for x (not in tau(x), and
leading to an x' of length one more than the length of x), but is type
II real for y'.  In this case there is a recursion formula something
like

P_{x,y'} = P_{x' , y'}.

(If s is type II imaginary for x, then there are two terms on the
right.)  In any case the right side is known by downward induction on
x, so finally we know all P_{x,y_j} for y_j in the Thicket.

3) P_{x,y'} = 0.

This computation is carried out in the member function fill().
    */
  class Thicket {

  public:

    struct Edge;
    typedef std::vector<Edge> EdgeList;

  private:

    typedef klsupport::PrimitiveRow::iterator PI;
    typedef KLRow::iterator KLI;

    /*!
\brief List of elements y_j in Thicket.
    */
    std::vector<size_t> d_vertices;

    /*!
\brief Entry j lists the edges ending at y_j.
    */
    std::vector<EdgeList> d_edges;

    /*!
\brief List of all x that are extremal with respect to some y in Thicket.
    */
    std::vector<size_t> d_xlist;

    /*!
\brief Entry j lists all x of length strictly less than l(y) that are
extremal with respect to y_j.

Extremal means that each descent for y_j is a descent for x.
    */
    std::vector<klsupport::PrimitiveRow> d_extr;

    /*!
\brief Entry j lists all x of length strictly less than l(y) that are
primitive with respect to y_j.

Primitive means that each descent for y_j is either a descent for x or
type II imaginary for x.
    */
    std::vector<klsupport::PrimitiveRow> d_prim;

    /*!
\brief
    */
    std::vector<KLRow> d_klr;

    /*!
\brief
    */
    std::vector<PI> d_firstPrim;

    /*!
\brief
    */
    std::vector<KLI> d_firstKL;
    Helper* d_helper;

  public:

    // constructors and destructors
    Thicket() {}

    Thicket(Helper&, size_t);

    ~Thicket() {}

    // accessors

    /*!
\brief Cross action of simple reflection s on block element y.
    */
    size_t cross(size_t s, size_t y) const {
      return d_helper->cross(s,y);
    }

    /*!
\brief List of RankMax unsigned chars; number s gives the descent
status 0-7 of simple root s for y.
    */
    const descents::DescentStatus& descent(size_t y) const {
      return d_helper->descent(y);
    }


    /*!
\brief Unsigned char between 0 and 7; gives the descent of
simple root s for y.
    */
    descents::DescentStatus::Value descentValue(size_t s, size_t y) const {
      return d_helper->descentValue(s,y);
    }

    /*!
\brief List of Edge's ending at y_j.
    */
    const EdgeList& edgeList(size_t j) const {
      return d_edges[j];
    }

    /*!
\brief KL polynomial P_{x,y_pos}, for any y_pos in Thicket.
    */
    KLPolRef klPol(size_t x, size_t pos) const {
      return d_helper->klPol(x,d_vertices[pos],d_firstKL[pos],d_firstPrim[pos],
			     d_prim[pos].end());
    }

    /*!
\brief List of all x of length strictly less than l(y) that are
primitive with respect to y_j.

Primitive means that each descent for y_j is either a descent for x or
type II imaginary for x.
    */
    const klsupport::PrimitiveRow& primitiveRow(size_t j) const {
      return d_prim[j];
    }

    size_t nonExtremal(size_t x) const;

    /*!
\brief Pointer to the KL polynomial one.
    */
    KLPtr one() const {
      return d_helper->d_one;
    }

    /*!
\brief Semisimple rank.
    */
    size_t rank() const {
      return d_helper->rank();
    }

    /*!
\brief Number of vertices in Thicket.
    */
    size_t size() const {
      return d_vertices.size();
    }

    /*!
\brief List of vertices in Thicket.
    */
    const std::vector<size_t>& vertices() const {
      return d_vertices;
    }

    /*!
\brief Pointer to the KL polynomial zero.
    */
    KLPtr zero() const {
      return d_helper->d_zero;
    }

    // manipulators
    bool ascentCompute(size_t x, size_t pos);

    void edgeCompute(size_t x, size_t pos, const Edge& e);

    void fill();

    void fillXList();

    KLRow& klRow(size_t y) {
      return d_helper->d_kl[y];
    }

    KLStore& store() {
      return d_helper->d_store;
    }

    }; // class Thicket

    /*!
\brief Pair (y , s x y) (with s type II real) in a Thicket.

    */
  struct Thicket::Edge {

    /*!
\brief Number of the block element at which the edge begins.
    */
    size_t source;

    /*!
\brief Number of the block element at which the edge ends.
    */
    size_t y;

    /*!
\brief Number of a simple root s, assumed to be type II real for y.
    */
    size_t s;

    /*!
\brief List of formulas for P_{x_i,y} + P_{x_i,source}, for i in a
list of elements primitive with respect to some y' in the Thicket.
    */
    std::vector<KLPol> recursion;
    // constructors and destructors
    Edge() {}
    Edge(size_t x, size_t y, size_t s, const std::vector<KLPol>& r)
      :source(x), y(y), s(s), recursion(r) {}
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

  }; // class ThicketIterator

  }
  }

/*****************************************************************************

        Chapter I -- The KLContext class.

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {
  using namespace atlas::kl::helper;

KLContext::KLContext(klsupport::KLSupport& kls)
  :d_support(&kls)

{
  d_zero = d_store.match(Zero);
  d_one = d_store.match(One);

  if (d_zero==d_one) std::cout << "Warning, One=Zero !\n";
}

/******** copy, assignment and swap ******************************************/


/*!
  \brief Copy constructor.

  Since we use indices instead of iterators, noting gets invalidated any more
*/
KLContext::KLContext(const KLContext& other)
  :d_state(other.d_state),
   d_support(other.d_support),
   d_prim(other.d_prim),
   d_kl(other.d_kl),
   d_mu(other.d_mu),
   d_store(other.d_store),
   d_zero(other.d_zero),
   d_one(other.d_one)
{}

KLContext& KLContext::operator= (const KLContext& other)

/*!
  \brief Assignment operator.

  Use copy constructor. This requires a check for self-assignment, or the
  source would be destroyed!
*/

{
  // handle self-assignment
  if (&other != this) {
    this->~KLContext();
    new(this) KLContext(other);
  }

  return *this;
}

void KLContext::swap(KLContext& other)

{
  d_state.swap(other.d_state);

  std::swap(d_support,other.d_support);

  d_prim.swap(other.d_prim);

  d_kl.swap(other.d_kl);
  d_mu.swap(other.d_mu);

  d_store.swap(other.d_store);
  std::swap(d_zero,other.d_zero);
  std::swap(d_one,other.d_one);

  return;
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/
void KLContext::fill()

/*!
  \brief Fills the kl- and mu-lists.

  Explanation: this is the main function in this module; all the work is
  deferred to the Helper class.
*/

{
#ifdef VERBOSE
  std::cerr << "computing kazhdan-lusztig polynomials ..." << std::endl;
#endif

  if (d_state.test(KLFilled))
    return;

  Helper help(*this);

  help.fill();
  swap(help);

  d_state.set(KLFilled);

#ifdef VERBOSE
  std::cerr << "done" << std::endl;
#endif

  return;
}

KLPolRef KLContext::klPol(size_t x, size_t y) const

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

{
  using namespace klsupport;

  typedef PrimitiveRow::const_iterator EI;

  const PrimitiveRow& pr = d_prim[y];
  const KLRow& klr = d_kl[y];

  if (d_support->primitivize(x,descentSet(y))) {
    EI xptr = std::lower_bound(pr.begin(),pr.end(),x);
    if (xptr != pr.end() and *xptr == x) { // result is nonzero
      size_t xpos = xptr - pr.begin();
      return d_store[klr[xpos]];
    } else {
      return d_store[d_zero];
    }
  }

  // if we get here, the primitivization hit a real compact ascent
  return d_store[d_zero];
}

MuCoeff KLContext::mu(size_t x, size_t y) const

/*!
  \brief Returns mu(x,y).

  Explanation: it is guaranteed that all the x'es such that mu(x,y) != 0
  occur in d_mu[y] (and in fact, that only those occur.) So it is a simple
  matter of looking up x.
*/

{
  const MuRow& mr = d_mu[y];

  if (not std::binary_search(mr.begin(),mr.end(),
	std::make_pair(x,static_cast<size_t>(0ul)),
			     MuCompare()))
    return 0;

  return std::lower_bound(mr.begin(),mr.end(),
	std::make_pair(x,static_cast<size_t>(0ul)),
			  MuCompare())->second;
}

}

/*****************************************************************************

        Chapter II -- The Helper class.

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {
  namespace helper {
Helper::Helper(const KLContext& kl)
  :KLContext(kl)

{}

/******** accessors **********************************************************/
MuCoeff Helper::ascentMu(size_t x, size_t y, size_t s) const

/*!
  \brief Computes mu(x,y) in a good ascent situation.

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
  case DescentStatus::RealNonparity: {
    return 0;
  }
  case DescentStatus::ImaginaryTypeI: {
    size_t x1 = cayley(s,x).first;
    return x1 == y ? 1 : 0;
  }
  case DescentStatus::ImaginaryTypeII: {
    BlockEltPair x1 = cayley(s,x);
    // mu(x,y) = mu(x1.first,y) + mu(x1.second,y)
    return (x1.first == y or x1.second == y) ? 1 : 0;
  }
  default: // this cannot happen
    assert(false);
    return MuCoeff(0); // keep compiler happy
  }
}

size_t Helper::firstDirectRecursion(size_t y) const

/*!
  \brief Returns the first descent generator that is not real type II

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

/*!
  \brief Gets mu(x,y) by a good descent recursion.

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
  case DescentStatus::RealTypeII: {
    size_t x1 = inverseCayley(s,x).first;
    return mu(x1,y1);
  }
  default: // this cannot happen
    assert(false);
    return MuCoeff(0); // keep compiler happy
  }
}

KLPolRef Helper::klPol(size_t x, size_t y,
			   KLRow::const_iterator klv,
			   klsupport::PrimitiveRow::const_iterator p_begin,
			   klsupport::PrimitiveRow::const_iterator p_end) const

/*!
  \brief Returns the Kazhdan-Lusztig polynomial for x corresponding to
  the given row.

  Precondition: klv holds the tail of the set of primitive Kazhdan-Lusztig
  polynomials for y, enough to find the required one by elementary lookup;
  [p_begin,p_end[ is the corresponding range of primitve elements.

  Algorithm: primitivize x w.r.t. the descents in y; if a real compact
  situation is encountered, return zero; otherwise look up the primitive
  x in the extremal range.
*/

{
  size_t xp = x;

  if (d_support->primitivize(xp,descentSet(y))) {
    if (std::binary_search(p_begin,p_end,xp)) { // result is nonzero
      size_t xpos = std::lower_bound(p_begin,p_end,xp) - p_begin;
      return d_store[klv[xpos]];
    } else {
      return d_store[d_zero];
    }
  } else {
    return d_store[d_zero];
  }
}

MuCoeff Helper::lengthOneMu(size_t x, size_t y) const

/*!
  \brief Computes mu(x,y) for l(y)-l(x) = 1.

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

void Helper::makeExtremalRow(klsupport::PrimitiveRow& e, size_t y) const

/*!
  \brief Puts in e the list of all x extremal w.r.t. y.

  Explanation: this means that either x = y, or length(x) < length(y),
  and every descent for y is a descent for x.
*/

{
  using namespace bitmap;

  BitMap b(size());
  size_t c = d_support->lengthLess(length(y));

  b.fill(c);
  b.insert(y);

  // extremalize
  d_support->extremalize(b,descentSet(y));
  // copy to list
  e.reserve(e.size()+b.size()); // ensure tight fit after copy
  std::copy(b.begin(),b.end(),back_inserter(e));

  return;
}

void Helper::makePrimitiveRow(klsupport::PrimitiveRow& e, size_t y) const

/*!
  \brief Puts in e the list of all x primitive w.r.t. y.

  Explanation: this means that either x = y, or length(x) < length(y),
  and every descent for y is either a descent, or an imaginary type II
  ascent for x.
*/

{
  using namespace bitmap;

  BitMap b(size());
  size_t c = d_support->lengthLess(length(y));

  b.fill(c);
  b.insert(y);

  // extremalize
  d_support->primitivize(b,descentSet(y));
  // copy to list
  e.reserve(e.size()+b.size()); // ensure tight fit after copy
  std::copy(b.begin(),b.end(),back_inserter(e));

  return;
}

MuCoeff Helper::recursiveMu(size_t x, size_t y) const

/*!
  \brief Gets mu(x,y) by direct recursion.

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

/*!
  \brief Gets mu(x,y) by type II recursion.

  Precondition: length(y) > 0; length(x) = length(y)-1; all previous
  mu-rows have been filled in; x is extremal w.r.t. y, and all descents
  of y are real type II;

  Explanation: in this situation, we have recursion formulas that will
  yield mu(x,y)+mu(x,s.y). It is known that proceeding in this way we must
  end up with a y' for which x is not extremal.

  Algorithm: we keep a stack of recursion formulas of the above kind, walking
  through the orbit of y under the simple descents, until we reach a y
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

/*!
  \brief Finishes the filling of the R-packet starting at y.

  Precondition: all the rows in the packet capable of a direct recursion
  are filled; at least one row is of this form;

  Explanation: what we do here, is completing the rows that can be completed
  from the already filled ones, using real type II recursions. Whatever is left
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
	std::vector<KLPol> klv;
	PrimitiveRow e;
	makeExtremalRow(e,y2);
	recursionRow(klv,e,y2,s);
	// klv[j] is P_{x,y2}+P_{x,y1}, for x = e[j]
	for (size_t j = 0; j < klv.size(); ++j) {
	  size_t x = e[j];
	  try {
	    klv[j].safeSubtract(klPol(x,y1));
	  }
	  catch (error::NumericUnderflow& e){
	    throw kl_error::KLError(x,y1,__LINE__,
				    static_cast<const KLContext&>(*this));
	  }
	}
	// write out row
	writeRow(klv,e,y2);
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

/*!
  \brief Fills in the row for y using a direct recursion.

  Precondition: s is either a complex, or a real type I descent generator for
  y.

  The real work is done by the recursionRow function, that will be used also
  in the real type II descents.
*/

{
  using namespace klsupport;

  std::vector<KLPol> klv;
  PrimitiveRow e;

  // put result of recursion formula in klv
  makeExtremalRow(e,y);
  recursionRow(klv,e,y,s);

  // write result
  writeRow(klv,e,y);

  return;
}

void Helper::fill()

/*!
  \brief Dispatches the work of filling the kl- and mu-lists.
*/

{
  // make sure the support is filled
  d_support->fill();

  // resize the lists
  d_prim.resize(d_support->size());
  d_kl.resize(d_support->size());
  d_mu.resize(d_support->size());

  // fill the lists
  size_t y = 0;
  size_t minLength = length(0);

  // do the minimal length cases; they come first in the enumeration
  for (; y < d_kl.size() and length(y) == minLength; ++y) {
    d_prim[y].push_back(y);
    // the k-l polynomial is 1
    d_kl[y].push_back(d_one);
    // there are no mu-coefficients
  }

  // do the other cases
  for (; y < d_kl.size(); ++y) {
#ifdef VERBOSE
    std::cerr << y << "\r";
#endif
    try {
      fillKLRow(y);
    } catch (kl_error::KLError& e) {
      e("error: negative coefficient in k-l construction");
    }
    fillMuRow(y);
  }

#ifdef VERBOSE
  std::cerr << std::endl;
#endif

  return;
}

void Helper::fillKLRow(size_t y)

/*!
  \brief Fills in the row for y in the kl-table.

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

/*!
  \brief Fills in the row for y in the mu-table.

  Precondition: the row for y in the kl-table has been filled; length(y) > 0;

  Explanation: for the elements of length < length(y) - 1, mu(x,y) can
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

  const PrimitiveRow& e = primitiveRow(y);

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
    if (d_store[klp].degree() < d)
      continue;
    // if we get here, we found a mu-coefficient for x
    d_mu[y].push_back(std::make_pair( x , d_store[klp][d] ));
  }

  // do cases of length ly-1
  size_t x_begin = d_support->lengthLess(ly-1);
  size_t x_end = d_support->lengthLess(ly);

  for (size_t x = x_begin; x < x_end; ++x) {
    MuCoeff mu = lengthOneMu(x,y);
    if (mu!=MuCoeff(0))
      d_mu[y].push_back(std::make_pair(x,mu));
  }

  return;
}

void Helper::fillThickets(size_t y)

/*!
  \brief Finishes the filling of the R-packet starting at y.

  Precondition: all the rows in the packet that are not part of a "thicket"
  have already been filled.

  Explanation: a thicket is an element that has only real type II descents,
  and that moreover is such that all elements in its conected component for
  the relation of being connected by a sequence of real type II cross-actions,
  is of the same form (for instance, it happens quite often that the principal
  series representations are a thicket.) Then we must use the "structural fact"
  in Lemma 6.2 of David's Park City notes: for each x lower than y, there is
  an element y' of the thicket for which x has an ascent (or if there is no
  such y', the k-l polynomial is zero.)

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

void Helper::muCorrection(std::vector<KLPol>& klv,
			  const klsupport::PrimitiveRow& e, size_t y, size_t s)

/*!
  \brief Subtracts from klv the correcting terms in the k-l recursion.

  Precondtion: klp contains the terms corresponding to c_s.c_y, for the x that
  are extremal w.r.t. y; the mu-table and kl-table has been filled in for
  elements of length <= y.

  Explanation: the recursion formula is of the form:

    lhs = c_s.c_{y1} - sum_{z} mu(z,y1)c_z

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
      // subtract x^d.mu.P_{x,z} from klv[j], where d = 1/2(l(y)-l(z))
      KLPolRef pol = klPol(x,z);
      Degree d = (l_y-l_z)/2;
      try {
	klv[j].safeSubtract(pol,d,mu);
      }
      catch (error::NumericUnderflow& e){
	throw kl_error::KLError(x,y,__LINE__,
				static_cast<const KLContext&>(*this));
      }
    }

  }

  return;
}

void Helper::recursionRow(std::vector<KLPol>& klv,
			  const klsupport::PrimitiveRow& e,size_t y, size_t s)

/*!
  \brief Puts in klv the right-hand side of the recursion formula for y
  corresponding to the descent s.

  Precondition: s is either a complex, or a real type I or type II descent
  generator for y.

  Explanation: the shape of the formula is:

    P_{x,y} = (c_s.c_{y1})-part - correction term

  where y1 = cross(s,y) when s is complex for y, one of the two elements in
  inverseCayley(s,y) when s is real. The (c_s.c_{y1})-part depends on the
  status of x w.r.t. s (we look only at extremal x, so we know it is a
  descent); the correction term, coming from sum_z mu(z,y1)c_z, depends only
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

  klv.resize(e.size());

  for (size_t j = 0; j < klv.size()-1; ++j) {
    size_t x = e[j];
    switch (descentValue(s,x)) {
    case DescentStatus::ImaginaryCompact: {
      // (q+1)P_{x,y1}
      klv[j] = klPol(x,y1).freeze();
      klv[j].safeAdd(klv[j],1);
    }
      break;
    case DescentStatus::ComplexDescent: {
      size_t x1 = cross(s,x);
      // P_{x1,y1}+q.P_{x,y1}
      klv[j] = klPol(x1,y1).freeze();
      klv[j].safeAdd(klPol(x,y1),1);
    }
      break;
    case DescentStatus::RealTypeI: {
      BlockEltPair x1 = inverseCayley(s,x);
      // P_{x1.first,y1}+P_{x1.second,y1}+(q-1)P_{x,y1}
      klv[j] = klPol(x1.first,y1).freeze();
      klv[j].safeAdd(klPol(x1.second,y1));
      klv[j].safeAdd(klPol(x,y1),1);
      try {
	klv[j].safeSubtract(klPol(x,y1));
      }
      catch (error::NumericUnderflow& e){
	throw kl_error::KLError(x,y,__LINE__,
				static_cast<const KLContext&>(*this));
      }
    }
      break;
    case DescentStatus::RealTypeII: {
      size_t x1 = inverseCayley(s,x).first;
      // P_{x_1,y_1}+qP_{x,y1}-P_{s.x,y1}
      klv[j] = klPol(x1,y1).freeze();
      klv[j].safeAdd(klPol(x,y1),1);
      try {
	klv[j].safeSubtract(klPol(cross(s,x),y1));
      }
      catch (error::NumericUnderflow& e){
	throw kl_error::KLError(x,y,__LINE__,
				static_cast<const KLContext&>(*this));
      }
    }
      break;
    default: // this cannot happen
      assert(false);
      break;
    }
  }

  // last k-l polynomial is 1
  klv.back() = One;

  // do mu-correction
  muCorrection(klv,e,y,s);

  return;
}

void Helper::writeRow(const std::vector<KLPol>& klv,
		      const klsupport::PrimitiveRow& er, size_t y)

/*!
  \brief Writes down row y in d_kl and d_prim.

  Precondition: klv contains the polynomials corresponding to the
  extremal values in the row; er contains the corresponding (extremal
  with respect to y) block elements.

  Explanation: the difficulty is that we want to write down the _primitive_
  elements, i.e., those x for which all descents for y are either descents
  or imaginary type II ascents; and moreover, we write down only those values
  for which the polynomial is nonzero. The values for which there is an
  imaginary type II ascent are computed "on the fly", from higher up
  values in the row.
*/

{
  using namespace blocks;
  using namespace klsupport;

  PrimitiveRow pr;
  makePrimitiveRow(pr,y);

  KLRow klr(pr.size());
  PrimitiveRow nzpr(pr.size());
  KLRow::iterator new_pol = klr.end();
  PrimitiveRow::iterator new_extr = nzpr.end();
  PrimitiveRow::iterator nzpr_end = nzpr.end();

  // set stops in pr at extremal elements
  std::vector<size_t> stop(er.size()+1);
  stop[0] = 0;

  for (size_t j = 0; j < er.size(); ++j) {
    size_t epos = std::lower_bound(pr.begin(),pr.end(),er[j]) - pr.begin() + 1;
    stop[j+1] = epos;
  }

  // write polynomials, beginning with largest x (which is y).
  for (size_t j = er.size(); j;) {
    --j;
    // insert extremal polynomial
    if (not klv[j].isZero()) {
      *--new_extr = er[j];
      *--new_pol  = d_store.match(klv[j]);
    }
    // do the others
    for (size_t i = stop[j+1]-1; i > stop[j];) {
      --i;
      size_t s = firstAscent(descent(pr[i]),descent(y),rank());
      BlockEltPair x1 = cayley(s,pr[i]);
      KLPol pol = klPol(x1.first,y,new_pol,new_extr,nzpr_end).freeze();
      pol.safeAdd(klPol(x1.second,y,new_pol,new_extr,nzpr_end));
      if (not pol.isZero()) {
	*--new_extr = pr[i];
	*--new_pol  = d_store.match(pol);
      }
    }
  }

  // commit
  d_prim[y].reserve(d_prim[y].size()+(nzpr.end()-new_extr));
  d_kl[y].reserve(d_kl[y].size()+(klr.end()-new_pol)); // ensure tight fit

  copy(new_extr,nzpr.end(),back_inserter(d_prim[y]));
  copy(new_pol,klr.end(),back_inserter(d_kl[y]));

  return;
}

  }
}

/*****************************************************************************

        Chapter III -- The Thicket class

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {
  namespace helper {

Thicket::Thicket(Helper& h, size_t y)
  :d_helper(&h)

/*!
  \brief Constructs a thicket from y.

  Explanation: the thicket of y is the connected component of y in the graph
  whose edges are the real type II cross-actions. It is therefore contained
  in the R-packet of y. The edges of the thicket are a spanning tree; each
  edge goes from y1 to y2, where y2 is obtained from y1 through real cross
  action by s (the position of y2 in the vertex list, and s, are recorded
  in the edge). The recursion recorded is the recursion for s as a descent of
  y2. In this way, we can move along a path of edges and have the required
  recursion formulas available for computing the polynomial for y2 from
  that of y1 (all edges are two-sided, i.e., we also make the corresponding
  edge from y2 to y1.)
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
	e1.push_back(Edge(y1,y2,s,std::vector<KLPol>()));
	e2.push_back(Edge(y2,y1,s,std::vector<KLPol>()));
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
    for (size_t i = 0; i < el.size(); ++i) {
      // replace el[i].y and el[i].source by its position in d_vertices
      el[i].source = std::lower_bound(d_vertices.begin(),d_vertices.end(),
				      el[i].source) - d_vertices.begin();
      el[i].y = std::lower_bound(d_vertices.begin(),d_vertices.end(),el[i].y)
	- d_vertices.begin();
    }
  }

  // make extremal lists
  d_extr.resize(thicket.size());

  for (size_t j = 0; j < d_extr.size(); ++j)
    d_helper->makeExtremalRow(d_extr[j],d_vertices[j]);

  // make primitive lists
  d_prim.resize(thicket.size());

  for (size_t j = 0; j < d_prim.size(); ++j)
    d_helper->makePrimitiveRow(d_prim[j],d_vertices[j]);

  // fill in recursions
  for (size_t j = 0; j < size(); ++j) {
    EdgeList& el = d_edges[j];
    for (size_t i = 0; i < el.size(); ++i) {
      size_t y2 = d_vertices[el[i].y];
      d_helper->recursionRow(el[i].recursion,d_extr[el[i].y],y2,el[i].s);
    }
  }

  // resize kl lists
  d_klr.resize(size());
  d_firstKL.resize(size());
  d_firstPrim.resize(size());
}

/******** accessors **********************************************************/
size_t Thicket::nonExtremal(size_t x) const

/*!
  \brief Returns the position in d_vertices of the first y that is not
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

bool Thicket::ascentCompute(size_t x, size_t pos)

/*!
  \brief Checks if x has an ascent in row y_pos, and computes the K-L pol in
  that case.

  Precondition: x is primitive in the row;

  Explanation: If the ascent exists, it will be imaginary type II. if
  x1 = (x1.first,x1.second) is the corresponding Cayley transform, the
  formula is P_x = P_x1.first + P_x1.second, both of which can be read off
  from the known part of the row.
*/

{
  using namespace blocks;

  size_t y = d_vertices[pos];
  size_t s = firstAscent(descent(x),descent(y),rank());

  if (s == rank()) // no ascent
    return false;

  BlockEltPair x1 = d_helper->cayley(s,x);
  KLPolRef t=klPol(x1.first,pos);
  KLPol pol = t.freeze();
  pol.safeAdd(klPol(x1.second,pos));

  if (not pol.isZero()) { // write pol
    --d_firstKL[pos];
    *d_firstKL[pos] = d_helper->d_store.match(pol);
    --d_firstPrim[pos];
    *d_firstPrim[pos] = x;
  }

  return true;
}

void Thicket::edgeCompute(size_t x, size_t pos, const Edge& e)

/*!
  \brief Computes the k-l polynomial for x in row pos, using the recurrence
  relation from e.

  Precondition: x is extremal in the row; e points towards pos; the polynomial
  for x for the source of e is known.

  Explanation: P_{x,pos} + P_{x,source} will be given by the recurrence
  relation corresponding to x.
*/

{
  using namespace klsupport;

  const PrimitiveRow& er = d_extr[pos];
  size_t xpos = std::lower_bound(er.begin(),er.end(),x) - er.begin();

  KLPol pol = e.recursion[xpos];

  try {
    pol.safeSubtract(klPol(x,e.source));
  }
  catch (error::NumericUnderflow& err){
    throw kl_error::KLError(x,d_vertices[pos],__LINE__,
			    static_cast<const KLContext&>(*d_helper));
  }

  if (not pol.isZero()) { // write pol
    *--d_firstKL[pos]   = d_helper->d_store.match(pol);
    *--d_firstPrim[pos] = x;
  }

  return;
}

void Thicket::fill()

/*!
  \brief Fills in the k-l polynomials for the elements of the thicket.

  Algorithm: for each element y_j of the thicket, we have the list of
  primitive elements pr[j], and a list of k-l polynomials klv[j]. We
  are going to fill in the polynomials downwards from the top (the top
  elements being all ones); at the end of the process, we will have
  iterators pointing into each pr[j] and each klv[j], where the
  extremal list from that iterator on will contain all x'es for which
  the klpol is non zero (in other words, the other ones have been
  "weeded out" from pr[j]), and klv[j] from the iterator on contains
  the corresponding polynomials. Then all that remains to do is copy
  these onto the corresponding lists of the context.  In order to
  achieve this, we put in d_xlist an ordered list of all elements that
  are primitive for _some_ y_j, not equal to y_j, and we move
  downwards in that list. It is guaranteed by Lemma 6.2 in David's
  notes, that for each x in d_xlist there is a y_j for which there is
  an easy reduction (i.e., x might be primitive for y_j, but not
  extremal.) So we can deduce P_{x,y_j} from the already known part of
  klv[j], and then traversing the tree structure of edge relations
  imposed on the thicket, we can compute all the other P_{x,y_i} that
  need to be (i.e., those for which x is in pr[j].) We record the
  result only if the corresponding polynomial is non-zero.
*/

{
  using namespace klsupport;

  fillXList();

  // initialize rows, fill in last element, initialize iterators
  for (size_t j = 0; j < size(); ++j) {
    const PrimitiveRow& pr = primitiveRow(j); // last element is y_j
    d_klr[j].resize(pr.size());
    d_klr[j].back() = d_helper->d_one;
    d_firstPrim[j] = d_prim[j].end()-1;
    d_firstKL[j] = d_klr[j].end()-1;
  }

  // fill in rows
  for (size_t j = d_xlist.size(); j;) {
    --j;
    size_t x = d_xlist[j];
    // find a vertex s.t. x is not extremal
    size_t iy = nonExtremal(x);
    if (iy == size()) // do nothing; zero polynomials are ignored
      continue;
    for (ThicketIterator i(*this,iy); i(); ++i) {
      size_t pos = *i;
      const PrimitiveRow& pr = primitiveRow(pos);
      if (not std::binary_search(pr.begin(),pr.end(),x))
	// x is not primitive
	continue;
      if (ascentCompute(x,pos)) // x is not extremal
	continue;
      // if we get here, x is extremal
      edgeCompute(x,pos,i.edge());
    }
  }

  // write result

  for (size_t j = 0; j < size(); ++j) {
    size_t y = d_vertices[j];

    d_helper->d_prim[y].reserve(d_helper->d_prim[y].size()+// ensure tight fit
				(d_prim[j].end()-d_firstPrim[j]));
    d_helper->d_kl[y].reserve(d_helper->d_kl[y].size()+    // ensure tight fit
			      (d_klr[j].end()-d_firstKL[j]));

    copy(d_firstPrim[j],d_prim[j].end(),back_inserter(d_helper->d_prim[y]));
    copy(d_firstKL[j],d_klr[j].end(),back_inserter(d_helper->d_kl[y]));
  }

  return;
}

void Thicket::fillXList()

/*!
  \brief Puts in d_xlist the union of the extremal lists for the elements
  of the thicket, top element excluded.
*/

{
  using namespace klsupport;

  std::set<size_t> xs;

  for (size_t j = 0; j < size(); ++j) {
    const PrimitiveRow& p = primitiveRow(j);
    for (size_t i = 0; i < p.size()-1; ++i)
      xs.insert(p[i]);
  }

  d_xlist.reserve(xs.size());
  copy(xs.begin(),xs.end(),std::back_inserter(d_xlist));

  return;
}

  }
}

/*****************************************************************************

        Chapter IV -- The ThicketIterator class

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {
  namespace helper {

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

/*!
  \brief Pre-increment operator.

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
	d_done.insert(i);
	return *this;
      }
    }
  }

  // if we get here, we have exhausted the thicket
  d_pos = d_thicket->size();

  return *this;
}

} //namespace helper
} //namespace kl


/*****************************************************************************

        Chapter V -- Functions declared in kl.h

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {

void wGraph(wgraph::WGraph& wg, const KLContext& klc)

/*!
  \brief Puts in wg the W-graph for this block.

  Explanation: the W-graph is a graph with one vertex for each element of the
  block; the corresponding descent set is the tau-invariant, i.e. the set of
  generators s that are either complex descents, real type I or II, or
  imaginary compact. Let x < y in the block such that mu(x,y) != 0, and
  descent(x) != descent(y). Then there is an edge from x to y unless descent(x)
  is contained in descent(y), and an edge from y to x unless descent(y) is
  contained in descent(x). Note that the latter always happens when the length
  difference is > 1, so that in that case there will only be an edge from
  x to y (the edge must be there because we already assumed that the descent
  sets were not equal.) In both cases, the coefficient corresponding to the
  edge is mu(x,y).

  NOTE: if I'm not mistaken, the edgelists come already out sorted.
*/

{
  using namespace bitset;
  using namespace graph;
  using namespace kl;
  using namespace wgraph;

  wg.reset();
  wg.resize(klc.size());

  // fill in descent sets
  for (size_t y = 0; y < klc.size(); ++y)
    wg.descent(y) = klc.descentSet(y);

  // fill in edges and coefficients
  for (size_t y = 0; y < klc.size(); ++y) {
    const RankFlags& d_y = wg.descent(y);
    const MuRow& mrow = klc.muRow(y);
    for (size_t j = 0; j < mrow.size(); ++j) {
      size_t x = mrow[j].first;
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
      RankFlags d = d_y;
      d &= d_x;
      if (d != d_y) { // d_y is not contained in d_x
	              // add edge from y to x
	wg.edgeList(y).push_back(x);
	wg.coeffList(y).push_back(mu);
      }
      if (d != d_x) { // d_x is not contained in d_y
	              // add edge from x to y
	wg.edgeList(x).push_back(y);
	wg.coeffList(x).push_back(mu);
      }
    }
  }

  return;
}

}

/*****************************************************************************

        Chapter VI -- Functions local to this module

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {
  namespace helper {

size_t firstAscent(const descents::DescentStatus& d1,
		   const descents::DescentStatus& d2, size_t rank)

/*!
  \brief Returns the first s that is an ascent for d1, and a descent for d2;
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

/*!
  \brief Returns the first s that is a non-ImaginaryTypeII ascent for
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
}
