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
#endif

#include <cassert>
#include <map>
#include <stack>
#include <stdexcept>

#include "bitmap.h"
#include "basic_io.h"
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

  namespace kl {
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

  };


  namespace helper {

  using namespace kl;

  size_t firstAscent(const descents::DescentStatus&,
		     const descents::DescentStatus&, size_t);

  size_t goodAscent(const descents::DescentStatus&,
		    const descents::DescentStatus&, size_t);

  class Thicket;


  typedef hashtable::HashTable<KLPolEntry,KLIndex> KLHashStore;

  class Helper
  : public KLContext
  {

    friend class Thicket;

    KLHashStore d_hashtable;

  public:

    // constructors and destructors
    Helper(const KLContext&);

    ~Helper() {}

    //accessors

    /*!
 \brief Cayley transform of block element y through simple root s.
    */
    blocks::BlockEltPair cayley(size_t s, BlockElt y) const {
      return d_support->block().cayley(s,y);
    }

    /*!
 \brief Cross action of simple root s on block element y.
    */
    BlockElt cross(size_t s, BlockElt y) const {
      return d_support->block().cross(s,y);
    }

    /*!
 \brief Returns vector of RANK_MAX unsigned char; entry s gives
 descent status of simple root s for block element y.
    */
    const descents::DescentStatus& descent(BlockElt y) const {
      return d_support->block().descent(y);
    }

    /*!
\brief Unsigned char whose value gives the descent status of
simple root s for block element y.
    */
    descents::DescentStatus::Value descentValue(size_t s, BlockElt y) const {
      return d_support->descentValue(s,y);
    }

    /*!
\brief Second coordinate (corresponding to K^vee orbit on G^vee/B^vee)
of pair of integers specifying block element y.
    */
    size_t dualOrbit(BlockElt y) const {
      return d_support->block().y(y);
    }

    size_t firstDirectRecursion(BlockElt y) const;

    blocks::BlockEltPair inverseCayley(size_t s, BlockElt y) const {
      return d_support->block().inverseCayley(s,y);
    }

    /* Declare that member function klPol from base class KLContext is also
       considered. This is necessary since it would otherwise be shadowed
       rather than overloaded by a new definition with different signature
       in the Helper class.
    */
    using KLContext::klPol;

    const KLPol& klPol(BlockElt x, BlockElt y, KLRow::const_iterator klv,
		       klsupport::PrimitiveRow::const_iterator p_begin,
		       klsupport::PrimitiveRow::const_iterator p_end) const;
    inline bool ascentMu(BlockElt x, BlockElt y, size_t s) const;

    inline MuCoeff goodDescentMu(BlockElt x, BlockElt y, size_t s) const;

    MuCoeff lengthOneMu(BlockElt x, BlockElt y) const;

    /*!
\brief First coordinate (corresponding to K orbit on G/B) of pair of
integers specifying block element y.
    */
    size_t orbit(BlockElt y) const {
      return d_support->block().x(y);
    }

    MuCoeff type2Mu(BlockElt x, BlockElt y) const;

    // manipulators
    void completePacket(BlockElt y);

    void directRecursion(BlockElt y, size_t s);

    void fill();

    void fillKLRow(BlockElt y);

    void fillMuRow(BlockElt y);

    void fillThickets(BlockElt y);

    void muCorrection(std::vector<KLPol>& klv,
		      const klsupport::PrimitiveRow& e,
		      BlockElt y, size_t s);

    void recursionRow(std::vector<KLPol> & klv,
		      const klsupport::PrimitiveRow& e, BlockElt y, size_t s);

    void writeRow(const std::vector<KLPol>& klv,
		  const klsupport::PrimitiveRow& e, BlockElt y);
  }; // class Helper




    // class Thicket

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
    std::vector<BlockElt> d_vertices;

    /*!
\brief Entry j lists the edges ending at y_j.
    */
    std::vector<EdgeList> d_edges;

    /*!
\brief List of all x that are extremal with respect to some y in Thicket.
    */
    std::vector<BlockElt> d_xlist;

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

    Thicket(Helper&, BlockElt);

    ~Thicket() {}

    // accessors

    /*!
\brief Cross action of simple reflection s on block element y.
    */
    BlockElt cross(size_t s, BlockElt y) const {
      return d_helper->cross(s,y);
    }

    /*!
\brief List of RankMax unsigned chars; number s gives the descent
status 0-7 of simple root s for y.
    */
    const descents::DescentStatus& descent(BlockElt y) const {
      return d_helper->descent(y);
    }


    /*!
\brief Unsigned char between 0 and 7; gives the descent of
simple root s for y.
    */
    descents::DescentStatus::Value descentValue(size_t s, BlockElt y) const {
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
    const KLPol& klPol(BlockElt x, size_t pos) const {
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

    size_t nonExtremal(BlockElt x) const;

    /*!
\brief Pointer to the KL polynomial one.
    */
    KLIndex one() const {
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
    const std::vector<BlockElt>& vertices() const {
      return d_vertices;
    }

    /*!
\brief Pointer to the KL polynomial zero.
    */
    KLIndex zero() const {
      return d_helper->d_zero;
    }

    // manipulators
    bool ascentCompute(BlockElt x, size_t pos);

    void edgeCompute(BlockElt x, size_t pos, const Edge& e);

    void fill();

    void fillXList();

    KLRow& klRow(BlockElt y) {
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
\brief Index in d_vertices of the block element at which the edge begins.
    */
    size_t source;

    /*!
\brief Index in d_vertices of the block element at which the edge ends.
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
    Edge(BlockElt x, BlockElt y, size_t s, const std::vector<KLPol>& r)
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

    size_t operator* () const { // current vertex index, or size() at end
      return d_pos;
    }

    const Thicket::Edge& edge() const { // current edge
      return *d_current;
    }

    // manipulators
    ThicketIterator& operator++();

  }; // class ThicketIterator

  }
  }

/*****************************************************************************

        Chapter I -- Methods of the KLContext and KLPolEntry classes.

 *****************************************************************************/

/* methods of KLContext */


namespace kl {
  using namespace atlas::kl::helper;

KLContext::KLContext(klsupport::KLSupport& kls)
  :d_state(),
   d_support(&kls),
   d_prim(),
   d_kl(),
   d_mu(),
   d_store(2),
   d_zero(0),
   d_one(1)
{
  d_store[d_zero]=Zero;
  d_store[d_one]=One;
}

/******** copy, assignment and swap ******************************************/


/*!
  \brief Copy constructor.

  Since we use indices instead of iterators, nothing gets invalidated any more
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

  d_store.swap(other.d_store);  // this puts the Helper store into base object
  std::swap(d_zero,other.d_zero);
  std::swap(d_one,other.d_one);
}

/******** accessors **********************************************************/

const KLPol& KLContext::klPol(BlockElt x, BlockElt y) const

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

  const PrimitiveRow& pr = d_prim[y];
  const KLRow& klr = d_kl[y];

  x=d_support->primitivize(x,descentSet(y));
  if (x>y) return d_store[d_zero]; // includes case x==blocks::UndefBlock
  PrimitiveRow::const_iterator xptr = std::lower_bound(pr.begin(),pr.end(),x);
  if (xptr == pr.end() or *xptr != x) // not found
    return d_store[d_zero];
  return d_store[klr[xptr - pr.begin()]];
}

MuCoeff KLContext::mu(BlockElt x, BlockElt y) const

/*!
  \brief Returns mu(x,y).

  Explanation: it is guaranteed that all the x'es such that mu(x,y) != 0
  occur in d_mu[y] (and in fact, that only those occur.) So it is a simple
  matter of looking up x.
*/

{
  const MuRow& mr = d_mu[y];

  std::vector<BlockElt>::const_iterator xloc=
    std::lower_bound(mr.first.begin(),mr.first.end(),x);

  if (xloc==mr.first.end() or *xloc!=x)
    return 0; // x not found in mr

  return mr.second[xloc-mr.first.begin()];
}

/* The following two methods were moved here form the Helper class, since
   they turn out to be useful even when no longer constructing the KLContext
*/
void KLContext::makeExtremalRow(klsupport::PrimitiveRow& e, BlockElt y) const

/*!
  \brief Puts in e the list of all x extremal w.r.t. y.

  Explanation: this means that either x = y, or length(x) < length(y),
  and every descent for y is a descent for x.
*/

{
  using namespace bitmap;

  BitMap b(size());
  size_t c = d_support->lengthLess(length(y));

  b.fill(c);     // start with all elements < y in length
  b.insert(y);   // and y itself

  // extremalize (filter out those that are not extremal)
  d_support->extremalize(b,descentSet(y));

  // copy from bitset b to list e
  e.reserve(e.size()+b.size()); // ensure tight fit after copy
  std::copy(b.begin(),b.end(),back_inserter(e));

}

void KLContext::makePrimitiveRow(klsupport::PrimitiveRow& e, BlockElt y) const

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

  b.fill(c);     // start with all elements < y in length
  b.insert(y);   // and y itself

  // primitivize (filter out those that are not primitive)
  d_support->primitivize(b,descentSet(y));
  // copy to list
  e.reserve(e.size()+b.size()); // ensure tight fit after copy
  std::copy(b.begin(),b.end(),back_inserter(e));

}


/******** manipulators *******************************************************/
void KLContext::fill()

/*!
  \brief Fills the KL- and mu-lists.

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
  swap(help); // swaps the base object, including the d_store

  d_state.set(KLFilled);

#ifdef VERBOSE
  std::cerr << "done" << std::endl;
#endif

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
{ const polynomials::Polynomial<KLCoeff>& P=*this;
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

    /* main constructor */

Helper::Helper(const KLContext& kl)
  :KLContext(kl)
  , d_hashtable(d_store) // hash table refers to our own d_store
{}

/******** accessors **********************************************************/

size_t Helper::firstDirectRecursion(BlockElt y) const

/*!
  \brief Returns the first descent generator that is not real type II

  Explanation: those are the ones that give a direct recursion formula for
  the K-L basis element.
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

const KLPol& Helper::klPol(BlockElt x, BlockElt y,
			   KLRow::const_iterator klv,
			   klsupport::PrimitiveRow::const_iterator p_begin,
			   klsupport::PrimitiveRow::const_iterator p_end) const

/*!
  \brief Returns the Kazhdan-Lusztig polynomial for x corresponding to
  the given row.

  Precondition: klv holds the tail of the set of primitive Kazhdan-Lusztig
  polynomials for y, enough to find the required one by elementary lookup;
  [p_begin,p_end[ is the corresponding range of primitive elements.

  Algorithm: primitivize x w.r.t. the descents in y; if a real nonparity
  situation is encountered, return zero; otherwise look up the primitive x in
  the range and return the corresponding element from klv
*/

{
  BlockElt xp = d_support->primitivize(x,descentSet(y));

  if (xp>y) return d_store[d_zero]; // includes case xp==blocks::UndefBlock
  klsupport::PrimitiveRow::const_iterator xptr =
    std::lower_bound(p_begin,p_end,xp);
  if (xptr == p_end or *xptr != xp) return d_store[d_zero];
  return d_store[klv[xptr-p_begin]];
}


inline bool Helper::ascentMu(BlockElt x, BlockElt y, size_t s) const

/*!
  \brief Computes whether mu(x,y)==1 in a good ascent situation (else it is 0)

  Preconditions: l(y) > 0; l(x) = l(y)-1; s is an ascent for x w.r.t. y;

  Explanation: this is the situation where mu(x,y) is directly expressible.
  The point is that x will ascend to an element of the same length as y,
  so that the corresponding K-L polynomial is zero unless x ascends to y.
*/

{
  using namespace blocks;
  using namespace descents;

  switch (descentValue(s,x)) {
  case DescentStatus::ComplexAscent: {
    BlockElt x1 = cross(s,x);
    return x1 == y;
  }
  case DescentStatus::RealNonparity: {
    return false;
  }
  case DescentStatus::ImaginaryTypeI: {
    BlockElt x1 = cayley(s,x).first;
    return x1 == y;
  }
  case DescentStatus::ImaginaryTypeII: {
    BlockEltPair x1 = cayley(s,x);
    // mu(x,y) = mu(x1.first,y) + mu(x1.second,y)
    return x1.first == y or x1.second == y;
  }
  default: // this cannot happen
    assert(false);
    return false; // keep compiler happy
  }
}

inline MuCoeff Helper::goodDescentMu(BlockElt x, BlockElt y, size_t s) const

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

  BlockElt y1;

  if (descentValue(s,y) == DescentStatus::ComplexDescent) {
    y1 = cross(s,y);
  } else { // s is real type I for y
    y1 = inverseCayley(s,y).first;
  }

  switch (descentValue(s,x)) {
  case DescentStatus::ImaginaryCompact: {
    return MuCoeff(0);
  }
  case DescentStatus::ComplexDescent: {
    BlockElt x1 = cross(s,x);
    return mu(x1,y1);
  }
  case DescentStatus::RealTypeI: {
    BlockEltPair x1 = inverseCayley(s,x);
    return mu(x1.first,y1)+mu(x1.second,y1);
  }
  case DescentStatus::RealTypeII: {
    BlockElt x1 = inverseCayley(s,x).first;
    return mu(x1,y1);
  }
  default: // this cannot happen
    assert(false);
    return MuCoeff(0); // keep compiler happy
  }
}


MuCoeff Helper::lengthOneMu(BlockElt x, BlockElt y) const

/*!
  \brief Computes mu(x,y) for l(y)-l(x) = 1.

  Preconditions: l(y) > 0; l(x) = l(y)-1; the mu-rows for y of smaller
  lengths have already been filled in.

  Explanation: these are the mu-values that can come up in possibly
  non-extremal situations. The value can be obtaind simply as the constant
  coefficient of klPol(x,y), but for efficiency reasons we handle some easy
  cases directly, avoiding the cost of calling klPol. In fact in all the cases
  these values of mu can be, and used to be, computed recursively by formulas
  using only other mu-values of the same kind.
*/

{
  // look if x has any ascents that are descents for y
  bitset::RankFlags ascents=descentSet(y); ascents.andnot(descentSet(x));
  if (ascents.any())
    return ascentMu(x,y,ascents.firstBit()) ? MuCoeff(1) : MuCoeff(0);

  /* Doing the following case separately costs more time than it gains

  // if we get here, x is extremal w.r.t. y
  s = firstDirectRecursion(y);
  if (s != rank()) // the answer is another mu
    return goodDescentMu(x,y,s);
  */

  // Now that the easy cases are gone, just lookup the KL polynomial
  KLPolRef p=klPol(x,y);
  return p.isZero() ? MuCoeff(0) : p[0];

  // the final case used to be: return type2Mu(x,y);

}

MuCoeff Helper::type2Mu(BlockElt x, BlockElt y) const

/*!
  \brief Gets mu(x,y) by type II recursion.

  This code is currently unused (see lengthOneMu above).

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

  std::map<BlockElt,size_t> rel;
  std::stack<BlockElt> toDo;

  rel.insert(std::make_pair(y,rank()));
  toDo.push(y);

  size_t s;
  BlockElt y1 = y;
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
	BlockElt y2 = cross(s,y1);
	if (rel.insert(std::make_pair(y2,s)).second) // new element
	  toDo.push(y2);
      }
    }
  }

 unwind:
  while (y1 != y) {
    std::map<BlockElt,size_t>::iterator i = rel.find(y1);
    s = i->second;
    y1 = cross(s,y1);
    // goodDescent() is mu + mu(x,y2)
    mu = goodDescentMu(x,y1,s) - mu;
  }

  return mu;
}

/******** manipulators *******************************************************/
void Helper::completePacket(BlockElt y)

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
  std::stack<BlockElt> filled;
  std::set<BlockElt> empty;

  for (BlockElt y1 = y; orbit(y1) == o; ++y1) {
    if (d_kl[y1].size() != 0)
      filled.push(y1);
    else
      empty.insert(y1);
  }

  while (not filled.empty()) {

    BlockElt y1 = filled.top();
    filled.pop();

    for (size_t s = 0; s < rank(); ++s) {
     if (descentValue(s,y1) != DescentStatus::RealTypeII)
	continue;
      BlockElt y2 = cross(s,y1);
      std::set<BlockElt>::iterator i = empty.find(y2);
      if (i != empty.end()) { // found a new row
 	// fill row y2
	std::vector<KLPol> klv;
	PrimitiveRow e;
	makeExtremalRow(e,y2);
	recursionRow(klv,e,y2,s);
	// klv[j] is P_{x,y2}+P_{x,y1}, for x = e[j]
	for (size_t j = 0; j < klv.size(); ++j) {
	  BlockElt x = e[j];
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

}

void Helper::directRecursion(BlockElt y, size_t s)

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

}

void Helper::fill()

/*!
  \brief Dispatches the work of filling the KL- and mu-lists.
*/

{
  // make sure the support is filled
  d_support->fill();

  // resize the outer lists to the block size
  d_prim.resize(d_support->size());
  d_kl.resize(d_support->size());
  d_mu.resize(d_support->size());

  // fill the lists
  BlockElt y = 0;
  size_t minLength = length(0);

  // do the minimal length cases; they come first in the enumeration
  for (; y < d_kl.size() and length(y) == minLength; ++y) {
    d_prim[y].push_back(y); // singleton list for this row
    // the K-L polynomial is 1
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

}

void Helper::fillKLRow(BlockElt y)

/*!
  \brief Fills in the row for y in the KL-table.

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
  for (BlockElt y1 = y; orbit(y1) == o; ++y1) {
    size_t s = firstDirectRecursion(y1);
    if (s != rank()) { // direct recursion
      directRecursion(y1,s);
    } else
      done = false;
  }

  // deal with the others
  if (not done)
    completePacket(y);

}

void Helper::fillMuRow(BlockElt y)

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

{
  using namespace klsupport;

  const PrimitiveRow& e = d_prim[y]; // list of nonzero polynomials

  size_t ly = length(y);
  if (ly==0) // we are in fact never called for |y| values of length 0
    return;  // but this is prudent, since next loop would fail for |ly==0|

  PrimitiveRow::const_iterator start= e.begin();
  // traverse lengths of opposite parity, up to ly-3
  for (size_t lx=(ly-1)%2,d = (ly-1)/2; d>0; --d,lx+=2) {// d=(ly-1-lx)/2

    PrimitiveRow::const_iterator stop =
      std::lower_bound(start,e.end(),d_support->lengthLess(lx+1));
    for (start= std::lower_bound(start,stop,d_support->lengthLess(lx));
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
  BlockElt x_begin = d_support->lengthLess(ly-1);
  BlockElt x_end = d_support->lengthLess(ly);

  for (BlockElt x = x_begin; x < x_end; ++x) {
    MuCoeff mu = lengthOneMu(x,y);
    if (mu!=MuCoeff(0)) {
      d_mu[y].first.push_back(x);
      d_mu[y].second.push_back(mu);
    }
  }

}

void Helper::fillThickets(BlockElt y)

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
  such y', the K-L polynomial is zero.)

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

  for (BlockElt y1 = y; orbit(y1) == o; ++y1)
    if (d_kl[y1].size() == 0) {
      Thicket thicket(*this,y1);
      thicket.fill();
    }

}

void Helper::muCorrection(std::vector<KLPol>& klv,
			  const klsupport::PrimitiveRow& e,
			  BlockElt y, size_t s)

/*!
  \brief Subtracts from klv the correcting terms in the K-L recursion.

  Precondtion: klp contains the terms corresponding to c_s.c_y, for the x that
  are extremal w.r.t. y; the mu-table and KL-table has been filled in for
  elements of length <= y.

  Explanation: the recursion formula is of the form:

    lhs = c_s.c_{y1} - sum_{z} mu(z,y1)c_z

  where z runs over the elements < y such that s is a descent for z,
  y1 is s.y, and lhs is c_y when s is a complex descent or real type I for y,
  and c_{y}+c_{s.y} when s is real type II.
*/

{
  using namespace descents;
  using namespace klsupport;
  using namespace polynomials;

  BlockElt y1;

  if (descentValue(s,y) == DescentStatus::ComplexDescent) {
    y1 = cross(s,y);
  } else { // s is real for y
    y1 = inverseCayley(s,y).first;
  }

  const MuRow& mrow = d_mu[y1];
  size_t l_y = length(y);

  for (size_t i = 0; i < mrow.first.size(); ++i) {

    BlockElt z = mrow.first[i];
    size_t l_z = length(z);

    DescentStatus::Value v = descentValue(s,z);
    if (not DescentStatus::isDescent(v))
      continue;

    MuCoeff mu = mrow.second[i]; // mu!=MuCoeff(0)

    Degree d = (l_y-l_z)/2; // power of q used in the loops below

    if (mu==MuCoeff(1)) // avoid useless multiplication by 1 if possible

      for (size_t j = 0; j < e.size(); ++j) {
	BlockElt x = e[j];
	if (length(x) > l_z) break;

	KLPolRef pol = klPol(x,z);
	try {
	  klv[j].safeSubtract(pol,d); // subtract q^d.P_{x,z} from klv[j]
	}
	catch (error::NumericUnderflow& e){
	  throw kl_error::KLError(x,y,__LINE__,
				  static_cast<const KLContext&>(*this));
	}
      } // for (j)

    else // mu!=MuCoeff(1)

      for (size_t j = 0; j < e.size(); ++j) {
	BlockElt x = e[j];
	if (length(x) > l_z) break;

	KLPolRef pol = klPol(x,z);
	try{
	  klv[j].safeSubtract(pol,d,mu); // subtract q^d.mu.P_{x,z} from klv[j]
	}
	catch (error::NumericUnderflow& e){
	  throw kl_error::KLError(x,y,__LINE__,
				  static_cast<const KLContext&>(*this));
	}
      } // for {j)

  } // for (i)

}

void Helper::recursionRow(std::vector<KLPol>& klv,
			  const klsupport::PrimitiveRow& e,BlockElt y, size_t s)

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

  BlockElt y1;

  if (descentValue(s,y) == DescentStatus::ComplexDescent) {
    y1 = cross(s,y);
  } else { // s is real type I or type II for y
    y1 = inverseCayley(s,y).first;
  }

  klv.resize(e.size());

  for (size_t j = 0; j < klv.size()-1; ++j) {
    BlockElt x = e[j];
    switch (descentValue(s,x)) {
    case DescentStatus::ImaginaryCompact: {
      // (q+1)P_{x,y1}
      klv[j] = klPol(x,y1);
      klv[j].safeAdd(klv[j],1);
    }
      break;
    case DescentStatus::ComplexDescent: {
      BlockElt x1 = cross(s,x);
      // P_{x1,y1}+q.P_{x,y1}
      klv[j] = klPol(x1,y1);
      klv[j].safeAdd(klPol(x,y1),1);
    }
      break;
    case DescentStatus::RealTypeI: {
      BlockEltPair x1 = inverseCayley(s,x);
      // P_{x1.first,y1}+P_{x1.second,y1}+(q-1)P_{x,y1}
      klv[j] = klPol(x1.first,y1);
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
      BlockElt x1 = inverseCayley(s,x).first;
      // P_{x_1,y_1}+qP_{x,y1}-P_{s.x,y1}
      klv[j] = klPol(x1,y1);
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

  // last K-L polynomial is 1
  klv.back() = One;

  // do mu-correction
  muCorrection(klv,e,y,s);

}

void Helper::writeRow(const std::vector<KLPol>& klv,
		      const klsupport::PrimitiveRow& er, BlockElt y)

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
      *--new_pol  = d_hashtable.match(klv[j]);
    }
    // do the others
    for (size_t i = stop[j+1]-1; i > stop[j];) {
      --i;
      size_t s = firstAscent(descent(pr[i]),descent(y),rank());
      BlockEltPair x1 = cayley(s,pr[i]);
      KLPol pol = klPol(x1.first,y,new_pol,new_extr,nzpr_end);
      pol.safeAdd(klPol(x1.second,y,new_pol,new_extr,nzpr_end));

      if (not pol.isZero()) {
	*--new_extr = pr[i];
	*--new_pol  = d_hashtable.match(pol);
      }
    }
  }

  // commit
  d_prim[y].reserve(klr.end() - new_pol);
  d_kl[y].reserve(klr.end() - new_pol);

  copy(new_extr,nzpr.end(),back_inserter(d_prim[y]));
  copy(new_pol,klr.end(),back_inserter(d_kl[y]));
}

  }
}

/*****************************************************************************

        Chapter III -- The Thicket class

  ... explain here when it is stable ...

 *****************************************************************************/

namespace kl {
  namespace helper {

Thicket::Thicket(Helper& h, BlockElt y)
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

  typedef std::map<BlockElt,EdgeList> M;
  typedef M::iterator MI;

  std::stack<BlockElt> toDo; // used as a bag (a queue would do just as well)
  M thicket;

  toDo.push(y);
  thicket.insert(std::make_pair(y,EdgeList()));

  while (not toDo.empty()) {
    BlockElt y1 = toDo.top();
    toDo.pop();
    EdgeList& e1 = thicket.find(y1)->second; // y1 is known to exist in thicket
    for (size_t s = 0; s < rank(); ++s) {
      if (descentValue(s,y1) != DescentStatus::RealTypeII)
	continue;
      BlockElt y2 = cross(s,y1);
      std::pair<MI,bool> ins = thicket.insert(std::make_pair(y2,EdgeList()));
      if (ins.second) { // y2 was new (otherwise just ignore edge y1 -> y2)
	EdgeList& e2 = ins.first->second; // locate empty list just inserted
	e1.push_back(Edge(y1,y2,s,std::vector<KLPol>()));
	e2.push_back(Edge(y2,y1,s,std::vector<KLPol>()));
	toDo.push(y2);
      }
    }
  }

  // write out result

  // first flatten map to pair of vectors d_vertices, d_edges
  d_vertices.resize(thicket.size());
  d_edges.resize(thicket.size());
  size_t j = 0;

  for (MI i = thicket.begin(); i != thicket.end(); ++i) {
    d_vertices[j] = i->first;
    d_edges[j].swap(i->second);
    ++j;
  }

  // now modify d_egdes so that the source and y fields of its edges
  // are no longer interpreted as BlockElt, but as index into d_vertices
  for (size_t j = 0; j < d_edges.size(); ++j) {
    EdgeList& el = d_edges[j];
    for (size_t i = 0; i < el.size(); ++i) {
      // replace el[i].y and el[i].source by their positions in d_vertices
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
      BlockElt y2 = d_vertices[el[i].y];
      d_helper->recursionRow(el[i].recursion,d_extr[el[i].y],y2,el[i].s);
    }
  }

  // resize kl lists
  d_klr.resize(size());
  d_firstKL.resize(size());
  d_firstPrim.resize(size());
}

/******** accessors **********************************************************/
size_t Thicket::nonExtremal(BlockElt x) const

/*!
  \brief Returns the position in d_vertices of the first y that is not
  extremal w.r.t. x, size() if there is no such y.
*/

{
  for (size_t j = 0; j < d_vertices.size(); ++j) {
    BlockElt y = d_vertices[j];
    size_t s = firstAscent(descent(x),descent(y),rank());
    if (s != rank()) // y was found
      return j;
  }

  return size();
}

/******** manipulators *******************************************************/

bool Thicket::ascentCompute(BlockElt x, size_t pos)

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

  BlockElt y = d_vertices[pos];
  size_t s = firstAscent(descent(x),descent(y),rank());

  if (s == rank()) // no ascent
    return false;

  BlockEltPair x1 = d_helper->cayley(s,x);
  KLPol pol = klPol(x1.first,pos);
  pol.safeAdd(klPol(x1.second,pos));

  if (not pol.isZero()) { // write pol
    *--d_firstKL[pos] = d_helper->d_hashtable.match(pol);
    *--d_firstPrim[pos] = x;
  }

  return true;
}

void Thicket::edgeCompute(BlockElt x, size_t pos, const Edge& e)

/*!
  \brief Computes the K-L polynomial for x in row pos, using the recurrence
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
    *--d_firstKL[pos]   = d_helper->d_hashtable.match(pol);
    *--d_firstPrim[pos] = x;
  }

}

void Thicket::fill()

/*!
  \brief Fills in the K-L polynomials for the elements of the thicket.

  Algorithm: for each element y_j of the thicket, we have the list of
  primitive elements pr[j], and a list of K-L polynomials klv[j]. We
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
  for (size_t j = d_xlist.size(); j>0;) {
    --j;
    BlockElt x = d_xlist[j];
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
    BlockElt y = d_vertices[j];

    d_helper->d_prim[y].reserve(d_prim[j].end() - d_firstPrim[j]);
    d_helper->d_kl[y].reserve(d_prim[j].end() - d_firstPrim[j]);

    copy(d_firstPrim[j],d_prim[j].end(),back_inserter(d_helper->d_prim[y]));
    copy(d_firstKL[j],d_klr[j].end(),back_inserter(d_helper->d_kl[y]));
  }

}

void Thicket::fillXList()

/*!
  \brief Puts in d_xlist the union of the extremal lists for the elements
  of the thicket, top element excluded.
*/

{
  using namespace klsupport;

  std::set<BlockElt> xs;

  for (size_t j = 0; j < size(); ++j) {
    const PrimitiveRow& p = primitiveRow(j);
    for (size_t i = 0; i < p.size()-1; ++i)
      xs.insert(p[i]);
  }

  d_xlist.reserve(xs.size());
  copy(xs.begin(),xs.end(),std::back_inserter(d_xlist));

}

  }
}

/*****************************************************************************

        Chapter IV -- The ThicketIterator class implementation

 *****************************************************************************/

namespace kl {
  namespace helper {

ThicketIterator::ThicketIterator(const Thicket& th, size_t j)
  : d_thicket(&th)

{
  d_stack.push(j);
  d_done.resize(th.size()); // all bits unset initially
  d_done.insert(j);
  d_current = th.edgeList(j).begin();
  d_currentEnd = th.edgeList(j).end();
  d_pos = j;
}

/******** manipulators *******************************************************/
ThicketIterator& ThicketIterator::operator++ ()

/*!
  \brief Pre-increment operator.

  Incrementing the iterator means going to the next new element, obtained as
  destination vertex of one of the remaining edges in the current edgelist, or
  in an edgelist for a vertex popped from the stack if necessary; here 'new'
  meaning the vertex has not yet been recorded in d_done. If no new
  destination vertex can be obtained even after trying all lists on the stack
  we return the thicket size as indication that no vertices are left. In the
  contrary case the vertex returned is recorded in d_done and pushed onto the
  stack for later processing of its neighbors.
*/

{
  for (; d_current != d_currentEnd; ++d_current) {
    // get index of destination of current edge in vertices
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

    // make edge list of element just popped the current one
    d_current = d_thicket->edgeList(j).begin();
    d_currentEnd = d_thicket->edgeList(j).end();

    // and do exactly what was done
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

 *****************************************************************************/

namespace kl {

void wGraph(wgraph::WGraph& wg, const KLContext& klc)

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

{
  using namespace bitset;
  using namespace graph;
  using namespace kl;
  using namespace wgraph;

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

}

} // namespace kl

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
