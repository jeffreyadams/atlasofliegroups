/*!
\file
\brief Implementation of the class KGB representing orbits of K on G/B.

  This module contains code for the construction of a block in the
  one-sided parameter set (in other words, the subset of the one-sided
  parameter set corresponding to a single real form.) As explained in
  my Palo Alto III notes, this is equivalent to parametrizing the set
  K\\G/B of (K,B)-orbits in G; hence the provocative title.
*/
/*
  This is kgb.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "kgb.h"

#include <cassert>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <list>

#include "bruhat.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "error.h"
#include "gradings.h"
#include "hashtable.h"
#include "realredgp.h"
#include "rootdata.h"
#include "tits.h"
#include "tori.h"
#include "weyl.h"

/*
  This module contains code for the construction of a block in the
  one-sided parameter set (in other words, the subset of the one-sided
  parameter set corresponding to a single real form.) As explained in
  my Palo Alto III notes, this is equivalent to parametrizing the set
  K\\G/B of (K,B)-orbits in G; hence the provocative title.

*/

namespace atlas {

  // namespace {
  namespace kgb {
    namespace {

/*
	                Local functions
*/

gradings::Grading baseGrading(const realredgp::RealReductiveGroup&);

void makeHasse(std::vector<set::SetEltList>&, const KGB&);



/*
	           Some small additional classes
*/

struct TI_Entry // To allow hash tables of TwistedInvolution values
  : public weyl::TwistedInvolution
{
  TI_Entry(const weyl::TwistedInvolution& tw): weyl::TwistedInvolution(tw) {}

  // members required for an Entry parameter to the HashTable template
  typedef std::vector<TI_Entry> Pooltype; // associated storage type
  size_t hashCode(size_t modulus) const;  // hash function
}; // class TI_Entry

struct TE_Entry // To allow hash tables of TitsElt values
  : public tits::TitsElt
{
  TE_Entry(const tits::TitsElt& t) : tits::TitsElt(t) {}

  // members required for an Entry parameter to the HashTable template
  typedef std::vector<TE_Entry> Pooltype; // associated storage type
  size_t hashCode(size_t modulus) const;  // hash function
}; // class TI_Entry

  /*!
\brief A |FiberData| object associates to each twisted involution a subspace
describing how corresponding Tits elements should be normalized.
  */
class FiberData
{
  TI_Entry::Pooltype pool;
  hashtable::HashTable<TI_Entry,unsigned int> hash_table;
  std::vector<latticetypes::SmallSubspace> data;
public:
  FiberData(const realredgp::RealReductiveGroup& GR);

  const latticetypes::SmallSubspace& mod_space(const tits::TitsElt& a) const
  {
    return data[hash_table.find(a.tw())];
  }
  void reduce(tits::TitsElt& a) const { mod_space(a).mod_reduce(a.t()); }

}; // class FiberData




/*
                 The |KGBHelp| class, a helper class for |KGB|

 */


class KGBHelp
{
  const size_t d_rank;
  const tits::TitsGroup& d_titsGroup;

  // data that will be exported (in permuted form) to the KGB class
  std::vector<KGBEltList> d_cross;
  std::vector<KGBEltList> d_cayley;
  std::vector<KGBInfo> d_info;

  // other data

  /*!
\brief List of Tits elements parametrizing KGB orbits.

  Access is usually via the hash table |d_tits|
  */
  TE_Entry::Pooltype d_pool;
  hashtable::HashTable<TE_Entry,KGBElt> d_tits;

  /*!
\brief Flags the noncompact imaginary roots among the simple roots for
G in the base theta-stable Borel.
  */
  gradings::Grading d_baseGrading;

  //! Permits reducing each Tits group element modulo its fiber denominator
  FiberData d_fiberData;

 public:

// constructors and destructors
  KGBHelp(realredgp::RealReductiveGroup&);

  ~KGBHelp() {};

// public accessor and manipulator:

//! deliver values to fields of a |KGB| object under construction.
  size_t export_tables(std::vector<KGBEltList>& cross,
		       std::vector<KGBEltList>& cayley,
		       weyl::TwistedInvolutionList& twisted,
		       std::vector<KGBInfo>& info) const;

//! fill the KGB set and return it
KGBHelp& fill();

// accessors (private)
private:
  void cayleyTransform(tits::TitsElt& a, size_t s) const {
    titsGroup().leftMult(a,s);
  }

  const tits::TitsGroup& titsGroup() const {
    return d_titsGroup;
  }

  const weyl::WeylGroup& weylGroup() const {
    return d_titsGroup.weylGroup();
  }

  size_t twist(size_t s) const {
    return d_titsGroup.twist(s);
  }

  void basedTwistedConjugate(tits::TitsElt& a, size_t s) const {
    titsGroup().twistedConjugate(a,s);
    if (d_baseGrading.test(s))
      a += d_titsGroup.simpleCoroot(s);
  }

// private manipulators
  void crossExtend(KGBElt parent);

  void cayleyExtend(KGBElt parent);

  void setStatus(KGBElt parent, KGBElt child, size_t s);

}; // class KGBHelp




/* the following auxiliary class serves for calling standard search and
   sorting routines for twisted involutions; it is used by |tauPacket|
*/
class InvolutionCompare {
private:
  const weyl::WeylGroup& W;
public:
  explicit InvolutionCompare(const weyl::WeylGroup& w) : W(w) {}

  // one should have a < b iff
  // (a) involutionLength(a) < involutionLength(b) or
  // (b) involutionLengths are equal and length(a) < length (b) or
  // (c) both lengths are equal and a < b
  bool operator()
   (const weyl::TwistedInvolution& a, const weyl::TwistedInvolution& b) const
  {
    if      (W.involutionLength(a) != W.involutionLength(b))
      return W.involutionLength(a) <  W.involutionLength(b) ;
    else if (W.length(a.w()) != W.length(b.w()))
      return W.length(a.w()) <  W.length(b.w());
    else
      return a < b;
  }
}; // class InvolutionCompare

/* This class is used when calling |std::sort| to sort the KGB elements by
   their associated twisted involutions, according to |InvolutionCompare|.
   However what will be actualy sorted is a list of indices, so we add as
   data member a reference to the vector that stores the Tits elements.
*/
class IndexCompare
  : public InvolutionCompare
{
private:
  const TE_Entry::Pooltype& pool;
public:
  explicit IndexCompare(const weyl::WeylGroup& W,const TE_Entry::Pooltype& p)
    : InvolutionCompare(W),pool(p) {}

  // here |unsigned long| arguments, as |setutils::Partition| contains such
  bool operator() (unsigned long i, unsigned long j) const {
    return InvolutionCompare::operator()(pool[i].tw(),pool[j].tw());
  }
}; // class Compare

} // namespace
} // namespace kgb


/*****************************************************************************

        Chapter I -- The KGB class, public methods

******************************************************************************/

namespace kgb {

/*!
  \brief Construct the KGB data structure for the given real form.
*/
KGB::KGB(realredgp::RealReductiveGroup& G)
  : d_rank(G.semisimpleRank())
  , d_cross()
  , d_cayley()
  , d_inverseCayley()
  , d_bruhat(NULL)
  , d_weylGroup(&G.weylGroup())
{
  size_t size =
    KGBHelp(G).fill().export_tables(d_cross,d_cayley,d_involution,d_info);

  // check that the size is right
  assert(size == G.kgbSize());

  // compute inverse Cayley tables

  d_inverseCayley.resize
    (d_rank,KGBEltPairList(size,std::make_pair(UndefKGB,UndefKGB)));

  for (KGBElt x=0; x<size; ++x)
    for (weyl::Generator s=0; s < d_rank; ++s)
      if (cayley(s,x)!=UndefKGB)
      {
	KGBEltPair& target=d_inverseCayley[s][cayley(s,x)];
	if (target.first==UndefKGB) target.first=x;
	else target.second=x;
      }
}

KGB::~KGB()

{
  delete d_bruhat;
}

/******** copy, assignment and swap ******************************************/


/******** accessors **********************************************************/


/*!
  \brief Returns the range of |KGBElt| values |z| with |involution(z)==w|.

  It is possible to return a range of |KGBElt|, since these elements have
  been sorted by considering the associated twisted involution first.

  [MvL: It was tempting to apply |std::find| with arguments of type
  |ctr_iterator::CounterIterator<KGBElt>| bounding the indices to
  |d_involution|, to find the range of such indices that compare equal to |w|
  under |KGB::compare|. I even believe this is why Fokko defined
  |CounterIterator| in the first place. The idea is doomed however, by the
  requirement that we supply to |std::find| as third argument a value of that
  type that itself compares equal to |w|, which means we would first have to
  find an appropriate index, which is as hard as our initial task (besides,
  such an index need not exist). Hence |InvolutionCompare| is needed instead.]
*/
KGBEltPair KGB::tauPacket(const weyl::TwistedInvolution& w) const
{
  using namespace weyl;

  typedef std::vector<TwistedInvolution>::const_iterator TI_it;

  InvolutionCompare comp(weylGroup());

  std::pair<TI_it,TI_it> range =
    std::equal_range(d_involution.begin(),d_involution.end(),w,comp);

  KGBElt first = range.first - d_involution.begin();
  KGBElt last = range.second - d_involution.begin();

  return std::make_pair(first,last);
}

/*!
  \brief Returns the length of the twisted involution |involution(z)|
   as a Weyl group element.

  NOTE: this is not inlined to avoid a dependency on weyl.h
*/
size_t KGB::weylLength(KGBElt z) const
{
  return weylGroup().length(d_involution[z].w());
}

/******** manipulators *******************************************************/

/*!
  \brief Constructs the BruhatOrder.

  NOTE: may throw a MemoryOverflow error. Commit-or-rollback is guaranteed.
*/
void KGB::fillBruhat()

{
  using namespace bruhat;
  using namespace set;

  if (d_state.test(BruhatConstructed)) // work was already done
    return;

  try {
    std::vector<SetEltList> hd; makeHasse(hd,*this);
    BruhatOrder* bp = new BruhatOrder(hd); // may throw here

    // commit
    delete d_bruhat; // this is overly careful: it must be NULL
    d_bruhat = bp;
    d_state.set(BruhatConstructed);
  }
  catch (std::bad_alloc) { // transform failed allocation into MemoryOverflow
    throw error::MemoryOverflow();
  }
}

} // namespace kgb





/*****************************************************************************

            Chapter II -- The auxiliary classes, methods

******************************************************************************/

// namespace {
namespace kgb {
  namespace {

size_t TI_Entry::hashCode(size_t modulus) const
{
  unsigned int hash=0;
  for (size_t i=constants::RANK_MAX; i-->0; )
    hash = 13*(hash+(*this)[i]);
  return hash & modulus-1;
}

size_t TE_Entry::hashCode(size_t modulus) const
{
  unsigned int hash=0;
  for (size_t i=constants::RANK_MAX; i-->0; )
    hash = 13*(hash+w()[i]);
  return hash+t().data().to_ulong() & modulus-1;
}


/*

  The |FiberData| constructor computes a subspace for each twisted involution
  $tw$, representing an involution $\tau$ of $H$ and an involution
  $q=tw.\delta$ of the weight lattice $X$. The space stores for it is the
  result of acting by $w^{-1}$ on the image $I$ in $X^\vee/2X^\vee$ of the
  $-1$ eigenspace $X^\vee_-$ of $q^t$. Since that eigenspace is stable by
  $(q^t)^{-1}=\delta.w^{-1}$, the subspace is also equal to $\delta.I$.

  The subspace $I$ is what one divides by to get the fiber group of the real
  Cartan associated to $H$ and $\tau$. The information of $w$ together with a
  point of the fiber can be combined in a Tits group element of the form $t.w$
  where $t\in X^\vee/2X^\vee$ lies in the image of $X^\vee_+ + X^\vee_-$, and
  only its coset modulo $I$ matters. However such an element is in fact stored
  in the form $w.t'$ where $t'=w^{-1}(t)$, and this torus part $t'$ is only
  relevant modulo the space stored in the |FiberData| object. We shall work
  with Tits elements whose torus parts are distinguished representatives of
  these cosets. The method |reduce| transforms Tits elements to such a form.

  $I$ also corresponds to the subset of elements of order 2 in $H$ that
  are of the form $h*\tau(h)^{-1}$.
*/
FiberData::FiberData(const realredgp::RealReductiveGroup& GR)
  : pool()
  , hash_table(pool)
  , data()
{
  const complexredgp::ComplexReductiveGroup& G=GR.complexGroup();
  const weyl::WeylGroup& W=G.weylGroup();
  const rootdata::RootDatum& rd=G.rootDatum();
  latticetypes::BinaryMap delta(G.distinguished().transposed());

  std::vector<latticetypes::BinaryMap> refl(G.semisimpleRank());
  for (weyl::Generator s=0; s<refl.size(); ++s)
  {
   // get endomorphism of weight lattice $X$ given by generator $twist(s)$
    latticetypes::LatticeMatrix r =
      rd.rootReflection(rd.simpleRootNbr(G.titsGroup().twist(s)));

    // reflection map is induced vector space endomorphism of $X^* / 2X^*$
    refl[s] = latticetypes::BinaryMap(r.transposed());
  }

  data.reserve(G.numInvolutions());
  for (bitmap::BitMap::iterator it=GR.cartanSet().begin(); it(); ++it)
  {
    cartanclass::CartanClass cc=G.cartan(*it);
    size_t i =
      hash_table.match
      (weyl::TwistedInvolution(weyl::WeylElt
	(rd.word_of_inverse_matrix(G.distinguished()*cc.involution()),W)
      ));

    assert(i==data.size()); // this twisted involution should be new
    {
      using namespace latticetypes;
      LatticeMatrix qtr= cc.involution().transposed();
      data.push_back(SmallSubspace(SmallBitVectorList(tori::minusBasis(qtr)),
				   G.rank()));
      data[i].apply(delta);
    }

    for ( ; i<data.size(); ++i) // |data.size()|  increases during loop
      for (size_t s=0; s<G.semisimpleRank(); ++s)
      {
	weyl::TwistedInvolution stw=W.twistedConjugated(pool[i],s);
	if (hash_table.match(stw)==data.size()) // then |stw| is new
	{
	  data.push_back(data[i]);
	  data.back().apply(refl[s]);
	}
      }
  }
}


/*!
  \brief The helper constructor just initializes the lists for the first
  element. The actual filling process is handled by the |fill| method.

  Algorithm: from the datum of the real form (or more precisely, from the
  corresponding fundamental grading) we recover the basic cocycle that
  transforms the whole construction into a Tits group computation.
*/

KGBHelp::KGBHelp(realredgp::RealReductiveGroup& GR)
  : d_rank(GR.semisimpleRank())
  , d_titsGroup(GR.titsGroup())

  , d_cross(d_rank)
  , d_cayley(d_rank)
  , d_info()

  , d_pool(1,tits::TitsElt(GR.rank())) // single element, identity
  , d_tits(d_pool)

  // only the final two fields use the real form of |GR|
  , d_baseGrading(baseGrading(GR))
  , d_fiberData(GR)
{
  KGBElt size = GR.kgbSize();

  // initialize descent, length, and associated (twisted) involution
  d_info.reserve(size); d_info.push_back(0);

  // set up cross and cayley tables with undefined values
  for (size_t j = 0; j < d_rank; ++j) {
    d_cross[j].resize(size,UndefKGB);
    d_cayley[j].resize(size,0); // leave undefined (set by |cayleyExtend|)
  }
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/

/*!
  \brief Constructs the full orbit set.

  Precondition: the object is in the initial state, the one it is put in by
  the call to its constructor;

  Algorithm: The idea is just to start out from the fundamental orbit,
  and then to saturate through cross actions and Cayley transforms.

  It is important that for each element the cross actions are defined before
  the Cayley transforms is tried
*/
KGBHelp& KGBHelp::fill()
{
  for (KGBElt j = 0; j < d_pool.size(); ++j) {
    // these calls will usually enlarge |d_pool.size()|;
    crossExtend(j);
    cayleyExtend(j);
  }
  return *this;
}

/*!
  \brief Tries to enlarge the parameter set by cross-actions.

  Precondition: |parent| is the index into |d_tits| of the parameter we are
  extending from.

  We use the fact that the |KGBHelp| structure is always kept in the state
  where all downward links are already filled in.
*/
void KGBHelp::crossExtend(KGBElt parent)
{
  const tits::TitsElt current = d_tits[parent];

  // try to get new elements by cross-action
  for (size_t s = 0; s < d_rank; ++s) {
    // see if link was already filled
    if (d_cross[s][parent] != UndefKGB)
      continue;

    // twisted-conjugate (using base grading) |current| by |s|
    tits::TitsElt a = current; basedTwistedConjugate(a,s);

    /* Check that result is independent of representative:
    {
      latticetypes::SmallBitVectorList b =
	d_fiberData.mod_space(current).basis();
      for (size_t i=0; i<b.size(); ++i)
      {
	TitsElt ai=current; ai+=b[i]; basedTwistedConjugate(ai,s);
	assert(ai.tw()==a.tw());
	d_fiberData.reduce(ai);
	assert(ai.t()==a.t());
      }
    }
    */

    // find the Tits element
    d_fiberData.reduce(a);
    KGBElt x = d_tits.match(a);
    if (x==d_info.size()) // add a new Tits element
    {
      bool same_fiber=d_tits[x].tw()==d_tits[parent].tw();
      d_info.push_back(d_info[parent].length+(same_fiber ? 0 : 1));
    }

    // set descent direction of complex |s| for |x| (whether it was new or old)
    if (d_info[parent].length!=d_info[x].length)
      d_info[x].desc.set(s,d_info[parent].length < d_info[x].length);

    // add new cross-links (with harmless redundancy if |parent==x|)
    d_cross[s][parent] = x;
    d_cross[s][x] = parent;

    // update status info
    setStatus(parent,x,s);
  }
}

/*!
  \brief Fills in status information between |from| and |to==cross[s][from]|
  as well as |cayley[s][from]|.

  Precondition: cross[s][from] = to;

  Explanation: if |from| and |to| have different root datum involutions, then
  |s| is complex for both of them. Otherwise, |s| is imaginary for |from|, and
  real for |cayley[s][from]|. To find whether |s| is compact or noncompact, we
  note that the grading may be computed from the base grading and the torus
  part of the d_tits[from].
*/
void KGBHelp::setStatus(KGBElt from, KGBElt to, size_t s)
{
  if (d_info[from].status[s] == gradings::Status::Real)
    return;

  if (d_tits[from].tw() != d_tits[to].tw()) {
    d_info[from].status.set(s,gradings::Status::Complex);
    d_info[to].status.set(s,gradings::Status::Complex);
  }
  else
  { // now the root $\alpha_s$ is imaginary
    bool b = d_baseGrading.test(s);
    b ^= scalarProduct(d_tits[from].t(),titsGroup().simpleRoot(twist(s)));

    if (b)
    {
      d_info[from].status.set(s,gradings::Status::ImaginaryNoncompact);
      d_info[to].status.set(s,gradings::Status::ImaginaryNoncompact);
    }
    else // here we certainly have |from == to|
      d_info[from].status.set(s,gradings::Status::ImaginaryCompact);
  }
}

/*!
  \brief Tries to enlarge the parameter set by cayley transforms from |parent|.

  Precondition: |parent| is the |KGBElt| (index in |d_pool|) we are extending
  from.

  It is assumed that the it is an invariant of the |KGBHelp| structure that
  all downward links are already filled in.
*/
void KGBHelp::cayleyExtend(KGBElt parent)
{
  const tits::TitsElt current = d_tits[parent];

  // try to get new elements by Cayley transform
  for (size_t s = 0; s < d_rank; ++s) {

    gradings::Status::Value v = d_info[parent].status[s];
    if (v != gradings::Status::ImaginaryNoncompact) {
      d_cayley[s][parent] = UndefKGB;
      continue;
    }

    // cayley-transform |current| by $\sigma_s$
    tits::TitsElt a = current; cayleyTransform(a,s);
    assert(a.t()==current.t()); // involution part should go up


    // now look up the correspondingly reduced Tits element
    d_fiberData.reduce(a); // subspace has grown, so mod out new supspace
    KGBElt x = d_tits.match(a);
    if (x==d_info.size()) // add a new Tits element
      d_info.push_back(d_info[parent].length+1); // length goes always up

    // add new cayley link
    d_cayley[s][parent] = x;

    // update descent table
    d_info[x].desc.set(s);

    // add status value
    d_info[x].status.set(s,gradings::Status::Real);
  }
}

size_t KGBHelp::export_tables(std::vector<KGBEltList>& cross,
			      std::vector<KGBEltList>& cayley,
			      weyl::TwistedInvolutionList& twisted,
			      std::vector<KGBInfo>& info) const
{
  using namespace setutils;
  size_t size=d_info.size();

  // sort (partially)
  Permutation a(size,1);
  IndexCompare comp(weylGroup(),d_pool);
  std::sort(a.begin(),a.end(),comp);

  // export the cross and cayley maps, permuting each constituent list
  cross.clear(); cayley.clear();
  cross.reserve(d_rank); cayley.reserve(d_rank);
  Permutation ai(a,-1);
  for (size_t s = 0; s < d_rank; ++s) {
    cross.push_back(a.pull_back(ai.renumber(d_cross[s])));
    cayley.push_back(a.pull_back(ai.renumber(d_cayley[s],UndefKGB)));
  }

  // export the other lists
  twisted.clear(); twisted.reserve(size);
  for (KGBElt x=0; x<size; ++x)
    twisted.push_back(d_pool[a[x]].tw());

  a.pull_back(d_info).swap(info);

  return size;
}

} // namespace
} // namespace kgb

/*****************************************************************************

        Chapter III -- Functions local to kgb.cpp

******************************************************************************/

namespace kgb {
  namespace {

/*!
  \brief Returns the the simple roots which are noncompact
  imaginary for the base grading of the fundamental Cartan in G..

Algorithm: the variable rset is first made to flag the noncompact
imaginary roots for the base grading of the fundamental Cartan in G.
Each simple root for G is tested for membership in rset; if it passes,
the corresponding bit of gr is set.
*/
gradings::Grading baseGrading(const realredgp::RealReductiveGroup& G)
{
  const rootdata::RootDatum& rd = G.rootDatum();
  rootdata::RootSet rset= G.noncompactRoots();

  gradings::Grading result;

  for (size_t s = 0; s < G.semisimpleRank(); ++s) {
    if (rset.isMember(rd.simpleRootNbr(s)))
      result.set(s);
  }

  return result;
}

void makeHasse(std::vector<set::SetEltList>& hd, const KGB& kgb)

/*!
  \brief Puts in hd the hasse diagram data for the Bruhat ordering on KGB.

  Explanation: this is the closure ordering of orbits. We use the
  algorithm from Richardson and Springer.
*/

{
  using namespace gradings;
  using namespace poset;
  using namespace set;

  hd.resize(kgb.size());

  for (KGBElt x = 0; x < kgb.size(); ++x) {

    const Descent& d = kgb.descent(x);
    if (d.none()) // element is minimal in bruhat order
      continue;

    size_t s = d.firstBit();
    KGBElt sx;

    if (kgb.status(s,x) == Status::Complex) { // s is complex
      sx = kgb.cross(s,x);
      hd[x].push_back(sx);
    } else { // s is real for x
      KGBEltPair sxp = kgb.inverseCayley(s,x);
      hd[x].push_back(sxp.first);
      if (sxp.second != UndefKGB) // s is real type I for x
	hd[x].push_back(sxp.second);
      sx = sxp.first;
    }

    for (size_t j = 0; j < hd[sx].size(); ++j) {
      KGBElt z = hd[sx][j];
      if (kgb.isAscent(s,z)) {
	if (kgb.status(s,z) == Status::ImaginaryNoncompact)
	  hd[x].push_back(kgb.cayley(s,z));
	else // s is complex for z
	  hd[x].push_back(kgb.cross(s,z));
      }
    }

    std::sort(hd[x].begin(),hd[x].end());
  }

}

} // namespace
} // namespace kgb

} // namespace atlas
