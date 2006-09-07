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

#include "bruhat.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "ctr_iterator.h"
#include "error.h"
#include "gradings.h"
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
    namespace helper {

  void pause() {;}

using namespace kgb;

void makeHasse(std::vector<set::SetEltList>&, const KGB&);

  /*!
\brief Used with an involution tau of H to store the subspace of
elements of order 2 in H of the form h*tau(h)^{-1}.

This is what one divides by to get the fiber group of the real Cartan.  A
Tits element representing an H conjugacy classes of strong involutions
is defined only as a coset of this subgroup.  We choose our Tits
elements as distinguished representatives of these cosets.  The member
function normalize gets these distinguished representatives, and the
member function conjugate changes FiberData when tau is replaced by a
conjugate.
  */
class FiberData {

  /*!
\brief Intended to allow modification of FiberData by an external
function conjugate.

[There seems to exist no such external function; presumably it was
replaced by the member function conjugate.  DV 7/21/06]
  */ 
  friend void conjugate(FiberData&, size_t, const tits::TitsGroup&);

 private:

  /*!
\brief Subspace of elements of order 2 in the Cartan of the form
h*tau(h)^{-1}.

This is stored as a subspace of the Z/2Z vector space of all elements
of order 2 (which is naturally (Z/2Z)^n, since the coweight lattice is
Z^n), endowed with the unique basis in row-reduced form.
  */
  latticetypes::ComponentSubspace d_subspace;

 public:
// constructors and destructors
  FiberData(const latticetypes::ComponentSubspace& sp):d_subspace(sp) {}
  FiberData(const latticetypes::ComponentList& b, size_t r):d_subspace(b,r) {}
  
// accessors
  /*!  
\brief Replaces the Tits element by the normal representative in its H
conjugacy class.
  */
  void normalize(tits::TitsElt& a) const {
    d_subspace.representative(a.t(),a.t());
  }

// manipulators
  /*!
\brief changes FiberData when the defining involution tau is replaced
  by r tau r^{-1}.
  */
  void conjugate(const latticetypes::ComponentMap& r) {
    d_subspace.apply(r);
  }
};

      /*! 
\brief Derived class of KGB, to carry out the construction of KGB.
      */
class Helper:public KGB {

private:

  /*!
\brief Constant iterator to the pair containing the WeylElt, or to the
  end of the map if no such pair exists. 
  */
  typedef std::map<weyl::WeylElt,FiberData>::const_iterator WI;

  /*!
\brief Constant iterator to the pair containing the TitsElt, or to the
  end of the map if no such pair exists. 
  */
  typedef std::map<tits::TitsElt,size_t>::const_iterator TI;

  // big auxiliary data

  /*!
\brief List of Tits elements parametrizing KGB orbits.

More precisely, these are the distinguished representatives in the H
conjugacy classes. This ensures that distinct elements correspond to
distinct orbits.
  */
  tits::TitsEltList d_tits;

  /*!
\brief Inverse function of d_tits.

This map allows us to search for a Tits element to see if it
is new. 
  */
  std::map<tits::TitsElt,size_t> d_titsMap;

  /*!
\brief Associates to each root datum involution tau the subspace of
elements of order 2 in the Cartan of the form h*tau(h)^{-1}.

This is what one divides to get the fiber group of the real Cartan.  A
Tits element representing an H conjugacy classes of strong involutions
is defined only as a coset of this subgroup.  We choose our Tits
elements as distinguished representatives of these cosets.  The
function FiberData::normalize gets these distinguished
representatives, and FiberData::conjugate changes FiberData when tau
is replaced by a conjugate.
  */
  std::map<weyl::WeylElt,FiberData> d_fiberData;

  // small auxiliary data

  /*!
\brief Pointer to the real group whose K orbits on G/B are being
  computed.
  */
  const realredgp::RealReductiveGroup* d_realGroup;

  /*!
\brief Flags the noncompact imaginary roots among the simple roots for
G in the base theta-stable Borel.
  */
  gradings::Grading d_baseGrading;

  // these lists have rank entries

  /*!
\brief Entry \#i is the mod 2 matrix of simple reflection \#i.
  */
  std::vector<latticetypes::ComponentMap> d_reflectionMap;

  /*!
\brief Entry \#i is the matrix of simple reflection \#i.
  */
  std::vector<latticetypes::LatticeMatrix> d_reflectionMatrix;

  /*!
\brief The mod 2 matrix of the distinguished involution for the inner
  class (transposed).
  */
  latticetypes::ComponentMap d_twistMatrix;

public:
  
// constructors and destructors
  Helper(realredgp::RealReductiveGroup&);

  virtual ~Helper() {};

// accessors
  void cayleyTransform(tits::TitsElt& a, size_t s) const {
    titsGroup().leftProd(a,s);
  }

  WI find (const weyl::WeylElt& w) const {
    return d_fiberData.find(w);
  }

  TI find (const tits::TitsElt& a) const {
    return d_titsMap.find(a);
  }

  bool isNew(WI i) const {
    return i == d_fiberData.end();
  }

  bool isNew(TI i) const {
    return i == d_titsMap.end();
  }

  const realredgp::RealReductiveGroup& realGroup() const {
    return *d_realGroup;
  }

  const latticetypes::ComponentMap& reflectionMap(size_t s) const {
    return d_reflectionMap[s];
  }

  const rootdata::RootDatum& rootDatum() const {
    return d_realGroup->rootDatum();
  }

  size_t size() const {
    return d_tits.size();
  }

  const tits::TitsGroup& titsGroup() const {
    return d_realGroup->titsGroup();
  }

  size_t twist(size_t s) const {
    return titsGroup().twist(s);
  }

  void twistedConjugate(tits::TitsElt& a, size_t s) const {
    titsGroup().twistedConjugate(a,s);
    if (d_baseGrading.test(s))
      a += titsGroup().simpleCoroot(s);
  }

// manipulators
  WI cayleyAdd(weyl::WeylElt&, size_t);

  TI cayleyAdd(size_t, size_t);

  WI crossAdd(weyl::WeylElt&, size_t);

  TI crossAdd(size_t, size_t);

  void crossExtend(size_t);

  void cayleyExtend(size_t);

  void fill();

  void saturate(const weyl::WeylElt&);

  void setStatus(size_t, size_t, size_t);

  void sort();
};

class Compare{
private:
  const KGB* d_kgb;
public:
  explicit Compare(const KGB& kgb):d_kgb(&kgb) {}
  // one should have i < j iff (a) d_length[i] < d_length[j] or
  // (b) lengths are equal and d_involution[i] < d_involution[j]
  bool operator() (unsigned long i, unsigned long j) const {
    if (d_kgb->length(i) < d_kgb->length(j))
      return true;
    else if (d_kgb->length(j) < d_kgb->length(i))
      return false;
    else if (d_kgb->weylLength(i) < d_kgb->weylLength(j))
      return true;
    else if (d_kgb->weylLength(j) < d_kgb->weylLength(i))
      return false;
    else
      return d_kgb->involution(i) < d_kgb->involution(j);
  }
};

class InvolutionCompare {
private:
  const weyl::WeylGroup* d_W;
public:
  explicit InvolutionCompare(const weyl::WeylGroup& W):d_W(&W) {}
  // one should have i < j iff 
  // (a) involutionLength(v) < involutinLength(w) or
  // (b) involutionLengths are equal and length(v) < length (w) or
  // (c) both lengths are equal and v < w
  bool operator() (const weyl::WeylElt& v, const weyl::WeylElt& w) const {
    if (d_W->involutionLength(v) < d_W->involutionLength(w))
      return true;
    else if (d_W->involutionLength(w) < d_W->involutionLength(v))
      return false;
    else if (d_W->length(v) < d_W->length(w))
      return true;
    else if (d_W->length(w) < d_W->length(v))
      return false;
    else
      return v < w;
  }
};
  
    }
  }

  // namespace {
  namespace kgb {
    namespace helper {
  using namespace kgb;

  void initGrading(gradings::Grading&, const realredgp::RealReductiveGroup&);
    }
  }

/*****************************************************************************

        Chapter I -- The KGB class

  ... explain here when it is stable ...

******************************************************************************/

namespace kgb {

  using namespace atlas::kgb::helper;

KGB::KGB(size_t rank)   
  :d_rank(rank),
   d_cross(d_rank),
   d_cayley(d_rank),
   d_inverseCayley(d_rank),
   d_bruhat(0)

/*!
  \brief Does only a trivial initialization.

  Explanation: preliminary to the Helper constructor; rank is the semisimple 
  rank of the underlying group.
*/

{}


KGB::KGB(realredgp::RealReductiveGroup& G)
  : d_bruhat(0)

/*!
  \brief Constructs the KGB data structure for the given real form.
*/

{
  using namespace bruhat;

  pause();
  Helper kgbhelp(G);

  // construct the orbit set
  kgbhelp.fill();

  // sort
  kgbhelp.sort();

  // transfer data to *this
  swap(kgbhelp);

  // check that the size is right
  assert(size() == G.kgbSize());
}

KGB::~KGB()

{
  delete d_bruhat;
}

/******** copy, assignment and swap ******************************************/

void KGB::swap(KGB& other)

{  
  std::swap(d_rank,other.d_rank);

  d_cross.swap(other.d_cross);
  d_cayley.swap(other.d_cayley);
  d_inverseCayley.swap(other.d_inverseCayley);
  d_descent.swap(other.d_descent);
  d_length.swap(other.d_length);
  d_involution.swap(other.d_involution);
  d_status.swap(other.d_status);

  d_state.swap(other.d_state);

  std::swap(d_bruhat,other.d_bruhat);

  std::swap(d_weylGroup,other.d_weylGroup);

  return;
}

/******** accessors **********************************************************/

bool KGB::compare (KGBElt x, KGBElt y) const

/*!
  \brief Returns whether involution(x) < involution(y). 

  Explanation: the ordering is involution-length first, then weyl-length, 
  then order of Weyl group elements in whatever way that comes out of their 
  representation. This is not an ordering, only a preordering of course.
*/

{
  if (length(x) < length(y))
    return true;
  else if (length(y) < length(x))
    return false;
  else if (weylLength(x) < weylLength(y))
    return true;
  else if (weylLength(y) < weylLength(x))
    return false;
  else
    return involution(x) < involution(y);
}

const weyl::WeylElt& KGB::involution(KGBElt z) const

/*!
  \brief Returns the root datum involution corresponding to z.

  In fact, returns the corresponding Weyl group element, s.t. w.delta = tau.

  NOTE: this is not inlined to avoid a dependency on weyl.h
*/

{
  return d_involution[z];
}

bool KGB::isAscent(size_t s, KGBElt x) const

/*!
  \brief Tells whether simple root \#s is an ascent generator for KGB
  element \#x.

  Explanation: ascent generators are the ones for which are either a complex
  ascent, or are imaginary noncompact.

  The dual information is recorded in the d_descent list.

  NOTE: this is not inlined to avoid a compiler dependency on bitset.h
*/

{
  using namespace gradings;

  return not isDescent(s,x) and not (status(s,x) == Status::ImaginaryCompact);
}

bool KGB::isDescent(size_t s, KGBElt x) const

/*!
  \brief Tells whether simple root \# s is a descent generator for KGB
  element \#x.

  Explanation: descent generators are the ones for which the twisted involution
  tau descends; i.e. either s does not commute with tau and s.tau.s has shorter
  length in the Weyl group, or it commutes, and s.tau has shorter length.

  This information is recorded in the d_descent list.

  NOTE: this is not inlined to avoid a compiler dependency on bitset.h
*/

{
  return d_descent[x].test(s);
}

KGBEltPair KGB::tauPacket(const weyl::WeylElt& w) const

/*!
  \brief Returns the range in which involution(z) = w.

  Precondition: the enumeration has been sorted by values of tau.
*/

{
  using namespace weyl;

  typedef std::pair<std::vector<WeylElt>::const_iterator,
    std::vector<WeylElt>::const_iterator> IP;

  InvolutionCompare comp(weylGroup());

  IP range = std::equal_range(d_involution.begin(),d_involution.end(),w,comp);

  KGBElt first = range.first - d_involution.begin();
  KGBElt last = range.second - d_involution.begin();

  return std::make_pair(first,last);
}

size_t KGB::weylLength(KGBElt z) const

/*!
  \brief Returns the length of involution(z) as a Weyl group element.

  NOTE: this is not inlined to avoid a dependency on weyl.h
*/

{
  return d_weylGroup->length(d_involution[z]);
}

/******** manipulators *******************************************************/
void KGB::fillBruhat()

/*!
  \brief Constructs the BruhatOrder.

  NOTE: may throw an Overflow error. Commit-or-rollback is guaranteed.

  NOTE: we use an auto_ptr for exception-safety. My first!
*/

{
  using namespace bruhat;
  using namespace set;

  if (d_state.test(BruhatConstructed)) // work was already done
    return;

  std::auto_ptr<BruhatOrder> bp;
  std::vector<SetEltList> hd;

  try {
    makeHasse(hd,*this);
    bp = std::auto_ptr<BruhatOrder>(new BruhatOrder(hd));
  }
  catch (...) {
    throw error::MemoryOverflow();
  }

  // commit
  delete d_bruhat;
  d_bruhat = bp.release();
  d_state.set(BruhatConstructed);

  return;
}

}

/*****************************************************************************

        Chapter II -- The Helper class

  ... explain here when it is stable ...

******************************************************************************/

// namespace {
namespace kgb {
  namespace helper {
  /*!
\brief

  */
Helper::Helper(realredgp::RealReductiveGroup& G)
  :KGB(G.semisimpleRank()),
   d_reflectionMap(rank()),
   d_reflectionMatrix(rank())

/*!
  \brief The helper constructor just initializes the lists for the first
  element. The filling process is handled by fill().

  Algorithm: from the datum of the real form (or more precisely, from the 
  corresponding fundamental grading) we recover the basic cocycle that 
  transforms the whole construction in a Tits group computation. 
*/

{  
  using namespace bitvector;
  using namespace gradings;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tits;
  using namespace tori;
  using namespace weyl;

  d_realGroup = &G;
  d_weylGroup = &G.weylGroup();

  // put in d_twistMatrix the mod 2 matrix of the distinguished involution
  mod2(d_twistMatrix,G.distinguished());
  d_twistMatrix.transpose();

  // put the matrices of the simple coreflections in d_reflectionMatrix
  // and the corresponding mod 2 maps in d_reflectionMap
  for (size_t s = 0; s < rank(); ++s) {
    RootNbr j = rootDatum().simpleRootNbr(s);
    rootDatum().rootReflection(d_reflectionMatrix[s],j);
    mod2(d_reflectionMap[s],d_reflectionMatrix[s]);
    // need to transpose to get the matrix on the coroot side
    d_reflectionMap[s].transpose();
  }

  // d_tits is a list of Tits group elements, that are constructed as 
  // representatives for our orbits. It is initialized with a single element.
  d_tits.push_back(TitsElt(G.rank()));

  // the map d_titsMap allows us to search for a Tits element to see if it is 
  // new. The associated integer is the index of the Tits element in d_tits
  d_titsMap.insert(std::make_pair(d_tits[0],static_cast<size_t>(0ul)));

  // the map d_fiberData associates to each root datum involution the 
  // corresponding representatives and projection. We need this to properly 
  // normalize the Tits elements
  const RealTorus& T = G.cartan(0).dualFiber().torus();
  const ComponentSubspace& sp = T.topology().subspace();
  d_fiberData.insert(std::make_pair(WeylElt(), FiberData(sp)));

  // write down the initial grading
  initGrading(d_baseGrading,G);

  // d_status[0] will be filled in while looking at cross-actions
  d_status.push_back(Status());

  // initialize descent, length and tau
  d_length.push_back(0);
  d_descent.push_back(Descent());
  d_involution.push_back(WeylElt());

  // initialize cross and cayley tables
  // the undefined value for Cayley is 0 (nothing Cayley-tranforms to 0)
  for (size_t j = 0; j < d_rank; ++j) {
    d_cross[j].push_back(UndefKGB);
    d_cayley[j].push_back(0);
    d_inverseCayley[j].push_back(std::make_pair(UndefKGB,UndefKGB));
  }
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/

Helper::WI Helper::cayleyAdd(weyl::WeylElt& w, size_t s)

/*!
  \brief Adds a new root datum involution by Cayley transform
  from the first argument, through simple root \#(second argument).

  Precondition: w corresponds to a known root datum involution, and s 
  twisted-commutes with w. 

  Algorithm: the new reflection is s_alpha.tau; the new normalizing
  subspace is computed from the definition. [DV 7/21/06: the old
  comment said "computed by adding m_alpha to the previous one."  This
  procedure would work for type I Cayley transforms but not for type
  II, and in any case it is not what the code does.]

  NOTE: this is the most expensive function in the construction, so we
  would like to avoid using it as much as possible
*/

{
  using namespace lattice;
  using namespace latticetypes;
  using namespace tits;
  using namespace tori;
  using namespace weyl;

  WeylElt sw = w;
  weylGroup().leftProd(sw,s);

  // make the lattice involution corresponding to sw
  LatticeMatrix q;
  identityMatrix(q,realGroup().rank());

  WeylWord ww;
  weylGroup().out(ww,sw);

  for (size_t j = 0; j < ww.size(); ++j) {
    q *= d_reflectionMatrix[ww[j]];
  }

  q *= realGroup().distinguished();

  // go over to [NEGATIVE, not inverse DV] inverse transpose
  q.transpose();
  q.negate();

  // make subspace
  WeightList minus;
  plusBasis(minus,q);
  ComponentList minus2;
  mod2(minus2,minus);

  FiberData fd(minus2,realGroup().rank());
  fd.conjugate(d_twistMatrix);

  d_fiberData.insert(std::make_pair(sw,fd));  
  saturate(sw);

  return d_fiberData.find(sw);
}

Helper::TI Helper::cayleyAdd(size_t from, size_t s)

/*!
  \brief Adds a new parameter by Cayley transform from d_tits[first
  argument] through simple root \#(second argument).
*/

{    
  using namespace gradings;
  using namespace tits;

  TitsElt a = d_tits[from];
  cayleyTransform(a,s);
  WI fdi = find(a.w());
  fdi->second.normalize(a);

  size_t to = d_tits.size(); // index of the new element
  d_tits.push_back(a);

  // increment the lists
  d_descent.push_back(Descent());
  d_length.push_back(d_length[from]+1); // length goes always up
  d_involution.push_back(a.w());
  d_status.push_back(Status());
  for (size_t j = 0; j < d_rank; ++j) {
    d_cross[j].push_back(UndefKGB);
    d_cayley[j].push_back(0);
    d_inverseCayley[j].push_back(std::make_pair(UndefKGB,UndefKGB));
  }

  return d_titsMap.insert(std::make_pair(a,to)).first;
}

Helper::WI Helper::crossAdd(weyl::WeylElt& w, size_t s)

/*!
  \brief Adds a new root datum involution through cross action.

  Precondition: w is a Weyl group element representing a known root
  datum involution; s is (the number of) a generator such that s.w.s
  is not in d_titsMap;

  Algorithm: we add a new element to d_titsMap, by conjugating the necessary
  data from those for w.
*/

{
  using namespace weyl;

  WeylElt sws = w;
  weylGroup().twistedConjugate(sws,s);

  FiberData fd = find(w)->second;
  fd.conjugate(reflectionMap(twist(s)));

  return d_fiberData.insert(std::make_pair(sws,fd)).first;
}

Helper::TI Helper::crossAdd(size_t from, size_t s)

/*!
  \brief Adds a new parameter.

  Precondition: from is the index of a known parameter; s is a generator
  such that s x from is new;
*/

{
  using namespace gradings;
  using namespace tits;

  TitsElt a = d_tits[from];
  twistedConjugate(a,s);
  WI fdi = find(a.w());
  fdi->second.normalize(a);

  size_t to = d_tits.size(); // index of the new element
  d_tits.push_back(a);

  // increment the lists
  d_descent.push_back(Descent());
  d_length.push_back(d_length[from]);
  d_involution.push_back(a.w());
  d_status.push_back(Status());
  for (size_t j = 0; j < d_rank; ++j) {
    d_cross[j].push_back(UndefKGB);
    d_cayley[j].push_back(0);
    d_inverseCayley[j].push_back(std::make_pair(UndefKGB,UndefKGB));
  }

  // adjust length
  if (a.w() != d_involution[from]) // length goes up
    ++d_length[to];

  return d_titsMap.insert(std::make_pair(a,to)).first;
}

void Helper::cayleyExtend(size_t from)

/*!
  \brief Tries to enlarge the parameter set by cayley transforms.

  Precondition: from is the index in d_tits of the parameter we are extending
  from.

  It is assumed that the Helper structure is always kept in the state where
  all downward links are already filled in.
*/

{    
  using namespace gradings;
  using namespace tits;
  using namespace weyl;

  TitsElt current = d_tits[from];

  // try to get new elements by cayley-action
  for (size_t s = 0; s < d_rank; ++s) {

    Status::Value v = d_status[from][s];
    if (v != Status::ImaginaryNoncompact) {
      d_cayley[s][from] = UndefKGB;
      continue;
    }

    // cayley-transform current by sigma_s
    TitsElt a = current;
    cayleyTransform(a,s);

    // find the root datum involution
    WI fdi = find(a.w());
    if (isNew(fdi)) // add a new twisted involution
      fdi = cayleyAdd(current.w(),s);

    // find the Tits element
    const FiberData& fd = fdi->second;
    fd.normalize(a);
    TI ti = find(a);
    if (isNew(ti)) // add a new Tits element
      ti = cayleyAdd(from,s);

    // add new cayley link
    d_cayley[s][from] = ti->second;

    // fill inverse cayley link
    KGBEltPair& ic = d_inverseCayley[s][ti->second];
    if (from < ic.first) {
      ic.second = ic.first;
      ic.first = from;
    } else {
      ic.second = from;
    }

    // update descent table
    d_descent[ti->second].set(s);

    // add status value
    d_status[ti->second].set(s,Status::Real);
  }

    return;
}

void Helper::crossExtend(size_t from)

/*!
  \brief Tries to enlarge the parameter set by cross-actions.

  Precondition: from is the index in d_tits of the parameter we are extending
  from.

  It is assumed that the Helper structure is always kept in the state where
  all downward links are already filled in.
*/

{    
  using namespace gradings;
  using namespace tits;
  using namespace weyl;

  TitsElt current = d_tits[from];

  // try to get new elements by cross-action
  for (size_t s = 0; s < d_rank; ++s) {
    // see if link was already filled
    if (d_cross[s][from] != UndefKGB)
      continue;

    // twist-conjugate current by sigma_s
    TitsElt a = current;
    twistedConjugate(a,s);

    // find the root datum involution
    WI fdi = find(a.w());
    if (isNew(fdi)) // add a new twisted involution
      fdi = crossAdd(current.w(),s);

    // find the Tits element
    const FiberData& fd = fdi->second;
    fd.normalize(a);
    TI ti = find(a);
    if (isNew(ti)) // add a new Tits element
      ti = crossAdd(from,s);

    // add new cross-links
    d_cross[s][from] = ti->second;
    d_cross[s][ti->second] = from;

    // update descent table
    if (d_length[from] < d_length[ti->second])
      d_descent[ti->second].set(s);

    // update status info
    setStatus(from,ti->second,s);
  }

    return;
}

void Helper::fill()

/*!
  \brief Constructs the full orbit set.

  Precondition: the object is in the initial state, the one it is put in by
  the call to its constructor;

  Algorithm: The idea is just to start out from the fundamental orbit,
  and then to saturate through Cayley transforms and cross actions.
*/

{
  for (size_t j = 0; j < d_tits.size(); ++j) {
    // these calls may enlarge d_tits;
    crossExtend(j);
    cayleyExtend(j);
  }

  return;

}

void Helper::saturate(const weyl::WeylElt& w)

/*!
  \brief Add all cross-transforms from from in the fiberDatum map.

  Explanation: the point is to avoid direct Cayley constructions, which are
  fairly expensive.
*/

{
  using namespace weyl;

  std::set<WeylElt> found;
  std::stack<WeylElt> toDo;

  found.insert(w);
  toDo.push(w);

  while (!toDo.empty()) {

    WeylElt v = toDo.top();
    toDo.pop();

    for (Generator s = 0; s < d_rank; ++s) {
      WeylElt v1 = v;
      weylGroup().twistedConjugate(v1,s);
      if (found.insert(v1).second)  { // found a new element
	toDo.push(v1);
	FiberData fd = find(v)->second;
	fd.conjugate(reflectionMap(twist(s)));
	d_fiberData.insert(std::make_pair(v1,fd));
      }
    }
  }

  return;
}

void Helper::setStatus(size_t from, size_t to, size_t s)

/*!
  \brief Fills in status information.

  Precondition: cross[s][from] = to;

  Explanation: if from and to have different root datum involutions, then
  s is complex both for from and to. Otherwise, s is imaginary for from,
  and real for cayley[s][from]. To find whether s is compact or noncompact,
  we note that the grading may be computed from the base grading and the
  torus part of the d_tits[from].
*/

{
  using namespace bitvector;
  using namespace gradings;

  if (d_status[from][s] == Status::Real)
    return;

  if (d_involution[from] != d_involution[to]) {
    d_status[from].set(s,Status::Complex);
    d_status[to].set(s,Status::Complex);
    return;
  }

  // now the root alpha_s is imaginary
  bool b = d_baseGrading.test(s);
  b ^= scalarProduct(d_tits[from].t(),titsGroup().simpleRoot(twist(s)));

  if (b) {
    d_status[from].set(s,Status::ImaginaryNoncompact);
    d_status[to].set(s,Status::ImaginaryNoncompact);
  }
  else // from == to
    d_status[from].set(s,Status::ImaginaryCompact);

  return;
}

void Helper::sort()

/*!
  \brief Sorts the entries in order of increasing length, and for a given
  length, in order of increasing root datum involution.

  Algorithm: make the sorting permutation, then apply it to the various tables.
*/

{
  using namespace setutils;

  Permutation a;
  identity(a,size());

  Compare comp(*this);
  std::sort(a.begin(),a.end(),comp);
  Permutation ai;
  invert(ai,a);

  // permute the cross and cayley maps
  for (size_t s = 0; s < d_rank; ++s) {
    for (size_t j = 0; j < size(); ++j) {
      d_cross[s][j] = ai[d_cross[s][j]];
      if (d_cayley[s][j] != UndefKGB)
	d_cayley[s][j] = ai[d_cayley[s][j]];
      if (d_inverseCayley[s][j].first != UndefKGB)
	d_inverseCayley[s][j].first = ai[d_inverseCayley[s][j].first];
      if (d_inverseCayley[s][j].second != UndefKGB)
	d_inverseCayley[s][j].second = ai[d_inverseCayley[s][j].second];
    }
    permute(ai,d_cross[s]);
    permute(ai,d_cayley[s]);
    permute(ai,d_inverseCayley[s]);
  }

  // permute the other lists
  permute(ai,d_descent);
  permute(ai,d_length);
  permute(ai,d_involution);
  permute(ai,d_status);

  return;
}

}
}

/*****************************************************************************

        Chapter III -- Functions local to kgb.cpp

  ... explain here when it is stable ...

******************************************************************************/

// namespace {

namespace kgb {
  namespace helper {

void initGrading(gradings::Grading& gr, const realredgp::RealReductiveGroup& G)

/*!
  \brief Flags in gr the the simple roots which are noncompact
  imaginary for the base grading of the fundamental Cartan in G..

Algorithm: the variable rset is first made to flag the noncompact
imaginary roots for the base grading of the fundamental Cartan in G.
Each simple root for G is tested for membership in rset; if it passes,
the corresponding bit of gr is set.
*/

{  
  using namespace rootdata;

  const RootDatum& rd = G.rootDatum();

  RootSet rset;
  G.grading(rset);

  for (size_t s = 0; s < G.semisimpleRank(); ++s) {
    RootNbr r = rd.simpleRootNbr(s);
    if (rset.isMember(r))
      gr.set(s);
  }

  return;
}

void makeHasse(std::vector<set::SetEltList>& hd, const KGB& kgb)

/*!
  \brief Puts in hd the hasse diagram data for the Bruhat ordering
  on kgb.
  
  Explanation: this is the closure ordering of orbits. We use the
  algorithm from Richardson and Springer.
*/

{
  using namespace gradings;
  using namespace poset;
  using namespace set;

  hd.resize(kgb.size());

  for (size_t x = 0; x < kgb.size(); ++x) {

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

  return;

}

}
}

}
