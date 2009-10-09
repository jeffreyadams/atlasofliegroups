/*!
\file
\brief Implementation of the class KGB representing orbits of K on G/B.

  This module contains code for the construction of a block in the
  one-sided parameter set (in other words, the subset of the one-sided
  parameter set corresponding to a single real form.) As explained in
  my Palo Alto III notes, this is equivalent to parameterizing the set
  K\\G/B of (K,B)-orbits in G; hence the provocative title.
*/
/*
  This is kgb.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007-2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "kgb.h"

#include <cassert>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>

#include <iostream>

#include "bitmap.h"
#include "bruhat.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "error.h"
#include "gradings.h"
#include "hashtable.h"
#include "latticetypes.h"
#include "realredgp.h"
#include "rootdata.h"
#include "tits.h"
#include "tori.h"
#include "weyl.h"

/*
  This module contains code for the construction of a block in the
  one-sided parameter set (in other words, the subset of the one-sided
  parameter set corresponding to a single real form.) As explained in
  my Palo Alto III notes, this is equivalent to parameterizing the set
  K\\G/B of (K,B)-orbits in G; hence the provocative title.

*/

namespace atlas {

  // namespace {
  namespace kgb {
    namespace {

/*
	                Local functions
*/

gradings::Grading
square_class_grading_offset(const cartanclass::Fiber& f,
			    cartanclass::square_class csc,
			    const rootdata::RootDatum& rd);
gradings::Grading
grading_offset_for(const realredgp::RealReductiveGroup& GR);

void makeHasse(std::vector<set::SetEltList>&, const KGB&);



/*
	           Some small additional classes
*/

  /*!
\brief A |FiberData| object associates to each twisted involution a subspace
describing how corresponding Tits elements should be normalized.

It also records the Cartan class that each twisted involution belongs to.
  */
class FiberData
{
  const tits::TitsGroup& Tits;
  weyl::TI_Entry::Pooltype pool;
  hashtable::HashTable<weyl::TI_Entry,unsigned int> hash_table;
  std::vector<latticetypes::SmallSubspace> data;
  std::vector<unsigned int> Cartan_class;
public:
  FiberData(complexredgp::ComplexReductiveGroup& G,
	    const bitmap::BitMap& Cartan_classes);

  void reduce(tits::TitsElt& a) const;

  size_t cartanClass(const weyl::TwistedInvolution& tw) const
  {
    return Cartan_class[hash_table.find(tw)];
  }

private: // the space actually stored need not be exposed
  const latticetypes::SmallSubspace& mod_space(const tits::TitsElt& a) const
  {
    size_t i=hash_table.find(a.tw());
    assert(i!=hash_table.empty);
    return data[i];
  }
}; // class FiberData




/*
                 The |KGBHelp| class, a helper class for |KGB|

 */

/* The helper class stores the basic data for KGB elements that will be
   exported to the |KGB| class, but also additional data that is used during
   generation but is not retained in the final representation; notably tables
   for and reducing Tits elements to standard form an looking them up. Much is
   stored in a table |d_fiberData| recording per-twisted-involution data. Upon
   exportation the numbering of the elements will be changed, and all
   internally used indices modified to reflect the renumbering.
 */


class KGBHelp
{
  const size_t d_rank;

  // data that will be exported (in permuted form) to the KGB class
  std::vector<KGBEltList> d_cross;
  std::vector<KGBEltList> d_cayley;
  std::vector<KGBInfo> d_info;

  // other data

  /*!\brief List of Tits elements parameterizing KGB orbits.

  Accessed usually via the hash table |hash|, using a |TE_Entry| as key
  */
  tits::TE_Entry::Pooltype Tits;
  hashtable::HashTable<tits::TE_Entry,KGBElt> hash;

  tits::BasedTitsGroup basePoint; // base point choice resides here

  //! Permits reducing each Tits group element modulo its fiber denominator
  FiberData d_fiberData;

 public:

// constructors and destructors
  KGBHelp(realredgp::RealReductiveGroup&);
  KGBHelp(complexredgp::ComplexReductiveGroup& G,
	  const bitmap::BitMap& Cartan_classes,
	  KGBElt size,
	  const tits::BasedTitsGroup& base,
	  const std::vector<tits::TitsElt>& seeds); // auxiliary constructor

  ~KGBHelp() {};

// public manipulators and accessor:

//!\brief fill the KGB set and return it
  KGBHelp& fill();

//!\brief deliver values to fields of a |KGB| object under construction.
  size_t export_tables(std::vector<KGBEltList>& cross,
		       std::vector<KGBEltList>& cayley,
		       tits::TitsEltList& tits,
		       std::vector<KGBInfo>& info,
		       tits::BasedTitsGroup*& base) const;

  // though used only from within |export_tables|, this method must be public
  bool comp(KGBElt x,KGBElt y) const
  {
    if (d_info[x].length!=d_info[y].length)
      return d_info[x].length<d_info[y].length;
    unsigned long lx=weylGroup().length(Tits[x].w());
    unsigned long ly=weylGroup().length(Tits[y].w());
    return lx!=ly ? lx<ly : Tits[x].w()<Tits[y].w();
  }

// accessors (private)
private:

  const tits::TitsGroup& titsGroup() const {
    return basePoint.titsGroup();
  }

  const weyl::WeylGroup& weylGroup() const {
    return titsGroup().weylGroup();
  }


// private manipulators
  void cross_extend(KGBElt parent);
  void cayleyExtend(KGBElt parent);

}; // class KGBHelp

// A non-method construction function
KGBHelp refined_helper(realredgp::RealReductiveGroup& GR,
		       const bitmap::BitMap& Cartan_classes);


/* The following auxiliary class provides a comparison object for calling
   standard search and sorting routines for twisted involutions. It serves as
   specification of the comparison that is used in sorting the KGB structure,
   but in fact it is so expensive for repeated use (since |W.involutionLength|
   has to recompute its result each time) that its use has been discontinued.
   In Fokko's code it was used by |tauPacket|, and formed the main bottleneck
   for the block construction; now the hash table |d_tau| together with the
   table |first_of_tau| provide a much faster way to implement |tauPacket|.
*/
class InvolutionCompare {
private:
  const weyl::TwistedWeylGroup& tW;
public:
  explicit InvolutionCompare(const weyl::TwistedWeylGroup& W) : tW(W) {}

  // one should have a < b iff
  // (a) involutionLength(a) < involutionLength(b) or
  // (b) involutionLengths are equal and length(a) < length (b) or
  // (c) both lengths are equal and a < b
  bool operator()
   (const weyl::TwistedInvolution& a, const weyl::TwistedInvolution& b) const
  {
    const weyl::WeylGroup& W=tW.weylGroup();
    if      (tW.involutionLength(a) != tW.involutionLength(b))
      return tW.involutionLength(a) <  tW.involutionLength(b) ;
    else if (W.length(a.w()) != W.length(b.w()))
      return W.length(a.w()) <  W.length(b.w());
    else
      return a < b;
  }
}; // class InvolutionCompare

/* For sorting the KGB elements on exportation from the helper class, we use
   the |comp| method of the helper class, which effectively applies the
   comparison of |InvolutionCompare| to the Weyl group parts of the stored
   Tits group elements, but using their stored length value rather than
   calling |WeylGroup::involutionLength|, for much improved speed. The class
   below serves to wrap that |comp| method into a proper comparison object.
*/
class IndexCompare
{
private:
  const KGBHelp& kgb;
public:
  explicit IndexCompare(const KGBHelp& h) : kgb(h) {}

  // here |unsigned long| arguments, as |setutils::Permutation| contains such
  bool operator() (unsigned long i, unsigned long j) const {
    return kgb.comp(i,j);
  }
}; // class IndexCompare

} // namespace
} // namespace kgb


/*****************************************************************************

        Chapter I -- The KGB class, public methods

******************************************************************************/

namespace kgb {

/*!
  \brief Construct the KGB data structure for the given real form,
  but if any Cartan classes are specified, only for those classes

  This really handles two cases (the standard construction and a specialized
  one for a small set of Cartan classes) in one. The reason that they are
  combined into a single constructor is that this allows a single constructor
  of a containing class to choose between the two methods (a constructor
  method must use a fixed constructor for each of its sub-objects, so having
  multiple constructors here would force the same for every containing class).
*/
KGB::KGB(realredgp::RealReductiveGroup& GR,
	 const bitmap::BitMap& Cartan_classes)
  : d_rank(GR.semisimpleRank())
  , d_torus_rank(GR.rank())
  , d_cross()
  , d_cayley()
  , d_inverseCayley()
  , d_info()
  , d_tits()
  , first_of_tau()
  , d_pool(), d_tau(d_pool)
  , d_bruhat(NULL)
  , d_base(NULL)
{
  bool traditional= Cartan_classes.size()==0; // whether traditional generation
  size_t size =
    (traditional ? KGBHelp(GR) : refined_helper(GR,Cartan_classes))
    .fill()
    .export_tables(d_cross,d_cayley,d_tits,d_info,d_base);

  // check that the size is right
  assert(size == (traditional ? GR.KGB_size()
		  : GR.complexGroup().KGB_size(GR.realForm(),Cartan_classes)
		  ));

  // reserve one more than the number of twisted involutions present
  first_of_tau.reserve((traditional ? GR.numInvolutions()
			:GR.complexGroup().numInvolutions(Cartan_classes)
			)+1);

  /* Since elements are sorted by twisted involution, collecting those is easy.
     By using the hash table |d_tau|, it gets initialized in the proper order
  */
  { // insert twisted involutions in order into hash table |d_tau|
    for (KGBElt x=0; x<size; ++x)
    {
      unsigned int old = d_tau.size();
      if (d_tau.match(involution(x))==old) // a new twisted involution
	first_of_tau.push_back(x); // record first |x| for each |tau|
    }
    first_of_tau.push_back(size); // put final sentinel
  }

  // Finally compute inverse Cayley tables

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

/******** copy, assignment and swap ******************************************/


/******** accessors **********************************************************/


/*!
  \brief Returns the range of |KGBElt| values |x| with |involution(x)==w|.

  It is possible to return a range of |KGBElt|, since these elements have
  been sorted by considering the associated twisted involution first.

  The fields |d_tau| and |first_of_tau| were added to |KGB| in order to make
  this method, central to the block construction, as fast as possible.
*/
KGBEltPair KGB::tauPacket(const weyl::TwistedInvolution& w) const
{
  unsigned int i=d_tau.find(w);
  if (i==d_tau.empty)
    return KGBEltPair(0,0);
  return KGBEltPair(first_of_tau[i],first_of_tau[i+1]);
}


/******** manipulators *******************************************************/

/*!
  \brief Constructs the BruhatOrder.

  NOTE: may throw a MemoryOverflow error. Commit-or-rollback is guaranteed.
*/
void KGB::fillBruhat()

{
  if (d_state.test(BruhatConstructed)) // work was already done
    return;

  try {
    std::vector<set::SetEltList> hd; makeHasse(hd,*this);
    bruhat::BruhatOrder* bp = new bruhat::BruhatOrder(hd); // may throw here

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

namespace kgb {

/*    II a. |FiberData|  */

namespace { // |FiberData| and |KGBHelp| are in anonymous namespace

/*
  The |FiberData| constructor computes a subspace for each twisted involution
  $tw$ (representing an involution $\tau$ of $H$ and an involution
  $q=tw.\delta$ of the weight lattice $X$) that occurs for $GR$. The subspace
  of $X^\vee/2X^\vee$ that will be stored for $tw$ is the image $I$ of the
  $-1$ eigenspace $X^\vee_-$ of $q^t$. When a Tits group element of the form
  $x.\sigma_w$ occurs in the KGB construction, only the coset of the left torus
  part $x$ modulo $I$ matters, and it will after computation be systematically
  normalized by reducing the left torus part $x$ modulo $I$.

  The image $I$ is what one divides by to get the fiber group of the real
  Cartan associated to $H$ and $\tau$, in which case the "numerator" is the
  image of $X^\vee_+ +X^\vee_-$. The possible values of the torus part $x$ are
  such that there are as much values modulo $I$ as there are in the fiber
  group, but they are not naturally in bijection, and $0$ need no be a valid
  value for $x$: the possible values are such that the Tits group element $a$
  given by $(x,w)$ satisfies $a*twisted(a)=e$, where the twisting is done
  with respect to the based Tits group (determined by the base grading).

  This constructor depends on the real form only via the set |Cartan_classes|
  that determines (limits) the set of twisted involutions to be considered.
*/
FiberData::FiberData(complexredgp::ComplexReductiveGroup& G,
		     const bitmap::BitMap& Cartan_classes)
  : Tits(G.titsGroup())
  , pool()
  , hash_table(pool)
  , data()
  , Cartan_class()
{
  { // dimension |data| and |Cartan_class|
    size_t n_inv=G.numInvolutions(Cartan_classes);
    data.reserve(n_inv); Cartan_class.reserve(n_inv);
  }

  const rootdata::RootDatum& rd=G.rootDatum();

   std::vector<latticetypes::BinaryMap> refl(G.semisimpleRank());
   for (weyl::Generator s=0; s<refl.size(); ++s)
    // get endomorphism of weight lattice $X$ given by generator $s$
    // reflection map is induced vector space endomorphism of $X^* / 2X^*$
     refl[s] = latticetypes::BinaryMap(rd.simple_reflection(s).transposed());

   for (bitmap::BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
   {
     size_t cn=*it;
     size_t i = hash_table.match(G.twistedInvolution(cn));
     assert(i==data.size()); // this twisted involution should be new

     { // store data for canonical twisted involution |i| of Cartan class |cn|
       latticetypes::LatticeMatrix q= G.cartan(cn).involution();
       data.push_back(tits::fiber_denom(q)); // compute subspace $I$
       Cartan_class.push_back(cn); // record number of Cartan class
     }


     // now generate all non-canonical twisted involutions for Cartan class
     for ( ; i<data.size(); ++i) // |data.size()|  increases during the loop
       for (size_t s=0; s<G.semisimpleRank(); ++s)
       {
	 weyl::TwistedInvolution stw =
	   G.twistedWeylGroup().twistedConjugated(pool[i],s);
	if (hash_table.match(stw)==data.size()) // then |stw| is new
	{
	  data.push_back(data[i]);     // start with copy of source subspace
	  data.back().apply(refl[s]);  // modify according to cross action used
	  Cartan_class.push_back(cn);  // record number of Cartan class
 	}
      }
  }

  // check that the number of generated elements is as predicted
  assert(data.size()==G.numInvolutions(Cartan_classes));
  assert(Cartan_class.size()==G.numInvolutions(Cartan_classes));
}


void FiberData::reduce(tits::TitsElt& a) const
{
  Tits.left_torus_reduce(a,mod_space(a));
}




/*    II b. The main helper class |KGBHelp|  */

/*
   The actual KGB construction takes place below. During the construction, the
   elements are represented as Tits group elements. The links in the KGB
   structure are realized by |BasedTitsGroup::basedTwistedConjugate| for the
   cross actions and by |basedTitsGroup::Cayley_transform| for Cayley
   transforms. After each of these, the result is subject to
   |d_fiberdata.reduce| to normalize the representation of a KGB element.
*/

/*! \brief
  This constructor sets |basePoint| for |GR|, and a trivial initial value.

  So we choose a different grading offset for each real form, but always start
  at the same Tits element. For an approach where different real forms can
  share a grading offset and thus make their KGB sets mesh together, see the
  next constructor. About the current (more or less) constructor Fokko said:
  ``from the datum of the real form (or more precisely, from the corresponding
  fundamental grading) we recover the basic cocycle that transforms the whole
  construction into a Tits group computation''. By the basic (1-)cocycle he
  seems to have meant the map from \f$W\f$ to the \f$W\f$-module \f$Y/2Y\f$
  defined by the mentioned adjoint fiber element as image of the identity, and
  extended by the natural grading shifts for the simple reflections (only the
  image of the identity was actually computed).
*/

KGBHelp::KGBHelp(realredgp::RealReductiveGroup& GR)
  : d_rank(GR.semisimpleRank())

  , d_cross(d_rank)
  , d_cayley(d_rank)
  , d_info()

  , Tits(), hash(Tits)

  // only the final fields depend on the real form of |GR|
  , basePoint(GR.complexGroup(),GR.grading_offset())
  , d_fiberData(GR.complexGroup(),GR.Cartan_set())
{
  KGBElt size = GR.KGB_size();

  Tits.reserve(size);
  d_info.reserve(size);

  // set up cross and Cayley tables with undefined values
  for (size_t j = 0; j < d_rank; ++j)
  {
    d_cross[j].resize(size,UndefKGB);
    d_cayley[j].resize(size,0); // leave unset (set by |cayleyExtend|)
  }

  // set identity Tits element as seed of the KGB construction
  hash.match(tits::TitsElt(titsGroup()));

  // store its length and Cartan class (the latter is in fact always 0)
  d_info.push_back(KGBInfo(0,d_fiberData.cartanClass(hash[0].tw())));
}

/* The following constructor serves the quasi-constructor |refined_helper|.
   Its main purpose is to allow pre-computation of the |basePoint| field as
   the base object of an |EnrichedTitsGroup| object, and a set of |seeds|, one
   for each minimal Cartan class. The other arguments serve mainly to transmit
   values that are computed in |refined_helper| anyway, allowing them to be
   used in the construction without recomputation.
*/
KGBHelp::KGBHelp(complexredgp::ComplexReductiveGroup& G,
		 const bitmap::BitMap& Cartan_classes,
		 KGBElt size, // size of KGB for selected |Cartan_classes|
		 const tits::BasedTitsGroup& base, // base point information)
		 const std::vector<tits::TitsElt>& seeds)
  : d_rank(G.semisimpleRank())

  , d_cross(d_rank)
  , d_cayley(d_rank)
  , d_info()

  , Tits() // no elements a priori, those from |seeds| are added below
  , hash(Tits)

  // only the final fields depend on the real form of |GR|
  , basePoint(base)
  , d_fiberData(G,Cartan_classes)
{
  Tits.reserve(size);
  d_info.reserve(size);

  // set up cross and Cayley tables with undefined values
  for (size_t j = 0; j < d_rank; ++j)
  {
    d_cross[j].resize(size,UndefKGB);
    d_cayley[j].resize(size,0); // leave undefined (set by |cayleyExtend|)
  }

  set::SetEltList m=G.Cartan_ordering().minima(Cartan_classes);
  assert(m.size()==seeds.size()); // one seed for every minimal Cartan class

  for (size_t i=0; i<seeds.size(); ++i)
  {
    tits::TitsElt a= seeds[i];
    d_fiberData.reduce(a);

    // now add KGB element for the reduced Tits group element
#ifndef NDEBUG
    size_t k=hash.match(a);
    assert(k==d_info.size()); // this KGB element should be new
#else
    hash.match(a); // enter new KGB element into the tables
#endif

    // add additional information (length,Cartan class) for this KGB element
    d_info.push_back(KGBInfo(G.twistedWeylGroup().involutionLength(a.tw()),
			     m[i]));
  }
}


/*!
  \brief
  This quasi-constructor builds a |KGBHelp| object initialized with an
  element for each minimal Cartan class, and base point for the square class.

  The base grading is set up to correspond to (the chosen adjoint fiber
  element in the fundamental Cartan for) the central square class of this
  real form, which is done by the call to |square_class_grading_offset|.

  The initial element then represents the place within its central square class
  of one chosen strong real form |srf| lying over this weak real form; it has
  the identity twisted involution, and a torus factor obtained by lifting the
  representative fiber group element |srf.first| via the |fromBasis| method of
  the fiber group back to the ``coweight lattice modulo 2'' \f$Y/2Y\f$.

  Here we actually look up the strong real form in order to get a proper
  initial Tits group element associated to this Cartan class
*/
KGBHelp refined_helper(realredgp::RealReductiveGroup& GR,
		       const bitmap::BitMap& Cartan_classes)
{
  tits::EnrichedTitsGroup
    square_class_base(tits::EnrichedTitsGroup::for_square_class(GR));

  complexredgp::ComplexReductiveGroup& G=GR.complexGroup();
  realform::RealForm rf=GR.realForm();
  assert(square_class_base.square()==
	 G.fundamental().central_square_class(rf));

  KGBElt size =  G.KGB_size(rf,Cartan_classes);

  set::SetEltList m=G.Cartan_ordering().minima(Cartan_classes);

  std::vector<tits::TitsElt> seeds; seeds.reserve(m.size());
  for (size_t i=0; i<m.size(); ++i)
    seeds.push_back
      (m.size()==1 // use backtrack only in case of multiple minimal classes
       ? square_class_base.grading_seed(G,rf,m[i])
       : square_class_base.backtrack_seed(G,rf,m[i])
      );

  return KGBHelp(G,Cartan_classes,size,square_class_base,seeds);
}

/******** public methods ****************************************************/

/*!
  \brief Constructs the full orbit set.

  Precondition: the object is in the initial state, the one it is put in by
  the call to its constructor (at least one seed element is present).

  Algorithm: The idea is just to start out from the given element,
  and then to saturate through cross actions and Cayley transforms.

  It is important that for each element the cross actions are defined before
  the Cayley transforms is tried, because the status information set by the
  former is used by the latter.
*/
KGBHelp& KGBHelp::fill()
{
  for (KGBElt j = 0; j < Tits.size(); ++j) {
    // these calls will usually enlarge |Tits.size()|;
    cross_extend(j);
    cayleyExtend(j);
  }
  return *this;
}

size_t KGBHelp::export_tables(std::vector<KGBEltList>& cross,
			      std::vector<KGBEltList>& cayley,
			      tits::TitsEltList& Tits_out,
			      std::vector<KGBInfo>& info,
			      tits::BasedTitsGroup*& base) const
{
  size_t size=d_info.size();

  // sort (partially)
  setutils::Permutation a(size,1);
  IndexCompare comp(*this);
  std::stable_sort(a.begin(),a.end(),comp); // better and faster than |sort|

  // export the cross and Cayley maps, permuting each constituent list
  cross.clear(); cayley.clear();
  cross.reserve(d_rank); cayley.reserve(d_rank);
  setutils::Permutation ai(a,-1); // compute inverse of |a|
  for (size_t s = 0; s < d_rank; ++s) {
    cross.push_back(a.pull_back(ai.renumbering(d_cross[s])));
    cayley.push_back(a.pull_back(ai.renumbering(d_cayley[s],UndefKGB)));
  }

  // export the other lists
  /* We cannot do |a.pull_back(Tits).swap(Tits_out);|, as we need to convert
     a vector of |TE_Entry| to a vector of |TitsElt|; a loop is required */
  Tits_out.clear(); Tits_out.reserve(size);
  for (KGBElt x=0; x<size; ++x)
    Tits_out.push_back(Tits[a[x]]);

  // but here there is no difficulty
  a.pull_back(d_info).swap(info); // pull back |d_info| through |a|, and export

  base=new tits::BasedTitsGroup(basePoint); // export base point information

  return size;
}

/******** manipulators *******************************************************/

/*!
  \brief Tries to enlarge the parameter set by cross-actions, without
  supposing that elements are generated by increasing order of twisted
  involution length

  Precondition: |parent| is the index into |Tits| of the parameter we are
  extending from.
*/
void KGBHelp::cross_extend(KGBElt parent)
{
  const tits::TitsElt current = Tits[parent];

  // try to get new elements by cross-action
  for (size_t s = 0; s < d_rank; ++s) {
    // see if link was already filled
    if (d_cross[s][parent] != UndefKGB)
      continue;

    // twisted-conjugate (using base grading) |current| by |s|
    tits::TitsElt a = current; basePoint.basedTwistedConjugate(a,s);

    /* Check that this operation is its own inverse
    {
      tits::TitsElt b=a; basePoint.basedTwistedConjugate(b,s);
      d_fiberData.reduce(b);
      assert(b==current);
    }
    */

    /* Check that result is independent of representative:
    {
      latticetypes::SmallBitVectorList b =
	d_fiberData.mod_space(current).basis();
      for (size_t i=0; i<b.size(); ++i)
      {
	TitsElt ai=current; ai+=b[i]; basePoint.basedTwistedConjugate(ai,s);
	assert(ai.tw()==a.tw());
	d_fiberData.reduce(ai);
	assert(ai.t()==a.t());
      }
    }
    */

    // find the Tits element
    d_fiberData.reduce(a);

    /* test whether action was (more or less) as computed in cartanclass.cpp
    if (a.tw()!=current.tw())
    { weyl::Generator t = titsGroup().twisted(s);
      tits::TorusPart x= current.t();
      titsGroup().reflect(x,t);
      if (not grading_offset.test(s))
	x += titsGroup().simpleCoroot(t);
      tits::TitsElt b(weylGroup().twistedConjugated(current.w(),s),x);
      d_fiberData.reduce(b);
      assert(a==b);
    }
    */

    // involution length changes by 0 for twisted commutation, or else by +/-1
    int lc=
      a.tw()==current.tw() ? 0 : weylGroup().length_change(s,current.w());

    // now find the Tits element
    KGBElt child = hash.match(a);
    if (child==d_info.size()) // add a new Tits element
      d_info.push_back(KGBInfo(d_info[parent].length+lc,
			       d_fiberData.cartanClass(a.tw())
			       ));

    // set cross links for |parent| and for |child| (whether new or old)
    d_cross[s][parent] = child;
    d_cross[s][child] = parent;

    if (lc!=0) // then set complex status, and whether ascent or descent
    {
      d_info[parent].status.set(s,gradings::Status::Complex);
      d_info[child].status.set(s,gradings::Status::Complex);
      d_info[parent].desc.set(s,lc<0);
      d_info[child].desc.set(s,lc>0);
    }
    else if (weylGroup().hasDescent(s,current.w())) // real
    {
      d_info[parent].status.set(s,gradings::Status::Real); // |child==parent|
      d_info[parent].desc.set(s,true); // real roots are always descents
    }
    else// imaginary
    {
      d_info[parent].status.set_imaginary
	(s,basePoint.simple_grading(current,s));
      d_info[parent].desc.set(s,false); // imaginary roots are never descents
      if (child!=parent) // give child same status (which will be noncompact)
      {
	d_info[child].status.set(s,d_info[parent].status[s]);
	d_info[child].desc.set(s,false);
      }
    }
  }
}


/*!
  \brief Tries to enlarge the parameter set by Cayley transforms from |parent|.

  Precondition: |parent| is a |KGBElt| whose |status| field in |d_info| has
  been set to the proper value.

  Calling this function also sets the fields |d_cayley[s][parent]| either to
  the appropriate value or to |UndefKGB| (it was 0, which is neither of those)

  It is assumed that the it is an invariant of the |KGBHelp| structure that
  all downward links are already filled in.
*/
void KGBHelp::cayleyExtend(KGBElt parent)
{
  const tits::TitsElt current = Tits[parent];

  // try to get new elements by Cayley transform
  for (size_t s = 0; s < d_rank; ++s) {

    gradings::Status::Value v = d_info[parent].status[s];
    if (v != gradings::Status::ImaginaryNoncompact) {
      d_cayley[s][parent] = UndefKGB;
      continue;
    }

    // Cayley-transform |current| by $\sigma_s$
    tits::TitsElt a = current; basePoint.Cayley_transform(a,s);
    assert(titsGroup().length(a)>titsGroup().length(current)); // should go up


    // now look up the correspondingly reduced Tits element
    d_fiberData.reduce(a); // subspace has grown, so mod out new subspace
    KGBElt x = hash.match(a);
    if (x==d_info.size()) // add a new Tits element
      d_info.push_back(KGBInfo(d_info[parent].length+1, // length goes up
			       d_fiberData.cartanClass(a.tw())
			       ));
    // add new Cayley link
    d_cayley[s][parent] = x;
  }
}



} // namespace
} // namespace kgb

/*****************************************************************************

        Chapter III -- Functions local to kgb.cpp

******************************************************************************/

namespace kgb {
  namespace {

/*!
  \brief Puts in |Hasse| the Hasse diagram of the Bruhat ordering on |kgb|.

  Explanation: this is the closure ordering of orbits. We use the algorithm
  from Richardson and Springer.
*/
void makeHasse(std::vector<set::SetEltList>& Hasse, const KGB& kgb)
{
  Hasse.resize(kgb.size());

  for (KGBElt x = 0; x < kgb.size(); ++x)
  {
    std::set<KGBElt> h_x;
    const kgb::Descent& d = kgb.descent(x);
    if (d.none()) // element is minimal in Bruhat order
      continue;

    size_t s = d.firstBit();
    KGBElt sx;

    if (kgb.status(s,x) == gradings::Status::Complex)
      sx = kgb.cross(s,x);
    else
    { // |s| is real for |x|
      KGBEltPair sxp = kgb.inverseCayley(s,x);
      sx = sxp.first; // and will be inserted below
      if (sxp.second != UndefKGB) // |s| is real type I for |x|
	h_x.insert(sxp.second);
    }
    h_x.insert(sx);

    for (size_t j = 0; j < Hasse[sx].size(); ++j) {
      KGBElt z = Hasse[sx][j];
      if (kgb.isAscent(s,z)) {
	if (kgb.status(s,z) == gradings::Status::ImaginaryNoncompact)
	  h_x.insert(kgb.cayley(s,z));
	else // s is complex for z
	  h_x.insert(kgb.cross(s,z));
      }
    }

    std::copy(h_x.begin(),h_x.end(),std::back_inserter(Hasse[x]));
  } // for |x|

}

} // namespace
} // namespace kgb

} // namespace atlas
