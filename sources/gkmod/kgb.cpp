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
  Copyright (C) 2007 Marc van Leeuwen
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
#include "cartanset.h"
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
  FiberData(const complexredgp::ComplexReductiveGroup& G,
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
   generation but is not retained in the final representation; notably for
   each KGB element the Tits group element with which it is associated (the
   correspondence between the two depends on fairly arbitrary choices, whence
   retaining the Tits group element after generation is not really useful).
   Also stored is a table |d_fiberData| recording per-twisted-involution data.
   Upon exportation the numbering of the elements will be changed, and all
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

  /*!\brief List of Tits elements parametrizing KGB orbits.

  Accessed usually via the hash table |d_tits| (and |d_tits[i]| is |d_pool[i]|)
  */
  tits::TE_Entry::Pooltype d_pool;
  hashtable::HashTable<tits::TE_Entry,KGBElt> d_tits;

  BasedTitsGroup basePoint; // owned reference, base point choice resides here

  //! Permits reducing each Tits group element modulo its fiber denominator
  FiberData d_fiberData;

 public:

// constructors and destructors
  KGBHelp(realredgp::RealReductiveGroup&);
  KGBHelp(const complexredgp::ComplexReductiveGroup& G,
	  const bitmap::BitMap& Cartan_classes,
	  KGBElt size,
	  const BasedTitsGroup& base,
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
		       BasedTitsGroup*& base,
		       bool traditional) const;

  // though used only from within |export_tables|, this method must be public
  bool comp(KGBElt x,KGBElt y) const
  {
    if (d_info[x].length!=d_info[y].length)
      return d_info[x].length<d_info[y].length;
    unsigned long lx=weylGroup().length(d_pool[x].w());
    unsigned long ly=weylGroup().length(d_pool[y].w());
    return lx!=ly ? lx<ly : d_pool[x].w()<d_pool[y].w();
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
   for the block construction; now the hash table |d_tau| togther with the
   table |first_of_tau| provide a much faster way to implement |tauPacket|.
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

  This really handles two cases (the standard construction and a specialised
  on for a small set of Cartan classes) in one. The reason that they are
  combined into a single constructor is that this allows a single constructor
  of a containing class to choose between the two methods (a constructor
  method must use a fixed constructor for each of its subobjects, so having
  multiple constructors here would force the same for every containing class).
*/
KGB::KGB(realredgp::RealReductiveGroup& GR,
	 const bitmap::BitMap& Cartan_classes)
  : d_rank(GR.semisimpleRank())
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
    .export_tables(d_cross,d_cayley,d_tits,d_info,d_base,traditional);

  // check that the size is right
  assert(size == (traditional ? GR.kgbSize()
		  : GR.complexGroup().cartanClasses().KGB_size
		     (GR.realForm(),Cartan_classes)
		  ));

  // reserve one more that the number of twisted involutions present
  first_of_tau.reserve((traditional ? GR.numInvolutions()
			:GR.complexGroup().numInvolutions(Cartan_classes)
			)+1);

  /* Since elements are sorted by twisted involution, collecting those is easy.
     By using the hash table |d_tau|, it gets initialised in the proper order
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
  \brief Returns the range of |KGBElt| values |z| with |involution(z)==w|.

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

namespace kgb {
  namespace {


/*    II a. |FiberData|  */

/*
  The |FiberData| constructor computes a subspace for each twisted involution
  $tw$ (representing an involution $\tau$ of $H$ and an involution
  $q=tw.\delta$ of the weight lattice $X$) that occurs for $GR$. The subspace
  of $X^\vee/2X^\vee$ that will be stored for $tw$ is the image $I$ of the
  $-1$ eigenspace $X^\vee_-$ of $q^t$. When a Tits group element of the form
  $x.\sigma_w$ occurs in the KGB construction, only the coset the left torus
  part $x$ modulo $I$ matters, and it will after computation be systematically
  normalised by reducing the left torus part $x$ modulo $I$.

  The image $I$ corresponds to the subset of elements of order 2 in $H$ that
  are of the form $h*\tau(h)^{-1}$. It is also what one divides by to get the
  fiber group of the real Cartan associated to $H$ and $\tau$. Note however
  that it is not true that such $t\in X^\vee/2X^\vee$ will always lie in the
  image of $X^\vee_+ + X^\vee_-$ that is used to form the fiber group.

  The Cartan class associated to $tw$ is also recorded.

  This constructor depends on the real form only via the set |Cartan_classes|
  that determines (limits) the set of twisted involutions to be considered.
*/
FiberData::FiberData(const complexredgp::ComplexReductiveGroup& G,
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
  latticetypes::BinaryMap delta(G.distinguished().transposed());

  std::vector<latticetypes::BinaryMap> refl(G.semisimpleRank());
  for (weyl::Generator s=0; s<refl.size(); ++s)
  {
   // get endomorphism of weight lattice $X$ given by generator $s$
    latticetypes::LatticeMatrix r = rd.rootReflection(rd.simpleRootNbr(s));

    // reflection map is induced vector space endomorphism of $X^* / 2X^*$
    refl[s] = latticetypes::BinaryMap(r.transposed());
  }

  for (bitmap::BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
  {
    size_t cn=*it;
    size_t i = hash_table.match(G.twistedInvolution(cn));
    assert(i==data.size()); // this twisted involution should be new

    { // store data for canonical twisted involution |i| of Cartan class |cn|
      using namespace latticetypes;
      LatticeMatrix qtr= G.cartan(cn).involution().transposed();
      data.push_back(SmallSubspace(SmallBitVectorList(tori::minusBasis(qtr)),
				   G.rank())); // compute subspace $I$
      Cartan_class.push_back(cn); // record number of Cartan class
    }


    // now generate all non-canonical twisted involutions for this Cartan class
    for ( ; i<data.size(); ++i) // |data.size()|  increases during the loop
      for (size_t s=0; s<G.semisimpleRank(); ++s)
      {
	weyl::TwistedInvolution stw=G.weylGroup().twistedConjugated(pool[i],s);
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
  a=tits::TitsElt(Tits,mod_space(a).mod_image(Tits.left_torus_part(a)),a.w());
}



/*    II b. The classes |BasedTitsGroup| and |EnrichedTitsGroup| */

BasedTitsGroup::BasedTitsGroup(realredgp::RealReductiveGroup& GR)
  : Tg(GR.titsGroup())
  , rd(GR.complexGroup().rootDatum())
  , grading_offset(grading_offset_for(GR))
{}

tits::TitsElt BasedTitsGroup::twisted(const tits::TitsElt& a) const
{
  tits::TitsElt result(Tg,Tg.twisted(Tg.left_torus_part(a)));
  weyl::WeylGroup W=Tg.weylGroup();
  weyl::WeylWord ww=W.word(a.w());
  for (size_t i=0; i<ww.size(); ++i)
  {
    weyl::Generator s=W.twisted(ww[i]);
    Tg.mult_sigma(result,s);
    if (grading_offset[s]) // $s$ noncompact imaginary: add $m_\alpha$
      Tg.right_add(result,Tg.simpleCoroot(s));
  }
  return result;
}

EnrichedTitsGroup::EnrichedTitsGroup(realredgp::RealReductiveGroup& GR,
				     const cartanclass::Fiber& fund)
  : BasedTitsGroup(GR)
  , srf(fund.strongRepresentative(GR.realForm()))
  , base_compact(fund.compactRoots(fund.class_base(srf.second)))
{
  // in this derived class, |grading_offset| depends only on the square class
  grading_offset=square_class_grading_offset(fund,srf.second,rd);
}

/* The reasoning given for |simple_grading| fails for non-simple roots, so we
   cannot easily compute the grading at arbitrary imaginary roots even in the
   fundamental fiber, using only |grading_offset|. However, |d_base_compact|
   knows the full grading by $x_0.\delta$; we can still use a scalar product.
*/
bool EnrichedTitsGroup::is_compact(const tits::TorusPart& x,
				   rootdata::RootNbr n) const
{
  latticetypes::SmallBitVector rn(rd.root(n));
  return base_compact.isMember(n) ^ bitvector::scalarProduct(x,rn);
}

/* To compute the grading at arbitrary imaginary roots at an arbitrary fiber,
   no simple closed formula seems avaiable. So we revert to appying a set of
   cross actions to convert the question into one about a simple root.
 */
bool EnrichedTitsGroup::grading(tits::TitsElt a, rootdata::RootNbr n) const
{
  if (not rd.isPosRoot(n))
    n=rd.rootMinus(n);

  assert(rd.isPosRoot(n));
  size_t i; // declare outside loop to allow inspection of final value
  do
    for (i=0; i<rd.semisimpleRank(); ++i)
      if (rd.scalarProduct(rd.simpleRootNbr(i),n)>0)
	if (n==rd.simpleRootNbr(i))
	  return simple_grading(a,i);
	else
	{
	  rd.rootReflect(n,i);
	  basedTwistedConjugate(a,i);
	  break;
	}
  while (i!=rd.semisimpleRank()); // i.e., until no change occurs any more

  assert(false); // root |n| cannot become negative without becoming simple
  return false;
}

/* The methods below can be used to compute an initial Tits element for a
   given real form |rf| at a given Cartan class |cn|. They employ various
   methods, and have various degrees of success, as explained below.
*/

/* The method |naive_seed| attempts to get an initial Tits group element by
   extracting the necessary information fom the |Fiber| object associated to
   the Cartan class |cn|, and lifting that information from the level of the
   fiber group back to the level of torus parts. But as the name indicates,
   the result is not always good; the main case where it works reliably is for
   the fundamental Cartan (|cn==0|). The circumstance that makes it useless in
   other fibers is that it fails to account for the different grading choices
   involved in identifying the (weak) real form in the fiber and in the KGB
   construction, so that there is no guarantee that the lifted element will
   appear to belong to the correct real form. In fact the computed Tits
   element |a| might not even give a strong involution $a.x_0.\delta$ at all.

   The main reason for leaving this (unused) method in the code is that it
   provides an alternative method (if it were called from the second |KGBHelp|
   constructor) to choosing the identity Tits element as starting point in the
   fundamental fiber, and that it illustrates at least how lifting of
   information from the fiber group should be done.
 */
tits::TitsElt EnrichedTitsGroup::naive_seed
 (const complexredgp::ComplexReductiveGroup& G,
  realform::RealForm rf, size_t cn) const
{
  // locate fiber, weak and strong real forms, and check central square class
  const cartanclass::Fiber& f=G.cartan(cn).fiber();
  const cartanset::CartanClassSet& ccs=G.cartanClasses();
  cartanclass::adjoint_fiber_orbit wrf = ccs.real_form_part(rf,cn);
  cartanclass::StrongRealFormRep srf=f.strongRepresentative(wrf);
  assert(srf.second==f.central_square_class(wrf));

  // now lift strong real form from fiber group to a torus part in |result|
  latticetypes::SmallBitVector v(bitset::RankFlags(srf.first),f.fiberRank());
  tits::TorusPart x = f.fiberGroup().fromBasis(v);

  // right-multiply this torus part by canonical twisted involution for |cn|
  tits::TitsElt result(titsGroup(),x,ccs.twistedInvolution(cn));

  return result; // result should be reduced immediatly by caller
}

/* The method |grading_seed| attempts to correct the shortcomings of
   |naive_seed| by insisting on obtaining an element exhibiting a grading that
   corresponds to the real form |rf|. Thus no element is actually recovered
   from any fiber group, but rather a set of equations for the torus part is
   set up and solved. Giving the right grading, one hopes that the Tits
   element automatically defines a strong involution in the square class; we
   test this by an |assert| statement

   Although the system of equations is highly underdetermined, it might
   suffice to fix a correct torus part modulo the subspace by which these
   torus parts are systematically reduced, and modulo torus parts that lie in
   $Z(G)$ (which cannot be detected by any grading of roots, but which for the
   same reason have no effect on the construction); in any case the seed
   produced here seems to work for partial KGB constructions that are useful
   or the construction of (small) blocks. In the more general case that there
   can be more than one minimal Cartan class in the set demanded, this method
   would risk giving non-coherent seeds in different classes, since nothing
   guarantees that the choices made will belong to the same strong real form;
   if this happens, too much elements will be generated Cartan classes that
   lie above more than one minimal Cartan class.
 */
tits::TitsElt EnrichedTitsGroup::grading_seed
  (const complexredgp::ComplexReductiveGroup& G,
   realform::RealForm rf, size_t cn) const
{
  // locate fiber and weak real form
  const cartanclass::Fiber& f=G.cartan(cn).fiber();
  const cartanset::CartanClassSet& ccs=G.cartanClasses();
  cartanclass::adjoint_fiber_orbit wrf = ccs.real_form_part(rf,cn);

  // get an element lying over the canonical twisted involution for |cn|
  tits::TitsElt a(titsGroup(),ccs.twistedInvolution(cn)); // trial element

  // get the grading of the imaginary root system given by the element |a|
  gradings::Grading base_grading;
  for (size_t i=0; i<f.imaginaryRank(); ++i)
    base_grading.set(i,grading(a,f.simpleImaginary(i)));

  // get the grading of the same system given by chosen representative of |wrf|
  gradings::Grading form_grading = f.grading(f.weakReal().classRep(wrf));

  /* set up equations with simple imaginary roots as left hand sides, and with
     the the bits of |base_grading-form_grading| as right hand side.
   */
  latticetypes::BinaryEquationList eqns(f.imaginaryRank());
  for (size_t i = 0; i < eqns.size(); ++i)
  {
    latticetypes::BinaryEquation& equation = eqns[i];
    equation = latticetypes::BinaryEquation(rd.root(f.simpleImaginary(i)));
    equation.pushBack((base_grading^form_grading)[i]);
  }

  // solve, and tack a solution |x| to the left of |a|.
  tits::TorusPart x(G.rank());
#ifdef DEBUG
  bool success=bitvector::firstSolution(x,eqns);
  assert(success);
#else
  bitvector::firstSolution(x,eqns);
#endif

  titsGroup().left_add(x,a); // form element $x.\sigma_w$

#ifdef DEBUG
  // double-check that we have found an element that gives the desired grading
  for (size_t i=0; i<f.imaginaryRank(); ++i)
    assert(grading(a,f.simpleImaginary(i))==form_grading[i]);

  tits::TitsElt check=a; Tg.mult(check,twisted(check));

  // maybe we should reduce here, but we don't have |FiberData| available...
  assert(check==tits::TitsElt(Tg));
#endif

  return a;  // result should be reduced immediatly by caller
}

/* In this final and most elaborate seeding function, which is also the most
   reliable one, we stoop down to simulating the KGB construction back from
   the fundamental fiber to the one for which we try to find a seed, and to
   try all the representatives in the fundamental fiber of the strong real
   form, until finding one that, along the chosen path of cross actions and
   Cayley transforms, proves to be suited for every necessary Cayley transform
   (making the simple root involved noncompact).
 */
tits::TitsElt EnrichedTitsGroup::backtrack_seed
 (const complexredgp::ComplexReductiveGroup& G,
  realform::RealForm rf, size_t cn) const
{
  const weyl::WeylGroup& W= titsGroup().weylGroup();

  const cartanset::CartanClassSet& ccs=G.cartanClasses();
  const weyl::TwistedInvolution& tw=ccs.twistedInvolution(cn);

  rootdata::RootList Cayley;
  weyl::WeylWord cross;
  cartanset::cayley_and_cross_part(Cayley,cross,tw,rd,W);

  /* at this point we can get from the fundamental fiber to |tw| by first
     applying cross actions according to |cross|, and then applying Cayley
     transforms in the strongly orthogonal set |Cayley|.
  */

  // transform strong orthogonal set |Cayley| back to distinguished involution
  for (size_t i=0; i<Cayley.size(); ++i)
    for (size_t j=cross.size(); j-->0; )
      rd.rootReflect(Cayley[i],cross[j]);

  /* at this point we can get from the fundamental fiber to |tw| by first
     applying Cayley transforms in the strongly orthogonal set |Cayley|, and
     then applying cross actions according to |cross|
  */

  /* Now find an element in the chosen strong real form at the fundamental
     fiber, that has noncompact grading on all the roots of |Cayley| (which
     are imaginary for $\delta$)
   */
  tits::TitsElt result(titsGroup());

  const cartanclass::Fiber& fund=G.fundamental();
  const partition::Partition& srp = fund.strongReal(square());
  for (unsigned long x=0; x<srp.size(); ++x)
    if (srp(x)==srp(f_orbit()))
    {
      latticetypes::SmallBitVector v
	(static_cast<bitset::RankFlags>(x),fund.fiberRank());
      tits::TorusPart t = fund.fiberGroup().fromBasis(v);
      for (size_t i=0; i<Cayley.size(); ++i)
	if (is_compact(t,Cayley[i]))
	  goto again; // none of the |Cayley[i]| should be compact

      // if we get here, |t| is OK as torus part
      result = tits::TitsElt(titsGroup(),t); // pure torus part
      goto found;
    again: {}
    }
  assert(false); // getting here means none of the orbit elements is in order

found:

  /* Now we must apply the Cayley transforms and cross actions to |result|.
     However, Cayley transforms by non-simple roots are not implemented, and
     so we reorder the operations as in |W.involution_expr(tw)|, which gives
     the same cross actions, but interspersed with simple Cayley transforms.
   */

  // transform |result| via Cayley transforms and cross actions
  std::vector<signed char> dec=W.involution_expr(tw);
  for (size_t j=dec.size(); j-->0; )
    if (dec[j]>=0)
    {
      assert(simple_grading(result,dec[j])); // must be noncompact
      Cayley_transform(result,dec[j]);
    }
    else
      basedTwistedConjugate(result,~dec[j]);

  assert(result.tw()==tw);

  return result;  // result should be reduced immediatly by caller
}


/*    II c. The main helper class |KGBHelp|  */

/*
   The actual KGB contruction takes place below. During the construction, the
   elements are represented as Tits group elements. The links in the KGB
   structure are realised by |BasedTitsGroup::basedTwistedConjugate| for the
   cross actions and by |basedTitsGroup::Cayley_transform| for Cayley
   transforms. After each of these, the result is subject to
   |d_fiberdata.reduce| to normalise the representation of a KGB element.
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

  Currently the fields |d_srf| and |d_base_compact| are only (possibly) used
  by the other constructor, and so are left just default-constructed here.
*/

KGBHelp::KGBHelp(realredgp::RealReductiveGroup& GR)
  : d_rank(GR.semisimpleRank())

  , d_cross(d_rank)
  , d_cayley(d_rank)
  , d_info()

  , d_pool(), d_tits(d_pool)

  // only the final fields depend on the real form of |GR|
  , basePoint(GR)
  , d_fiberData(GR.complexGroup(),GR.cartanSet())
{
  KGBElt size = GR.kgbSize();

  d_pool.reserve(size);
  d_info.reserve(size);

  // set up cross and cayley tables with undefined values
  for (size_t j = 0; j < d_rank; ++j)
  {
    d_cross[j].resize(size,UndefKGB);
    d_cayley[j].resize(size,0); // leave unset (set by |cayleyExtend|)
  }

  // set identity Tits element as seed of the KGB construction
  d_tits.match(tits::TitsElt(titsGroup()));

  // store its length and Cartan class (the latter is in fact always 0)
  d_info.push_back(KGBInfo(0,d_fiberData.cartanClass(d_tits[0].tw())));
}

/* The following constructor serves the pseudo-constructor |refined_helper|.
   Its main purpose is to allow precomputation of the |basePoint| field as
   the base object of a |EnrichedTitsGroup| object, and a set of |seeds|, one
   for each minimal Cartan class. The other arguments serve mainly to transmit
   values that are computed in |refined_helper| anyway, allowing them to be
   used in the construction without recompuation.
*/
KGBHelp::KGBHelp(const complexredgp::ComplexReductiveGroup& G,
		 const bitmap::BitMap& Cartan_classes,
		 KGBElt size, // size of KGB for selected |Cartan_classes|
		 const BasedTitsGroup& base, // base point information)
		 const std::vector<tits::TitsElt>& seeds)
  : d_rank(G.semisimpleRank())

  , d_cross(d_rank)
  , d_cayley(d_rank)
  , d_info()

  , d_pool() // no elements a priori, those from |seeds| are added below
  , d_tits(d_pool)

  // only the final fields depend on the real form of |GR|
  , basePoint(base)
  , d_fiberData(G,Cartan_classes)
{
  d_pool.reserve(size);
  d_info.reserve(size);

  // set up cross and cayley tables with undefined values
  for (size_t j = 0; j < d_rank; ++j)
  {
    d_cross[j].resize(size,UndefKGB);
    d_cayley[j].resize(size,0); // leave undefined (set by |cayleyExtend|)
  }

  set::SetEltList m=G.cartanClasses().ordering().minima(Cartan_classes);
  assert(m.size()==seeds.size()); // one seed for every minimal Cartan class

  for (size_t i=0; i<seeds.size(); ++i)
  {
    tits::TitsElt a= seeds[i];
    d_fiberData.reduce(a);

    // now add KGB element for the reduced Tits group element
#ifdef DEBUG
    size_t k=d_tits.match(a);
    assert(k==d_info.size()); // this KGB element should be new
#else
    d_tits.match(a); // enter new KGB element into the tables
#endif

    // add additional infomation (length,Cartan class) for this KGB element
    d_info.push_back(KGBInfo(G.weylGroup().involutionLength(a.tw()),m[i]));
  }
}


/*!
  \brief
  This pseudo-constructor builds a |KGBHelp| object initialized with an
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
  const complexredgp::ComplexReductiveGroup& G=GR.complexGroup();
  EnrichedTitsGroup square_class_base(GR,G.fundamental());

  realform::RealForm rf=GR.realForm();
  assert(square_class_base.square()==
	 G.fundamental().central_square_class(rf));

  const cartanset::CartanClassSet& ccs=G.cartanClasses();
  KGBElt size =  ccs.KGB_size(rf,Cartan_classes);

  set::SetEltList m=ccs.ordering().minima(Cartan_classes);

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
  for (KGBElt j = 0; j < d_pool.size(); ++j) {
    // these calls will usually enlarge |d_pool.size()|;
    cross_extend(j);
    cayleyExtend(j);
  }
  return *this;
}

size_t KGBHelp::export_tables(std::vector<KGBEltList>& cross,
			      std::vector<KGBEltList>& cayley,
			      tits::TitsEltList& tits,
			      std::vector<KGBInfo>& info,
			      BasedTitsGroup*& base,
			      bool traditional) const
{
  size_t size=d_info.size();

  // sort (partially)
  setutils::Permutation a(size,1);
  IndexCompare comp(*this);
  if (traditional)
    std::sort(a.begin(),a.end(),comp);
  else
    std::stable_sort(a.begin(),a.end(),comp); // better, faster, more reliable

  // export the cross and cayley maps, permuting each constituent list
  cross.clear(); cayley.clear();
  cross.reserve(d_rank); cayley.reserve(d_rank);
  setutils::Permutation ai(a,-1); // compute inverse of |a|
  for (size_t s = 0; s < d_rank; ++s) {
    cross.push_back(a.pull_back(ai.renumber(d_cross[s])));
    cayley.push_back(a.pull_back(ai.renumber(d_cayley[s],UndefKGB)));
  }

  // export the other lists
  /* We cannot do |a.pull_back(d_pool).swap(tits);|, as we need to convert
     a vector of |TE_Entry| to a vector of |TitsElt|; a loop is required */
  tits.clear(); tits.reserve(size);
  for (KGBElt x=0; x<size; ++x)
    tits.push_back(d_pool[a[x]]);

  // but here there is no difficulty
  a.pull_back(d_info).swap(info); // pull back |d_info| through |a|, and export

  base=new BasedTitsGroup(basePoint); // export base point infomation

  return size;
}

/******** manipulators *******************************************************/

/*!
  \brief Tries to enlarge the parameter set by cross-actions, without
  supposing that elements are generated by increasing order of twisted
  involution length

  Precondition: |parent| is the index into |d_tits| of the parameter we are
  extending from.
*/
void KGBHelp::cross_extend(KGBElt parent)
{
  const tits::TitsElt current = d_tits[parent];

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
    KGBElt child = d_tits.match(a);
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
  const tits::TitsElt current = d_tits[parent];

  // try to get new elements by Cayley transform
  for (size_t s = 0; s < d_rank; ++s) {

    gradings::Status::Value v = d_info[parent].status[s];
    if (v != gradings::Status::ImaginaryNoncompact) {
      d_cayley[s][parent] = UndefKGB;
      continue;
    }

    // cayley-transform |current| by $\sigma_s$
    tits::TitsElt a = current; basePoint.Cayley_transform(a,s);
    assert(titsGroup().length(a)>titsGroup().length(current)); // should go up


    // now look up the correspondingly reduced Tits element
    d_fiberData.reduce(a); // subspace has grown, so mod out new supspace
    KGBElt x = d_tits.match(a);
    if (x==d_info.size()) // add a new Tits element
      d_info.push_back(KGBInfo(d_info[parent].length+1, // length goes up
			       d_fiberData.cartanClass(a.tw())
			       ));
    // add new cayley link
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

/*! \brief
  Returns the grading offset (on simple roots) adapted to |G|. This flags the
  simple roots that are noncompact imaginary at the fundamental Cartan in G.

Algorithm: the variable |rset| is first made to flag, among the imaginary
roots of the fundamental Cartan, those that are noncompact for the chosen
representative (in the adjoint fiber) of the real form of |G|. The result is
formed by extracting only the information concerning the presence of the
\emph{simple} roots in |rset|.
*/
gradings::Grading grading_offset_for(const realredgp::RealReductiveGroup& G)
{
  const rootdata::RootDatum& rd = G.rootDatum();
  rootdata::RootSet rset= G.noncompactRoots();

  return cartanclass::restrictGrading(rset,rd.simpleRootList());
}

/*!
 \brief Returns the grading offset for the base real form of the square class
 (coset in adjoint fiber group) |csc|; |f| and |rd| are corresponding values
 */
gradings::Grading
square_class_grading_offset(const cartanclass::Fiber& f,
			    cartanclass::square_class csc,
			    const rootdata::RootDatum& rd)
{
  rootdata::RootSet rset = f.noncompactRoots(f.class_base(csc));
  return cartanclass::restrictGrading(rset,rd.simpleRootList());
}


/*!
  \brief Puts in |Hasse| the Hasse diagram of the Bruhat ordering on |kgb|.

  Explanation: this is the closure ordering of orbits. We use the algorithm
  from Richardson and Springer.
*/
void makeHasse(std::vector<set::SetEltList>& Hasse, const KGB& kgb)
{
  using gradings::Status;

  Hasse.resize(kgb.size());

  for (KGBElt x = 0; x < kgb.size(); ++x)
  {
    std::set<KGBElt> h_x;
    const kgb::Descent& d = kgb.descent(x);
    if (d.none()) // element is minimal in Bruhat order
      continue;

    size_t s = d.firstBit();
    KGBElt sx;

    if (kgb.status(s,x) == Status::Complex)
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
	if (kgb.status(s,z) == Status::ImaginaryNoncompact)
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
