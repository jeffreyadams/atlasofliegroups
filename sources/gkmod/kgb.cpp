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
namespace kgb {

/*
	                Local functions
*/

namespace {

void makeHasse(std::vector<set::SetEltList>&, const KGB&);

} // |namespace|

/*
	           Some small additional classes
*/


// a |GlobalTitsElement| can be hashed as (fingerprint,twisted_inv) pair
struct KGB_elt_entry
{
  latticetypes::RatLatticeElt fingerprint;
  weyl::TI_Entry ti;

  // obligatory fields for hashable entry
  typedef std::vector<KGB_elt_entry> Pooltype;
  size_t hashCode(size_t modulus) const; // hash function
  bool operator !=(const KGB_elt_entry& x) const
    { return ti!=x.ti or fingerprint!=x.fingerprint; }

  KGB_elt_entry (const latticetypes::RatLatticeElt f,
		 const weyl::TwistedInvolution& tw)
    : fingerprint(f), ti(tw) {}
}; //  |struct KGB_elt_entry|

size_t KGB_elt_entry::hashCode(size_t modulus) const
{
  unsigned long d= fingerprint.denominator();
  const latticetypes::LatticeElt& num=fingerprint.numerator();
  size_t h=ti.hashCode(modulus);
  for (size_t i=0; i<num.size(); ++i)
    h=((h*d)+num[i])&(modulus-1);
  return h;
}

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
}; // |class FiberData|


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


/*
                 The |KGBHelp| class declaration, a helper class for |KGB|

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

}; // |class KGBHelp|

// A non-method construction function (could also have been static method)
KGBHelp refined_helper(realredgp::RealReductiveGroup& GR,
		       const bitmap::BitMap& Cartan_classes);


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
}; // |class IndexCompare|



/*

        The |KGB_base| class implementation, public methods

*/

KGBEltPair KGB_base::tauPacket(const weyl::TwistedInvolution& w) const
{
  unsigned int i=inv_hash.find(w);
  if (i==inv_hash.empty)
    return KGBEltPair(0,0);
  return KGBEltPair(first_of_tau[i],first_of_tau[i+1]);
}

size_t KGB_base::packet_size(const weyl::TwistedInvolution& w) const
{
  unsigned int i=inv_hash.find(w);
  if (i==inv_hash.empty)
    return 0;
  return first_of_tau[i+1]-first_of_tau[i];
}

/*

        The |KGB_base| class implementation, protected methods

*/

// resserve space in |data| and |info| arrays
void KGB_base::reserve (size_t n)
{
  info.reserve(n);
  for (size_t s=0; s<data.size(); ++s)
    data[s].reserve(n);
}

// add a slot to the |info| array, and reserve space in the |data| arrays
void KGB_base::add_element
  (unsigned int length, unsigned int Cartan_class, weyl::TwistedInvolution tw)
{
  info.push_back(KGBEltInfo(length,Cartan_class,tw));
  for (size_t s=0; s<data.size(); ++s)
    data[s].push_back(KGBfields());
}


/*

        The |global_KGB| class implementation, public methods

*/

global_KGB::global_KGB(complexredgp::ComplexReductiveGroup& G_C,
		       const tits::GlobalTitsGroup& Tits_group)
  : KGB_base(G_C.twistedWeylGroup())
  , G(G_C)
  , Tg(Tits_group)
  , fiber_data(G,(generate_involutions(G_C.numInvolutions()),inv_hash))
  , elt()
{
  size_t size = G.global_KGB_size();
  elt.reserve(size);
  KGB_base::reserve(size);

  { // get elements at the fundamental fiber
    const latticetypes::SmallSubquotient& fg = G.fundamental().fiberGroup();
    first_of_tau.push_back(0); // start of fundamental fiber

    unsigned long n= 1ul << Tg.square_class_generators().size();
    for (unsigned long c=0; c<n; ++c)
    {
      gradings::Grading gr=bitvector::combination
	(Tg.square_class_generators(),bitset::RankFlags(c));
      latticetypes::RatWeight rw (G.rank());
      for (gradings::Grading::iterator it=gr.begin(); it(); ++it)
	rw += G.rootDatum().fundamental_coweight(*it);

      tits::TorusElement t(rw/=2); // halve: image of $x\mapsto\exp(\pi\ii x)$

      for (unsigned long i=0; i<fg.size(); ++i)
      {
	tits::TorusElement s=t;
	s += fg.fromBasis // add |TorusPart| from fiber group, use bits from |i|
	  (latticetypes::SmallBitVector(bitset::RankFlags(i),fg.dimension()));
	elt.push_back(tits::GlobalTitsElement(s));
      } // |for (i)|
    } // |for (c)|
    first_of_tau.push_back(elt.size()); // end of fundamental fiber
  }

  generate_elements(size);

}


bool global_KGB::compact(rootdata::RootNbr n, // assumed imaginary at |a|
			 const tits::GlobalTitsElement& a) const
{
  const latticetypes::Weight& alpha= G.rootDatum().root(n);
  return a.torus_part().negative_at(alpha) == // question was: whether compact
    fiber_data.at_rho_imaginary(alpha,a.tw());
}


/*

        The |global_KGB| class implementation, private methods

*/

// after construction of the |KGB_base| we precompute and hash all involutions
void global_KGB::generate_involutions(size_t n)
{
  inv_pool.reserve(n);  // filled below
  first_of_tau.reserve(n+1); // to be filled by constructor of derived object
  const weyl::TwistedWeylGroup& W=twistedWeylGroup();

  inv_hash.match(weyl::TwistedInvolution()); // set initial element
  for (size_t i=0; i<inv_pool.size(); ++i) // pool grows from 1 to |n|
  {
    const weyl::TI_Entry& parent=inv_pool[i];
    for (weyl::Generator s=0; s<W.rank(); ++s)
      if (W.hasTwistedCommutation(s,parent))
	inv_hash.match(W.prod(s,parent));
      else
        inv_hash.match(W.twistedConjugated(parent,s));
  }
  assert(inv_pool.size()==n);
}

void global_KGB::generate_elements(size_t n_elts)
{
  const weyl::TwistedWeylGroup& W = Tg;
  const rootdata::RootDatum& rd = G.rootDatum();

  KGB_elt_entry::Pooltype elt_pool; elt_pool.reserve(n_elts);
  hashtable::HashTable<KGB_elt_entry,unsigned long> elt_hash(elt_pool);

  {
    weyl::TwistedInvolution e; // identity
    for (size_t i=0; i<elt.size(); ++i)
    {
      add_element(0,0,e); // blank Cartan, length, involution fields
      elt_hash.match(KGB_elt_entry(fiber_data.fingerprint(elt[i]),e));
    }
    assert(elt_hash.size()==elt.size()); // all distinct; in fact an invariant
  }

  for (size_t inv_nr=0; inv_nr<inv_pool.size(); ++inv_nr)
  {
    const weyl::TwistedInvolution& tw = inv_pool[inv_nr];
    for (weyl::Generator s=0; s<W.rank(); ++s)
    {
      weyl::TwistedInvolution new_tw =W.twistedConjugated(tw,s);
      size_t new_nr = inv_hash.find(new_tw);
      bool is_new = new_nr+1 >= first_of_tau.size();
      if (is_new)
	assert(new_nr+1==first_of_tau.size()); // since we mimick generation
      bool imaginary = new_nr==inv_nr and not W.hasDescent(s,tw);

      // generate cross links
      for (KGBElt x=first_of_tau[inv_nr]; x<first_of_tau[inv_nr+1]; ++x)
      {
 	tits::GlobalTitsElement child=elt[x]; //start out with a copy
	int d = Tg.cross(s,child); // cross act and record length difference
	assert(child.tw()==new_tw);
	KGB_elt_entry ee(fiber_data.fingerprint(child),new_tw);
	KGBElt k = elt_hash.match(ee);
	if (k==elt.size()) // then new
	{
	  assert(is_new);
	  elt.push_back(child);
	  assert(fiber_data.Cartan_class(new_tw)==Cartan_class(x));
	  add_element(length(x)+d,Cartan_class(x),new_tw);
	}
	data[s][x].cross_image=k;

	if (d!=0) // just made complex cross action
	{
	  info[x].status.set(s,gradings::Status::Complex);
	  info[x].desc.set(s,d<0);
	}
	else if (imaginary)
	{
	  info[x].status.set_imaginary // always true (noncompact) at identity
	    (s,not elt[x].torus_part().negative_at(rd.simpleRoot(s)));
	  info[x].desc.set(s,false); // imaginary roots are never descents
	}
	else // real
	{
	  assert(x==k);
	  info[x].status.set(s,gradings::Status::Real);
	  info[x].desc.set(s,true); // real roots are always descents
	}
      } // |for(x)|

      // generate Cayley links
      if (imaginary)
      {
	new_tw = W.prod(s,tw);
	new_nr = inv_hash.find(new_tw);
	is_new = new_nr+1 >= first_of_tau.size();
	if (is_new)
	  assert(new_nr+1==first_of_tau.size()); // since we mimick generation

	for (KGBElt x=first_of_tau[inv_nr]; x<first_of_tau[inv_nr+1]; ++x)
	  if (info[x].status[s]==gradings::Status::ImaginaryNoncompact)
	  {
	    tits::GlobalTitsElement child=Tg.Cayley(s,elt[x]);
	    assert(child.tw()==new_tw);
	    KGBElt k=elt_hash.match
	      (KGB_elt_entry(fiber_data.fingerprint(child),new_tw));
	    if (k==elt.size()) // then new
	    {
	      elt.push_back(child);
	      add_element(length(x)+1,fiber_data.Cartan_class(new_tw),new_tw);
	    }
	    data[s][x].Cayley_image = k;
	    if (data[s][k].inverse_Cayley_image.first==UndefKGB)
	      data[s][k].inverse_Cayley_image.first=x;
	    else
	      data[s][k].inverse_Cayley_image.second=x;

	  } // |for (x)|
      } // |if (imaginary)|

      if (is_new) // in upwards cross action or Cayley transform, never both
	first_of_tau.push_back(elt.size()); // close the tau packet


    } // |for(s)|
  } // |for(inv_nr)|

  assert(elt.size()==n_elts);
}



/*****************************************************************************

        The KGB class, public methods

******************************************************************************/

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
  : KGB_base(GR.twistedWeylGroup())
  , d_tits()
  , d_state()
  , d_bruhat(NULL)
  , d_base(NULL)
{
  bool traditional= Cartan_classes.size()==0; // whether traditional generation

  std::vector<KGBEltList> d_cross;
  std::vector<KGBEltList> d_cayley;
  std::vector<KGBInfo> d_info;

  size_t size =
    (traditional ? KGBHelp(GR) : refined_helper(GR,Cartan_classes))
    .fill()
    .export_tables(d_cross,d_cayley,d_tits,d_info,d_base);

  // check that the size is right
  assert(size == (traditional ? GR.KGB_size()
		  : GR.complexGroup().KGB_size(GR.realForm(),Cartan_classes)
		  ));

  KGB_base::reserve(size);
  first_of_tau.reserve(traditional
		       ? GR.numInvolutions()+1
		       : GR.complexGroup().numInvolutions(Cartan_classes)+1);

  for (KGBElt x=0; x<size; ++x)
  {
    KGB_base::add_element(d_info[x].length,d_info[x].cartan,d_tits[x].tw());
    info[x].status=d_info[x].status; info[x].desc=d_info[x].desc;
    unsigned int h=inv_hash.match(d_tits[x].tw());
    if (h==first_of_tau.size()) // then twisted involution was not seen before
      first_of_tau.push_back(x); // push KGB element as first for involution
    else // check that all non-new involutions are equal to most recent one
      assert(h==first_of_tau.size()-1);
    for (weyl::Generator s=0; s<rank(); ++s)
    {
      data[s][x].cross_image=d_cross[s][x];
      data[s][x].Cayley_image = d_cayley[s][x];
    }
  }
  first_of_tau.push_back(size); // push terminating number

  for (KGBElt x=0; x<size; ++x)
    for (weyl::Generator s=0; s<rank(); ++s)
    {
      KGBElt c=data[s][x].Cayley_image;
      if (c!=UndefKGB)
      {
	KGBEltPair& target=data[s][c].inverse_Cayley_image;
	if (target.first==UndefKGB) target.first=x;
	else target.second=x;
      }
    }
}

/******** copy, assignment and swap ******************************************/


/******** accessors **********************************************************/


// Looks up a |tits::TitsElt| value and returns its KGB number, or |size()|
// Since KGB does not have mod space handy, must assume |a| already reduced
KGBElt KGB::lookup(const tits::TitsElt& a, const tits::TitsGroup& Tg) const
{
  KGBEltPair p = tauPacket(a.tw());
  tits::TorusPart t = Tg.left_torus_part(a);
  for (KGBElt x=p.first; x<p.second; ++x)
    if (Tg.left_torus_part(titsElt(x))==t)
      return x;
  return size(); // report failure
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



/*****************************************************************************

            Chapter II -- The auxiliary classes, methods

******************************************************************************/

/*    II a. |GlobalFiberData|  */

// precompute data for all involutions, given a hash table |h| enumerating them
// data: Cartan class, projector of torus parts and $\check\rho_{im}$
GlobalFiberData::GlobalFiberData
  (complexredgp::ComplexReductiveGroup& G,
   const hashtable::HashTable<weyl::TI_Entry,unsigned int>& h)
  : hash_table(h)
  , Cartans(h.size())
  , proj(h.size())
  , check_2rho_imag(h.size())
{
  const rootdata::RootDatum& rd=G.rootDatum();
  const weyl::TwistedWeylGroup& W = G.twistedWeylGroup();

  std::vector<latticetypes::LatticeMatrix> refl(G.semisimpleRank());
  for (weyl::Generator s=0; s<refl.size(); ++s)
    // get automorphism of lattice $X_*$ given by generator $s$
    // reflection map will be used as automorphism of $X_*\tensor\Q/X_*$
    refl[s] = rd.simple_reflection(s).transposed();

  std::vector<size_t> to_do(h.size());
  bitmap::BitMap seen(h.size()); // flags involution numbers in |to_do| array

  for (size_t cn=0; cn<G.numCartanClasses(); ++cn)
  {
    const cartanclass::CartanClass& cc = G.cartan(cn);
    size_t first = h.find(G.twistedInvolution(cn));

    { // store data for canonical twisted involution |i| of Cartan class |cn|
      Cartans[first]=cn;

      latticetypes::LatticeMatrix A =cc.involution().transposed();
      for (size_t i=0; i<G.rank(); ++i)
	A(i,i) += 1;
      // now $A=\theta_x^t+1$ is a matrix whose kernel is $(X_*)^{-\theta_x^t}$
      proj[first]=lattice::row_saturate(A);

      check_2rho_imag[first]=rd.twoRho(cc.dualFiber().imaginaryRootSet());
    }

    to_do.assign(1,first);
    seen.insert(first);

    for (size_t i=0; i<to_do.size(); ++i) // |to_do| grows during the loop
      for (size_t s=0; s<G.semisimpleRank(); ++s)
      {
	weyl::TwistedInvolution tw = W.twistedConjugated(h[to_do[i]],s);
	size_t k = h.find(tw);
	assert (k!=h.empty);
	if (seen.isMember(k))
	  continue;
	seen.insert(k);
	to_do.push_back(k);
	Cartans[k]=cn;                  // record number of Cartan class
	proj[k]=proj[to_do[i]]*refl[s]; // apply $refl[s]^{-1}$ to old kernel
	check_2rho_imag[k]=refl[s].apply(check_2rho_imag[to_do[i]]);
      }

  } // |for(cn)|
  // check that the number of generated elements is as predicted
}

//  this demonstrates what the |proj matrices can be used for
bool GlobalFiberData::equivalent(const tits::GlobalTitsElement& x,
				 const tits::GlobalTitsElement& y) const
{
  unsigned int k = hash_table.find(x.tw());
  assert (hash_table.find(y.tw())==k);
  latticetypes::RatWeight t= (x.torus_part()-y.torus_part()).as_rational();
  latticetypes::LatticeElt p = proj[k].apply(t.numerator());

  for (size_t i=0; i<p.size(); ++i)
    if (p[i]%t.denominator()!=0)
      return false;

  return true;
}

// computing "fingerprints" allows direct comparison without using |equivalent|
latticetypes::RatLatticeElt
  GlobalFiberData::fingerprint(const tits::GlobalTitsElement& x) const
{
  unsigned int k = hash_table.find(x.tw());
  latticetypes::RatLatticeElt t = x.torus_part().as_rational();
  latticetypes::LatticeElt p = proj[k].apply(t.numerator());

  // reduce modulo integers and return
  for (size_t i=0; i<p.size(); ++i)
    p[i]= intutils::remainder(p[i],t.denominator());
  return latticetypes::RatLatticeElt(p,t.denominator()).normalize();
}


/*    II b. |FiberData|  */


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
    // reflection map is induced vector space endomorphism of $X_* / 2X_*$
    refl[s] = latticetypes::BinaryMap(rd.simple_reflection(s).transposed());

  for (bitmap::BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
  {
    size_t cn=*it;
    size_t i = hash_table.match(G.twistedInvolution(cn));
    assert(i==data.size()); // this twisted involution should be new

    { // store data for canonical twisted involution |i| of Cartan class |cn|
      const latticetypes::LatticeMatrix &q = G.cartan(cn).involution();
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




/*    II c. The main helper class |KGBHelp|, implementation  */

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
  for (size_t x = 0; x < d_rank; ++x)
  {
    d_cross[x].resize(size,UndefKGB);
    d_cayley[x].resize(size,0); // leave unset (set by |cayleyExtend|)
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
  for (size_t x = 0; x < d_rank; ++x)
  {
    d_cross[x].resize(size,UndefKGB);
    d_cayley[x].resize(size,0); // leave undefined (set by |cayleyExtend|)
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

  Here we actually look up the strong real form in order to get a proper
  initial Tits group element associated to this Cartan class

  The base grading is set up to correspond to (the chosen adjoint fiber
  element in the fundamental Cartan for) the central square class of this
  real form; all this is done inside the |EnrichedTitsGroup| constructor.

  The initial element then represents the place within its central square class
  of one chosen strong real form |srf| lying over this weak real form; it has
  the identity twisted involution, and a torus factor obtained by lifting the
  representative fiber group element |srf.first| via the |fromBasis| method of
  the fiber group back to the ``coweight lattice modulo 2'' \f$Y/2Y\f$.
*/
KGBHelp refined_helper(realredgp::RealReductiveGroup& GR,
		       const bitmap::BitMap& Cartan_classes)
{
  tits::EnrichedTitsGroup square_class_base(GR);

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
  for (KGBElt x = 0; x < Tits.size(); ++x) {
    // these calls will usually enlarge |Tits.size()|;
    cross_extend(x);
    cayleyExtend(x);
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
      assert(child==parent);
      d_info[parent].status.set(s,gradings::Status::Real); // |child==parent|
      d_info[parent].desc.set(s,true); // real roots are always descents
    }
    else // imaginary
    {
      d_info[parent].status.set_imaginary
	(s,basePoint.simple_grading(current,s));
      d_info[parent].desc.set(s,false); // imaginary roots are never descents
      if (child!=parent) // give child same status (which will be noncompact)
      {
	assert(d_info[parent].status[s]==gradings::Status::ImaginaryNoncompact);
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

  It is assumed that it is an invariant of the |KGBHelp| structure that
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


/*****************************************************************************

        Chapter III -- Functions local to kgb.cpp

******************************************************************************/

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
