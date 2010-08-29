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
#include "subdatum.h"
#include "set.h"
#include "tits.h"
#include "tori.h"
#include "weyl.h"

#include <iostream>
#include <iomanip>
#include "basic_io.h"
#include "prettyprint.h"
#include "ioutils.h"

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

void makeHasse(std::vector<set::SetEltList>&, const KGB_base&);

} // |namespace|



size_t KGB_elt_entry::hashCode(size_t modulus) const
{
  unsigned long d= fingerprint.denominator();
  const latticetypes::LatticeElt& num=fingerprint.numerator();
  size_t h=tw.hashCode(modulus);
  for (size_t i=0; i<num.size(); ++i)
    h=((h*d)+num[i])&(modulus-1);
  return h;
}




/*

        The |KGB_base| class implementation, public methods

*/

KGBElt KGB_base::cross_act(const weyl::WeylWord ww, KGBElt x) const
{
  for (size_t i=ww.size(); i-->0;)
    x = cross(ww[i],x);
  return x;
}

KGBElt KGB_base::cross_act(KGBElt x, const weyl::WeylWord ww) const
{
  for (size_t i=0; i<ww.size(); ++i)
    x = cross(ww[i],x);
  return x;
}

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
  (unsigned int length, weyl::TwistedInvolution tw)
{
  info.push_back(KGBEltInfo(length,tw));
  for (size_t s=0; s<data.size(); ++s)
    data[s].push_back(KGBfields());
}


/*

        The |global_KGB| class implementation, public methods

*/

// create structure incorporating all KGB structures for a given inner class
global_KGB::global_KGB(complexredgp::ComplexReductiveGroup& G)
  : KGB_base(G.twistedWeylGroup())
  , Tg(G) // construct global Tits group as subobject
  , fiber_data(G,(generate_involutions(G.numInvolutions()),inv_hash))
  , elt()
{
  size_t size = G.global_KGB_size();
  elt.reserve(size);
  KGB_base::reserve(size);

  { // get elements at the fundamental fiber, for "all" square classes
    const latticetypes::SmallSubquotient& fg = G.fundamental().fiberGroup();
    first_of_tau.push_back(0); // start of fundamental fiber

    // get number of subsets of generators, and run through them
    unsigned long n= 1ul << Tg.square_class_generators().size();
    for (unsigned long c=0; c<n; ++c)
    {
      gradings::Grading gr=bitvector::combination
	(Tg.square_class_generators(),bitset::RankFlags(c));
      latticetypes::RatWeight rw (G.rank());
      for (gradings::Grading::iterator // for flagged (imaginary) simple roots
	     it=gr.begin(); it(); ++it)
	rw += G.rootDatum().fundamental_coweight(*it); // a sum of f. coweights

      tits::TorusElement t(rw/=2); // halve: image of $x\mapsto\exp(\pi\ii x)$

      weyl::TwistedInvolution e; // identity
      for (unsigned long i=0; i<fg.size(); ++i)
      {
	tits::TorusElement s=t;
	s += fg.fromBasis // add |TorusPart| from fiber group; bits from |i|
	  (latticetypes::SmallBitVector(bitset::RankFlags(i),fg.dimension()));
	elt.push_back(tits::GlobalTitsElement(s));
	add_element(0,e); // create in base; length, tw both trivial
      } // |for (i)|
    } // |for (c)|
    first_of_tau.push_back(elt.size()); // end of fundamental fiber
  }

  generate(G.rootDatum(),size);

}

global_KGB::global_KGB(complexredgp::ComplexReductiveGroup& G,
		       const tits::GlobalTitsElement& x)
  : KGB_base(G.twistedWeylGroup())
  , Tg(G) // construct global Tits group as subobject
  , fiber_data(G,(generate_involutions(G.numInvolutions()),inv_hash))
  , elt()
{
  tits::GlobalTitsElement a=x; // start at an element that we certainly want
  weyl::Generator s;
  while ((s=Tg.weylGroup().leftDescent(a.tw()))<rank())
    if (Tg.hasTwistedCommutation(s,a.tw()))
      Tg.inverse_Cayley(s,a);
    else
      Tg.cross(s,a);

  assert (a.tw()==weyl::TwistedInvolution()); // we are at fundamental fiber
  { // get elements at the fundamental fiber
    first_of_tau.push_back(0); // start of fundamental fiber
    KGB_elt_entry::Pooltype elt_pool;
    hashtable::HashTable<KGB_elt_entry,unsigned long> elt_hash(elt_pool);
    add_element(0,a.tw()); // create in base; length, tw
    elt_hash.match(fiber_data.pack(a));

    /* Generate fundamental fiber. We use that cross actions by imaginary,
       non-compact, simple roots suffice, but without knowing a nice proof. */
    for (size_t i=0; i<elt_hash.size(); ++i) // |elt_hash| grows during loop
      for (weyl::Generator s=0; s<rank(); ++s)
	if (Tg.twisted(s)==s // then imaginary, since fiber is fundamental
	    and not Tg.compact(s,a=elt_hash[i].repr())) // noncompact
	{
	  Tg.cross(s,a);
	  size_t old_size = elt_hash.size();
	  if (elt_hash.match(fiber_data.pack(a))==old_size) // then it's new
 	    add_element(0,a.tw()); // create in base; length, tw
	}
    first_of_tau.push_back(elt_hash.size()); // end of fundamental fiber

    elt.reserve(elt_hash.size()); // now copy elements from hash table to |elt|
    for (size_t i=0; i<elt_hash.size(); ++i)
      elt.push_back(elt_hash[i].repr());
  }
  generate(G.rootDatum(),0); // complete element generation, no predicted size
} // |global_KGB::global_KGB|

bool global_KGB::compact(const rootdata::RootDatum& rd,
			 rootdata::RootNbr n, // assumed imaginary at |a|
			 const tits::GlobalTitsElement& a) const
{
  const latticetypes::Weight& alpha= rd.root(n);
  return a.torus_part().negative_at(alpha) == // question was: whether compact
    (fiber_data.at_rho_imaginary(alpha,a.tw())%2!=0); // CHECK ME!
}

kgb::KGBElt global_KGB::lookup(const tits::GlobalTitsElement& a) const
{
  KGBEltPair p = tauPacket(a.tw());
  for (KGBElt x=p.first; x<p.second; ++x)
    if (fiber_data.equivalent(element(x),a))
      return x;
  return size(); // report failure
}

std::ostream& global_KGB::print(std::ostream& strm, KGBElt x) const
{
  latticetypes::RatWeight t = torus_part(x).as_rational();
  return strm << std::setw(3*t.size()+3)<< t;
}

/*

        The |global_KGB| class implementation, private methods

*/

// after construction of the |KGB_base| we precompute and hash all involutions
void global_KGB::generate_involutions(size_t n)
{
  inv_pool.reserve(n);  // filled below
  first_of_tau.reserve(n+1); // to be filled by constructor of derived object

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
} // |global_KGB::generate_involutions|

void global_KGB::generate(const rootdata::RootDatum& rd, size_t predicted_size)
{
  const weyl::TwistedWeylGroup& W = Tg; // for when |GlobalTitsGroup| not used

  KGB_elt_entry::Pooltype elt_pool; elt_pool.reserve(predicted_size);
  hashtable::HashTable<KGB_elt_entry,KGBElt> elt_hash(elt_pool);

  {
    weyl::TwistedInvolution e; // identity
    for (size_t i=0; i<elt.size(); ++i)
      elt_hash.match(fiber_data.pack(elt[i]));

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
	KGB_elt_entry ee = fiber_data.pack(child);
	KGBElt k = elt_hash.match(ee);
	if (k==elt.size()) // then new
	{
	  assert(is_new);
	  elt.push_back(child);
	  assert(fiber_data.Cartan_class(new_tw)==Cartan_class(x));
	  add_element(length(x)+d,/* Cartan_class(x), */new_tw);
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
	    KGBElt k=elt_hash.match(fiber_data.pack(child));
	    if (k==elt.size()) // then new
	    {
	      elt.push_back(child);
	      add_element(length(x)+1,/*fiber_data.Cartan_class(new_tw),*/
			  new_tw);
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

  assert(elt.size()==predicted_size or predicted_size==0);
} // |global_KGB::generate|

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
  std::vector<unsigned int> involution_length;
  std::vector<unsigned int> Cartan_class;
public:
  FiberData(complexredgp::ComplexReductiveGroup& G,
	    const bitmap::BitMap& Cartan_classes);

  size_t n_involutions() const { return pool.size(); }
  const weyl::TwistedInvolution& involution(size_t i) const { return pool[i]; }

  unsigned int seq_no(const weyl::TwistedInvolution& tw) const
  { return hash_table.find(tw); }

  unsigned int length(size_t i) const { return involution_length[i]; }

  unsigned int cartanClass(const weyl::TwistedInvolution& tw) const
  { return Cartan_class[seq_no(tw)]; }

  void reduce(tits::TitsElt& a) const;

private: // the space actually stored need not be exposed
  const latticetypes::SmallSubspace& mod_space(const tits::TitsElt& a) const
  {
    size_t i=hash_table.find(a.tw());
    assert(i!=hash_table.empty);
    return data[i];
  }
}; // |class FiberData|




/*

        The KGB class, public methods

*/

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
  , rd(GR.rootDatum())
  , Cartan()
  , left_torus_part()
  , d_state()
  , d_bruhat(NULL)
  , d_base(NULL)
{
  bool traditional= Cartan_classes.size()==0; // whether traditional generation
  size_t size = generate(GR,Cartan_classes);

  // check that the size is right
  assert(size == (traditional ? GR.KGB_size()
		  : GR.complexGroup().KGB_size(GR.realForm(),Cartan_classes)
		  ));

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

std::ostream& KGB::print(std::ostream& strm, KGBElt x) const
{
  return prettyprint::prettyPrint(strm,torus_part(x));
}

/*

     The KGB class, private manipulators

*/

// Auxiliary class whose declaration must have been seen in |KGB::generate|
class invol_compare
{
  const FiberData& fd;
  const weyl::WeylGroup& W;
public:
  invol_compare(const FiberData& f, const weyl::WeylGroup& WG)
    : fd(f), W(WG) {}
  bool operator() (unsigned int i, unsigned int j) const
  {
    if (fd.length(i)!=fd.length(j)) // compare involution lengths
      return fd.length(i)<fd.length(j);
    unsigned long li=W.length(fd.involution(i)); // then Weyl group lengths
    unsigned long lj=W.length(fd.involution(j));
    return li!=lj ? li<lj
      : fd.involution(i)<fd.involution(j); // if all else fails, internal |<|
  }
}; // |class invol_compare|


size_t KGB::generate
  (realredgp::RealReductiveGroup& GR, const bitmap::BitMap& Cartan_classes)
{
  complexredgp::ComplexReductiveGroup& G=GR.complexGroup();
  size_t rank = G.semisimpleRank();
  bool traditional=Cartan_classes.size()==0; // whether traditional generation

  tits::TE_Entry::Pooltype elt_pool; // of size |size()|
  hashtable::HashTable<tits::TE_Entry,KGBElt> elt_hash(elt_pool);

  //! Permits reducing each Tits group element modulo its fiber denominator
  FiberData fiber_data(G, traditional ? GR.Cartan_set() : Cartan_classes);

  size_t size =
    traditional ? GR.KGB_size() : G.KGB_size(GR.realForm(),Cartan_classes);

  elt_pool.reserve(size);
  KGB_base::reserve(size);

  if (traditional)
  {
    d_base = new tits::BasedTitsGroup(G,GR.grading_offset());

    elt_hash.match(tits::TitsElt(titsGroup())); // identity Tits element is seed
    // store its length and Cartan class (both are in fact always 0)
    const weyl::TwistedInvolution& tw = elt_hash[0].tw();
    KGB_base::add_element(0,tw);
  }
  else
  {
    tits::EnrichedTitsGroup square_class_base(GR);
    d_base = new tits::BasedTitsGroup(square_class_base);
    realform::RealForm rf=GR.realForm();
    assert(square_class_base.square()==
	   G.fundamental().central_square_class(rf));

    set::SetEltList m=G.Cartan_ordering().minima(Cartan_classes);

    for (size_t i=0; i<m.size(); ++i)
    {
      tits::TitsElt a=
	(m.size()==1 // use backtrack only in case of multiple minimal classes
	 ? square_class_base.grading_seed(G,rf,m[i])
	 : square_class_base.backtrack_seed(G,rf,m[i])
	 );

      fiber_data.reduce(a);

      size_t k=elt_hash.match(a);
      assert(k==info.size()); // this KGB element should be new
      const weyl::TwistedInvolution& tw = elt_hash[k].tw();
      unsigned int l = G.twistedWeylGroup().involutionLength(tw);

      // add additional information (length, involution) for this KGB element
      KGB_base::add_element(l,tw);
    } // |for (i)|
  } // |if(traditional)|

  for (KGBElt x=0; x<elt_hash.size(); ++x) // loop makes |elt_hash| grow
  {
    // extend by cross actions
    const tits::TitsElt& current = elt_pool[x];

    for (weyl::Generator s=0; s<rank; ++s)
    {
      tits::TitsElt a = current; d_base->basedTwistedConjugate(a,s);
      fiber_data.reduce(a);

      int lc=
	a.tw()==current.tw() ? 0 : weylGroup().length_change(s,current.w());

      // now find the Tits element in |elt_hash|, or add it if new
      KGBElt child = elt_hash.match(a);
      if (child==info.size()) // add a new Tits element
      {
	unsigned int l = info[x].length+lc;
	KGB_base::add_element(l,a.tw());
      }

      // set cross links for |x|
      data[s][x].cross_image = child;

      if (lc!=0) // then set complex status, and whether ascent or descent
      {
	info[x].status.set(s,gradings::Status::Complex);
	info[x].desc.set(s,lc<0);
      }
      else if (weylGroup().hasDescent(s,current.w())) // real
      {
	assert(child==x);
	info[x].status.set(s,gradings::Status::Real);
	info[x].desc.set(s); // real roots are always descents
      }
      else // imaginary
      {
	info[x].status.set_imaginary
	  (s,d_base->simple_grading(current,s));
	info[x].desc.reset(s); // imaginary roots are never descents

	if (info[x].status[s] == gradings::Status::ImaginaryNoncompact)
	{
	  // Cayley-transform |current| by $\sigma_s$
	  tits::TitsElt a = current; d_base->Cayley_transform(a,s);
	  assert(titsGroup().length(a)>titsGroup().length(current));
	  fiber_data.reduce(a); // subspace has grown, mod out new subspace

	  KGBElt child = elt_hash.match(a);
	  if (child==info.size()) // add a new Tits element
	  {
	    unsigned int l = info[x].length+1; // length goes up
	    KGB_base::add_element(l,a.tw());
	  }
	  // add new Cayley link
	  data[s][x].Cayley_image = child;

	} // if |ImaginaryNoncompact|

      } // complex/real/imaginary disjunction

    } // |for(s)|

  } // |for (x)|

  // now sort and export to |tits|
  setutils::Permutation a(size),a1(size); // reorderd indices, and inverse
  { // sort involutions (rather than KGB elements, for efficiency)
    setutils::Permutation p(fiber_data.n_involutions(),1);
    invol_compare i_cmp(fiber_data,G.weylGroup());
    std::sort(p.begin(),p.end(),i_cmp); // all entries are distinct
    setutils::Permutation p1(p,-1); // inverse of |p|: where involutions went

    // compute number of strict predecessors each |fiber_data.involution(i)|
    std::vector<size_t> pred(p1.size()+1,0);
    for (KGBElt x=0; x<size; ++x)
      ++pred[1+p1[fiber_data.seq_no(elt_pool[x].w())]];
    for (size_t s=0,i=1; i<pred.size(); ++i) // cumulate
      s = pred[i] += s;

    for (KGBElt x=0; x<size; ++x) // finally compute KGB ordering positions
    {
      size_t i = p1[fiber_data.seq_no(elt_pool[x].w())];
      a[a1[x] = pred[i]++]=x;
    }
  }

  for (weyl::Generator s=0; s<rank; ++s)
  {
    a.pull_back(data[s]).swap(data[s]);
    for (KGBElt x=0; x<size; ++x)
    {
      data[s][x].cross_image = a1[data[s][x].cross_image];
      if (data[s][x].Cayley_image != UndefKGB)
	data[s][x].Cayley_image = a1[data[s][x].Cayley_image];
    }
  }

  a.pull_back(info).swap(info); // only permute here, no data to renumber

  left_torus_part.reserve(size);
  for (KGBElt x=0; x<size; ++x) // pull back through |a| while constructing
    left_torus_part.push_back(titsGroup().left_torus_part(elt_pool[a[x]]));

  Cartan.reserve(traditional
		 ? GR.numInvolutions()
		 : G.numInvolutions(Cartan_classes));
  first_of_tau.reserve(Cartan.capacity()+1);

  for (KGBElt x=0; x<size; ++x)
  {
    unsigned int h=inv_hash.match(involution(x));
    if (h==first_of_tau.size()) // then twisted involution was not seen before
    {
      first_of_tau.push_back(x); // push KGB element as first for involution
      Cartan.push_back(fiber_data.cartanClass(involution(x)));
    }
    else // check that all non-new involutions are equal to most recent one
      assert(h==first_of_tau.size()-1);
  }
  first_of_tau.push_back(size); // push terminating number


  return size;
} // |KGB::generate|

// find Cartan class: locate |x| within |first_of_tau|, then index |Cartan|
size_t KGB::Cartan_class(KGBElt x) const
{
  size_t i = involution_index(x); // index of first larger entry, minus one
  assert(i<Cartan.size());
  return Cartan[i];
}

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



/*

        The subsys_KGB class, constructor

*/


subsys_KGB::subsys_KGB(const KGB_base& kgb,
		       const subdatum::SubSystem& sub,
		       const weyl::TwistedWeylGroup& sub_W,
		       kgb::KGBElt x)
: KGB_base(sub_W)
, in_parent()
, in_child(kgb.size(),UndefKGB)
{
  assert(sub_W.rank()==sub.rank());

  size_t sub_rank = sub.rank();
  const bitset::RankFlags full(constants::lMask[sub_rank]);

  { // modify |x|, descending to minimal element for |subsys|
    weyl::Generator s;
    do
    {
      for(s=0; s<sub_rank; ++s)
      {
	KGBElt xx=kgb.cross_act(sub.to_simple(s),x);
	if (kgb.isDescent(sub.simple(s),xx))
	{
	  if (kgb.status(sub.simple(s),xx)==gradings::Status::Real)
	    xx=kgb.inverseCayley(sub.simple(s),xx).first;
	  else // complex descent
	    xx=kgb.cross(sub.simple(s),xx);
	  x = kgb.cross_act(xx,sub.to_simple(s));
	  break;
	} // |if(isDescent)|
      } // |for(s)|
    } while(s<sub_rank); // loop until no descents found in |subsys|
  }

  // set up mapping tables in two directions
  in_parent.push_back(x); in_child[x]=0; // map element 0 to |x|
  // value in |x| is now no longer needed; variable is reused below

  { // get elements at the fundamental fiber
    weyl::TwistedInvolution triv;
    first_of_tau.push_back(0); // start of fundamental fiber
    inv_hash.match(triv); // set initial (trivial) twisted involution

    add_element(0,inv_pool[0]); // create in base; length, tw

    /* Generate fundamental fiber. We use that cross actions by imaginary,
       non-compact, simple roots suffice, but without knowing a nice proof. */

    for (KGBElt cur=0; cur<in_parent.size(); ++cur) // |in_parent| grows
    {
      x = in_parent[cur];
      for (weyl::Generator s=0; s<sub_rank; ++s)
 	if (sub_W.twisted(s)==s) // imaginary
 	{
	  KGBElt y=kgb::cross(kgb,x,sub.simple(s),sub.to_simple(s));
 	  assert(kgb.involution(y)==kgb.involution(x)); // stay in fiber
	  if (in_child[y]==UndefKGB)
	  {
	    in_child[y]=in_parent.size();
	    in_parent.push_back(y); // map new element to y
	    add_element(0,triv); // create in base; length, tw
	  }
	  // no need so define a link here, code below takes care of this
 	}
    } // |for (cur)|
    first_of_tau.push_back(in_parent.size()); // end of fundamental fiber
  }



  // now generate fibers at other twisted involutions
  for (size_t i=0; i<inv_pool.size(); ++i) // inv_pool grows
  {
    const weyl::TwistedInvolution source=inv_pool[i]; // keep a copy for speed
    unsigned int l=length(first_of_tau[i]); // length of any |x| below

    for (weyl::Generator s=0; s<sub_rank; ++s)
    {
      unsigned int k = inv_hash.match(sub_W.twistedConjugated(source,s));
      assert(k<first_of_tau.size()); // at most one new twisted involution
      weyl::TwistedInvolution tw=inv_pool[k];
      int d = i==k ? 0 : sub_W.weylGroup().length_change(s,source);

      // generate cross action neighbours if needed, and links to them
      for (KGBElt cur=first_of_tau[i]; cur<first_of_tau[i+1]; ++cur) // again
      {
	x = in_parent[cur];
	KGBElt xx = kgb.cross_act(sub.to_simple(s),x);
	KGBElt y = kgb.cross_act(kgb.cross(sub.simple(s),xx),sub.to_simple(s));
	if (in_child[y]==UndefKGB)
	{
	  in_child[y]=in_parent.size(); // number |y| as new element in child
	  in_parent.push_back(y); // and map this new element to |y| in parent
	  add_element(l+d,tw); // create; length, tw
	}
	data[s][cur].cross_image=in_child[y]; // set link from |x| to child
	info[cur].status.set(s,kgb.status(sub.simple(s),xx));
	info[cur].desc.set(s,kgb.isDescent(sub.simple(s),xx));

	if (status(s,cur)==gradings::Status::ImaginaryNoncompact) // do Cayley
	{
	  k=inv_hash.match(sub_W.prod(s,source));
	  assert(k<first_of_tau.size()); // at most one new twisted involution
	  tw=inv_pool[k]; // redefine

	  y= kgb.cross_act(kgb.cayley(sub.simple(s),xx),sub.to_simple(s));
	  if (in_child[y]==UndefKGB)
	  {
	    in_child[y]=in_parent.size(); // number |y| as new element in child
	    in_parent.push_back(y); // and map this new element to |y| in parent
	    add_element(l+1,tw); // create; length, tw
	  }

	  // make links
	  data[s][cur].Cayley_image=in_child[y];
	  if (data[s][in_child[y]].inverse_Cayley_image.first==UndefKGB)
	    data[s][in_child[y]].inverse_Cayley_image.first=cur;
	  else
	    data[s][in_child[y]].inverse_Cayley_image.second=cur;
	} // |if(ImaginaryNoncompact)|
      } // |for(x)|
      if (k==first_of_tau.size()-1) // a new twisted involution was visited
      {
	assert(in_parent.size()>first_of_tau.back()); // new fiber non empty
	first_of_tau.push_back(in_parent.size()); // then signal end of fiber
      }
    } // |for (s)|

  } // |for(i)| (loop over wisted involutions)


} // |subsys_KGB::subsys_KGB|

std::ostream& subsys_KGB::print(std::ostream& strm, KGBElt x) const
{
  int width = ioutils::digits(in_child.size()-1,10ul);
  return strm << "(=" << std::setw(width)<< in_parent[x] << ')';
}


/*****************************************************************************

            Chapter II -- The auxiliary classes, methods

******************************************************************************/

/*    II a. |GlobalFiberData|  */

// precompute data for all involutions, given a hash table |h| enumerating them
// data: Cartan class, projector of torus parts and $\check\rho_{im}$
GlobalFiberData::GlobalFiberData
  (complexredgp::ComplexReductiveGroup& G,
   hashtable::HashTable<weyl::TI_Entry,unsigned int>& h)
  : hash_table(h)
  , info(h.size())
  , refl(G.semisimpleRank())
{
  const rootdata::RootDatum& rd=G.rootDatum();
  const weyl::TwistedWeylGroup& W = G.twistedWeylGroup();

  for (weyl::Generator s=0; s<refl.size(); ++s)
    // get automorphism of lattice $X_*$ given by generator $s$
    // reflection map will be used as automorphism of $X_*\tensor\Q/X_*$
    refl[s] = rd.simple_reflection(s).transposed();

  std::vector<size_t> to_do(h.size()); // this really only reserves the space
  bitmap::BitMap seen(h.size()); // flags involution numbers in |to_do| array

  for (size_t cn=0; cn<G.numCartanClasses(); ++cn)
  {
    const cartanclass::CartanClass& cc = G.cartan(cn);
    size_t first = h.find(G.twistedInvolution(cn));

    { // store data for canonical twisted involution |i| of Cartan class |cn|
      latticetypes::LatticeMatrix A =cc.involution().transposed();
      for (size_t i=0; i<G.rank(); ++i)
	A(i,i) += 1;
      // now $A=\theta_x^t+1$, a matrix whose kernel is $(X_*)^{-\theta_x^t}$

      info[first]=inv_info
	(cn,
	 lattice::row_saturate(A),
	 rd.dual_twoRho(cc.imaginaryRootSet()),
	 rootdata::RootSet(rd.numRoots(),cc.simpleImaginary()),
	 rootdata::RootSet(rd.numRoots(),cc.simpleReal()));
    }

    to_do.assign(1,first); // reset list to singleton containing |first|
    seen.insert(first); // we don't bother to remove old members from |seen|

    for (size_t i=0; i<to_do.size(); ++i) // |to_do| grows during the loop
    {
      const inv_info& cur = info[to_do[i]];
      for (weyl::Generator s=0; s<G.semisimpleRank(); ++s)
      {
	weyl::TwistedInvolution tw = W.twistedConjugated(h[to_do[i]],s);
	size_t k = h.find(tw);
	assert (k!=h.empty); // since |h| is closed under twisted conjugation
	if (seen.isMember(k))
	  continue;
	seen.insert(k);
	to_do.push_back(k);
	info[k]=inv_info
	  (cn,
	   cur.proj*refl[s], // apply $refl[s]^{-1}$ to old kernel
	   refl[s].apply(cur.check_2rho_imag),
	   rd.simple_root_permutation(s).renumbering(cur.simple_imag),
	   rd.simple_root_permutation(s).renumbering(cur.simple_real));
      } // |for (s)|
    } // |for (i)|

  } // |for(cn)|
} // |GlobalFiberData::GlobalFiberData|

GlobalFiberData::GlobalFiberData
(const subdatum::SubSystem& sub,
   hashtable::HashTable<weyl::TI_Entry,unsigned int>& h)
  : hash_table(h)
  , info()
  , refl(sub.rank())
{

  for (weyl::Generator s=0; s<refl.size(); ++s)
    // get automorphism of lattice $X_*$ given by generator $s$
    // reflection map will be used as automorphism of $X_*\tensor\Q/X_*$
    refl[s] = sub.parent_datum().root_reflection(sub.parent_nr_simple(s));
}

void GlobalFiberData::add_class(const subdatum::SubSystem& sub,
				const tits::GlobalTitsGroup& Tg,
				const weyl::TwistedInvolution& tw)
{
  const rootdata::RootDatum pd = sub.parent_datum();
  const weyl::TwistedWeylGroup& W = Tg;

  const size_t prev_inv = hash_table.size();
  if (hash_table.match(tw)<prev_inv)
    return; // all involutions present in hash table are fully catered for

  std::vector<std::pair<size_t,weyl::Generator> > history; // record ancestors
  history.reserve(0x100); // avoid some initial reallocations for efficiency

  cartanclass::InvolutionData id = // involution data in subsystem, for |tw|
    cartanclass::InvolutionData::build(sub,Tg,tw);

  latticetypes::LatticeMatrix A = Tg.involution_matrix(tw); // parent side
  for (size_t i=0; i<A.numRows(); ++i)
    A(i,i) -= 1; // now $A=\theta-1$, a matrix whose kernel is $(X^*)^\theta$

  rootdata::RootSet imaginary(pd.numRoots()),
    simple_imaginary(pd.numRoots()), simple_real(pd.numRoots());

  // convert from numbering of subsystem to parent numbering
  for (rootdata::RootSet::iterator it=id.imaginary_roots().begin(); it(); ++it)
    imaginary.insert(sub.parent_nr(*it));
  for (size_t i=0; i<id.imaginary_rank(); ++i)
    simple_imaginary.insert(sub.parent_nr(id.imaginary_basis(i)));
  for (size_t i=0; i<id.real_rank(); ++i)
    simple_real.insert(sub.parent_nr(id.real_basis(i)));

  // saturate by conjugation before reserving space and filling |info|
  for (size_t i=prev_inv; i<hash_table.size(); ++i) // hash table grows
  {
    size_t last = hash_table.size()-1; // |hash_table| contains at least |tw|

    for (weyl::Generator s=0; s<W.rank(); ++s)
    {
      weyl::TwistedInvolution ti = W.twistedConjugated(hash_table[i],s);
      size_t k = hash_table.match(ti);
      if (k>last) // a new twisted involution
      {
	history.push_back(std::make_pair(i,s));
	last=k;
      }
    } // |for(s)|
  } // |for(i<hash_table.size())|

    // reserve table space for new Cartan class
  info.reserve(hash_table.size());

  // create initial element
  info.push_back(inv_info(info.empty() ? 0 : info.back().Cartan+1,
			  lattice::row_saturate(A),
			  pd.twoRho(imaginary), // |pd| is dual w.r.t. |Tg|
			  simple_imaginary, // but |id| built for |Tg|, no dual
			  simple_real));

  // do remaining elements using history
  for (size_t k=prev_inv + 1; k<hash_table.size(); ++k)
  {
    size_t i=history[k-(prev_inv + 1)].first;
    assert(i<k); // this should refer to earlier element
    const inv_info& cur = info[i];
    weyl::Generator s=history[k-(prev_inv + 1)].second;
    info.push_back(inv_info
		   (cur.Cartan,   // same Cartan class
		    cur.proj*refl[s], // apply $refl[s]^{-1}$ to old kernel
		    refl[s].apply(cur.check_2rho_imag),
		    pd.simple_root_permutation(s).renumbering(cur.simple_imag),
		    pd.simple_root_permutation(s).renumbering(cur.simple_real)
		    ));
  }
  assert(info.size()==hash_table.size());
}


//  this demonstrates what the |proj matrices can be used for
bool GlobalFiberData::equivalent(const tits::GlobalTitsElement& x,
				 const tits::GlobalTitsElement& y) const
{
  unsigned int k = hash_table.find(x.tw());
  assert (hash_table.find(y.tw())==k);
  latticetypes::RatWeight t= (x.torus_part()-y.torus_part()).as_rational();
  latticetypes::LatticeElt p = info[k].proj.apply(t.numerator());

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
  assert (k!=hash_table.empty);
  latticetypes::RatLatticeElt t = x.torus_part().as_rational();
  latticetypes::LatticeElt p = info[k].proj.apply(t.numerator());

  // reduce modulo integers and return
  for (size_t i=0; i<p.size(); ++i)
    p[i]= intutils::remainder(p[i],t.denominator());
  return latticetypes::RatLatticeElt(p,t.denominator()).normalize();
}

  tits::GlobalTitsElement GlobalFiberData::imaginary_cross
  (const rootdata::RootDatum& dual_rd, // dual for pragmatic reasons
     rootdata::RootNbr alpha, // any imaginary root
     tits::GlobalTitsElement a) const
{
  const latticetypes::Weight& v = dual_rd.coroot(alpha); // imaginary
  if (a.torus_part().negative_at(v) == (at_rho_imaginary(v,a.tw())%2==0))
    a.torus_part() += tits::TorusPart(dual_rd.root(alpha));
  return a;
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
  , involution_length()
  , Cartan_class()
{
  { // dimension |data|, |elememt_length|, and |Cartan_class|
    size_t n_inv=G.numInvolutions(Cartan_classes);
    data.reserve(n_inv); Cartan_class.reserve(n_inv);
  }

  const rootdata::RootDatum& rd = G.rootDatum();
  const weyl::TwistedWeylGroup& tW = G.twistedWeylGroup();

  std::vector<latticetypes::BinaryMap> refl(G.semisimpleRank());
  for (weyl::Generator s=0; s<refl.size(); ++s)
    // get endomorphism of weight lattice $X$ given by generator $s$
    // reflection map is induced vector space endomorphism of $X_* / 2X_*$
    refl[s] = latticetypes::BinaryMap(rd.simple_reflection(s).transposed());

  for (bitmap::BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
  {
    size_t cn=*it;
    const weyl::TwistedInvolution& canonical = G.twistedInvolution(cn);
    size_t i = hash_table.match(canonical);
    assert(i==data.size()); // this twisted involution should be new

    { // store data for canonical twisted involution of Cartan class |cn|
      const latticetypes::LatticeMatrix &q = G.cartan(cn).involution();
      data.push_back(tits::fiber_denom(q)); // compute subspace $I$
      involution_length.push_back(tW.involutionLength(canonical));
      Cartan_class.push_back(cn); // record number of Cartan class
    }


    // now generate all non-canonical twisted involutions for Cartan class
    for ( ; i<data.size(); ++i) // |data.size()|  increases during the loop
      for (size_t s=0; s<G.semisimpleRank(); ++s)
      {
	weyl::TwistedInvolution stw = tW.twistedConjugated(pool[i],s);
	if (hash_table.match(stw)==data.size()) // then |stw| is new
	{
	  data.push_back(data[i]);     // start with copy of source subspace
	  data.back().apply(refl[s]);  // modify according to cross action used
	  involution_length.push_back(involution_length[i] +
				      (tW.hasDescent(s,pool[i]) ? -1 : +1));
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



/*****************************************************************************

        Chapter III -- Functions declared in kgb.h

******************************************************************************/

// general cross action in (non simple) root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt cross(const KGB_base& kgb, KGBElt x,
	     weyl::Generator s, weyl::WeylWord ww)
{
  return kgb.cross_act(kgb.cross(s,kgb.cross_act(ww,x)),ww);
}

// general Cayley transform in (non simple) non-compact imaginary root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt Cayley (const KGB_base& kgb, KGBElt x,
	       weyl::Generator s, weyl::WeylWord ww)
{
  return kgb.cross_act(kgb.cayley(s,kgb.cross_act(ww,x)),ww);
}

// general inverse Cayley transform (choice) in (non simple) real root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt inverse_Cayley (const KGB_base& kgb, KGBElt x,
		       weyl::Generator s, weyl::WeylWord ww)
{
  return kgb.cross_act(kgb.inverseCayley(s,kgb.cross_act(ww,x)).first,ww);
}


/*****************************************************************************

        Chapter IV -- Functions local to kgb.cpp

******************************************************************************/

namespace {

/*!
  \brief Puts in |Hasse| the Hasse diagram of the Bruhat ordering on |kgb|.

  Explanation: this is the closure ordering of orbits. We use the algorithm
  from Richardson and Springer.
*/
void makeHasse(std::vector<set::SetEltList>& Hasse, const KGB_base& kgb)
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

    for (set::SetEltList::const_iterator
	   it=Hasse[sx].begin(); it!= Hasse[sx].end(); ++it)
    {
      KGBElt z = *it; // element below |sx| in Bruhat order
      switch (kgb.status(s,z))
      {
      case gradings::Status::ImaginaryNoncompact:
	h_x.insert(kgb.cayley(s,z));
	break;
      case gradings::Status::Complex:
	if (not kgb.isDescent(s,z))
	  h_x.insert(kgb.cross(s,z)); // complex ascent
	break;
      default: break;
      }
    }

    Hasse[x].assign(h_x.begin(),h_x.end());
  } // for |x|

}

} // namespace

} // namespace kgb
} // namespace atlas
