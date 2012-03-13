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
#include <iomanip>
#include <algorithm>

#include "arithmetic.h"	// functions
#include "bruhat.h"
#include "cartanclass.h"// |CartanClass| methods
#include "complexredgp.h"
#include "lattice.h"	// |row_saturate|
#include "realredgp.h"  // |KGB| constructor
#include "rootdata.h"
#include "subsystem.h"
#include "tits.h"
#include "weyl.h"
#include "involutions.h" // for |InvolutionTable|

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

void makeHasse(std::vector<set::EltList>&, const KGB_base&);

} // |namespace|



size_t KGB_elt_entry::hashCode(size_t modulus) const
{
  unsigned long d= fingerprint.denominator();
  const int_Vector& num=fingerprint.numerator();
  size_t h=tw.hashCode(modulus);
  for (size_t i=0; i<num.size(); ++i)
    h=((h*d)+num[i])&(modulus-1);
  return h;
}

bool KGB_elt_entry::operator !=(const KGB_elt_entry& x) const
{ return tw!=x.tw or fingerprint!=x.fingerprint; } // ignore repr

KGB_elt_entry::KGB_elt_entry (const RatWeight& f,
			      const GlobalTitsElement& y)
  : t_rep(y.torus_part()), tw(y.tw()), fingerprint(f) {}

GlobalTitsElement KGB_elt_entry::repr() const
{ return GlobalTitsElement(t_rep,tw); }



/*

        The |KGB_base| class implementation, public methods

*/

const RootDatum& KGB_base::rootDatum() const { return G.rootDatum(); }
const WeylGroup& KGB_base::weylGroup() const { return G.weylGroup(); }
const TwistedWeylGroup& KGB_base::twistedWeylGroup() const
  { return G.twistedWeylGroup(); }

KGBElt KGB_base::cross(const WeylWord& ww, KGBElt x) const
{
  for (size_t i=ww.size(); i-->0;)
    x = cross(ww[i],x);
  return x;
}

KGBElt KGB_base::cross(KGBElt x, const WeylWord& ww) const
{
  for (size_t i=0; i<ww.size(); ++i)
    x = cross(ww[i],x);
  return x;
}

KGBEltPair KGB_base::tauPacket(const TwistedInvolution& w) const
{
  unsigned int i=inv_hash.find(w);
  if (i==inv_hash.empty)
    return KGBEltPair(0,0);
  return KGBEltPair(first_of_tau[i],first_of_tau[i+1]);
}

size_t KGB_base::packet_size(const TwistedInvolution& w) const
{
  unsigned int i=inv_hash.find(w);
  if (i==inv_hash.empty)
    return 0;
  return first_of_tau[i+1]-first_of_tau[i];
}

// compute Cartan class using |G.involution_table()|
CartanNbr KGB_base::Cartan_class(KGBElt x) const
{
  const Cartan_orbits& i_tab = G.involution_table();
  // compute passing by involution index, |TwistedInvolution|, |InvolutionNbr|
  return i_tab.Cartan_class(involution(x));
}

unsigned int KGB_base::length(KGBElt x) const
{
  const Cartan_orbits& i_tab = G.involution_table();
  return i_tab.length(involution(x));
}

const WeightInvolution& KGB_base::involution_matrix(KGBElt x) const
{
  const Cartan_orbits& i_tab = G.involution_table();
  return i_tab.matrix(involution(x));
}

InvolutionNbr KGB_base::inv_nr(KGBElt x) const
{
  const Cartan_orbits& i_tab = G.involution_table();
  return i_tab.nr(involution(x));
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
void KGB_base::add_element()
{
  info.push_back(EltInfo());
  for (size_t s=0; s<data.size(); ++s)
    data[s].push_back(KGBfields());
}


/*

        The |global_KGB| class implementation, public methods

*/

// create structure incorporating all KGB structures for a given inner class
global_KGB::global_KGB(ComplexReductiveGroup& G)
  : KGB_base(G,G.semisimpleRank())
  , Tg(G) // construct global Tits group as subobject
  , elt()
{
  G.involution_table().add(G,~ BitMap(G.numCartanClasses()));
  generate_involutions(G.numInvolutions());

  size_t size = G.global_KGB_size();
  elt.reserve(size);
  KGB_base::reserve(size);

  { // get elements at the fundamental fiber, for "all" square classes
    const SmallSubquotient& fg = G.fundamental().fiberGroup();
    first_of_tau.push_back(0); // start of fundamental fiber

    // get number of subsets of generators, and run through them
    unsigned long n= 1ul << Tg.square_class_generators().size();
    for (unsigned long c=0; c<n; ++c)
    {
      Grading gr=bitvector::combination
	(Tg.square_class_generators(),RankFlags(c));
      RatWeight rw (G.rank());
      for (Grading::iterator // for flagged (imaginary) simple roots
	     it=gr.begin(); it(); ++it)
	rw += G.rootDatum().fundamental_coweight(*it); // a sum of f. coweights

      for (unsigned long i=0; i<fg.size(); ++i)
      {
	TorusElement t=y_values::exp_pi(rw);
	t += fg.fromBasis // add |TorusPart| from fiber group; bits from |i|
	  (SmallBitVector(RankFlags(i),fg.dimension()));
	elt.push_back(GlobalTitsElement(t));
	add_element(); // create in base
      } // |for (i)|
    } // |for (c)|
    first_of_tau.push_back(elt.size()); // end of fundamental fiber
  }

  generate(size);

}

global_KGB::global_KGB(ComplexReductiveGroup& G,
		       const GlobalTitsElement& x)
  : KGB_base(G,G.semisimpleRank())
  , Tg(G) // construct global Tits group as subobject
  , elt()
{
  generate_involutions(G.numInvolutions());
  assert(Tg.is_valid(x)); // unless this holds, we cannot hope to succeed

  Cartan_orbits& i_tab = G.involution_table();
  const RootDatum& rd = G.rootDatum();

  i_tab.add(G,~ BitMap(G.numCartanClasses())); // generate all

  GlobalTitsElement a=x; // start at an element that we certainly want
  weyl::Generator s;
  while ((s=Tg.weylGroup().leftDescent(a.tw()))<rank())
    if (Tg.hasTwistedCommutation(s,a.tw()))
      Tg.do_inverse_Cayley(s,a);
    else
      Tg.cross_act(s,a);

  assert(Tg.is_valid(a)); // check that we still have a valid element
  assert (a.tw()==TwistedInvolution()); // we are at fundamental fiber

  { // get elements at the fundamental fiber
    InvolutionNbr inv = i_tab.nr(a.tw());
    first_of_tau.push_back(0); // start of fundamental fiber
    KGB_elt_entry::Pooltype elt_pool;
    hashtable::HashTable<KGB_elt_entry,unsigned long> elt_hash(elt_pool);

    elt_hash.match(i_tab.x_pack(a));

    // generating reflections are for simple-imaginary roots at |inv|
    const RootNbrList& gen_root = i_tab.imaginary_basis(inv);
    for (size_t i=0; i<elt_hash.size(); ++i) // |elt_hash| grows during loop
    {
      add_element(); // create in base
      for (size_t k=0; k<gen_root.size(); ++k)
	elt_hash.match(i_tab.x_pack
	  (Tg.cross(rd.reflectionWord(gen_root[k]),elt_hash[i].repr())));
    }

    first_of_tau.push_back(elt_hash.size()); // end of fundamental fiber

    elt.reserve(elt_hash.size()); // now copy elements from hash table to |elt|
    for (size_t i=0; i<elt_hash.size(); ++i)
      elt.push_back(elt_hash[i].repr());
  } // got elements at the fundamental fiber

  generate(0); // complete element generation, no predicted size
} // |global_KGB::global_KGB|

bool global_KGB::compact(RootNbr n, // assumed imaginary at |a|
			 const GlobalTitsElement& a) const
{
  const Cartan_orbits& i_tab = G.involution_table();
  TorusElement t=a.torus_part();
  t += i_tab.check_rho_imaginary(i_tab.nr(a.tw()));

  return not t.negative_at(rootDatum().root(n)); //  whether compact
}

KGBElt global_KGB::lookup(const GlobalTitsElement& a) const
{
  const Cartan_orbits& i_tab = G.involution_table();

  KGBEltPair p = tauPacket(a.tw());
  for (KGBElt x=p.first; x<p.second; ++x)
    if (i_tab.x_equiv(element(x),a))
      return x;
  return size(); // report failure
}

std::ostream& global_KGB::print(std::ostream& strm, KGBElt x) const
{
  RatWeight t = torus_part(x).log_2pi();
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

  const TwistedWeylGroup& W = twistedWeylGroup();

  InvolutionNbr end_length=0; // end of the length-interval under construction
  inv_hash.match(TwistedInvolution()); // set initial element

  while (end_length<inv_hash.size())
  {
    const InvolutionNbr start_length=end_length; // start new interval here
    end_length = inv_pool.size(); // and run until current end

    for (weyl::Generator s=0; s<W.rank(); ++s)
      for (InvolutionNbr i=start_length; i<end_length; ++i)
      {
	const weyl::TI_Entry& parent=inv_pool[i];
	if (W.hasTwistedCommutation(s,parent))
	  inv_hash.match(W.prod(s,parent));
	else
	  inv_hash.match(W.twistedConjugated(parent,s));
      } // |for(i)|; |for(s)|;
  } // while length interval non-empty

  assert(inv_pool.size()==n);
} // |global_KGB::generate_involutions|

void global_KGB::generate(size_t predicted_size)
{
  const Cartan_orbits& i_tab = G.involution_table();
  const TwistedWeylGroup& W = Tg; // for when |GlobalTitsGroup| is not used

  KGB_elt_entry::Pooltype elt_pool; elt_pool.reserve(predicted_size);
  hashtable::HashTable<KGB_elt_entry,KGBElt> elt_hash(elt_pool);

  KGBElt end_length=0; // end of the length-interval under construction

  for (size_t i=0; i<elt.size(); ++i)
    elt_hash.match(i_tab.x_pack(elt[i]));

  assert(elt_hash.size()==elt.size()); // all distinct; in fact an invariant

  while (end_length<first_of_tau.size()-1)
  {
    const KGBElt start_length=end_length; // start new interval here
    end_length = first_of_tau.size()-1; // and run until current end

    for (weyl::Generator s=0; s<W.rank(); ++s)
      for (InvolutionNbr inv_nr=start_length; inv_nr<end_length; ++inv_nr)
      {
	const TwistedInvolution& tw = inv_pool[inv_nr];
	TwistedInvolution new_tw =W.twistedConjugated(tw,s);
	InvolutionNbr new_nr = inv_hash.find(new_tw);
	bool is_new = new_nr+1 >= first_of_tau.size();
	if (is_new)
	  assert(new_nr+1==first_of_tau.size()); // since we mimick generation
	bool imaginary = new_nr==inv_nr and not W.hasDescent(s,tw);

	// generate cross links
	for (KGBElt x=first_of_tau[inv_nr]; x<first_of_tau[inv_nr+1]; ++x)
	{
	  GlobalTitsElement child=elt[x]; //start out with a copy
	  int d = Tg.cross_act(s,child); // cross act; |d| is length difference
	  assert(child.tw()==new_tw);
	  KGB_elt_entry ee = i_tab.x_pack(child);
	  KGBElt k = elt_hash.match(ee);
	  if (k==elt.size()) // then new
	  {
	    assert(is_new);
	    elt.push_back(child);
	    assert(i_tab.Cartan_class(new_tw)==Cartan_class(x));
	    add_element();
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
	      (s,
	       not elt[x].torus_part().negative_at(rootDatum().simpleRoot(s)));
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
	      GlobalTitsElement child=Tg.Cayley(s,elt[x]);
	      assert(child.tw()==new_tw);
	      KGBElt k=elt_hash.match(i_tab.x_pack(child));
	      if (k==elt.size()) // then new
	      {
		elt.push_back(child);
		add_element();
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


      } // |for(inv_nr)|;  // |for(s)|
  } // while length interval non-empty

  assert(elt.size()==predicted_size or predicted_size==0);
} // |global_KGB::generate|



/*

        The KGB class, public methods

*/

KGB::KGB(RealReductiveGroup& GR,
	 const BitMap& Cartan_classes)
  : KGB_base(GR.complexGroup(),GR.complexGroup().semisimpleRank())
  , Cartan()
  , left_torus_part()
  , d_state()
  , d_bruhat(NULL)
  , d_base(NULL)
{
  ComplexReductiveGroup& G=GR.complexGroup();
  //const TitsGroup& Tg = G.titsGroup();
  size_t rank = G.semisimpleRank(); // |G.rank()| does not interest us here

  {// check |Cartan_classes|
    BitMap test_set = GR.Cartan_set();
    test_set.andnot(Cartan_classes);
    for (BitMap::iterator it=test_set.begin(); it(); ++it)
      assert(G.Cartan_ordering().below(*it).disjoint(Cartan_classes));
  }

  // make sure |G| has information about involutions for |Cartan_classes|
  Cartan_orbits& i_tab = G.involution_table();
  i_tab.add(G,Cartan_classes);

  tits::TE_Entry::Pooltype elt_pool; // of size |size()|
  hashtable::HashTable<tits::TE_Entry,KGBElt> elt_hash(elt_pool);

  size_t size = G.KGB_size(GR.realForm(),Cartan_classes);

  elt_pool.reserve(size);
  KGB_base::reserve(size);

  {
    tits::EnrichedTitsGroup square_class_base(GR);
    d_base = new TitsCoset(square_class_base);
    RealFormNbr rf=GR.realForm();
    assert(square_class_base.square()==
	   G.fundamental().central_square_class(rf));

    set::EltList mins=G.Cartan_ordering().minima(Cartan_classes);

    for (set::EltList::iterator it=mins.begin(); it!=mins.end(); ++it)
    {
      const CartanNbr cn = *it;
      TitsElt a=
	(mins.size()==1 // use backtrack only in (unused) multiple minima case
	 ? square_class_base.grading_seed(G,rf,cn)
	 : square_class_base.backtrack_seed(G,rf,cn)
	 );

      i_tab.reduce(a);
      size_t k=elt_hash.match(a);
      assert(k==info.size()); // this KGB element should be new
      KGB_base::add_element(); // add additional info for element |k|
    } // |for (it)|
  }

  for (KGBElt x=0; x<elt_hash.size(); ++x) // loop makes |elt_hash| grow
  {
    // extend by cross actions
    const TitsElt& current = elt_pool[x];
    EltInfo& my_info = info[x];

    for (weyl::Generator s=0; s<rank; ++s)
    {
      KGBfields& my_s = data[s][x];
      TitsElt a = current; d_base->basedTwistedConjugate(a,s);
      i_tab.reduce(a);

      int lc=
	a.tw()==current.tw() ? 0 : weylGroup().length_change(s,current.w());

      // now find the Tits element in |elt_hash|, or add it if new
      KGBElt child = elt_hash.match(a);
      if (child==info.size()) // add a new Tits element
	KGB_base::add_element();

      // set cross link for |x|
      my_s.cross_image = child;

      if (lc!=0) // then set complex status, and whether ascent or descent
      {
	my_info.status.set(s,gradings::Status::Complex);
	my_info.desc.set(s,lc<0);
      }
      else if (weylGroup().hasDescent(s,current.w())) // real
      {
	assert(child==x);
	my_info.status.set(s,gradings::Status::Real);
	my_info.desc.set(s); // real roots are always descents
      }
      else // imaginary
      {
	my_info.status.set_imaginary(s,d_base->simple_grading(current,s));
	my_info.desc.reset(s); // imaginary roots are never (KGB) descents

	if (my_info.status[s] == gradings::Status::ImaginaryNoncompact)
	{
	  // Cayley-transform |current| by $\sigma_s$
	  TitsElt a = current; d_base->Cayley_transform(a,s);
	  assert(titsGroup().length(a)>titsGroup().length(current));
	  i_tab.reduce(a); // subspace has grown, mod out new subspace

	  KGBElt child = elt_hash.match(a);
	  if (child==info.size()) // add a new Tits element
	    KGB_base::add_element();

	  // add new Cayley link
	  my_s.Cayley_image = child;

	} // if |ImaginaryNoncompact|

      } // complex/real/imaginary disjunction

    } // |for(s)|

  } // |for (x)|

  assert(KGB_base::size()==size);

  // now sort
  Permutation a1; // will be reordering assignment
  { // first sort involutions

    std::vector<InvolutionNbr> invs;
    invs.reserve(size); // first use needs less

    for (BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
      for (InvolutionNbr i=i_tab[*it].start; i<i_tab[*it].end(); ++i)
	invs.push_back(i);

    std::stable_sort(invs.begin(),invs.end(),i_tab.less());

    inv_pool.reserve(invs.size());
    std::transform(invs.begin(), invs.end(),back_inserter(inv_pool),
		   i_tab.as_map());
    inv_hash.reconstruct(); // now |inv_hash| numbers in sorted order

    invs.clear();
    for (KGBElt x=0; x<size; ++x) // list involution indices; cannot yet use
      invs.push_back(inv_hash.find(elt_pool[x].tw())); // |involution(x)| !

    a1 = permutations::standardize(invs,inv_pool.size(),&first_of_tau);
  }

  Permutation a(a1,-1); // |a[x]| locates KGB element that should become |x|

  // permute |data|, sending entry |x| to place |a1[x]}, renumbering internally
  for (weyl::Generator s=0; s<rank; ++s)
  {
    std::vector<KGBfields>& ds = data[s];
    a.pull_back(ds).swap(ds); // permute all the links for |s|
    for (KGBElt x=0; x<size; ++x) // now renumber internally through |a1|
    {
      KGBfields& dsx = ds[x];
      dsx.cross_image = a1[dsx.cross_image];
      if (dsx.Cayley_image != UndefKGB)
	dsx.Cayley_image = a1[dsx.Cayley_image];
    }
  }

  a.pull_back(info).swap(info); // similarly permute |info|, no renumbering

  left_torus_part.reserve(size);
  for (KGBElt x=0; x<size; ++x) // pull back through |a| while constructing
    left_torus_part.push_back(titsGroup().left_torus_part(elt_pool[a[x]]));

  Cartan.reserve(inv_pool.size());
  for (InvolutionNbr i=0; i<inv_pool.size(); ++i)
    Cartan.push_back(i_tab.Cartan_class(inv_pool[i]));

  // finally install inverse Cayley links
  for (KGBElt x=0; x<size; ++x)
    for (weyl::Generator s=0; s<rank; ++s)
    {
      KGBElt c=data[s][x].Cayley_image;
      if (c!=UndefKGB)
      {
	KGBEltPair& target=data[s][c].inverse_Cayley_image;
	if (target.first==UndefKGB) target.first=x;
	else target.second=x;
      }
    }
} // |KGB::KGB(GR,Cartan_classes,i_tab)|




KGB::~KGB() { delete d_bruhat; delete d_base; }

/******** copy, assignment and swap ******************************************/


/******** accessors **********************************************************/

TitsElt KGB::titsElt(KGBElt x) const
{ return TitsElt(titsGroup(),left_torus_part[x],involution(x)); }

size_t KGB::torus_rank() const { return titsGroup().rank(); }

// Looks up a |TitsElt| value and returns its KGB number, or |size()|
// Since KGB does not have mod space handy, must assume |a| already reduced
KGBElt KGB::lookup(const TitsElt& a, const TitsGroup& Tg) const
{
  KGBEltPair p = tauPacket(a.tw());
  tits::TorusPart t = Tg.left_torus_part(a);
  for (KGBElt x=p.first; x<p.second; ++x)
    if (Tg.left_torus_part(titsElt(x))==t)
      return x;
  return size(); // report failure
}

const poset::Poset& KGB::bruhatPoset() // this creates full poset on demand
{ return bruhatOrder().poset(); }

std::ostream& KGB::print(std::ostream& strm, KGBElt x) const
{
  return prettyprint::prettyPrint(strm,torus_part(x));
}

/*

     The KGB class, private manipulators

*/

// Construct the BruhatOrder

void KGB::fillBruhat()
{
  if (d_state.test(BruhatConstructed)) // work was already done
    return;

  std::vector<set::EltList> hd; makeHasse(hd,*this);
  BruhatOrder* bp = new BruhatOrder(hd); // pointer stored immediately: safe

  // commit
  assert(d_bruhat==NULL); // so no |delete| is needed
  d_bruhat = bp;
  d_state.set(BruhatConstructed);
}


/*!
\brief A |FiberData| object associates to each twisted involution a subspace
describing how corresponding Tits elements should be normalized.

It also records the Cartan class that each twisted involution belongs to.
*/
class FiberData
{
  const TitsGroup& Tits;
  weyl::TI_Entry::Pooltype pool;
  hashtable::HashTable<weyl::TI_Entry,unsigned int> hash_table;
  std::vector<BinaryMap> refl; // transposed simple reflections mod 2

  // the following three vectors are indexed by |hash_table| sequence numbers
  std::vector<SmallSubspace> data;
  std::vector<unsigned int> involution_length;
  std::vector<unsigned int> Cartan_class;
public:
  FiberData(ComplexReductiveGroup& G,
	    const BitMap& Cartan_classes); // for |KGB::generate|

  FiberData(const TitsGroup& Tg); // for |subsys_KGB|; start without classes

  InvolutionNbr n_involutions() const { return pool.size(); }
  const TwistedInvolution& involution(InvolutionNbr i) const { return pool[i]; }

  InvolutionNbr seq_no(const TwistedInvolution& tw) const
  { return hash_table.find(tw); }

  bool unseen(const TwistedInvolution& tw) const
  { return seq_no(tw)==hash_table.empty; }

  unsigned int length(InvolutionNbr i) const { return involution_length[i]; }

  CartanNbr cartanClass(const TwistedInvolution& tw) const
  { return Cartan_class[seq_no(tw)]; }

  void reduce(TitsElt& a) const;

  // manipulator
  void add_class(const TwistedInvolution& tw,
		 const WeightInvolution& theta);

private: // the space actually stored need not be exposed
  const SmallSubspace& mod_space(const TitsElt& a) const
  {
    size_t i=hash_table.find(a.tw());
    assert(i!=hash_table.empty);
    return data[i];
  }
  void complete_class(unsigned int initial); // extend by conjugations
}; // |class FiberData|




/*

        The subsys_KGB class, constructor

*/

// initially thought of as an extracting constructor, we currently
// do-it-ourselves, and only try to map back to the parent KGB at the end
subsys_KGB::subsys_KGB
  (const KGB& kgb,
   const subdatum::SubDatum& sub,       // subsystem is on dual side
   KGBElt x)
  : KGB_base(kgb.complexGroup(),sub.semisimple_rank())
  , in_parent()
  , Cartan()
  , torus_part()
  , parent_size(kgb.size())
{
  // no need to extend |G.involution_table()| since building |kgb| did that

  size_t sub_rank = sub.semisimple_rank();

  const TitsCoset Tc(sub,kgb.base_grading()); // constructor converts grading
  const TitsGroup& Tg = Tc.titsGroup(); // we need this several times
  const WeylGroup& W = Tc.weylGroup();

  TitsElt cur_x (Tg,kgb.torus_part(x),sub.init_twisted());

  { // modify |cur_x|, descending to minimal element for |subsys|
    weyl::Generator s;
    do
      for(s=0; s<sub_rank; ++s)
	if (W.hasDescent(s,cur_x.tw()))
	{
	  if (Tc.hasTwistedCommutation(s,cur_x.tw())) // real
	  {
	    SmallSubspace mod_space = tits::fiber_denom(sub.matrix(cur_x.tw()));
	    Tc.inverse_Cayley_transform(cur_x,s,mod_space);
	  }
	  else // complex descent
	    Tc.basedTwistedConjugate(cur_x,s);
	  break;
	} // |if(isDescent), |for(s)|
    while(s<sub_rank); // loop until no descents found in |sub|
  }

  tits::TE_Entry::Pooltype elt_pool(1,tits::TE_Entry(cur_x));
  hashtable::HashTable<tits::TE_Entry,KGBElt> elt_hash(elt_pool);

  FiberData fd(Tg);

  { // get elements at the fundamental fiber
    TwistedInvolution tw = cur_x.tw();
    fd.add_class(tw,sub.matrix(tw));
    first_of_tau.push_back(0); // start of fundamental fiber
    inv_hash.match(tw); // set initial twisted involution

    for (KGBElt cur=0; cur<elt_pool.size(); ++cur) // |elt_pool| grows
    {
      add_element(); // create |cur| in |info| array of base object
      for (weyl::Generator s=0; s<sub_rank; ++s)
 	if (not W.hasDescent(s,tw) // currently this is assured, could change
	    and Tc.hasTwistedCommutation(s,tw)) // imaginary simple reflection
 	{
	  TitsElt x = elt_pool[cur]; // take a copy
	  Tc.basedTwistedConjugate(x,s);
 	  assert(x.tw()==tw); // stay in fiber
	  fd.reduce(x);
	  elt_hash.match(x); // add to table if new
 	}
    } // |for (cur)|
    first_of_tau.push_back(info.size()); // end of fundamental fiber
  }

  // now generate fibers at other twisted involutions
  for (InvolutionNbr i=0; i<inv_pool.size(); ++i) // inv_pool grows
  {
    const TwistedInvolution source=inv_pool[i]; // keep a copy for speed

    for (weyl::Generator s=0; s<sub_rank; ++s)
    {
      InvolutionNbr k = inv_hash.match(Tg.twistedConjugated(source,s));
      const TwistedInvolution tw = inv_pool[k];
      int d = i==k ? 0 : W.length_change(s,source);

      // generate cross action neighbours if needed, and make links to them
      for (KGBElt cur=first_of_tau[i]; cur<first_of_tau[i+1]; ++cur)
      {
	TitsElt x = elt_pool[cur]; // take a copy
	Tc.basedTwistedConjugate(x,s);
	fd.reduce(x);
	KGBElt child = elt_hash.match(x);
	if (child==info.size()) // then newborn
	  add_element(); // create its |info| entry
	data[s][cur].cross_image=child; // set link from |x| to child
	if (d!=0)
	  info[cur].status.set(s,gradings::Status::Complex);
	else if (W.hasDescent(s,tw))
	  info[cur].status.set(s,gradings::Status::Real);
	else if (Tc.simple_grading(x,s))
	  info[cur].status.set(s,gradings::Status::ImaginaryNoncompact);
	else
	  info[cur].status.set(s,gradings::Status::ImaginaryCompact);

	info[cur].desc.set(s,W.hasDescent(s,tw)); // complex descent or real

	if (status(s,cur)==gradings::Status::ImaginaryNoncompact) // do Cayley
	{
	  k=inv_hash.match(W.prod(s,source));
	  assert(k<first_of_tau.size()); // at most one new twisted involution

	  const TwistedInvolution tw = inv_pool[k]; // mask previous
	  TitsElt x = elt_pool[cur];  // but do no modify their values
	  Tc.Cayley_transform(x,s);
	  if (fd.unseen(x.tw()))
	    fd.add_class(tw,sub.matrix(tw));
	  fd.reduce(x);
	  child = elt_hash.match(x);
	  if (child==info.size()) // then newborn
	    add_element(); // create its |info| entry, setting length, tw

	  // make links
	  data[s][cur].Cayley_image = child;
	  if (data[s][child].inverse_Cayley_image.first==UndefKGB)
	    data[s][child].inverse_Cayley_image.first=cur;
	  else
	    data[s][child].inverse_Cayley_image.second=cur;
	} // |if(ImaginaryNoncompact)|
      } // |for(cur)|
      if (k==first_of_tau.size()-1) // a new twisted involution was visited
      {
	assert(info.size()>first_of_tau.back()); // new fiber non empty
	first_of_tau.push_back(info.size());     // signal end of fiber |k|
      }
    } // |for (s)|

  } // |for(i)| (loop over twisted involutions)

  assert(elt_pool.size()==info.size());
  in_parent.reserve(info.size());
  Cartan.reserve(info.size()); torus_part.reserve(info.size());

  // convert twisted involutions to original form
  {
    const WeylGroup& pW= kgb.weylGroup(); // now work in parent Weyl group
    const WeylElt w_base = pW.element(sub.base_twisted_in_parent());

    for (InvolutionNbr i=0; i<inv_pool.size(); ++i)
    {
      WeylElt w = w_base;
      WeylWord ww = W.word(inv_pool[i]);
      for (size_t j=ww.size(); j-->0; ) // laboriously compute parent involution
	pW.leftMult(w,sub.reflection(ww[j]));

      inv_pool[i] = w; // overwrite by involution for parent group
      for (KGBElt x=first_of_tau[i]; x<first_of_tau[i+1]; ++x)
      {
	torus_part.push_back(Tg.left_torus_part(elt_pool[x]));
	Cartan.push_back(fd.cartanClass(elt_pool[x].tw()));
	in_parent.push_back // find the |KGBElt| value in parent for |x|
	  (kgb.lookup(TitsElt(Tg,torus_part.back(),w), // reconstruct with |w|
		      Tg)); // |Tg| is needed again to extract torus parts...
      }
    }
    inv_hash.reconstruct();
  }
} // |subsys_KGB::subsys_KGB|

std::ostream& subsys_KGB::print(std::ostream& strm, KGBElt x) const
{
  int width = ioutils::digits(parent_size-1,10ul);
  return prettyprint::prettyPrint
    (strm << '(' << std::setw(width)<< in_parent[x] << '=',torus_part[x])
    << ')';
}



/*****************************************************************************

            Chapter II -- The auxiliary classes, methods

******************************************************************************/

/*    II a. |GlobalFiberData|  */

// precompute data for all involutions, given a hash table |h| enumerating them
// data: Cartan class, projector of torus parts and $\check\rho_{im}$
GlobalFiberData::GlobalFiberData
  (ComplexReductiveGroup& G,
   hashtable::HashTable<weyl::TI_Entry,unsigned int>& h)
  : hash_table(h)
  , info(h.size())
  , refl(G.semisimpleRank())
{
  const RootDatum& rd=G.rootDatum();
  const TwistedWeylGroup& W = G.twistedWeylGroup();

  for (weyl::Generator s=0; s<refl.size(); ++s)
    // get automorphism of lattice $X_*$ given by generator $s$
    // reflection map will be used as automorphism of $X_*\tensor\Q/X_*$
    refl[s] = rd.simple_reflection(s).transposed();

  std::vector<size_t> to_do(h.size()); // this really only reserves the space
  BitMap seen(h.size()); // flags involution numbers in |to_do| array

  for (size_t cn=0; cn<G.numCartanClasses(); ++cn)
  {
    const CartanClass& cc = G.cartan(cn);
    size_t first = h.find(G.twistedInvolution(cn));

    { // store data for canonical twisted involution |i| of Cartan class |cn|
      CoweightInvolution A =cc.involution().transposed();
      for (size_t i=0; i<G.rank(); ++i)
	A(i,i) += 1;
      // now $A=\theta_x^t+1$, a matrix whose kernel is $(X_*)^{-\theta_x^t}$

      info[first]=inv_info
	(cn,
	 lattice::row_saturate(A),
	 rd.dual_twoRho(cc.imaginaryRootSet()),
	 cc.simpleImaginary(),
	 cc.simpleReal());
    }

    to_do.assign(1,first); // reset list to singleton containing |first|
    seen.insert(first); // we don't bother to remove old members from |seen|

    for (size_t i=0; i<to_do.size(); ++i) // |to_do| grows during the loop
    {
      const inv_info& cur = info[to_do[i]];
      for (weyl::Generator s=0; s<G.semisimpleRank(); ++s)
      {
	TwistedInvolution tw = W.twistedConjugated(h[to_do[i]],s);
	size_t k = h.find(tw);
	assert (k!=h.empty); // since |h| is closed under twisted conjugation
	if (seen.isMember(k))
	  continue;
	seen.insert(k);
	to_do.push_back(k);
	info[k]=inv_info
	  (cn,
	   cur.proj*refl[s], // apply $refl[s]^{-1}$ to old kernel
	   refl[s]*cur.check_2rho_imag,
	   rd.simple_root_permutation(s).renumbering(cur.simple_imag),
	   rd.simple_root_permutation(s).renumbering(cur.simple_real));
      } // |for (s)|
    } // |for (i)|

  } // |for(cn)|
} // |GlobalFiberData::GlobalFiberData|


GlobalFiberData::GlobalFiberData
  (const GlobalTitsGroup& Tg,
   hashtable::HashTable<weyl::TI_Entry,unsigned int>& h)
  : hash_table(h)
  , info()
  , refl(Tg.semisimple_rank())
{
  for (weyl::Generator s=0; s<refl.size(); ++s)
  {
    refl[s] = CoweightInvolution(Tg.rank()); // identity matrix
    const Weight& alpha     = Tg.parent_simple_root(s);
    const Coweight& alpha_v = Tg.parent_simple_coroot(s);
    for (size_t i=0; i<Tg.rank(); ++i)
      for (size_t j=0; j<Tg.rank(); ++j)
	refl[s](i,j) -= alpha[i]*alpha_v[j];
  }
}

GlobalFiberData::GlobalFiberData
  (const SubSystem& sub,
   hashtable::HashTable<weyl::TI_Entry,unsigned int>& h)
  : hash_table(h)
  , info()
  , refl(sub.rank())
{
  PreRootDatum parent = sub.pre_root_datum();
  size_t pr = parent.rank();
  const WeightList&   root   = parent.roots();
  const CoweightList& coroot = parent.coroots();
  for (weyl::Generator s=0; s<refl.size(); ++s)
  {
    refl[s] =  CoweightInvolution(pr); // identity matrix
    const Weight& alpha     = root[s];
    const Coweight& alpha_v = coroot[s];
    for (size_t i=0; i<pr; ++i)
      for (size_t j=0; j<pr; ++j)
	refl[s](i,j) -= alpha[i]*alpha_v[j];
  }
}

InvInfo::InvInfo(const SubSystem& subsys,
		 hashtable::HashTable<weyl::TI_Entry,unsigned int>& h)
  : GlobalFiberData(subsys,h), sub(subsys), n_Cartans(0)
{}

bool InvInfo::add_involution
  (const TwistedInvolution& tw, const GlobalTitsGroup& Tg)
{
  assert(info.size()==hash_table.size());
  if (hash_table.match(tw)<info.size())
    return false; // involution is already tabulated;

  const RootDatum& pd = sub.parent_datum(); // needed for 2rho_imaginary
  WeightInvolution A = Tg.involution_matrix(tw); // generate matrix $\theta$
  InvolutionData id = sub.involution_data(A); // roots of |sub| only
  RootNbrSet imaginary=id.real_roots(); // imaginary from |sub| side
  RootNbrList simple_imaginary=id.real_basis();
  RootNbrList simple_real=id.imaginary_basis();

  // subtract identity from |A|, and row-saturate it
  for (size_t i=0; i<A.numRows(); ++i)
    A(i,i) -= 1; // now $A=\theta-1$, a matrix whose kernel is $(X^*)^\theta$
  lattice::row_saturate(A).swap(A);

  info.push_back(inv_info(n_Cartans++,
			  A,
			  pd.twoRho(imaginary), // |pd| is dual w.r.t. |Tg|
			  simple_imaginary, // and |id| nomenclature for |Tg|
			  simple_real));

  assert(info.size()==hash_table.size()); // reestablished size match
  return true;
}

bool InvInfo::add_cross_neighbor
  (const TwistedInvolution& tw,
   unsigned int old_inv,
   weyl::Generator s)
{
  assert(info.size()==hash_table.size());
  if (hash_table.match(tw)<info.size())
    return false; // involution is already tabulated;

  const RootDatum& pd = sub.parent_datum(); // needed for 2rho_imaginary
  const Permutation& sigma = pd.root_permutation(sub.parent_nr_simple(s));

  assert(old_inv<hash_table.size()); // should refer to known element
  const inv_info& cur = info[old_inv];
  info.push_back(inv_info
     (cur.Cartan,   // same Cartan class
      cur.proj*refl[s], // apply $refl[s]^{-1}$ to old kernel
      refl[s]*cur.check_2rho_imag,
      sigma.renumbering(cur.simple_imag),
      sigma.renumbering(cur.simple_real)
      ));

  assert(info.size()==hash_table.size());
  return true;
}


//  this demonstrates what the |proj| matrices can be used for
bool GlobalFiberData::equivalent(const GlobalTitsElement& x,
				 const GlobalTitsElement& y) const
{
  unsigned int k = hash_table.find(x.tw());
  if (hash_table.find(y.tw())!=k)
    return false;

  RatWeight t= (x.torus_part()-y.torus_part()).log_2pi();
  int_Vector p = info[k].proj*t.numerator();

  for (size_t i=0; i<p.size(); ++i)
    if (p[i]%t.denominator()!=0)
      return false;

  return true;
}

// computing "fingerprints" allows direct comparison without using |equivalent|
RatWeight
  GlobalFiberData::fingerprint(const GlobalTitsElement& x) const
{
  unsigned int k = hash_table.find(x.tw());
  assert (k!=hash_table.empty);
  RatWeight t = x.torus_part().log_2pi();
  int_Vector p = info[k].proj*t.numerator();

  // reduce modulo integers and return
  for (size_t i=0; i<p.size(); ++i)
    p[i]= arithmetic::remainder(p[i],t.denominator());
  return RatWeight(p,t.denominator()).normalize();
}

GlobalTitsElement GlobalFiberData::imaginary_cross
  (const RootDatum& dual_rd, // dual for pragmatic reasons
     RootNbr alpha, // any imaginary root
     GlobalTitsElement a) const
{
  const Coweight& v = dual_rd.coroot(alpha); // imaginary
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
  such that there are as many values modulo $I$ as there are in the fiber
  group, but they are not naturally in bijection, and $0$ need no be a valid
  value for $x$: the possible values are such that the Tits group element $a$
  given by $(x,w)$ satisfies $a*twisted(a)=e$, where the twisting is done
  with respect to the based Tits group (determined by the base grading).

  This constructor depends on the real form only via the set |Cartan_classes|
  that determines (limits) the set of twisted involutions to be considered.
*/
FiberData::FiberData(ComplexReductiveGroup& G,
		     const BitMap& Cartan_classes)
  : Tits(G.titsGroup())
  , pool()
  , hash_table(pool)
  , refl(G.semisimpleRank())
  , data()
  , involution_length()
  , Cartan_class()
{
  { // dimension |data|, |element_length|, and |Cartan_class|
    size_t n_inv=G.numInvolutions(Cartan_classes);
    data.reserve(n_inv); Cartan_class.reserve(n_inv);
  }

  const RootDatum& rd = G.rootDatum();
  const TwistedWeylGroup& tW = G.twistedWeylGroup();

  for (weyl::Generator s=0; s<refl.size(); ++s)
    // get endomorphism of weight lattice $X$ given by generator $s$
    // reflection map is induced vector space endomorphism of $X_* / 2X_*$
    refl[s] = BinaryMap(rd.simple_reflection(s).transposed());

  for (BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
  {
    size_t cn=*it;
    const TwistedInvolution& canonical = G.twistedInvolution(cn);
    size_t i = hash_table.match(canonical);
    assert(i==data.size()); // this twisted involution should be new

    { // store data for canonical twisted involution of Cartan class |cn|
      const WeightInvolution &q = G.cartan(cn).involution();
      data.push_back(tits::fiber_denom(q)); // compute subspace $I$
      involution_length.push_back(tW.involutionLength(canonical));
      Cartan_class.push_back(cn); // record number of Cartan class
    }

    complete_class(i);
  } // |for (it)|

  // check that the number of generated elements is as predicted
  assert(data.size()==G.numInvolutions(Cartan_classes));
  assert(Cartan_class.size()==G.numInvolutions(Cartan_classes));
}


FiberData::FiberData(const TitsGroup& Tg)
  : Tits(Tg)
  , pool()
  , hash_table(pool)
  , refl(Tg.semisimple_rank(),BinaryMap(Tg.rank(),0))
  , data()
  , involution_length()
  , Cartan_class()
{
  for (weyl::Generator s=0; s<Tg.semisimple_rank(); ++s)
    for (size_t i=0; i<Tg.rank(); ++i)
    {
      SmallBitVector v(Tg.rank(),i); // $X_*$ basis vector $e_i$
      Tg.reflect(v,s); // apply simple reflection on this torus part
      refl[s].addColumn(v);
    }
}

void FiberData::add_class(const TwistedInvolution& tw,
			  const WeightInvolution& theta)
{
  size_t cn = Cartan_class.empty() ? 0 : Cartan_class.back()+1;
  if (hash_table.match(tw)==data.size()) // then this is a new class
  {
    data.push_back(tits::fiber_denom(theta)); // compute subspace $I$
    involution_length.push_back(Tits.involutionLength(tw));
    Cartan_class.push_back(cn); // record number of Cartan class

    complete_class(data.size()-1); // argument equals |hash_table.match(tw)|
  }
}

void FiberData::complete_class(unsigned int initial)
{
  size_t cn = Cartan_class[initial];
  const TwistedWeylGroup& tW = Tits;
  // now generate all non-canonical twisted involutions for Cartan class
  for (unsigned int i=initial; i<data.size(); ++i) // |data.size()|  increases
    for (weyl::Generator s=0; s<tW.rank(); ++s)
    {
      TwistedInvolution stw = tW.twistedConjugated(pool[i],s);
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



void FiberData::reduce(TitsElt& a) const
{
  a.reduce(mod_space(a));
}



/*****************************************************************************

        Chapter III -- Functions declared in kgb.h

******************************************************************************/

// general cross action in (non simple) root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt cross(const KGB_base& kgb, KGBElt x,
	     weyl::Generator s, const WeylWord& ww)
{
  return kgb.cross(kgb.cross(s,kgb.cross(ww,x)),ww);
}

// general Cayley transform in (non simple) non-compact imaginary root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt Cayley (const KGB_base& kgb, KGBElt x,
	       weyl::Generator s, const WeylWord& ww)
{
  return kgb.cross(kgb.cayley(s,kgb.cross(ww,x)),ww);
}

// general inverse Cayley transform (choice) in (non simple) real root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt inverse_Cayley (const KGB_base& kgb, KGBElt x,
		       weyl::Generator s, const WeylWord& ww)
{
  return kgb.cross(kgb.inverseCayley(s,kgb.cross(ww,x)).first,ww);
}

// status of |alpha| in |kgb|: conjugate to simple root follow cross actions
gradings::Status::Value status(const KGB_base& kgb, KGBElt x,
			       const RootSystem& rs, RootNbr alpha)
{
  weyl::Generator s;
  while (alpha!=rs.simpleRootNbr(s=rs.find_descent(alpha)))
  {
    rs.simple_reflect_root(alpha,s);
    x=kgb.cross(s,x);
  }
  return kgb.status(s,x);
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
void makeHasse(std::vector<set::EltList>& Hasse, const KGB_base& kgb)
{
  Hasse.resize(kgb.size());

  for (KGBElt x = 0; x < kgb.size(); ++x)
  {
    std::set<KGBElt> h_x;
    const DescentSet& d = kgb.descent(x);
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

    for (set::EltList::const_iterator
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
