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
  Copyright (C) 2007-2013 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

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
#include "innerclass.h"
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
  const Ratvec_Numer_t& num=fingerprint.numerator();
  size_t h=tw.hashCode(modulus);
  for (size_t i=0; i<num.size(); ++i)
    h=(h*d)+num[i];
  return h&(modulus-1);
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

const RootDatum& KGB_base::rootDatum() const { return ic.rootDatum(); }
const WeylGroup& KGB_base::weylGroup() const { return ic.weylGroup(); }
const TwistedWeylGroup& KGB_base::twistedWeylGroup() const
  { return ic.twistedWeylGroup(); }

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
  InvolutionNbr n = ic.involution_table().nr(w);
  if (n==HashTable<weyl::TI_Entry,InvolutionNbr>::empty)
    return KGBEltPair(0,0);  // involution not (yet) generated in |G|
  inv_index i=inv_loc[n];
  if (i==-1u)
    return KGBEltPair(0,0); // involution not known at this real form
  return KGBEltPair(first_of_tau[i],first_of_tau[i+1]);
}

size_t KGB_base::packet_size(const TwistedInvolution& w) const
{
  InvolutionNbr n = ic.involution_table().nr(w);
  if (n==HashTable<weyl::TI_Entry,InvolutionNbr>::empty)
    return 0;  // involution not (yet) generated in |G|, so empty packet here
  inv_index i=inv_loc[n];
  if (i==-1u) // or involution known in |G| but not for |G_R|, empty packet
    return 0;
  return first_of_tau[i+1]-first_of_tau[i];
}

// compute Cartan class with aid of |ic.involution_table()|
CartanNbr KGB_base::Cartan_class(KGBElt x) const
{ // compute passing by involution index, |InvolutionNbr|
  return ic.involution_table().Cartan_class(inv_nr(x));
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
global_KGB::global_KGB(InnerClass& G_C, bool dual_twist)
  : KGB_base(G_C,G_C.semisimpleRank())
  , Tg(G_C) // construct global Tits group as subobject
  , elt()
{
  // fill the entire involution table of |G_C|:
  for (CartanNbr i=0; i<G_C.numCartanClasses(); ++i)
    G_C.generate_Cartan_orbit(i);

  // generate all involutions in generation order into the list |inv_nrs|
  generate_involutions(G_C.numInvolutions());

  size_t size = G_C.global_KGB_size();
  elt.reserve(size);
  KGB_base::reserve(size);

  { // get elements at the fundamental fiber, for "all" square classes
    first_of_tau.push_back(0); // start of fundamental fiber

    // get number of subsets of generators, and run through them
    unsigned long n= 1ul << Tg.square_class_generators().size();
    for (unsigned long c=0; c<n; ++c)
    {
      Grading gr=bitvector::combination
	(Tg.square_class_generators(),RankFlags(c));
      RatWeight rcw (ic.rank());
      for (Grading::iterator // for flagged (imaginary) simple roots
	     it=gr.begin(); it(); ++it)
	rcw += ic.rootDatum().fundamental_coweight(*it); // a sum of f. coweights

      for (unsigned long i=0; i<ic.fundamental_fiber_size(); ++i)
      {
	TorusElement t=y_values::exp_pi(rcw);
	t += ic.lift_from_fundamental_fiber(i); // add |TorusPart| from for |i|
	elt.push_back(GlobalTitsElement(t));
	add_element(); // create in base
      } // |for (i)|
    } // |for (c)|
    first_of_tau.push_back(elt.size()); // end of fundamental fiber
  }

  generate(size,dual_twist);

}

global_KGB::global_KGB(InnerClass& G_C,
		       const GlobalTitsElement& x, bool dual_twist)
  : KGB_base(G_C,G_C.semisimpleRank())
  , Tg(G_C) // construct global Tits group as subobject
  , elt()
{
  generate_involutions(G_C.numInvolutions());
  assert(Tg.is_valid(x)); // unless this holds, we cannot hope to succeed

  const RootDatum& rd = ic.rootDatum();

  for (CartanNbr i=0; i<G_C.numCartanClasses(); ++i)
    G_C.generate_Cartan_orbit(i);

  const Cartan_orbits& i_tab = G_C.involution_table();

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
    HashTable<KGB_elt_entry,unsigned long> elt_hash(elt_pool);

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

  generate(0,dual_twist); // complete element generation, no predicted size
} // |global_KGB::global_KGB|

/*
  Simple-imaginary roots evaluate at torus part $t$ to $1$ iff they are compact.
  To get a grading valid for all imaginary roots, we must add the half sum of
  positive imaginary coroots to $t$, flipping values at the simple-imaginaries.
*/
bool global_KGB::compact(RootNbr n, // assumed imaginary at |a|
			 const GlobalTitsElement& a) const
{
  const auto& i_tab = ic.involution_table();
  const auto& rd = ic.rootDatum();
  RatCoweight grading_cwt= a.torus_part().as_Qmod2Z();
  grading_cwt +=
    RatCoweight(rd.dual_twoRho(i_tab.imaginary_roots(i_tab.nr(a.tw()))),2);

  return grading_cwt.dot(rd.root(n))%2==0; //  whether compact
}

KGBElt global_KGB::lookup(const GlobalTitsElement& a) const
{
  const Cartan_orbits& i_tab = ic.involution_table();

  KGBEltPair p = tauPacket(a.tw());
  for (KGBElt x=p.first; x<p.second; ++x)
    if (i_tab.x_equiv(element(x),a))
      return x;
  return UndefKGB; // report failure
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
void global_KGB::generate_involutions(size_t num_inv)
{
  inv_nrs.reserve(num_inv);
  inv_loc.assign(ic.numInvolutions(),-1);

  const TwistedWeylGroup& W = twistedWeylGroup();
  const Cartan_orbits& i_tab = ic.involution_table();

  inv_index end_length=0; // end of the length-interval under construction

  inv_nrs.push_back(i_tab.nr(TwistedInvolution())); // set initial element
  inv_loc[inv_nrs[0]]=0; // and record its index |0| in |inv_loc|

  while (end_length<inv_nrs.size())
  {
    const inv_index start_length=end_length; // start new interval here
    end_length = inv_nrs.size(); // and run until current end

    for (weyl::Generator s=0; s<W.rank(); ++s)
      for (inv_index i=start_length; i<end_length; ++i)
      {
	const TwistedInvolution& parent = i_tab.involution(inv_nrs[i]);
	const TwistedInvolution child = W.hasTwistedCommutation(s,parent)
	  ? W.prod(s,parent) : W.twistedConjugated(parent,s);

	InvolutionNbr inv = i_tab.nr(child);
	if (inv_loc[inv]==-1u) // involution not yet seen here
	{
	  inv_loc[inv] = inv_nrs.size(); // will index element pushed next:
	  inv_nrs.push_back(inv);
	}
      } // |for(i)|; |for(s)|;
  } // while length interval non-empty

  assert(inv_nrs.size()==num_inv);
} // |global_KGB::generate_involutions|

void global_KGB::generate(size_t predicted_size, bool dual_twist)
{
  const Cartan_orbits& i_tab = ic.involution_table();
  const TwistedWeylGroup& W = Tg; // for when |GlobalTitsGroup| is not used

  first_of_tau.reserve(inv_nrs.size()+1); // currently |first_of_tau.size()==2|

  KGB_elt_entry::Pooltype elt_pool; elt_pool.reserve(predicted_size);
  HashTable<KGB_elt_entry,KGBElt> elt_hash(elt_pool);

  inv_index end_length=0; // end of the length-interval under construction

  for (size_t i=0; i<elt.size(); ++i)
    elt_hash.match(i_tab.x_pack(elt[i]));

  assert(elt_hash.size()==elt.size()); // all distinct; in fact an invariant

  while (end_length<first_of_tau.size()-1)
  {
    const inv_index start_length=end_length; // start new interval here
    end_length = first_of_tau.size()-1; // and run until current end

    for (weyl::Generator s=0; s<W.rank(); ++s)
      for (inv_index index=start_length; index<end_length; ++index)
      {
	const TwistedInvolution& tw = nth_involution(index);
	TwistedInvolution new_tw =W.twistedConjugated(tw,s);
	inv_index new_nr = inv_loc[i_tab.nr(new_tw)];
	bool is_new = new_nr >= first_of_tau.size()-1;
	if (is_new) // check that involution indices come in predicted order
	  assert(new_nr==first_of_tau.size()-1); // since we mimick generation

	// record whether |s| is imaginary for current involution
	bool imaginary = new_nr==index // that is: twisted conjugation stable
	  and not W.hasDescent(s,tw); // and |s| is a ascent (otherwise: real)

	// generate cross links
	for (KGBElt x=first_of_tau[index]; x<first_of_tau[index+1]; ++x)
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
	  new_nr = inv_loc[i_tab.nr(new_tw)];
	  is_new = new_nr >= first_of_tau.size()-1;
	  if (is_new) // check that involution indices come in predicted order
	    assert(new_nr==first_of_tau.size()-1); // since we mimick generation

	  for (KGBElt x=first_of_tau[index]; x<first_of_tau[index+1]; ++x)
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


      } // |for(index)|;  // |for(s)|
  } // while length interval non-empty

  assert(elt.size()==predicted_size or predicted_size==0);

  // finally set the Hermitian dual links
  for (KGBElt i=0; i<elt.size(); ++i)
    info[i].dual = lookup
      (dual_twist ? Tg.dual_twisted(elt[i]) : Tg.twisted(elt[i]));
} // |global_KGB::generate|



/*

        The KGB class, public methods

*/

KGB::KGB(RealReductiveGroup& G,
	 const BitMap& Cartan_classes, bool dual_twist)
  : KGB_base(G.innerClass(),G.innerClass().semisimpleRank())
  , G(G)
  , Cartan()
  , left_torus_part()
  , d_state()
  , d_bruhat(NULL)
  , d_base(NULL)
{
  //const TitsGroup& Tg = ic.titsGroup();
  size_t rank = ic.semisimpleRank(); // |ic.rank()| does not interest us here

  {// check that |Cartan_classes| is upwards closed within |G.Cartan_set()|
    BitMap test_set = G.Cartan_set();
    test_set.andnot(Cartan_classes);
    for (BitMap::iterator it=test_set.begin(); it(); ++it)
      assert(Cartan_classes.disjoint(ic.Cartan_ordering().below(*it)));
  }

  // make sure |G| has information about involutions for |Cartan_classes|
  const Cartan_orbits& i_tab = ic.involution_table();

  for (auto it=Cartan_classes.begin(); it(); ++it)
    G.innerClass().generate_Cartan_orbit(*it);

  tits::TE_Entry::Pooltype elt_pool; // of size |size()|
  HashTable<tits::TE_Entry,KGBElt> elt_hash(elt_pool);

  size_t size = ic.KGB_size(G.realForm(),Cartan_classes);

  elt_pool.reserve(size);
  KGB_base::reserve(size);

  // check if we are being called to do a full or small KGB construction
  if (Cartan_classes.isMember(0)) // start from fundamental Cartan: do full KGB
  {
    const Grading gr = G.base_grading();
    d_base = new TitsCoset(ic,gr);
    TitsElt a (d_base->titsGroup(),G.x0_torus_part());
    i_tab.reduce(a); // should be unnecessary
    elt_hash.match(a); // plant the seed
    KGB_base::add_element(); // add additional info for initial element
  }
  else // partial KGB construction; needs more work to get initial element(s)
  {
    tits::EnrichedTitsGroup square_class_base(G);
    d_base = new TitsCoset(square_class_base); // copy construct from base
    RealFormNbr rf=G.realForm();
    assert(square_class_base.square()==ic.xi_square(rf));

    set::EltList mins=ic.Cartan_ordering().minima(Cartan_classes);
    // for (small) block there should be just one minimum, but we loop anyway
    for (set::EltList::iterator it=mins.begin(); it!=mins.end(); ++it)
    {
      TitsElt a= square_class_base.backtrack_seed(ic,rf,CartanNbr(*it));
      i_tab.reduce(a);

      size_t k=elt_hash.match(a);
      assert(k==info.size()); // this KGB element should be new
      ndebug_use(k);
      KGB_base::add_element(); // add additional info for element |k|
    } // |for (it)|
  }

  // now inductively fill the table |elt_pool|/|elt_hash|, and related arrays
  for (KGBElt x=0; x<elt_hash.size(); ++x) // loop makes |elt_hash| grow
  {
    // extend by cross actions
    const TitsElt& current = elt_pool[x];
    EltInfo& my_info = info[x];

    for (weyl::Generator s=0; s<rank; ++s)
    {
      KGBfields& my_s = data[s][x];
      TitsElt a = current; basedTitsGroup().basedTwistedConjugate(a,s);
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
	my_info.status.set_imaginary
	  (s,basedTitsGroup().simple_grading(current,s));
	my_info.desc.reset(s); // imaginary roots are never (KGB) descents

	if (my_info.status[s] == gradings::Status::ImaginaryNoncompact)
	{
	  // Cayley-transform |current| by $\sigma_s$
	  TitsElt a = current; basedTitsGroup().Cayley_transform(a,s);
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

    inv_nrs.reserve(ic.numInvolutions(Cartan_classes));
    for (BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
      for (InvolutionNbr i=i_tab[*it].start; i<i_tab[*it].end(); ++i)
	inv_nrs.push_back(i);
    inv_loc.assign(ic.numInvolutions(),-1);

    // sort involutions length: Weyl length; internal number
    std::stable_sort(inv_nrs.begin(),inv_nrs.end(),i_tab.less());
    for (inv_index i=0; i<inv_nrs.size(); ++i)
      inv_loc[inv_nrs[i]] = i;

    std::vector<inv_index> invs; // entries are indices into |inv_nrs|
    invs.reserve(size);
    for (KGBElt x=0; x<size; ++x) // list involution indices; cannot yet use
      invs.push_back(inv_loc[i_tab.nr(elt_pool[x].tw())]); // |involution(x)| !

    // record order |a1| of traversing the |x| by increasing |involution|, and
    // write to |first_of_tau| the boundaries of involution classes of |x|
    a1 = permutations::standardization(invs,inv_nrs.size(),&first_of_tau);
  }

  Permutation a(a1,-1); // |a[x]| locates KGB element that should become |x|

  // permute |data|, sending entry |x| to place |a1[x]}, renumbering internally
  for (weyl::Generator s=0; s<rank; ++s)
  {
    std::vector<KGBfields>& ds = data[s];
    ds=a.pull_back(ds); // permute all the links for |s| (mode assign |vector|)
    for (KGBElt x=0; x<size; ++x) // now renumber internally through |a1|
    {
      KGBfields& dsx = ds[x];
      dsx.cross_image = a1[dsx.cross_image];
      if (dsx.Cayley_image != UndefKGB)
	dsx.Cayley_image = a1[dsx.Cayley_image];
    }
  }

  info=a.pull_back(info); // similarly permute |info|, no renumbering

  left_torus_part.reserve(size);
  for (KGBElt x=0; x<size; ++x) // pull back through |a| while constructing
    left_torus_part.push_back(titsGroup().left_torus_part(elt_pool[a[x]]));

  Cartan.reserve(inv_nrs.size());
  for (auto it=inv_nrs.begin(); it!=inv_nrs.end(); ++it)
    Cartan.push_back(i_tab.Cartan_class(*it));

  TorusPart shift(ic.rank());
  bool do_dual_twist = dual_twist and is_dual_twist_stable(G,shift);
  if (do_dual_twist)
  {
    // see whether some element (our initial seed) maps into the block
    TitsElt test = titsGroup().dual_twisted(elt_pool[0],shift);
    i_tab.reduce(test);
    if (elt_hash.find(test)==elt_hash.empty)
      do_dual_twist = false; // twist stablises the square class, but not KGB
  }
  // finally install inverse Cayley and twist links
  for (KGBElt x=0; x<size; ++x)
  {
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
    if (not dual_twist)
      info[x].dual = lookup(titsGroup().twisted(titsElt(x)));
    else if (do_dual_twist)
      info[x].dual = lookup(titsGroup().dual_twisted(titsElt(x),shift));
    else
      info[x].dual = UndefKGB;
  }
} // |KGB::KGB(G,Cartan_classes,i_tab)|


bool KGB::is_dual_twist_stable
  (const RealReductiveGroup& G, TorusPart& shift) const
{
  // although |G| will be for "dualrealform", we use non-dualised nomenclature
  const RootDatum& rd = G.rootDatum();

  Grading base = basedTitsGroup().base_grading();
  RatWeight rw (titsGroup().rank());
  for (Grading::iterator it=base.begin(); it(); ++it) // flagged simple roots
  {
    rw -= rd.fundamental_coweight(*it); // because relative to implicit base
    rw += rd.fundamental_coweight(weylGroup().Chevalley_dual(*it));
  }
  // before requiring integrality, we need to mod out by equivalence
  int_Matrix A = G.innerClass().distinguished();
  A.transpose() += 1;
  int_Matrix projector = lattice::row_saturate(A);
  projector.apply_to(rw.numerator()); // this may change the size of |rw|
  rw.normalize();
  if (rw.denominator()!=1)
    return false;

  // now there is an integer vector with same projection as |rw|; find one
  BinaryEquationList eqns; eqns.reserve(projector.numRows());
  for (unsigned int i=0; i<projector.numRows(); ++i)
  {
    eqns.push_back(BinaryEquation(projector.row(i)));
    eqns.back().pushBack(rw.numerator()[i]%2!=0); // rhs of equation
  }

  bool success = bitvector::solvable(eqns,shift);
  assert(success);
  return success;
} // |KGB::is_dual_twist_stable|



KGB::~KGB() { delete d_bruhat; delete d_base; }

/******** copy, assignment and swap ******************************************/


/******** accessors **********************************************************/

TitsElt KGB::titsElt(KGBElt x) const
{ return TitsElt(titsGroup(),left_torus_part[x],involution(x)); }

size_t KGB::torus_rank() const { return titsGroup().rank(); }

RatCoweight KGB::base_grading_vector() const
{
  return G.g_rho_check();
}

// this rational coweight for |x| is central to synthetic real forms
RatCoweight KGB::torus_factor(KGBElt x) const
{
  RatCoweight tf = base_grading_vector();
  tf -= lift(torus_part(x));

  // finally ensure result is $\theta^t$-fixed
  return symmetrise(tf,involution_matrix(x));
}

// Looks up a |TitsElt| value and returns its KGB number, or |size()|
// Since KGB does not have mod space handy, must assume |a| already reduced
KGBElt KGB::lookup(TitsElt a) const
{
  ic.involution_table().reduce(a); // make sure |a| is reduced before searching
  const TitsGroup& Tg=titsGroup();
  KGBEltPair p = tauPacket(a.tw());
  TorusPart t = Tg.left_torus_part(a);
  for (KGBElt x=p.first; x<p.second; ++x)
    if (Tg.left_torus_part(titsElt(x))==t)
      return x;
  return UndefKGB; // report failure
}

// act by external twist |delta| on KGB element |x|
KGBElt KGB::twisted(KGBElt x,const WeightInvolution& delta) const
{
  auto a = titsElt(x);
  auto delta_twist = rootdata::twist(G.rootDatum(),delta);
  auto delta2 = BinaryMap(delta);
  RatCoweight diff = (base_grading_vector()-base_grading_vector()*delta)
    .normalize();
  if (diff.denominator()!=1)
    return UndefKGB;
  TorusPart corr(diff.numerator());
  TitsElt twisted_a
    (titsGroup(),
     titsGroup().left_torus_part(a)*delta2+corr,
     weylGroup().translation(a.w(),delta_twist)
     );
  return lookup(twisted_a);
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



/*****************************************************************************

            Chapter II -- The auxiliary classes, methods

******************************************************************************/


/*****************************************************************************

        Chapter III -- Functions declared in kgb.h

******************************************************************************/

// general cross action in root $\alpha$
KGBElt cross(const KGB_base& kgb, KGBElt x, RootNbr alpha)
{
  const RootSystem& rs = kgb.rootDatum();
  const WeylWord ww = conjugate_to_simple(rs,alpha);
  return kgb.cross(ww,kgb.cross(alpha,kgb.cross(x,ww)));
}

// general Cayley transform in root $\alpha$, non-compact imaginary
KGBElt any_Cayley (const KGB_base& kgb, KGBElt x, RootNbr alpha)
{
  const RootSystem& rs = kgb.rootDatum();
  const WeylWord ww = conjugate_to_simple(rs,alpha);
  x=kgb.cross(x,ww);
  KGBElt Cx = kgb.any_Cayley(alpha,x);
  if (Cx==UndefKGB)
  {
    std::string str (kgb.status(alpha,x)==gradings::Status::Complex
		     ? "complex" : "imaginary compact");
    throw std::runtime_error("Cayley transform by "+str+" root");
  }
  return kgb.cross(ww,Cx);
}



// status of |alpha| in |kgb|: conjugate to simple root follow cross actions
gradings::Status::Value status(const KGB_base& kgb, KGBElt x, RootNbr alpha)
{
  const RootSystem& rs = kgb.rootDatum();
  make_positive(rs,alpha);
  weyl::Generator s;

  while (alpha!=rs.simpleRootNbr(s=rs.find_descent(alpha)))
  {
    rs.simple_reflect_root(s,alpha);
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

} // |namespace|

} // |namespace kgb|
} // |namespace atlas|
