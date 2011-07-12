/*!
\file
\brief Class definition and function declarations for the class KGB
representing orbits of K on G/B.
*/
/*
  This is kgb.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KGB_H  /* guard against multiple inclusions */
#define KGB_H

#include "kgb_fwd.h"

#include "bitmap_fwd.h"
#include "realredgp_fwd.h"

#include "subdatum.h"

#include "weyl.h"
#include "tits.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "bruhat.h" // class definition needed for inlined KGB destructor
#include "gradings.h"
#include "hashtable.h"

#include <algorithm>
#include <iostream> // for virtual print method

namespace atlas {

namespace kgb {


/******** type definitions **************************************************/

//! \brief per KGB element information
struct KGBEltInfo
{
  gradings::Status status; ///< status of each simple root for this element
  unsigned int length; ///< dimension of the K orbit on G/B, minus minimal one
  Descent desc; ///< flags which simple reflections give a descent

  weyl::TwistedInvolution inv;

  KGBEltInfo(unsigned int l, const weyl::TwistedInvolution& tw)
  : status(), length(l), desc(), inv(tw) {}

}; // |struct KGBEltInfo|

/*
   The following base class follows a somewhat particular design: it is not an
   abstract base class in the sense that it has no virtual methods (and
   therefore a fortiori no purely virtual ones), yet it is not intended to
   have any independent instances. The point is that the base class provides
   all methods needed for basic \emph{use}, in for instance the block
   construction, but it does not provide sufficient methods for constructing
   the KGB set. Indeed the basic constructor is made |protected| to emphasise
   that it is up to derived classes to actually fill the tables in the
   structure, using concrete representations for the KGB elements (possibly
   using |tits::TitsElt|) that are mainly relevant during construction. Once
   constructed there is no need to slice off the base object (although it
   could be done) and the full derived object can be used for instance to
   print more information about block elements than available in |KGB_base|.
 */
class KGB_base
{
 protected: // available during construction from derived classes
  const weyl::TwistedWeylGroup& W; // hold a reference for convenience

  struct KGBfields // data parametrized by simple reflection and KGB element
  {
    KGBElt cross_image; // cross action image
    KGBElt Cayley_image;
    KGBEltPair inverse_Cayley_image;

    KGBfields()
    : cross_image(UndefKGB)
    , Cayley_image(UndefKGB)
      , inverse_Cayley_image(std::make_pair(UndefKGB,UndefKGB)) {}
  }; // | KGBfields|

  std::vector<std::vector<KGBfields> > data; // first index: simple reflection
  std::vector<KGBEltInfo> info; // per element information

  //!\brief tables to map twisted involutions to their sequence number
  weyl::TI_Entry::Pooltype inv_pool;
  hashtable::HashTable<weyl::TI_Entry,unsigned int> inv_hash;

  //!\brief to help find range of elements with fixed twisted involution
  std::vector<KGBElt> first_of_tau; // size: |numInvolutions()+1|



 protected: // constructor is only meant for use from derived classes
  explicit KGB_base(const weyl::TwistedWeylGroup& Wg)
    : W(Wg)
    , data(W.rank())
    , info()
    , inv_pool(), inv_hash(inv_pool)
    , first_of_tau()
  {}

 public:
  KGB_base (const KGB_base& org) // copy contructor
    : W(org.W) // share
    , data(org.data)
    , info(org.info)
    , inv_pool(org.inv_pool) // copy
    , inv_hash(inv_pool) // reconstruct, using our own |pool|
    , first_of_tau(org.first_of_tau)
  {}

  virtual ~KGB_base(){} // maybe overly cautious; no |KGB_base*| is intended
// accessors

  size_t rank() const { return data.size(); } // number of simple reflections
  size_t size() const { return info.size(); } // number of KGB elements
  unsigned int nr_involutions() const { return inv_pool.size(); }

  const weyl::TwistedWeylGroup& twistedWeylGroup() const { return W; }
  const weyl::WeylGroup& weylGroup() const { return W.weylGroup(); }

  KGBElt cross(weyl::Generator s, KGBElt x) const
    { return data[s][x].cross_image; }
  KGBElt cayley(weyl::Generator s, KGBElt x) const
    { return data[s][x].Cayley_image; }
  KGBEltPair inverseCayley(weyl::Generator s, KGBElt x) const
    { return data[s][x].inverse_Cayley_image; }

  KGBElt cross(const weyl::WeylWord& ww, KGBElt x) const;
  KGBElt cross(KGBElt x, const weyl::WeylWord& ww) const;

  size_t length(KGBElt x) const { return info[x].length; }
  weyl::TwistedInvolution involution(KGBElt x) const { return info[x].inv; }
  const Descent& descent(KGBElt x) const { return info[x].desc; }
  bool isDescent(weyl::Generator s, KGBElt x) const
    { return descent(x).test(s); }
  const gradings::Status& status(KGBElt x) const { return info[x].status; }
  gradings::Status::Value status(weyl::Generator s, KGBElt x) const
   { return status(x)[s]; }

  bool isAscent(weyl::Generator s, KGBElt x) const // not true for imag cpct!
    { return not isDescent(s,x)
      and not (status(s,x) == gradings::Status::ImaginaryCompact);
    }

  size_t weylLength(KGBElt x) const // needed in sorting, in case |length| ties
    { return weylGroup().length(involution(x).w()); }

  weyl::TwistedInvolution nth_involution(unsigned int n) const
    { return inv_pool[n]; }

  size_t involution_index(KGBElt x) const
  { return std::upper_bound(first_of_tau.begin(),first_of_tau.end(),x)
      -first_of_tau.begin() -1;
  }
  KGBEltPair packet (KGBElt x) const
  { size_t i = involution_index(x);
    return KGBEltPair(first_of_tau[i],first_of_tau[i+1]);
  }

  // range of KGB elements with given twisted involution, needed for block
  KGBEltPair tauPacket(const weyl::TwistedInvolution&) const;
  size_t packet_size(const weyl::TwistedInvolution&) const;

// virtual methods
  virtual size_t Cartan_class(KGBElt x) const { return ~0ul; }
  // print derived-class specific per-element information
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const
  { return strm; }

 protected:
  void reserve (size_t n); // prepare for generating |n| elements
  void add_element // create entry in |data| and |info| structures
   (unsigned int length, weyl::TwistedInvolution tw);


}; // |class KGB_base|




//		      |global_KGB| and subsidiary types



// a |GlobalTitsElement| can be hashed as (fingerprint,twisted_inv) pair
struct KGB_elt_entry
{
  tits::TorusElement t_rep; // a representative torus element, ignored in test
  weyl::TI_Entry tw;
  latticetypes::RatLatticeElt fingerprint; // charcterizes the torus element

  // obligatory fields for hashable entry
  typedef std::vector<KGB_elt_entry> Pooltype;
  size_t hashCode(size_t modulus) const; // hash function
  bool operator !=(const KGB_elt_entry& x) const
  { return tw!=x.tw or fingerprint!=x.fingerprint; } // ignore repr

  KGB_elt_entry (const latticetypes::RatLatticeElt& f,
		 const tits::GlobalTitsElement& y)
  : t_rep(y.torus_part()), tw(y.tw()), fingerprint(f) {}

  tits::GlobalTitsElement repr() const
  { return tits::GlobalTitsElement(t_rep,tw); }
}; //  |struct KGB_elt_entry|

/*!

A |GlobalFiberData| object associates to each twisted involution a reduced
description of its -1 eigenspace, which allows a quick test for equivalence of
|GlobalTitsElement| values at this twised involution. The object also records
the halfsum of the positive imaginary coroots, and the Cartan class.

Originally conceived for global use within an inner class (the constructor
uses its list of Cartan classes to generate all involutions), functionality
now also includes dynamically adding Cartan classes of involutions, by
specifying a new twisted involution with its involution matrix, from which a
whole conjugacy class will be added to the tables. This can even be done using
only the Weyl group of a (co)root subsystem, giving smaller Cartan classes.

The terminology remains attached to the original use, even though the use with
dynamical addition of Cartan classes is associated to a |GlobalTitsGroup|,
which will be employed from the dual side.
 */
class GlobalFiberData
{
  hashtable::HashTable<weyl::TI_Entry,unsigned int>& hash_table;

  struct inv_info
  {
    unsigned int Cartan;
    latticetypes::LatticeMatrix proj; // projectors for equivalence
    latticetypes::LatticeElt check_2rho_imag;
    rootdata::RootList simple_imag,simple_real;

  inv_info() : Cartan(~0), proj(), check_2rho_imag(), simple_imag()
    {} // allow uninitialized
  inv_info(unsigned int c,
	     const latticetypes::LatticeMatrix& p,
	     const latticetypes::LatticeElt& c2i,
	     const rootdata::RootList& si,
	     const rootdata::RootList& sr)
  : Cartan(c),proj(p),check_2rho_imag(c2i),simple_imag(si),simple_real(sr) {}
  };

  std::vector<inv_info> info;
  std::vector<latticetypes::LatticeMatrix> refl; // reflections at dual side

public:
  GlobalFiberData(complexredgp::ComplexReductiveGroup& G,
		  hashtable::HashTable<weyl::TI_Entry,unsigned int>& h);

  // contructor that does not install eny Cartan classes yet
  GlobalFiberData(const tits::GlobalTitsGroup& Tg,
		  hashtable::HashTable<weyl::TI_Entry,unsigned int>& h);

  GlobalFiberData(const GlobalFiberData& org) // copy contructor, handle ref
    : hash_table(org.hash_table) // share
    , info(org.info)
  {}

  //accessors
  size_t Cartan_class(const weyl::TwistedInvolution& tw) const
  {
    return info[hash_table.find(tw)].Cartan;
  }

  const rootdata::RootList& imaginary_basis(const weyl::TwistedInvolution& tw)
    const
  { return info[hash_table.find(tw)].simple_imag; }
  const rootdata::RootList& real_basis(const weyl::TwistedInvolution& tw)
    const
  { return info[hash_table.find(tw)].simple_real; }

  bool equivalent(const tits::GlobalTitsElement& x,
		  const tits::GlobalTitsElement& y) const;

  latticetypes::RatLatticeElt // a value characterizing the equivalence class
    fingerprint(const tits::GlobalTitsElement& x) const;

  int at_rho_imaginary(const latticetypes::Weight& alpha, // imaginary root
		       const weyl::TwistedInvolution& tw) const
    { return alpha.dot(info[hash_table.find(tw)].check_2rho_imag)/2; }

  tits::GlobalTitsElement
    imaginary_cross(const rootdata::RootDatum& dual_rd, // pragmatic reason
		    rootdata::RootNbr alpha, // any imaginary root
		    tits::GlobalTitsElement a) const; // |a| is by-value


  KGB_elt_entry pack(const tits::GlobalTitsElement& y) const
  { return KGB_elt_entry(fingerprint(y),y); }

  // manipulators
  void add_class(const subdatum::SubSystem& sub, // determines root system
		 const tits::GlobalTitsGroup& Tg, // interprets |tw| and values
		 const weyl::TwistedInvolution& tw); // derived from it

}; // |class GlobalFiberData|




class global_KGB : public KGB_base
{
  const tits::GlobalTitsGroup Tg;
  GlobalFiberData fiber_data;
  std::vector<tits::GlobalTitsElement> elt;

  global_KGB(const global_KGB& org); // forbid copying

 public:
  global_KGB(complexredgp::ComplexReductiveGroup& G);

  global_KGB(complexredgp::ComplexReductiveGroup& G,
	     const tits::GlobalTitsElement& x); // generate KGB containing |x|

// accessors
  const tits::GlobalTitsGroup& globalTitsGroup() const { return Tg; }

  tits::TorusElement torus_part(KGBElt x) const { return elt[x].torus_part(); }
  const tits::GlobalTitsElement& element(KGBElt x) const { return elt[x]; }

  bool compact(const rootdata::RootDatum& rd, rootdata::RootNbr alpha,
	       const tits::GlobalTitsElement& a) const;
  kgb::KGBElt lookup(const tits::GlobalTitsElement& x) const;

// virtual methods
  virtual size_t Cartan_class(KGBElt x) const
  { return fiber_data.Cartan_class(involution(x)); }
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const;

 private:
  void generate_involutions(size_t n);
  void generate(const rootdata::RootDatum& rd, size_t predicted_size);

}; // |class global_KGB|



//			     Fokko's |class KGB|


/*!
\brief Represents the orbits of K on G/B for a particular real form.

This class adds some information with respect to that kept in |KGB_base|, and
most importantly carries out the actual filling of the base object. As
additional data that are held in this derived class there is the
|TitsCoset| used during construction, and the torus parts (relative to
the base point) that distinguish elements in the same fiber. This class also
provides the possibility to generate and store the Bruhat order on the set.
*/

class KGB : public KGB_base
{

  enum State { BruhatConstructed, NumStates };

  const rootdata::RootDatum& rd; // needed for |half_rho| only
  std::vector<unsigned int> Cartan; ///< records Cartan classes of elements

  std::vector<tits::TorusPart> left_torus_part; // of size |size()|
  bitset::BitSet<NumStates> d_state;

/*! \brief Owned pointer to the Bruhat order on KGB (or NULL).

The class BruhatOrder contains a Poset describing the full partial order,
and in addition the Hasse diagram (set of all covering relations).
*/
  bruhat::BruhatOrder* d_bruhat;

  //! \brief Owned pointer to the based Tits group.
  tits::TitsCoset* d_base; // pointer, because contructed late by |generate|

 public:

// constructors and destructors
  explicit KGB(realredgp::RealReductiveGroup& GR,
	       const bitmap::BitMap& Cartan_classes =  bitmap::BitMap(0));

  ~KGB()
    { delete d_bruhat; delete d_base; } // these are is owned (if non NULL)

// copy, assignment and swap
// these are currently reserved; if defined, they shoud take care of |d_bruhat|
 private:
  KGB(const KGB&);
  KGB& operator=(const KGB&);
 public:

// accessors

  const rootdata::RootDatum& rootDatum() const { return rd; }
//! \brief The based Tits group.
  const tits::TitsCoset& basedTitsGroup() const { return *d_base; }
//! \brief The Tits group.
  const tits::TitsGroup& titsGroup() const { return d_base->titsGroup(); }
//! \brief The Weyl group.

  latticetypes::RatWeight half_rho() const
  { return latticetypes::RatWeight(rd.twoRho(),4); }

  virtual size_t Cartan_class(KGBElt x) const;

  tits::TorusPart torus_part(KGBElt x) const { return left_torus_part[x]; }
  tits::TitsElt titsElt(KGBElt x) const
    { return tits::TitsElt(titsGroup(),left_torus_part[x],involution(x)); }
//! \brief The (non-semisimple) rank of torus parts.
  size_t torus_rank() const { return titsGroup().rank(); }

  gradings::Grading base_grading() const { return d_base->base_grading(); }

  kgb::KGBElt lookup(const tits::TitsElt& a, const tits::TitsGroup& Tg) const;

#if 0
/*!
  \brief Method that used to return whether involution(x) < involution(y).

  Explanation: the ordering is involution-length first, then weyl-length, then
  order by representative Weyl elements (weyl::TwistedInvolution::operator<).
  This is only a partial ordering, that does not distinguish elements of a
  fiber over one same twisted involution.

  A similar function is used to sort the elements of |KGB| upon construction,
  so this method should hold for any |x<y| unless their involutions are the
  same. As the method is actually never used, I have excluded the code.  MvL.
*/
  bool compare(KGBElt x, KGBElt y) const
  {
    if      (length(x) != length(y))
      return length(x) < length(y);
    else if (weylLength(x) != weylLength(y))
      return weylLength(x) < weylLength(y);
    else
      return involution(x) < involution(y);
  }
#endif
// manipulators

// Creates Hasse diagram for Bruhat order on KGB and returns reference to it
  bruhat::BruhatOrder& bruhatOrder()
  { fillBruhat(); return *d_bruhat; }

  const poset::Poset& bruhatPoset() // this creates full poset on demand
  { return bruhatOrder().poset(); }

// virtual method
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const;

// private methods
private:
  size_t generate
    (realredgp::RealReductiveGroup& GR,const bitmap::BitMap& Cartan_classes);

  void fillBruhat();

}; // |class KGB|

// extract the KGB for a root subsystem; DV says this is not always possible
class subsys_KGB : public KGB_base
{
  gradings::Grading base_grading;

 public:
  std::vector<KGBElt> in_parent; // map elements back to parent, if possible
  std::vector<unsigned int> Cartan; ///< records Cartan classes of elements
  std::vector<tits::TorusPart> torus_part; // of size |size()|
  KGBElt parent_size;

  subsys_KGB(const kgb::KGB& kgb,
	     const subdatum::SubDatum& sub, // must be constructed before
	     kgb::KGBElt x);

// virtual methods
  virtual size_t Cartan_class(KGBElt x) const;
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const;

}; // |struct subsys_KGB|

/* ****************** function definitions **************************** */

// general cross action in (non simple) root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt cross(const KGB_base& kgb, KGBElt x,
	     weyl::Generator s, const weyl::WeylWord& ww);

// general Cayley transform in (non simple) non-compact imaginary root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt Cayley (const KGB_base& kgb, KGBElt x,
	       weyl::Generator s, const weyl::WeylWord& ww);

// general inverse Cayley transform (choice) in (non simple) real root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt inverse_Cayley (const KGB_base& kgb, KGBElt x,
		       weyl::Generator s, const weyl::WeylWord& ww);


} // |namespace kgb|

} // |namespace atlas|


#endif
