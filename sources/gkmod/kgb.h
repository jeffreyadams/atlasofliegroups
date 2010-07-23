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

#include "weyl.h"
#include "tits.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "bruhat.h" // class definition needed for inlined KGB destructor
#include "gradings.h"
#include "hashtable.h"

namespace atlas {

namespace kgb {


/******** type definitions **************************************************/

//! \brief per KGB element information
struct KGBEltInfo
{
  gradings::Status status; ///< status of each simple root for this element
  unsigned int length; ///< dimension of the K orbit on G/B, minus minimal one
  unsigned int cartan; ///< records to which Cartan class this element belongs
  Descent desc; ///< flags which simple reflections give a descent

  weyl::TwistedInvolution inv;

  KGBEltInfo(unsigned int l, unsigned int c, weyl::TwistedInvolution tw)
  : status(), length(l), cartan(c), desc(), inv(tw) {}

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
  const rootdata::RootDatum& rd; // needed to handle non-simple roots
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
  explicit KGB_base(const weyl::TwistedWeylGroup& Wg,
		    const rootdata::RootDatum& datum)
    : rd(datum)
    , W(Wg)
    , data(W.rank())
    , info()
    , inv_pool(), inv_hash(inv_pool)
    , first_of_tau()
  {}

 public:
  KGB_base (const KGB_base& org) // copy contructor
    : rd(org.rd)  // share
    , W(org.W) // share
    , data(org.data)
    , info(org.info)
    , inv_pool(org.inv_pool) // copy
    , inv_hash(inv_pool) // reconstruct, using our own |pool|
    , first_of_tau(org.first_of_tau)
  {}
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

  size_t length(KGBElt x) const { return info[x].length; }
  size_t Cartan_class(KGBElt x) const { return info[x].cartan; }
  weyl::TwistedInvolution involution(KGBElt x) const { return info[x].inv; }
  const Descent& descent(KGBElt x) const { return info[x].desc; }
  bool isDescent(weyl::Generator s, KGBElt x) const
    { return descent(x).test(s); }
  const gradings::Status& status(KGBElt x) const { return info[x].status; }
  gradings::Status::Value status(weyl::Generator s, KGBElt x) const
   { return status(x)[s]; }

  bool root_is_descent(rootdata::RootNbr alpha, KGBElt x) const;
  gradings::Status::Value root_status(rootdata::RootNbr alpha, KGBElt x) const;

  bool isAscent(weyl::Generator s, KGBElt x) const // not true for imag cpct!
    { return not isDescent(s,x)
      and not (status(s,x) == gradings::Status::ImaginaryCompact);
    }

  size_t weylLength(KGBElt x) const // needed in sorting, in case |length| ties
    { return weylGroup().length(involution(x).w()); }

  weyl::TwistedInvolution nth_involution(unsigned int n) const
    { return inv_pool[n]; }

  // range of KGB elements with given twisted involution, needed for block
  KGBEltPair tauPacket(const weyl::TwistedInvolution&) const;
  size_t packet_size(const weyl::TwistedInvolution&) const;

 protected:
  void reserve (size_t n); // prepare for generating |n| elements
  void add_element // create entry in |data| and |info| structures
   (unsigned int length, unsigned int Cartan_class, weyl::TwistedInvolution tw);


}; // |class KGB_base|

/*!

A |GlobalFiberData| object associates to each twisted involution a reduced
description of its -1 eigenspace, which allows a quick test for equivalence of
|GlobalTitsElement| values at this twised involution. The object also records
the halfsum of the positive imaginary coroots, and the Cartan class.
*/
class GlobalFiberData
{
  const hashtable::HashTable<weyl::TI_Entry,unsigned int>& hash_table;
  std::vector<unsigned int> Cartans;
  std::vector<latticetypes::LatticeMatrix> proj; // projectors for equivalence
  std::vector<latticetypes::LatticeElt> check_2rho_imag;
public:
  GlobalFiberData(complexredgp::ComplexReductiveGroup& G,
		  const hashtable::HashTable<weyl::TI_Entry,unsigned int>& h);

  GlobalFiberData(const GlobalFiberData& org) // copy contructor, handle ref
    : hash_table(org.hash_table) // share
    , Cartans(org.Cartans)
    , proj(org.proj)
    , check_2rho_imag(org.check_2rho_imag)
  {}

  size_t Cartan_class(const weyl::TwistedInvolution& tw) const
  {
    return Cartans[hash_table.find(tw)];
  }

  bool equivalent(const tits::GlobalTitsElement& x,
		  const tits::GlobalTitsElement& y) const;

  latticetypes::RatLatticeElt // a value characterizing the equivalence class
    fingerprint(const tits::GlobalTitsElement& x) const;

  int at_rho_imaginary(const latticetypes::Weight& alpha, // imaginary root
		       const weyl::TwistedInvolution& tw) const
    { return alpha.dot(check_2rho_imag[hash_table.find(tw)])/2; }
}; // |class GlobalFiberData|


class global_KGB : public KGB_base
{
  complexredgp::ComplexReductiveGroup& G;
  const tits::GlobalTitsGroup& Tg;
  GlobalFiberData fiber_data;
  std::vector<tits::GlobalTitsElement> elt;

 public:
  global_KGB(complexredgp::ComplexReductiveGroup& G,
	     const tits::GlobalTitsGroup& Tg);

  global_KGB(complexredgp::ComplexReductiveGroup& G,
	     const tits::GlobalTitsElement& x); // generate KGB containing |x|

  global_KGB(const global_KGB& org) // copy contructor, handle references
    : KGB_base(org)
    , G(org.G) // share
    , Tg(org.Tg) // share
    , fiber_data(org.fiber_data)
    , elt(org.elt)
  {}

// accessors
  const tits::GlobalTitsGroup& globalTitsGroup() const { return Tg; }
  const rootdata::RootDatum& rootDatum() const { return G.rootDatum(); }
  const complexredgp::ComplexReductiveGroup& complexGroup() const
    { return G; }

  tits::TorusElement torus_part(KGBElt x) const { return elt[x].torus_part(); }
  const tits::GlobalTitsElement& element(KGBElt x) const { return elt[x]; }

  bool compact(rootdata::RootNbr alpha, const tits::GlobalTitsElement& a) const;
  kgb::KGBElt lookup(const tits::GlobalTitsElement& x) const;

 private:
  void generate_involutions(size_t n);
  void generate(size_t predicted_size);

}; // |class global_KGB|



/*!
\brief Represents the orbits of K on G/B for a particular real form.

This class adds some information with respect to that kept in |KGB_base|, and
most importantly carries out the actual filling of the bas object. As
additional data that are held in this derived class there is the
|BasedTitsGroup| used during construction, and the torus parts (relative to
the base point) that distinguish elements in the same fiber. This class also
provides the possibility to generate and store the Bruhat order on the set.
*/

class KGB : public KGB_base
{

  enum State { BruhatConstructed, NumStates };

  std::vector<tits::TorusPart> left_torus_part; // of size |size()|
  bitset::BitSet<NumStates> d_state;

/*! \brief Owned pointer to the Bruhat order on KGB (or NULL).

The class BruhatOrder contains a Poset describing the full partial order,
and in addition the Hasse diagram (set of all covering relations).
*/
  bruhat::BruhatOrder* d_bruhat;

  //! \brief Owned pointer to the based Tits group.
  tits::BasedTitsGroup* d_base; // pointer, because constructed in helper class

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

//! \brief The based Tits group.
  const tits::BasedTitsGroup& basedTitsGroup() const { return *d_base; }
//! \brief The Tits group.
  const tits::TitsGroup& titsGroup() const { return d_base->titsGroup(); }
//! \brief The Weyl group.

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
  so this method should hold whenever |x<y| and |involution(x)!=involution(y)|
  and it turns out the method is never used, so I have excluded the code. MvL.
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

// private methods
private:
  size_t generate
    (realredgp::RealReductiveGroup& GR,const bitmap::BitMap& Cartan_classes);

  void fillBruhat();

}; // |class KGB|

// extract the KGB for a root subsystem
struct subsys_KGB : public KGB_base
{
  subsys_KGB(const KGB_base& kgb,
	     const rootdata::RootDatum& rd,
	     const rootdata::RootList& subsys,
	     const rootdata::RootDatum& sub_datum, // must be constructed before
	     const weyl::TwistedWeylGroup& sub_W,  // call for lifetime reasons
	     kgb::KGBElt x);
}; // |struct subsys_KGB|

/* ****************** function definitions **************************** */

// general cross action in (non simple) root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt cross(const KGB_base& kgb, KGBElt x,
	     weyl::Generator s, weyl::WeylWord ww);

// general Cayley transform in (non simple) non-compact imaginary root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt Cayley (const KGB_base& kgb, KGBElt x,
	       weyl::Generator s, weyl::WeylWord ww);

// general inverse Cayley transform (choice) in (non simple) real root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt inverse_Cayley (const KGB_base& kgb, KGBElt x,
		       weyl::Generator s, weyl::WeylWord ww);


} // |namespace kgb|

} // |namespace atlas|


#endif
