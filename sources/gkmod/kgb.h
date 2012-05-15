/*!
\file
\brief Class definition and function declarations for the class KGB
representing orbits of K on G/B.
*/
/*
  This is kgb.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2011 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KGB_H  /* guard against multiple inclusions */
#define KGB_H

#include "subdatum.h"

#include "gradings.h"	// containment in |KGBEltInfo|
#include "hashtable.h"	// containment in |KGB_base|

#include <algorithm>
#include <iostream> // for virtual print method

namespace atlas {

namespace kgb {


/******** type definitions **************************************************/


/*
   The following base class follows a somewhat particular design: it is not an
   abstract base class in the sense that it has no purely virtual methods (in
   the original design there were no virtual methods, but a few have been
   added since for printing purposes), yet it is not intended to have any
   independent instances. The point is that the base class provides all
   methods needed for basic \emph{use}, in for instance the block
   construction, but it does not provide sufficient methods for constructing
   the KGB set. Indeed the basic constructor is made |protected| to emphasise
   that it is up to derived classes to actually fill the tables in the
   structure, using concrete representations for the KGB elements (possibly
   using |TitsElt|) that are mainly relevant during construction. Once
   constructed there is no need to slice off the base object (although it
   could be done) and the full derived object can be used for instance to
   print more information about block elements than available in |KGB_base|.
 */
class KGB_base
{
 protected: // available during construction from derived classes
  const ComplexReductiveGroup& G; // hold a reference for convenience

  // per KGB element information
  struct EltInfo
  {
    gradings::Status status; ///< status of each simple root for this element
    DescentSet desc; ///<  which simple reflections are complex descent or real

  EltInfo() : status(), desc() {}

  }; // |struct EltInfo|

  struct KGBfields // per KGB element and simple reflection data
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
  std::vector<EltInfo> info; // per element information

  //!\brief tables to map twisted involutions to their sequence number
  weyl::TI_Entry::Pooltype inv_pool;
  hashtable::HashTable<weyl::TI_Entry,unsigned int> inv_hash;

  //!\brief to help find range of elements with fixed twisted involution
  std::vector<KGBElt> first_of_tau; // size: |numInvolutions()+1|



 protected: // constructor is only meant for use from derived classes
  explicit KGB_base(const ComplexReductiveGroup& GC, unsigned int ss_rank)
  : G(GC)
  , data(ss_rank)
  , info()
  , inv_pool(), inv_hash(inv_pool)
  , first_of_tau()
  {}


 public:
  KGB_base (const KGB_base& org) // copy contructor
    : G(org.G) // share
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

  const ComplexReductiveGroup& complexGroup() const { return G; }
  const RootDatum& rootDatum() const;
  const WeylGroup& weylGroup() const;
  const TwistedWeylGroup& twistedWeylGroup() const;

  KGBElt cross(weyl::Generator s, KGBElt x) const
    { return data[s][x].cross_image; }
  KGBElt cayley(weyl::Generator s, KGBElt x) const
    { return data[s][x].Cayley_image; }
  KGBEltPair inverseCayley(weyl::Generator s, KGBElt x) const
    { return data[s][x].inverse_Cayley_image; }

  KGBElt cross(const WeylWord& ww, KGBElt x) const;
  KGBElt cross(KGBElt x, const WeylWord& ww) const;

  unsigned int length(KGBElt x) const;

  TwistedInvolution nth_involution(unsigned int n) const
    { return inv_pool[n]; }

  InvolutionNbr involution_index(KGBElt x) const // internal index of involution
  { return std::upper_bound(first_of_tau.begin(),first_of_tau.end(),x)
      -first_of_tau.begin() -1;
  }
  const TwistedInvolution& involution(KGBElt x) const // after construction only
  { return inv_pool[involution_index(x)]; }

  const WeightInvolution & involution_matrix(KGBElt x) const;
  InvolutionNbr inv_nr(KGBElt x) const; // external number (within inner class)

  const DescentSet& descent(KGBElt x) const { return info[x].desc; }
  bool isDescent(weyl::Generator s, KGBElt x) const
    { return descent(x).test(s); }
  const gradings::Status& status(KGBElt x) const { return info[x].status; }
  gradings::Status::Value status(weyl::Generator s, KGBElt x) const
   { return status(x)[s]; }

  bool isComplexDescent(weyl::Generator s, KGBElt x) const
  { return status(x).isComplex(s) and isDescent(s,x); }

  bool isDoubleCayleyImage(weyl::Generator s, KGBElt x) const
  { return inverseCayley(s,x).second!=UndefKGB; }

  bool isAscent(weyl::Generator s, KGBElt x) const // not true for imag cpct!
    { return not isDescent(s,x)
      and not (status(s,x) == gradings::Status::ImaginaryCompact);
    }

  size_t weylLength(KGBElt x) const // needed in sorting, in case |length| ties
    { return weylGroup().length(involution(x).w()); }

  KGBEltPair packet (KGBElt x) const
  { InvolutionNbr i = involution_index(x);
    return KGBEltPair(first_of_tau[i],first_of_tau[i+1]);
  }

  // range of KGB elements with given twisted involution, needed for block
  KGBEltPair tauPacket(const TwistedInvolution&) const;
  size_t packet_size(const TwistedInvolution&) const;

// virtual methods
  virtual CartanNbr Cartan_class(KGBElt x) const; // default, uses |G| tables
  // print derived-class specific per-element information
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const
  { return strm; }

 protected:
  void reserve (size_t n); // prepare for generating |n| elements
  void add_element();  // create entry in |data| and |info|


}; // |class KGB_base|




//		      |global_KGB| and subsidiary types



// a |GlobalTitsElement| can be hashed as (fingerprint,twisted_inv) pair
struct KGB_elt_entry
{
  TorusElement t_rep; // a representative, ignored in test
  weyl::TI_Entry tw;
  RatWeight fingerprint; // charcterizes the torus element

  // obligatory fields for hashable entry
  typedef std::vector<KGB_elt_entry> Pooltype;
  size_t hashCode(size_t modulus) const; // hash function
  bool operator !=(const KGB_elt_entry& x) const; // ignores repr

  KGB_elt_entry (const RatWeight& f,
		 const GlobalTitsElement& y);
  GlobalTitsElement repr() const;

}; //  |struct KGB_elt_entry|

/*

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
 protected: // data will also be maintained by derived class |InvInfo|
  hashtable::HashTable<weyl::TI_Entry,unsigned int>& hash_table;

  struct inv_info
  {
    unsigned int Cartan;
    int_Matrix proj; // projectors for equivalence
    Weight check_2rho_imag;
    RootNbrList simple_imag,simple_real;

  inv_info()
  : Cartan(~0), proj(), check_2rho_imag(), simple_imag(), simple_real()
    {} // allow uninitialized construction
  inv_info(unsigned int c,
	   const int_Matrix& p,
	   const Weight& c2ri,
	   const RootNbrList& si,
	   const RootNbrList& sr)
  : Cartan(c),proj(p),check_2rho_imag(c2ri),simple_imag(si),simple_real(sr) {}
  };

  std::vector<inv_info> info;
  std::vector<WeightInvolution> refl; // simple reflections at dual side

public:
  GlobalFiberData(ComplexReductiveGroup& G,
		  hashtable::HashTable<weyl::TI_Entry,unsigned int>& h);

  // contructor that does not install eny Cartan classes yet
  GlobalFiberData(const GlobalTitsGroup& Tg,
		  hashtable::HashTable<weyl::TI_Entry,unsigned int>& h);

 protected: // this one is for use by |InvInfo||
  GlobalFiberData(const SubSystem& sub,
		  hashtable::HashTable<weyl::TI_Entry,unsigned int>& h);
 public:

  GlobalFiberData(const GlobalFiberData& org) // copy contructor, handle ref
    : hash_table(org.hash_table) // share
    , info(org.info)             // copy
  {}

  //accessors
  unsigned int find(const TwistedInvolution& tw) const
  { return hash_table.find(tw); }

  CartanNbr Cartan_class(const TwistedInvolution& tw) const
  { return info[find(tw)].Cartan;}

  const RootNbrList& imaginary_basis(const TwistedInvolution& tw) const
  { return info[find(tw)].simple_imag; }
  const RootNbrList& real_basis(const TwistedInvolution& tw) const
  { return info[find(tw)].simple_real; }

  bool equivalent(const GlobalTitsElement& x,
		  const GlobalTitsElement& y) const;

  RatWeight // a value characterizing the equivalence class
    fingerprint(const GlobalTitsElement& x) const;

  int at_rho_imaginary(const rootdata::Root& alpha, // imaginary root
		       const TwistedInvolution& tw) const
    { return alpha.dot(info[hash_table.find(tw)].check_2rho_imag)/2; }

  GlobalTitsElement
    imaginary_cross(const RootDatum& dual_rd, // pragmatic reason
		    RootNbr alpha, // any imaginary root
		    GlobalTitsElement a) const; // |a| is by-value


  KGB_elt_entry pack(const GlobalTitsElement& y) const
  { return KGB_elt_entry(fingerprint(y),y); }

  // manipulators: none here, but |InvInfo| provides two of them

}; // |class GlobalFiberData|

struct InvInfo : public GlobalFiberData
{
  const SubSystem& sub;
  CartanNbr n_Cartans; // number of Cartan classes generated

  InvInfo(const SubSystem& subsys,
	  hashtable::HashTable<weyl::TI_Entry,unsigned int>& h);

//manipulators
  // add involution |tw| with |Tg.involution_matrix(tw)|; report whether new
  bool add_involution(const TwistedInvolution& tw, const GlobalTitsGroup& Tg);

  // add |tw|, a neighbor of |info[old_inv]| by cross action by |s| of |sub|
  bool add_cross_neighbor(const TwistedInvolution& tw,
			  unsigned int old_inv, weyl::Generator s);
}; // |InvInfo|



class global_KGB : public KGB_base
{
  const GlobalTitsGroup Tg;
  std::vector<GlobalTitsElement> elt;

  global_KGB(const global_KGB& org); // forbid copying

 public:
  global_KGB(ComplexReductiveGroup& G);

  global_KGB(ComplexReductiveGroup& G,
	     const GlobalTitsElement& x); // generate KGB containing |x|

// accessors
  const GlobalTitsGroup& globalTitsGroup() const { return Tg; }

  TorusElement torus_part(KGBElt x) const
  { return elt[x].torus_part(); }
  const GlobalTitsElement& element(KGBElt x) const { return elt[x]; }

  bool compact(RootNbr alpha, const GlobalTitsElement& a) const;
  KGBElt lookup(const GlobalTitsElement& x) const;

// virtual methods
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const;

 private:
  void generate_involutions(size_t n);
  void generate(size_t predicted_size);

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

  std::vector<unsigned int> Cartan; ///< records Cartan classes of elements

  std::vector<tits::TorusPart> left_torus_part; // of size |size()|
  BitSet<NumStates> d_state;

/*! \brief Owned pointer to the Bruhat order on KGB (or NULL).

The class BruhatOrder contains a Poset describing the full partial order,
and in addition the Hasse diagram (set of all covering relations).
*/
  BruhatOrder* d_bruhat;

  //! \brief Owned pointer to the based Tits group.
  TitsCoset* d_base; // pointer, because contructed late by |generate|

 public:

// constructors and destructors
  explicit KGB(RealReductiveGroup& GR,
	       const BitMap& Cartan_classes);

  ~KGB(); // { delete d_bruhat; delete d_base; } // these are owned (or NULL)

// copy, assignment and swap
// these are currently reserved; if defined, they shoud take care of |d_bruhat|
 private:
  KGB(const KGB&);
  KGB& operator=(const KGB&);
 public:

// accessors

//! \brief The based Tits group.
  const TitsCoset& basedTitsGroup() const { return *d_base; }
//! \brief The Tits group.
  const TitsGroup& titsGroup() const { return d_base->titsGroup(); }

  RatWeight half_rho() const { return RatWeight(rootDatum().twoRho(),4); }

  tits::TorusPart torus_part(KGBElt x) const { return left_torus_part[x]; }
  TorusElement torus_part_global // |torus_part| but coded as in |global_KGB|
    (const RootDatum&rd, KGBElt x) const; // needs root datum (for base grading)

  TitsElt titsElt(KGBElt x) const;
  size_t torus_rank() const; // the (non-semisimple) rank of torus parts.

  Grading base_grading() const { return d_base->base_grading(); }

  // as the name suggests, the following assumes |alpha| simple-imaginary
  bool simple_imaginary_grading(KGBElt x,RootNbr alpha) const
  { return d_base->simple_imaginary_grading(torus_part(x),alpha); }

  KGBElt lookup(const TitsElt& a, const TitsGroup& Tg) const;


// manipulators

// Creates Hasse diagram for Bruhat order on KGB and returns reference to it
  BruhatOrder& bruhatOrder() { fillBruhat(); return *d_bruhat; }

  const poset::Poset& bruhatPoset(); // this creates full poset on demand

// virtual method
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const;

// private methods
private:
  size_t generate (RealReductiveGroup& GR,const BitMap& Cartan_classes);

  void fillBruhat();

}; // |class KGB|

// extract the KGB for a root subsystem; DV says this is not always possible
class subsys_KGB : public KGB_base
{
  Grading base_grading;

 public:
  std::vector<KGBElt> in_parent; // map elements back to parent, if possible
  std::vector<unsigned int> Cartan; ///< records Cartan classes of elements
  std::vector<tits::TorusPart> torus_part; // of size |size()|
  KGBElt parent_size;

  subsys_KGB(const KGB& kgb,
	     const subdatum::SubDatum& sub, // must be constructed before
	     KGBElt x);

// virtual methods
  virtual CartanNbr Cartan_class(KGBElt x) const // override with subsys Cartan
  { return Cartan[x];}
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const;

}; // |struct subsys_KGB|

/* ****************** function definitions **************************** */

// general cross action in (non simple) root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt cross(const KGB_base& kgb, KGBElt x,
	     weyl::Generator s, const WeylWord& ww);

// general Cayley transform in (non simple) non-compact imaginary root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt Cayley (const KGB_base& kgb, KGBElt x,
	       weyl::Generator s, const WeylWord& ww);

// general inverse Cayley transform (choice) in (non simple) real root
// root is given as simple root + conjugating Weyl word to simple root
KGBElt inverse_Cayley (const KGB_base& kgb, KGBElt x,
		       weyl::Generator s, const WeylWord& ww);

gradings::Status::Value status(const KGB_base& kgb, KGBElt x,
			       const RootSystem& rs, RootNbr alpha);

} // |namespace kgb|

} // |namespace atlas|


#endif
