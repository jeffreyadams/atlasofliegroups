/*!
\file
\brief Class definition and function declarations for the class KGB
representing orbits of K on G/B.
*/
/*
  This is kgb.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2013 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef KGB_H  /* guard against multiple inclusions */
#define KGB_H

#include "atlas_types.h"

#include "gradings.h"	// containment in |KGBEltInfo|
#include "hashtable.h"	// containment in |KGB_base|
#include "weyl.h"       // |weyl::TI_Entry::Pooltype|
#include "tits.h"       // containment |GlobalTitsGroup|
#include "y_values.h"   // containment |TorusElement|

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
    gradings::Status status; // status of each simple root for this element
    DescentSet desc; //  which simple reflections are complex descent or real
    KGBElt dual; // number of Hermitian dual of this element

  EltInfo() : status(), desc(), dual(UndefKGB) {}

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
  HashTable<weyl::TI_Entry,unsigned int> inv_hash;

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

  KGBElt Hermitian_dual(KGBElt x) const { return info[x].dual; }

  const TwistedInvolution& nth_involution(unsigned int n) const
  { return inv_pool[n]; } // useful mostly in traversing all our involutions

  const TwistedInvolution& involution(KGBElt x) const // after construction only
  { return inv_pool[involution_index(x)]; } // the one associated to |x|

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

 private: // this internal index of involution is only visible to methods above
  InvolutionNbr involution_index(KGBElt x) const
  { return std::upper_bound(first_of_tau.begin(),first_of_tau.end(),x)
      -first_of_tau.begin() -1;
  }

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
  RatWeight fingerprint; // characterizes the torus element, modulo equivalence

  // obligatory fields for hashable entry
  typedef std::vector<KGB_elt_entry> Pooltype;
  size_t hashCode(size_t modulus) const; // hash function
  bool operator !=(const KGB_elt_entry& x) const; // ignores repr

  KGB_elt_entry (const RatWeight& f,
		 const GlobalTitsElement& y);
  GlobalTitsElement repr() const;

}; //  |struct KGB_elt_entry|

class global_KGB : public KGB_base
{
  const GlobalTitsGroup Tg;
  std::vector<GlobalTitsElement> elt;

  global_KGB(const global_KGB& org); // forbid copying

 public:
  global_KGB(ComplexReductiveGroup& G, bool dual_twist=false);

  global_KGB(ComplexReductiveGroup& G,
	     const GlobalTitsElement& x,
	     bool dual_twist=false); // generate KGB containing |x|

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
  void generate(size_t predicted_size, bool dual_twist);

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

  std::vector<TorusPart> left_torus_part; // of size |size()|
  BitSet<NumStates> d_state;

/*! \brief Owned pointer to the Bruhat order on KGB (or NULL).

The class BruhatOrder contains a Poset describing the full partial order,
and in addition the Hasse diagram (set of all covering relations).
*/
  BruhatOrder* d_bruhat;

  //! \brief Owned pointer to the based Tits group.
  TitsCoset* d_base; // pointer, because constructed late by |generate|

 public:

// constructors and destructors
  explicit KGB(RealReductiveGroup& GR,
	       const BitMap& Cartan_classes, bool dual_twist=false);

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

  RatWeight half_rho() const;

  TorusPart torus_part(KGBElt x) const { return left_torus_part[x]; }
  // reconstruct from |torus_part| a |TorusElement| as in |global_KGB|
  RatCoweight base_grading_vector() const; // offset for |torus_part_global|
  RatCoweight torus_part_global(KGBElt x) const; // will be $\theta^t$-fixed

  TitsElt titsElt(KGBElt x) const; // get KGB element |x| as a |TitsElt|
  size_t torus_rank() const; // the (non-semisimple) rank of torus parts.

  Grading base_grading() const { return d_base->base_grading(); }

  // as the name suggests, the following assumes |alpha| simple-imaginary
  bool simple_imaginary_grading(KGBElt x,RootNbr alpha) const
  { return d_base->simple_imaginary_grading(torus_part(x),alpha); }

  KGBElt lookup(TitsElt a) const; // by value


// manipulators

// Creates Hasse diagram for Bruhat order on KGB and returns reference to it
  BruhatOrder& bruhatOrder() { fillBruhat(); return *d_bruhat; }

  const poset::Poset& bruhatPoset(); // this creates full poset on demand

// virtual method
  virtual std::ostream& print(std::ostream& strm, KGBElt x) const;

// private methods
private:
  bool is_dual_twist_stable(const RealReductiveGroup& GR, TorusPart& shift)
    const; // auxiliary to see if this dual KGB can be twist-stable

  void fillBruhat();

}; // |class KGB|

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
