/*!
\file
\brief Class definition and function declarations for the class KGB
representing orbits of K on G/B.
*/
/*
  This is kgb.h

  Copyright (C) 2004,2005 Fokko du Cloux
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

/******** constant declarations *********************************************/

const KGBElt UndefKGB = ~0ul;


/******** type definitions **************************************************/

//! \brief per KGB element information
struct KGBInfo
{
  gradings::Status status; ///< status of each simple root for this element
  unsigned int length; ///< dimension of the K orbit on G/B, minus minimal one
  unsigned int cartan; ///< records to which Cartan class this element belongs
  Descent desc; ///< flags which simple reflections give a descent

  KGBInfo(unsigned int l, unsigned int c) : status(),length(l),cartan(c),desc()
  {} // set length explicitly, both other fields set to their default value
};



/*!
\brief Represents the orbits of K on G/B for a particular real form.

Each orbit x defines an involution \f$\theta_x\f$ of H, represented by a
twisted involution in the Weyl group. In fact each orbit is characterized by
an element of the extended Weyl group lying over that twisted involution; thus
collection of orbits defining the same involution is parametrized using (the
equivalent of) elements of the fiber group. The elements occuring in the full
KGB structure for one involution form an orbit of the imaginary Weyl group.

These orbits are needed first of all for the parametrization of
irreducible representations of the real form (see the class Block).

In this class, an orbit is represented by KGBElt, which is a number specifying
the position of the orbit on a list. For each number, the involution
\f$\theta_x\f$ is retained (as |d_involution[KGBElt]|), but not the fiber
information distinguishing different orbits with the same involution. Instead,
the class retains the cross action of (each simple reflection in) W on orbits,
and the Cayley transform. Each of these is stored as a vector (indexed by
orbit numbers) of lists of KGBElt's, one for each simple reflection (often the
KGBElt UndefKGB = ~0 in the case of the Cayley transform). The twisted
involution and |KGBInfo| record associated to each element are stored as well,
as are the inverse Cayley transforms deduced from the forward ones. This
information is all that is used by the Kazhdan-Lusztig algorithm.

The actual construction of the orbit is carried out by the (no longer
derived!) class |KGBHelper| defined in the implementation module kgb.cpp. It
is that class that works with actual elements of the Tits group (which encode
both a twisted involution and an element of the associlated fiber group).
*/

class KGB {

  enum State { BruhatConstructed, NumStates };

  size_t d_rank; //!< semisimple rank
  size_t d_torus_rank; //!< full rank used in torus parts

/*!
 \brief The list d_cross[j] lists the images of the cross action
of simple root \#j on the KGB elements.
*/
  std::vector<KGBEltList> d_cross;

/*!
\brief The list |d_cayley[j]| lists the images of the Cayley transform
action of the simple root |\#j| on KGB elements, or |UndefKGB| if not defined
*/
  std::vector<KGBEltList> d_cayley;

/*! \brief The list d_inverseCayley[j] lists the images of the inverse Cayley
transform of simple root \#j on KGB elements.

For fixed |j| each KGB element can have up to two inverse Cayley images. For
those for which simple root \#j is not real, both components will be
|UndefKGB|, and otherwise the second component might still be |UndefKGB|.
*/
  std::vector<KGBEltPairList> d_inverseCayley;

/*!
\brief Information about individual KGB elements
*/
  std::vector<KGBInfo> d_info;

/*!
\brief Twisted Weyl group element defining Cartan involution for KGB
element \#j.
*/
  tits::TitsEltList d_tits; // of size size()

  //!\brief to help find range of elements with fixed twisted involution
  std::vector<KGBElt> first_of_tau; // size: #distinct twisted involutions +1

  //!\brief tables to map twisted involutions to their sequence number
  weyl::TI_Entry::Pooltype d_pool;
  hashtable::HashTable<weyl::TI_Entry,unsigned int> d_tau;

/*!
\brief Bit 0 flags whether or not the Bruhat order on KGB has already
been constructed.
*/
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

/*! \brief Takes the Cayley transform of KGB element \#x in the direction of
simple root \#s; returns |UndefKGB| unless that root was noncompact imaginary.
*/
  KGBElt cayley(size_t s, KGBElt x) const {
    return d_cayley[s][x];
  }

/*!
  \brief Method that used to return whether involution(x) < involution(y).

  Explanation: the ordering is involution-length first, then weyl-length, then
  order by representative Weyl elements (weyl::TwistedInvolution::operator<).
  This is only a partial ordering, that does not distinguish elements of a
  fiber over one same twisted involution.

  A similar function is used to sort the elements of |KGB| upon construction,
  so this method should hold whenever |x<y| and |involution(x)!=involution(y)|
  and it turns out the method is never used, so I have commented it out. MvL.

  bool compare(KGBElt x, KGBElt y) const {
  if      (length(x) != length(y))
    return length(x) < length(y);
  else if (weylLength(x) != weylLength(y))
    return weylLength(x) < weylLength(y);
  else
    return involution(x) < involution(y);
}
*/

/*!
 \brief Applies the cross action of simple root \#s to KGB element \#x.
*/
  KGBElt cross(size_t s, KGBElt x) const {
    return d_cross[s][x];
  }

/*!
\brief Applies the inverse Cayley transform in simple root \#s to
KGB element \#x.

If simple root \#s is not real, returns a pair of UndefKGB.
*/
  KGBEltPair inverseCayley(size_t s, KGBElt x) const {
    return d_inverseCayley[s][x];
  }

/*!
\brief Semisimple rank of the underlying group.
*/
  size_t rank() const {
    return d_rank;
  }


/*!
\brief Number of K orbits on G/B.
*/
  size_t size() const {
    return d_info.size();
  }


/*! \brief Returns a Bitset flagging the simple roots which are descents for
KGB element \#x.
*/
  const Descent& descent(KGBElt x) const {
    return d_info[x].desc;
  }
/*!
\brief Length (dimension of the K orbit on G/B, minus minimal orbit dimension)
for KGB element \#x.
*/
  size_t length(KGBElt x) const {
    return d_info[x].length;
  }

/*!
\brief Cartan class associated to KGB element \#x.
*/
  size_t Cartan_class(KGBElt x) const {
    return d_info[x].cartan;
  }



/*!
\brief Flag telling whether each simple root is real, complex, imaginary
compact, or imaginary noncompact for KGB element \#x.
 */
  const gradings::Status& status(KGBElt x) const {
    return d_info[x].status;
  }


/*!
\brief Tells whether simple root \#s is real, complex, imaginary compact,
or imaginary noncompact for KGB element \#x.
 */
  gradings::Status::Value status(size_t s, KGBElt x) const {
    return status(x)[s];
  }

/*!
  \brief Tells whether simple root \# s is a descent generator for KGB
  element \#x.

  This only depends on the twisted involution \f$\tau\f$ associated to |x|: it
  holds whenever the simple reflection for |s| either makes \f$\tau\f$ shorter
  by twisted conjugation (|s| is a complex descent for |tau|), or if it
  twisted commutes by left multiplication (|s| is real for |tau|).
*/
bool isDescent(size_t s, KGBElt x) const
{
  return d_info[x].desc.test(s);
}

/*!
  \brief Tells whether simple root \#s is an ascent generator for KGB
  element \#x.

  This does not depend only on the twisted involution: it holds for all
  complex roots that are not descents, but for imaginary roots only if they
  noncompact.
*/
bool isAscent(size_t s, KGBElt x) const
{
  return not isDescent(s,x)
     and not (status(s,x) == gradings::Status::ImaginaryCompact);
}

//! \brief The based Tits group.
  const tits::BasedTitsGroup& basedTitsGroup() const { return *d_base; }
//! \brief The Tits group.
  const tits::TitsGroup& titsGroup() const { return d_base->titsGroup(); }
//! \brief The Weyl group.
  const weyl::WeylGroup& weylGroup() const
    { return d_base->titsGroup().weylGroup(); }
//! \brief The twisted Weyl group.
  const weyl::TwistedWeylGroup& twistedWeylGroup() const
    { return d_base->titsGroup(); } // in fact its base object

/*!
  \brief Returns the root datum involution corresponding to x.

  In fact, returns the twisted involution $w$ with \f$w.\delta = \tau\f$.
*/
  weyl::TwistedInvolution involution(KGBElt x) const { return d_tits[x].tw(); }
  tits::TorusPart torus_part(KGBElt x) const
   { return d_base->titsGroup().left_torus_part(d_tits[x]); }
  tits::TitsElt titsElt(KGBElt x) const { return d_tits[x]; }
//! \brief The (non-semisimple) rank of torus parts.
  size_t torus_rank() const { return d_torus_rank; }


  KGBEltPair tauPacket(const weyl::TwistedInvolution&) const;


  gradings::Grading base_grading() const { return d_base->base_grading(); }

  size_t weylLength(KGBElt x) const {
    return weylGroup().length(involution(x).w());
  }

// manipulators

// Creates Hasse diagram for Bruhat order on KGB and returns reference to it
  bruhat::BruhatOrder& bruhatOrder()
  { fillBruhat(); return *d_bruhat; }

  const poset::Poset& bruhatPoset() // this creates full poset on demand
  { return bruhatOrder().poset(); }

// private methods
private:
  void fillBruhat();

}; // class KGB

} // namespace kgb

} // namespace atlas

#endif
