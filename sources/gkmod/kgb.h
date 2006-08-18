/*!
\file
\brief Class definition and function declarations for the class KGB
representing orbits of K on G/B. 
*/
/*
  This is kgb.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef KGB_H  /* guard against multiple inclusions */
#define KGB_H

#include "kgb_fwd.h"

#include "bruhat_fwd.h"
#include "realredgp_fwd.h"
#include "weyl_fwd.h"

#include "gradings.h"

namespace atlas {

/******** function declarations *********************************************/

/******** constant declarations *********************************************/

namespace kgb {

const KGBElt UndefKGB = ~0ul;

}

/******** type definitions **************************************************/

namespace kgb {

  /*!
\brief Represents the orbits of K on G/B for a particular real form.

Each orbit x defines an involution theta_x of H, coming from the
extended Weyl group.  The collection of orbits defining the same
involution is parametrized using the Fiber class: in the end it is a
particular orbit of the imaginary Weyl group on the fiber group.

These orbits are needed first of all for the parametrization of
irreducible representations of the real form (see the class Block).

In the class, an orbit is represented by KGBElt, which is a number
specifying the position of the orbit on a list.  For each number, the
involution theta_x is retained (as d_involution[KGBElt]), but not the fiber information
distinguishing different orbits with the same involution. Instead, the
class retains the cross action of (each simple reflection in) W on
orbits, and the Cayley transform.  Each of these is stored as a vector
(indexed by orbit numbers) of lists of KGBElt's, one for each simple
reflection (often the KGBElt UndefKGB ~0 in the case of the Cayley
transform)

The actual construction of the orbit is carried out by the derived
class Helper.  It is the Helper class that works with actual elements
of the Tits group (by means of fiber groups).  After the construction,
the fiber information (which uses a lot of memory) is discarded.  What
is retained by the class kgb is only the Cayley transforms (always
going from more compact to less compact involutions) and cross
actions.  This information is all that is used by the Kazhdan-Lusztig
algorithm.

*/

class KGB {

 protected:

  enum State { BruhatConstructed, NumStates };

/*! 
\brief Semisimple rank of the underlying group.
*/
  size_t d_rank; 

/*!
 \brief The list d_cross[j] lists the images of the cross action
of simple root \#j on the KGB elements.
*/
  std::vector<KGBEltList> d_cross;

/*! 
\brief The list d_cayley[j] lists the images of the Cayley transform
action of the simple root \#j on KGB elements. 

If simple root \#j is not noncompact imaginary for KGB element \#x,
then element \#x of the list d_cayley[j] is 0, which
means that it is undefined: because Cayley transforms always move to
more split Cartans, kgb element \#0 (which is on the fundamental
Cartan) can never be a Cayley transform.  
*/ 
std::vector<KGBEltList> d_cayley;

/*!  
\brief The list d_inverseCayley[j] lists the images of the
inverse Cayley transform of simple root \#j on KGB elements.

If simple root \#j is not real for KGB element \#x, then element \#x of the
list d_inverseCayley[j] is UndefKGB.
*/
  std::vector<KGBEltPairList> d_inverseCayley;

/*! 
\brief The Bitset d_descent[j] flags the simple roots which are
descents for KGB element \#j.
*/
  DescentList d_descent;

/*! 
\brief Length (dimension of the K orbit on G/B, minus minimal orbit
dimension) for KGB element \#j.  
*/ 
  std::vector<size_t> d_length;

/*! 
\brief Twisted Weyl group element defining Cartan involution for KGB
element \#j.
*/
  weyl::WeylEltList d_involution;

/*! 
\brief The Status in d_status[j] tells whether each simple root is
real, complex, imaginary compact, or imaginary noncompact for KGB element \#j. 
 */
  gradings::StatusList d_status;

/*! 
\brief Bit 0 flags whether or not the Bruhat order on KGB has already
been constructed.
*/
  bitset::BitSet<NumStates> d_state;
  
/*! 
\brief Pointer to the Bruhat order on KGB.

The class BruhatOrder provides the symmetrized bitmatrix of the
partial order, and the Hasse diagram of covering relations.
*/
bruhat::BruhatOrder* d_bruhat;

/*! 
\brief Pointer to the (twisted) Weyl group. 
*/
  const weyl::WeylGroup* d_weylGroup;

 public:

// constructors and destructors
  explicit KGB(size_t);

  explicit KGB(realredgp::RealReductiveGroup&);

  virtual ~KGB();

// copy, assignment and swap
  void swap(KGB&);

// accessors

/*!
\brief Pointer to the Bruhat order on KGB.

The class BruhatOrder provides the symmetrized bitmatrix of the
partial order, and the Hasse diagram of covering relations.
*/
  const bruhat::BruhatOrder& bruhatOrder() const {
    return *d_bruhat;
  }

/*! 
\brief Takes the Cayley transform of KGB element \#x in the direction of simple root \#s.

If root \#s is not noncompact imaginary, then the returned value is 0,
which means that it is undefined: because Cayley transforms always
move to more split Cartans, kgb element \#0 (which is on the
fundamental Cartan) can never be a Cayley transform.
*/
  KGBElt cayley(size_t s, KGBElt x) const {
    return d_cayley[s][x];
  }

  bool compare(KGBElt, KGBElt) const;

/*!
 \brief Applies the cross action of simple root \#s to KGB element \#x.
*/
  KGBElt cross(size_t s, KGBElt x) const {
    return d_cross[s][x];
  }

/*! 
\brief Returns a Bitset flagging the simple roots which are
descents for KGB element \#x.
*/
  const Descent& descent(KGBElt x) const {
    return d_descent[x];
  }

/*! 
\brief Applies the inverse Cayley transform in simple root \#s to 
KGB element \#x.  

If simple root \#s is not real, returns a pair of UndefKGB.
*/
  KGBEltPair inverseCayley(size_t s, KGBElt x) const {
    return d_inverseCayley[s][x];
  }

  bool isAscent(size_t s, KGBElt x) const;

  bool isDescent(size_t s, KGBElt x) const;

/*! 
\brief Length (dimension of the K orbit on G/B, minus minimal orbit
dimension) for KGB element \#x.  
*/ 
  size_t length(KGBElt x) const {
    return d_length[x];
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
    return d_length.size();
  }

/*!
\brief Flag telling whether each simple root is
real, complex, imaginary compact, or imaginary noncompact for KGB element \#x. 
 */
  const gradings::Status& status(KGBElt x) const {
    return d_status[x];
  }


/*!
\brief Tells whether simple root \#s is
real, complex, imaginary compact, or imaginary noncompact for KGB element \#x. 
 */
  gradings::Status::Value status(size_t s, KGBElt x) const {
    return d_status[x][s];
  }

  const weyl::WeylElt& involution(KGBElt) const;

  KGBEltPair tauPacket(const weyl::WeylElt&) const;

/*! 
\brief The (twisted) Weyl group. 
*/
  const weyl::WeylGroup& weylGroup() const {
    return *d_weylGroup;
  }

  size_t weylLength(KGBElt) const;

// manipulators
  void fillBruhat(); // might throw

};

}

}

#endif
