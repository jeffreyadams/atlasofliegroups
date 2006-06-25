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
class Helper (in the unnamed namespace of kgb, and therefore currently
undocumented by doxygen).  It is the Helper class that works with
actual elements of the Tits group (by means of fiber groups).  After
the construction, the fiber information is discarded.  What is
retained is only the Cayley transforms (always going from more compact
to less compact involutions) and cross actions.
  */
class KGB {

 protected:

  enum State { BruhatConstructed, NumStates };

  size_t d_rank;

  std::vector<KGBEltList> d_cross;
  std::vector<KGBEltList> d_cayley;
  std::vector<KGBEltPairList> d_inverseCayley;
  DescentList d_descent;
  std::vector<size_t> d_length;
  weyl::WeylEltList d_involution;
  gradings::StatusList d_status;

  bitset::BitSet<NumStates> d_state;

  bruhat::BruhatOrder* d_bruhat;

  const weyl::WeylGroup* d_weylGroup;

 public:

// constructors and destructors
  explicit KGB(size_t);

  explicit KGB(realredgp::RealReductiveGroup&);

  virtual ~KGB();

// copy, assignment and swap
  void swap(KGB&);

// accessors
  const bruhat::BruhatOrder& bruhatOrder() const {
    return *d_bruhat;
  }

  KGBElt cayley(size_t s, KGBElt x) const {
    return d_cayley[s][x];
  }

  bool compare(KGBElt, KGBElt) const;

  KGBElt cross(size_t s, KGBElt x) const {
    return d_cross[s][x];
  }

  const Descent& descent(KGBElt x) const {
    return d_descent[x];
  }

  KGBEltPair inverseCayley(size_t s, KGBElt x) const {
    return d_inverseCayley[s][x];
  }

  bool isAscent(size_t s, KGBElt x) const;

  bool isDescent(size_t s, KGBElt x) const;

  size_t length(KGBElt x) const {
    return d_length[x];
  }

  size_t rank() const {
    return d_rank;
  }

  size_t size() const {
    return d_length.size();
  }

  const gradings::Status& status(KGBElt x) const {
    return d_status[x];
  }

  gradings::Status::Value status(size_t s, KGBElt x) const {
    return d_status[x][s];
  }

  const weyl::WeylElt& involution(KGBElt) const;

  KGBEltPair tauPacket(const weyl::WeylElt&) const;

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
