/*!
\file
\brief Class definitions and function declarations for the classes
TitsGroup and TitsElt.
*/
/*
  This is tits.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef TITS_H  /* guard against multiple inclusions */
#define TITS_H

#include "tits_fwd.h"

#include "rootdata_fwd.h"

#include "bitvector.h"
#include "constants.h"
#include "latticetypes.h"
#include "weyl.h"

namespace atlas {

/******** type declarations *************************************************/

namespace tits {

  /*!
\brief Element of (Z/2Z)^rank, representing an element of T(2).

  The cocharacter lattice X_*(T) of T is always represented by Z^rank; the map
  X_*(T)-> T given by lambda^vee \mapsto exp(pi i lambda^vee) has 2X_*(T) as
  kernel, and induces a bijection of X_*(T)/2X_*(T) with the group T(2) of
  elements of order 2; hence T(2) is represented by (Z/2Z)^rank
  */
  typedef latticetypes::SmallBitVector TorusPart;

}

/******** function declarations *********************************************/

namespace tits {

  weyl::Twist makeTwist(const latticetypes::LatticeMatrix&,
			const rootdata::RootDatum&);

}

/******** type definitions **************************************************/

namespace tits {

/* We define two main classes, |TitsElt| and |TitsGroup|, as for Weyl groups.
   A |TitsElt| value stores both a |weyl::WeylElt| value and a |TorusPart|,
   which together specify the value. To compute with such elements however one
   needs additional element-independent data stored in the associated
   |TitsGroup| object, so many operations like multiplication are in fact
   methods of the latter class. The Tits group contains a normal subgroup
   isomorphic to $H(2)$ (the pure torus parts) for which the quotient is
   canonically the Weyl group $W$. The group is not a semidirect product,
   since $W$ does not lift to a subgroup, but each Weyl group element does
   have a canonical lift (multiplication of such elements stays within the set
   of canonical lifts if and only if the lengths add up). The values stored
   represent an element in the quotient and an element in the coset of its
   canonical lift. Computation is quite similar to what would be done for a
   semidirect product, computing separately in the quotient and the coset, but
   adjusting the computation in the coset in a way determined by the values in
   the quotient.

   An important design decision that was made for |TitsElt| is that the class
   should hide whether internally we represent a left or a right coset
   (actually this decision was retrofitted onto an implementation that
   initially did expose the choice it used, so that that the consequences of
   this choice were hard to weed out). This means that we cannot expose any
   form of reference to a "torus part", although we can export values computed
   as left or right coset components. For the Weyl group part we can export a
   constant reference, but not a variable one since changing it cannot leave
   both left and right coset components unchanged. Also all computations with
   |TitsElt| values need some access to the |TitsGroup| (most often by being
   methods of that class), even though some would not need to use it depending
   on the implementation. In order to make this possible, we are obliged to
   declare the |TitsGroup| class a |friend| of |TitsElt|.
*/

/* We define |TitsElt| first, so |TitsGroup| inlines can call its methods; as
   the opposite occurs also, we postpone those inlines of |titsElt|
*/

  /*!
\brief Represents an element of a Tits group.

An element is always written $w.t$, with $w$ the canonical representative
in the Tits group of a Weyl group element $w$ and \f$t \in T(2)\f$.
  */
class TitsElt {

  friend class TitsGroup; // necessary, as |TitsGroup| manipulates |TitsElt|s

 private:
  /*!
\brief Factor in T(2) of the Tits group element.
  */
  TorusPart d_t; // in fact thought of as to the left of |d_w|

  /*!
\brief Canonical Weyl part of the Tits group element.

If the Weyl group element $w$ has a reduced decomposition $s_1,...,s_r$, then
the canonical representative is the product of the corresponding Tits group
elements \f$\sigma_1,...,\sigma_r\f$. Because the \f$\sigma_i\f$ satisfy the
braid relations, this canonical representative is independent of the choice of
reduced decomposition.
  */
  weyl::TI_Entry d_w; // use |TI_Entry| rather than |WeylElt| for |TE_Entry|


 public:

// constructors and destructors

// All constructors take |TitsGroup| argument, though some do not need it

/*!\brief Constructs the identity element in the group */
  explicit TitsElt(const TitsGroup& Tits);

  /*!\brief The canonical representative of |w| in |Tits|. */
  TitsElt(const TitsGroup& Tits,const weyl::WeylElt& w); // group defines rank

  TitsElt(const TitsGroup&, TorusPart t) // pure torus part
  :  d_t(t),d_w(weyl::WeylElt())
  {}

  TitsElt(const TitsGroup&, const weyl::WeylElt& w, TorusPart t);

  TitsElt(const TitsGroup&, TorusPart t, const weyl::WeylElt& w)
  : d_t(t),d_w(w)
  {}

// copy and assignment can be left to their defaults

// accessors

// only the Weyl group component is exposed as constant reference.

/*!\brief Canonical Weyl part of the Tits group element. */
  const weyl::WeylElt& w() const {
    return d_w;
  }

// the same componenent under another name (to make it smell sweeter)

/*!\brief twisted involution represented by canonical Weyl part */
  const weyl::TwistedInvolution tw() const {
    return weyl::TwistedInvolution(d_w);
  }


// for the rest, only equality tests can bypass any use of the |TitsGroup|
  bool operator== (const TitsElt& a) const {
    return d_w == a.d_w and d_t == a.d_t;
  }

  bool operator!= (const TitsElt& a) const {
    return d_w != a.d_w or d_t != a.d_t;
  }

/* exceptionally we expose the raw torus part to derived classes; this should
   only be used for non-mathematical purposes like computing a hash code */
 protected:
  TorusPart t() const { return d_t; }

// to use |hashCode| method, we also give access to |d_w| as a |TI_Entry|
  const weyl::TI_Entry& ti() const { return d_w; }

 public:
/* no public manipulators: any operation defined without using the Tits group
   object would expose the implementation in an inacceptable way */
}; // class TitsElt

/*!\brief
  Represents a finite subgroup of the normalizer in $G$ of the Cartan $H$.

  We use a slight variant of the Tits group (also called extended Weyl
  group) as defined in J. Tits, J. of Algebra 4 (1966), pp. 96-116.

  The slight variant is that we include all elements of order two in the
  torus, instead of just the subgroup generated by the \f$m_\alpha\f$ (denoted
  \f$h_\alpha\f$ in Tits' paper.) Tits' original group may be defined by
  generators \f$\sigma_\alpha\f$ for \f$\alpha\f$ simple, subject to the braid
  relations and to \f$\sigma_\alpha^2= m_\alpha\f$; to get our group we just
  add a basis of elements of $H(2)$ as additional generators, and express the
  \f$m_\alpha\f$ in this basis. This makes for a simpler implementation, where
  totus parts are just elements of the $\Z/2\Z$-vector space $H(2)$.

  We have not tried to be optimally efficient here, as it is not expected that
  Tits computations will be significant computationally (dixit Fokko).

  Note on independence of choices: given a root \f$\alpha\f$ f $T$ in $G$,
  the corresponding homomorphism \f$\phi_\alpha\f$ from $SL(2)$ to $G$ is
  defined only up to conjugation by $T$. This means that the generator
  \f$\sigma_\alpha\f$ of the Tits group appears to be defined only up to
  multiplication by the image of the coroot \f$\alpha^\vee\f$. But we are
  fixing a pinning, which means in particular that the maps \f$\phi_\alpha\f$
  for \f$\alpha\f$ simple are canonically defined (by the requirement that the
  standard pinning of $SL(2)$ be carried to the pinning for $G$). This means
  that the generator \f$\sigma_\alpha\f$ (still for \f$\alpha\f$ simple) is
  canonically defined.
*/
class TitsGroup {

/*! \brief Diagram automorphism given by \f$\delta\f$  */
  weyl::Twist d_twist; // must come before |d_weyl|, needed in its construction

/*! \brief Owned reference to the underlying Weyl group. */
  weyl::WeylGroup& d_weyl;

  /*!
\brief Dimension of the Cartan T. This is the size of torus parts of elements.
  */
  size_t d_rank;


  /*!
\brief List of the images in character lattice mod 2 of the simple roots.

Regarded as elements of order two in the dual torus \f$T^\vee\f$, these are
the elements \f$m_\alpha^\vee\f$.
  */
std::vector<TorusPart> d_simpleRoot;

  /*!
\brief List of the elements \f$m_\alpha\f$ (for \f$\alpha\f$ simple) in $T(2)$.
  */
  std::vector<TorusPart> d_simpleCoroot;

  /*!
\brief Transpose of the reduction mod 2 of the matrix of the defining
involution of the inner class.

Gives the action of the involution \f$\delta\f$ on $T(2)$, for computing in
the \f$\delta\f$ coset of the Tits group.

This data member is currently unused (|d_twist| suffices), but the correct
value is available for future use in methods [MvL 19 June 2007]
  */
  latticetypes::BinaryMap d_involution;


// copy and assignment
// reserve and implement when necessary
  TitsGroup(const TitsGroup&);
  TitsGroup& operator= (const TitsGroup&);

 public:

// constructors and destructors
  TitsGroup(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~TitsGroup() { delete &d_weyl; } // Weyl group was owned

// All methods being accessors, we classify by behaviour w.r.t. TitsElt objects

// methods not involving |TitsElt|

  /*!\brief Rank of the torus. */
  const size_t rank() const { return d_rank; }

  /*!\brief Element m_\alpha of T(2) for simple coroot \#j. */
  TorusPart simpleCoroot(size_t j) const {
    return d_simpleCoroot[j];
  }

  /*!\brief Image in the character lattice mod 2 of simple root \#j. */
  TorusPart simpleRoot(size_t j) const {
    return d_simpleRoot[j];
  }

  /*!\brief Image under inner class diagram involution of node \#j. */
  size_t twist(size_t j) const {
    return d_twist[j];
  }

// methods only involving a |TorusPart|

// convert between torus parts $x$ and $y$ for which $x.w=w.y$ in Tits group
  TorusPart push_across(TorusPart x, const weyl::WeylElt& w) const;
  TorusPart pull_across(const weyl::WeylElt& w, TorusPart y) const;

// methods that only access some |TitsElt|

  /*!\brief Length of the underlying Weyl group element. */
  unsigned long length(const TitsElt& a) const {
    return d_weyl.length(a.w());
  }

  TorusPart left_torus_part(const TitsElt& a) const { return a.d_t; }

  TorusPart right_torus_part(const TitsElt& a) const
  { return push_across(a.d_t,a.d_w); }


// methods that manipulate a |TitsElt|

// left and right multiplication by $\sigma_s$ and by its inverse
  void sigma_mult(weyl::Generator s,TitsElt& a) const;
  void sigma_inv_mult(weyl::Generator s,TitsElt& a) const;
  void mult_sigma(TitsElt&, weyl::Generator) const;
  void mult_sigma_inv(TitsElt&, weyl::Generator) const;

  TitsElt prod(const TitsElt& a, TitsElt b) const;               // "a * b"
  void mult(TitsElt& a, const TitsElt& b) const { a=prod(a,b); } // "a *= b"

  void right_add(TitsElt& a,TorusPart t) const
  { a.d_t += pull_across(a.d_w,t); }
  void left_add (TorusPart t,TitsElt& a) const
  { a.d_t += t; }

/*!\brief Conjugate |a| by \f$\sigma_\alpha\f$, where \f$\alpha=\alpha_s\f$. */
  void sigma_conjugate(TitsElt& a, weyl::Generator s) const
  {
    sigma_mult(s,a); mult_sigma_inv(a,s);
  }
// the inverse operation of |sigma_conjugate(a,s)|
  void sigma_inv_conjugate(TitsElt& a, weyl::Generator s) const
  {
    sigma_inv_mult(s,a); mult_sigma(a,s);
  }




  /*!
\brief Twisted conjugates the TitsElt |a| by |sigma_alpha| where |alpha| is
  simple root \#s.

This corresponds to conjugation of the coset a.delta, with delta the defining
involution of the inner class. Note that the inverse of the generator
\sigma_\alpha is \sigma_\alpha.m_\alpha. Therefore this operation is
\emph{not} in general an involution with respect to |a| for fixed |s|.
However, if for elements whose Weyl group part is a twisted involution the
torus parts are projected to the fiber group over that twisted involution, as
is done in the KGB construction, it induces an involution on the quotient set.
  */
  void twistedConjugate(TitsElt& a, weyl::Generator s) const {
    sigma_mult(s,a); mult_sigma_inv(a,twist(s));
  }

  // the inverse operation: twisted conjugation by $\sigma_s^{-1}$
  void inverseTwistedConjugate(TitsElt& a, weyl::Generator s) const {
    sigma_inv_mult(s,a); mult_sigma(a,twist(s));
  }


  const weyl::WeylGroup& weylGroup() const {
    return d_weyl;
  }

// manipulators (none)

 private:
// private methods

  /*!
\brief Applies to the element x of T(2) simple reflection s.

This is the same thing as conjugating |x|, viewed as embedded in the Tits
group, by \f$\sigma_s\f$, or equivalently (since $\sigma_s^2$ commute with
$x$) by its inverse. This is used internally to compute the proper commutation
relations between Weyl group and torus parts of a Tits element.
  */
  void reflect(TorusPart& x, weyl::Generator s) const {
    if (bitvector::scalarProduct(x,d_simpleRoot[s]))
      x += d_simpleCoroot[s];
  }


}; // class TitsGroup

// postponed inlines of |TitsElt|
inline TitsElt::TitsElt(const TitsGroup& Tits)
: d_t(Tits.rank()), d_w(weyl::WeylElt())
{}

inline TitsElt::TitsElt(const TitsGroup& Tits,const weyl::WeylElt& w)
: d_t(Tits.rank()), d_w(w)
{}

inline TitsElt::TitsElt
  (const TitsGroup& Tits, const weyl::WeylElt& w, TorusPart t)
: d_t(Tits.pull_across(w,t)), d_w(w)
{}

struct TE_Entry // To allow hash tables of TitsElt values
  : public TitsElt
{
  TE_Entry(const TitsElt& t) : TitsElt(t) {}

  // members required for an Entry parameter to the HashTable template
  typedef std::vector<TE_Entry> Pooltype; // associated storage type
  size_t hashCode(size_t modulus) const // hash function
    {
      return ti().hashCode(modulus)+t().data().to_ulong() & modulus-1;
    }
}; // class TE_Entry


} // namespace tits

} // namespace atlas

#endif
