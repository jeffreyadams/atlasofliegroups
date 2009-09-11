/*!
\file
\brief Class definitions and function declarations for the classes
TitsGroup and TitsElt.
*/
/*
  This is tits.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef TITS_H  /* guard against multiple inclusions */
#define TITS_H

#include "tits_fwd.h"

#include "gradings_fwd.h"
#include "complexredgp_fwd.h"
#include "realredgp_fwd.h"
#include "tags.h"

#include "rootdata.h"
#include "bitvector.h"
#include "constants.h"
#include "latticetypes.h"
#include "weyl.h"
#include "realform.h"

#include "cartanclass.h"
#include "tori.h"

namespace atlas {

namespace tits {

/******** function declarations *********************************************/

  // 2-subgroup by which each |TorusPart| at involution |inv| will be reduced
  inline
  latticetypes::SmallSubspace fiber_denom(latticetypes::LatticeMatrix inv)
  {
    latticetypes::SmallBitVectorList b(tori::minusBasis(inv.transposed()));
    return latticetypes::SmallSubspace(b,inv.numRows());
  }

  gradings::Grading
  square_class_grading_offset(const cartanclass::Fiber& f,
			      cartanclass::square_class csc,
			      const rootdata::RootSystem& rs);


/******** type definitions **************************************************/

  /*!
\brief Element of (Z/2Z)^rank, representing an element of T(2).

  The cocharacter lattice X_*(T) of T is always represented by Z^rank; the map
  X_*(T)-> T given by lambda^vee \mapsto exp(pi i lambda^vee) has 2X_*(T) as
  kernel, and induces a bijection of X_*(T)/2X_*(T) with the group T(2) of
  elements of order 2; hence T(2) is represented by (Z/2Z)^rank
  */
  typedef latticetypes::SmallBitVector TorusPart;


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
  const weyl::WeylElt& w() const { return d_w; }

/* the same componenent under another name (to make it smell sweeter); note
   however that this returns a value, not a reference (we would have none if
   |weyl::TwistedInvolutio| were a distinct type from |weyl::WeylElt|)
*/
/*!\brief twisted involution represented by canonical Weyl part */
  const weyl::TwistedInvolution tw() const
    { return weyl::TwistedInvolution(d_w); }


// for the rest, only equality tests can bypass any use of the |TitsGroup|
  bool operator== (const TitsElt& a) const {
    return d_w == a.d_w and d_t == a.d_t;
  }

  bool operator!= (const TitsElt& a) const {
    return d_w != a.d_w or d_t != a.d_t;
  }

  bool operator< (const TitsElt& a) const // comparison for STL use
  { // test torus part first, which is easier
    return d_t != a.d_t ? d_t < a.d_t : d_w < a.d_w;
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
  torus parts are just elements of the $\Z/2\Z$-vector space $H(2)$.

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
class TitsGroup : public weyl::TwistedWeylGroup
{
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
  */
  latticetypes::BinaryMap d_involution;


// copy and assignment
// reserve and implement when necessary
  TitsGroup(const TitsGroup&);
  TitsGroup& operator= (const TitsGroup&);

 public:

// constructors and destructors

  //!\brief Ordinary constructor for inner class (with precomputed twist)
  TitsGroup(const rootdata::RootDatum&,
	    const weyl::WeylGroup& W,
	    const latticetypes::LatticeMatrix& d,
	    const weyl::Twist&);

  //!\brief Constructor for semisimple adjoint group
  TitsGroup(const latticetypes::LatticeMatrix& Cartan_matrix,
	    const weyl::WeylGroup& W,
	    const weyl::Twist& twist);

  //\brief Like a copy contructor, but reference |W| rather than share or copy
  TitsGroup(const TitsGroup& Tg, const weyl::WeylGroup& W)
    : TwistedWeylGroup(W,Tg.twist())
    , d_rank(Tg.d_rank)
    , d_simpleRoot(Tg.d_simpleRoot)
    , d_simpleCoroot(Tg.d_simpleCoroot)
    , d_involution(Tg.d_involution)
  {}

// All methods being accessors, we classify by behaviour w.r.t. TitsElt objects

// methods not involving |TitsElt|

  /*!\brief Rank of the torus. */
  const size_t rank() const { return d_rank; }

  //!\brief Element m_\alpha of T(2) for simple coroot \#j.
  TorusPart simpleCoroot(size_t j) const { return d_simpleCoroot[j]; }

  //!\brief Image in the character lattice mod 2 of simple root \#j.
  TorusPart simpleRoot(size_t j) const { return d_simpleRoot[j]; }

// methods only involving a |TorusPart|

// convert between torus parts $x$ and $y$ for which $x.w=w.y$ in Tits group
  TorusPart push_across(TorusPart x, const weyl::WeylElt& w) const;
  TorusPart pull_across(const weyl::WeylElt& w, TorusPart y) const;

  using TwistedWeylGroup::twisted; // overloaded in this class
  using TwistedWeylGroup::twist;   // idem

  //!\brief Binary matrix*vector product to compute twist on torus part
  TorusPart twisted(const TorusPart& x) const { return d_involution.apply(x); }

  //!\brief Reflection of |TorusPart|s defined by a twisted involution
  latticetypes::BinaryMap involutionMatrix(const weyl::WeylWord& tw) const;

  //!\brief In-place imperative version of |twisted(TorusPart x)|
  void twist(TorusPart& x) const { d_involution.apply(x,x); }

// methods that only access some |TitsElt|

  /*!\brief Length of the underlying Weyl group element. */
  unsigned long length(const TitsElt& a) const
    { return weylGroup().length(a.w()); }

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
  void twistedConjugate(TitsElt& a, weyl::Generator s) const
  { sigma_mult(s,a); mult_sigma_inv(a,twisted(s)); }

  // the inverse operation: twisted conjugation by $\sigma_s^{-1}$
  void inverseTwistedConjugate(TitsElt& a, weyl::Generator s) const
  { sigma_inv_mult(s,a); mult_sigma(a,twisted(s)); }

  void left_torus_reduce(tits::TitsElt& a, latticetypes::SmallSubspace V) const
    { a.d_t=V.mod_image(a.d_t); }

  void right_torus_reduce(tits::TitsElt& a, latticetypes::SmallSubspace V)
    const
  { a=TitsElt(*this,a.w(),V.mod_image(right_torus_part(a))); }

// manipulators (none)

 private:
// private methods

  /*!
\brief Applies to the element x of T(2) simple reflection s.

This is the same thing as conjugating |x|, viewed as embedded in the Tits
group, by \f$\sigma_s\f$, or equivalently by its inverse (as \f$\sigma_s^2\f$
commutes with $x$). This is used internally to compute the proper commutation
relations between Weyl group and torus parts of a Tits element.
  */
  void reflect(TorusPart& x, weyl::Generator s) const
  { if (d_simpleRoot[s].dot(x))
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
      return (ti().hashCode(modulus)+t().data().to_ulong()) & (modulus-1);
    }
}; // class TE_Entry


/*! \brief
Augments the |TitsGroup| class with the choice of a base point, which makes
possible the additional methods |simple_grading|, |basedTwistedConjugate| and
|Cayley_transform| (but the last does not actually use the base point choice).
*/
class BasedTitsGroup
{
  tits::TitsGroup* my_Tits_group; // pointer indicating ownership
  const tits::TitsGroup& Tg;

  /*! \brief
    Flags the noncompact imaginary roots for the basic strong involution,
    among the simple roots for G (only).

    This is an important parameter in the KGB generation. It depends on an
    implicitly chosen "basic strong involution" in the fundamental fiber,
    i.e., one of the form \f$x_0\delta\f$ with \f$x_0\in H\f$. The element
    $x_0$ is determined, up a factor in $Z(G)$, by the evaluations of simple
    roots at it; these are given by \f$\alpha_i(x_0)=(-1)^{g_i}\f$ for all
    $i$, where $g_i$ denotes |gradingOffset[i]|. When \f$\alpha_i\f$ is
    imaginary for \f$\delta\f$, the value $g_i$ is fixed by the choice of a
    $W_{im}$ orbit representative in the real form; otherwise (\f$\alpha_i\f$
    is complex for \f$\delta\f$) we choose $g_i=0$, which choice can be
    accommodated by $H$-conjugation of \f$x_0\delta\f$. Thanks to the latter
    convention, we are able to compute, using |grading_offset|, gradings at
    simple roots that are imaginary in the fiber given by \emph{any} twisted
    involution, even if those roots were complex in the fundamental fiber.
    (Although by continuity, the same gradings would remain valid if we would
    replace \f$x_0\delta\f$ by any $H$-conjugate.)
  */
  gradings::Grading grading_offset;

  const rootdata::RootSystem& rs;

 public:
  BasedTitsGroup(const complexredgp::ComplexReductiveGroup& G,
		 gradings::Grading base_grading);

  BasedTitsGroup(const complexredgp::ComplexReductiveGroup& G);// adjoint case

  BasedTitsGroup(const complexredgp::ComplexReductiveGroup& G,
		 tags::DualTag);// dual adjoint case

  ~BasedTitsGroup() { delete(my_Tits_group); }

  /* accessors */

  const tits::TitsGroup& titsGroup() const { return Tg; }

  gradings::Grading base_grading() const { return grading_offset; }

  // grading of KGB element represented by |a| at simple root |s|
  inline bool simple_grading(const tits::TitsElt& a, size_t s) const;

  // operation defining cross action of simple roots
  inline void basedTwistedConjugate(tits::TitsElt& a, size_t s) const;

  void basedTwistedConjugate(tits::TitsElt& a, const weyl::WeylWord& w) const;
  void basedTwistedConjugate(const weyl::WeylWord& w, tits::TitsElt& a) const;

  // operation defining Cayley transform in non-compact imaginary simple roots
  void Cayley_transform(tits::TitsElt& a, size_t s) const
    { Tg.sigma_mult(s,a); } // set |a| to $\sigma_s.a$

  // inverse Cayley transform for real simple roots
  // this requires knowing the subspace by which torus part has been reduced
  void inverse_Cayley_transform(tits::TitsElt& a, size_t s,
				const latticetypes::SmallSubspace& mod_space)
    const;

  void left_torus_reduce(tits::TitsElt& a, latticetypes::SmallSubspace V) const
  { titsGroup().left_torus_reduce(a,V); }

  void right_torus_reduce(tits::TitsElt& a, latticetypes::SmallSubspace V)
  { titsGroup().right_torus_reduce(a,V); }

  // conjugate Tits group element by $\delta_1$
  tits::TitsElt twisted(const tits::TitsElt& a) const;
  // also keep vversion for torus part only
  tits::TorusPart twisted(const tits::TorusPart& t) const
  { return Tg.twisted(t);}

  // general case of grading of a KGB element at an imaginary root
  bool grading(tits::TitsElt a, rootdata::RootNbr n) const;

  // various methods that provide a starting KGB element for any Cartan class
  tits::TitsElt naive_seed (complexredgp::ComplexReductiveGroup& G,
			    realform::RealForm rf, size_t cn) const;
  tits::TitsElt grading_seed (complexredgp::ComplexReductiveGroup& G,
			      realform::RealForm rf, size_t cn) const;

}; // BasedTitsGroup

// A richer version of |BasedTitsGroup|, with more methods
class EnrichedTitsGroup : public BasedTitsGroup
{
  // strong real form representative
  const cartanclass::StrongRealFormRep srf;

  // unique constructor is private; use the |for_square_class| method
  EnrichedTitsGroup(const realredgp::RealReductiveGroup&,
		    const cartanclass::Fiber&);


 public:
  static EnrichedTitsGroup
    for_square_class(const realredgp::RealReductiveGroup& GR);

  cartanclass::fiber_orbit f_orbit() const { return srf.first; }
  cartanclass::square_class square() const { return srf.second; }

  // whether arbitrary root |n| is compact for |x| in fundamental fiber
  bool is_compact(const tits::TorusPart& x, rootdata::RootNbr n) const
    { return not grading(TitsElt(titsGroup(),x),n); }

  tits::TitsElt backtrack_seed (const complexredgp::ComplexReductiveGroup& G,
				realform::RealForm rf, size_t cn) const;

}; // EnrichedTitsGroup

// inline definitions

/* The amazingly simple way to compute the grading at simple roots, for any
   twisted involution for which they are imaginary. Let $w$ be the twisted
   involution associated to the Tits element |a| (i.e., its Weyl group part)
   and let $\alpha=\alpha_s$ be a simple root, with image $\beta=\alpha_t$
   under the distinguished involution $\delta$ (so $t=twisted(s)$, and $s=t$
   if $\alpha$ is imaginary for $\delta$). The precondition is that $s$ is
   imaginary for $w$ which means that $s.w.t=w$ in $W$ and $l(s.w)>l(w)$, so
   that $s.w=w.t$ are both reduced. This means that this equation lifts to
   canonical Weyl elements in the Tits group, so in that group conjugation by
   $\sigma_w$ sends $\sigma_t$ to $\sigma_s$ (rather than its inverse).
   Therefore in $G$, the action of $Ad(\sigma_w.\delta)$ is the identity on
   the $SL(2)$ or $PSL(2)$ subgroup containing $\sigma_s$ whose Lie algebra
   contains $X_\alpha$ and $X_{-\alpha}$, and so those root spaces are fixed
   under the adjoint action of $\sigma_w.\delta$. Then the adjoint action of
   the strong involution $x.\sigma_w.x_0.\delta$ (which |a| represents) on the
   root vector $X_\alpha$ gives $(Ad(x).Ad(\sigma_w).Ad(x_0))(X_\beta)
   =\alpha(x)\beta(x_0)X_\alpha$. So the grading $g$ to be computed here
   should satisfy $(-1)^g=\alpha(x)\beta(x_0)$, or $g=\<\alpha,x>+\<beta,x_0>$
   in $\Z/2\Z$. The first term is computed by |scalarProduct|, while the
   second term is the grading of $x_0.\delta$ on the root $X_\beta=X_\alpha$
   if it is imaginary for $\delta$ (i.e., $s=t$), or $0$ if $X_\alpha$ and
   $X_beta$ are complex (thanks to the choice of $x_0$); either way,
   $\<beta,x_0>$ equals |grading_offset[s]|.

   Note that by using the a left torus factor |x| rather than a right factor,
   we need not refer to $\beta$, i.e., |simpleRoot(twisted(s))| below.
*/
bool BasedTitsGroup::simple_grading(const tits::TitsElt& a, size_t s) const
{
  return grading_offset[s] ^ Tg.simpleRoot(s).dot(Tg.left_torus_part(a));
}

/* Apart from the grading methods, this is another place where we use
   |grading_offset|. The task is to compute the conjugate of the strong
   involution $a.x_0.\delta$ by $\sigma_s$ or by $\sigma_s^{-1}=\sigma_s^3$.
   Note that since conjugation by $\sigma_s^2=m_s$ is not a trivial
   operation at the Tits group level, there is a difference between these
   operations, and they do not define an involution, but the difference can
   be seen to disappear into the modular reduction that is systematically
   applied to torus parts of Tits group elements; therefore this operation
   will define an involution at the level of the KGB structure. In fact we
   shall prefer to implement conjugation by $\sigma_s^{-1}$.

   For the twisted involution (the Weyl group part of |a|) this operation is
   just twisted conjugation by $s^{-1}=s$. For the entire Tits group part |a|,
   twisted conjugation by $\sigma_s^{-1}$ gives an element $a'$ such that
   $\sigma_s^{-1}a=a'\sigma_t^{-1}$ where $t=twisted(s)$, but this is not yet
   the whole story; in addition there is a contribution from interchanging
   $x_0$ and $\sigma_t$. Since $\sigma_t$ represents $t$ in the Weyl group,
   one has $Ad(\sigma_t)(x_0)=t(x_0)$, so we should add $\<\alpha_t,x_0>m_t$
   at the right to the torus part of |a|. By our choice of $x_0$, this scalar
   product can be nonzero only if $s=t$, and is equal to |grading_offset[s]|;
   moreover the value of $m_s$ in the coordinates of the torus part is
   available as |Tg.simpleCoroot(s)|. By the remark above about conjugating by
   $m_s$, the contribution may be added into the left torus part instead of
   the right torus part, which is what is done in the code below.
*/
void BasedTitsGroup::basedTwistedConjugate(tits::TitsElt& a, size_t s) const
{
  Tg.inverseTwistedConjugate(a,s);
  if (grading_offset[s]) // this implies |Tg.twisted(s)==s|
    Tg.left_add(Tg.simpleCoroot(s),a);
}


} // namespace tits

} // namespace atlas

#endif
