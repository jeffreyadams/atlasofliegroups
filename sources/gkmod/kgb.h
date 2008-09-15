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

/******** function declarations *********************************************/


/******** constant declarations *********************************************/

namespace kgb {

const KGBElt UndefKGB = ~0ul;

} // namespace kgb

/******** type definitions **************************************************/

namespace kgb {

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

/*! \brief
Augments the |TitsGroup| class with the choice of a base point, which makes
possible the additional methods |simple_grading|, |basedTwistedConjugate| and
|Cayley_transform| (but the last does not actually use the base point choice).
*/

class BasedTitsGroup
{
  const tits::TitsGroup& Tg;
 protected:
  const rootdata::RootDatum& rd;

  /*! \brief
    Flags the noncompact imaginary roots for the basic strong involution,
    among the simple roots for G.

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
  */
  gradings::Grading grading_offset;

 public:
  BasedTitsGroup(const realredgp::RealReductiveGroup&);

  /* accessors */

  const tits::TitsGroup& titsGroup() const { return Tg; }
  const rootdata::RootDatum& rootDatum() const { return rd; }
  gradings::Grading base_grading() const { return grading_offset; }

  // grading of KGB element represented by |a| at simple root |s|
  inline bool simple_grading(const tits::TitsElt& a, size_t s) const;

  // operation defining cross action of simple roots
  inline void basedTwistedConjugate(tits::TitsElt& a, size_t s) const;

  // operation defining Cayley transform in non-compact imaginary simple roots
  void Cayley_transform(tits::TitsElt& a, size_t s) const {
    Tg.sigma_mult(s,a); // set |a| to $\sigma_s.a$
  }

  // inverse Cayley transform for real simple roots
  void inverse_Cayley_transform(tits::TitsElt& a, size_t s) const {
    Tg.sigma_inv_mult(s,a); // set |a| to $\sigma_s.a$
  }

  // conjugate Tits group element by $\delta_1$
  tits::TitsElt twisted(const tits::TitsElt& a) const;

}; // BasedTitsGroup

// A richer version of |BasedTitsGroup|, with more methods
class EnrichedTitsGroup : public BasedTitsGroup
{
  // strong real form representative
  const cartanclass::StrongRealFormRep srf;
  // all compact imaginary roots for basic form
  const rootdata::RootSet base_compact;

  // unique constructor is private; use the |for_square_class| method
  EnrichedTitsGroup(const realredgp::RealReductiveGroup&,
		    const cartanclass::Fiber&);

 public:
  static EnrichedTitsGroup
    for_square_class(const realredgp::RealReductiveGroup& GR);

  cartanclass::fiber_orbit f_orbit() const { return srf.first; }
  cartanclass::square_class square() const { return srf.second; }

  // whether arbitrary root |n| is compact for |x| in fundamental fiber
  bool is_compact(const tits::TorusPart& x, rootdata::RootNbr n) const;

  // general case of grading of a KGB element at an imaginary root
  bool grading(tits::TitsElt a, rootdata::RootNbr n) const;

  // various methods that provide a starting KGB element for any Cartan class
  tits::TitsElt naive_seed (const complexredgp::ComplexReductiveGroup& G,
			    realform::RealForm rf, size_t cn) const;
  tits::TitsElt grading_seed (const complexredgp::ComplexReductiveGroup& G,
			      realform::RealForm rf, size_t cn) const;
  tits::TitsElt backtrack_seed (const complexredgp::ComplexReductiveGroup& G,
				realform::RealForm rf, size_t cn) const;

}; // EnrichedTitsGroup


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
  BasedTitsGroup* d_base; // pointer, because constructed in helper class

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
  const BasedTitsGroup& basedTitsGroup() const {
    return *d_base;
  }
//! \brief The Tits group.
  const tits::TitsGroup& titsGroup() const {
    return d_base->titsGroup();
  }
//! \brief The (twisted) Weyl group.
  const weyl::WeylGroup& weylGroup() const {
    return d_base->titsGroup().weylGroup();
  }
/*!
  \brief Returns the root datum involution corresponding to z.

  In fact, returns the twisted involution $w$ with \f$w.\delta = \tau\f$.
*/
  weyl::TwistedInvolution involution(KGBElt x) const {
    return d_tits[x].tw();
  }
  tits::TorusPart torus_part(KGBElt x) const {
    return d_base->titsGroup().left_torus_part(d_tits[x]);
  }
  tits::TitsElt titsElt(KGBElt x) const {
    return d_tits[x];
  }
//! \brief The (non-semisimple) rank of torus parts.
  size_t torus_rank() const {
    return d_base->rootDatum().rank();
  }


  KGBEltPair tauPacket(const weyl::TwistedInvolution&) const;


  gradings::Grading base_grading() const { return d_base->base_grading(); }

  size_t weylLength(KGBElt z) const {
    return weylGroup().length(involution(z).w());
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

// inline definitions

/* The amazingly simple way to compute the grading at simple roots, for any
   twisted involution for which they are imaginary. Let $w$ be the twisted
   involution associated to the Tits element |a| (i.e., its Weyl group part)
   and let $\alpha=\alpha_s$ be a simple root, with image $\beta=\alpha_t$
   under the distinguished involution $\delta$ (so $t=twist(s)$, and $s=t$ if
   $\alpha$ is imaginary for $\delta$). The precondition is that $s$ is
   imaginary for $w$ which means that $s.w.t=w$ in $W$ and $l(s.w)>l(w)$, so
   that $s.w=w.t$ are both reduced. This means that this equation lifts to
   canonical Weyl elements in the Tits group, so in that group conjugation by
   $\sigma_w$ sends $\sigma_t$ to $\sigma_s$ (rather than its inverse).
   Therefore in $G$, the action of $Ad(\sigma_w.\delta)$ is the identity on
   the $SL(2)$ or $PSL(2)$ subgroup containing $\sigma_s$ whose Lie algebra
   contains $X_\alpha$ and $X_alpha$, and so those roots are fixed under the
   adjoint action of $\sigma_w.\delta$. Then the adjoint action of the strong
   involution $x.\sigma_w.x_0.\delta$ (which |a| represents) on the root
   vector $X_\alpha$ gives $(Ad(x).Ad(\sigma_w).Ad(x_0))(X_\beta)
   =\alpha(x)\beta(x_0)X_\alpha$. So the grading $g$ to be computed here
   should satisfy $(-1)^g=\alpha(x)\beta(x_0)$, or $g=\<\alpha,x>+\<beta,x_0>$
   in $\Z/2\Z$. The first term is computed by |scalarProduct|, while the
   second term is the grading of $x_0.\delta$ on the root $X_\beta=X_\alpha$
   if it is imaginary for $\delta$ (i.e., $s=t$), or $0$ if $X_\alpha$ and
   $X_beta$ are complex (thanks to the choice of $x_0$); either way,
   $\<beta,x_0>$ equals |grading_offset[s]|.

   Note that by using the a left torus factor |x| rather than a right
   factor, we need not refer to $\beta$, i.e., |simpleRoot(twist(s))| below.
*/
bool BasedTitsGroup::simple_grading(const tits::TitsElt& a, size_t s) const
{
  return grading_offset[s]
    ^bitvector::scalarProduct(Tg.left_torus_part(a),Tg.simpleRoot(s));
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
   $\sigma_s^{-1}a=a'\sigma_t^{-1}$ where $t=twist(s)$, but this is not yet
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
  if (grading_offset[s]) // this implies |Tg.twist(s)==s|
    Tg.left_add(Tg.simpleCoroot(s),a);
}

} // namespace kgb

} // namespace atlas

#endif
