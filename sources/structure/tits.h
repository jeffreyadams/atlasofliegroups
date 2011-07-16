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

#include "tags.h"
#include "ratvec.h"   // |RatWeight| contained in |TorusElement|
#include "arithmetic_fwd.h"

#include "atlas_types.h"

#include "bitvector.h" // contained in |TitsGrooup|
#include "prerootdata.h" // contained in |GlobalTitsGrooup|
#include "weyl.h"      // contained in |TitsElt|

#include "complexredgp_fwd.h"
#include "realredgp_fwd.h"

#include "cartanclass.h"

#include "subdatum_fwd.h" // not subdatum.h, which includes us

namespace atlas {

namespace tits {

/******** function declarations *********************************************/

  inline TorusElement exp_pi(const RatWeight r);
  inline TorusElement exp_2pi(const RatWeight r);

  // 2-subgroup by which each |TorusPart| at involution |inv| will be reduced
  subquotient::SmallSubspace fiber_denom(const WeightInvolution& inv);

  Grading
  square_class_grading_offset(const Fiber& f,
			      cartanclass::square_class csc,
			      const RootSystem& rs);


/******** type definitions **************************************************/


  /*!
\brief Element of (Z/2Z)^rank, representing an element of H(2).

  The cocharacter lattice X_*(T) of T is always represented by Z^rank; the map
  X_*(T)-> T given by lambda^vee \mapsto exp(pi i lambda^vee) has 2X_*(T) as
  kernel, and induces a bijection of X_*(T)/2X_*(T) with the group H(2) of
  elements of order 2; hence H(2) is represented by (Z/2Z)^rank
  */
  typedef SmallBitVector TorusPart;



/* The following classes are newer than |TitsElement|, |TitsGroup| and
   |BasedTitsGroup| below, and less efficient, but provide a framework in
   which te sets modelled by the latter classes can be fitted.
*/

/* An element of finite order in $H$. Externally this is like a rational vector
   modulo integers in the coordinates. Internally we store twice the numerator
   of the rational vector, which facilitates adding halves to its coordinates.
*/
class TorusElement
{
  RatWeight repr; // represents $\exp(\pi i repr)$, no factor 2!

  TorusElement(const RatWeight r,tags::UnnormalizedTag)
   : repr(r) {} // raw constructor, like exp(pi i r), but no modular reduction
 public:
  explicit TorusElement(size_t rank) : repr(rank) {} // identity torus element

  // $\exp(\pi\ii r )$ or $\exp(2\pi\ii r )$
  TorusElement(const RatWeight& r,bool two);

  // accessors

  // provide access to value, but convert to "mod Z^rank" not "mod 2Z^rank"
  RatWeight log_pi(bool normalize) const;
  RatWeight log_2pi() const;

  // more often it will be practical to have acces to that "mod 2Z^rank" form
  const RatWeight& as_Qmod2Z() const { return repr; }

  bool operator== (const TorusElement& a) const { return repr==a.repr; }
  bool operator!= (const TorusElement& a) const { return repr!=a.repr; }
  bool operator<  (const TorusElement& a) const { return repr<a.repr; }

  TorusElement operator +(const TorusElement& t) const;
  TorusElement operator -(const TorusElement& t) const;

  // this method is to be used only at weights |alpha| taking value +1 or -1
  bool negative_at(const Coweight& alpha) const
    { return repr.scalarProduct(alpha)%2!=0; }

  // evaluation giving rational number modulo 2 (|negative_at| iff equals 1)
  arithmetic::Rational evaluate_at(const Coweight& alpha) const;

  // manipulators

  TorusElement& operator+=(TorusPart v); // arg by value since it is small

  // the following method assumes |prd| is on dual side with respect to torus
  void simple_reflect(const PreRootDatum& prd, weyl::Generator s)
  { prd.simpleReflect(repr.numerator(),s); } // numerator is weight for |prd|
}; // |class TorusElement|

inline TorusElement exp_pi(const RatWeight r)
{ return TorusElement(r,false); }
inline TorusElement exp_2pi(const RatWeight r)
{ return TorusElement(r,true); }





class GlobalTitsElement
{
  friend class GlobalTitsGroup; // which does almost all of the work for us

  // element stored as $(t,w)$ is interpreted as $t * \sigma_w * \delta_1$
  TorusElement t;
  TwistedInvolution w;

 public:

/*!\brief Constructs the identity element in the group */
  explicit GlobalTitsElement(size_t rank) : t(rank),w() {}

  /*!\brief The canonical representative $\sigma_w$ of |w| in |Tits|. */
  GlobalTitsElement(const WeylElt& we,size_t rank) : t(rank),w(we) {}

  explicit GlobalTitsElement(const TorusElement& te) : t(te),w() {}

  GlobalTitsElement(const TorusElement& te, const WeylElt& we)
  : t(te),w(we) {}

  // we compute modulo $t.\sigma_w\delta\equiv\sigma_w\delta.t$
  // therefore no constructor to form $\sigma_w\delta.t$ is needed

// copy and assignment can be left to their defaults

// accessors

// both components are exposed as constant references
/*!\brief twisted involution whose fiber we are in */
  const TorusElement& torus_part() const { return t; }
  const TwistedInvolution& tw() const { return w; }

// a method defined without help of |GlobalTitsGroup| (and which ignores |w|)
  GlobalTitsElement simple_imaginary_cross
  (const RootDatum& dual_rd, // dual for pragmatic reasons
   RootNbr alpha) const; // any simple-imaginary root


// manipulators
  TorusElement& torus_part() { return t; } // one may modify just torus part


// STL obligatories
  bool operator== (const GlobalTitsElement& a) const
  { return w == a.w and t == a.t; }

  bool operator!= (const GlobalTitsElement& a) const
  { return w != a.w or t != a.t; }

  bool operator< (const GlobalTitsElement& a) const // comparison, for STL use
  { return w!=a.w ? w<a.w : t<a.t ; }

}; // |class GlobalTitsElement|



//			   |class GlobalTitsGroup|

/*
  |GlobalTitsGroup| is a support class for working with |GlobalTitsElement|
   values. The main purpose is to allow defining cross actions and (inverse)
   Cayley transforms of such elements, producing new ones. The modular
   reduction of the |TorusElement| component that could be applied after
   Cayley transforms (and should be for deciding equality) is not implemented
   here, since it depends on the twisted involution in a way that is more
   efficiently handled by tabulation than by on-the-fly computation (cf.
   |kgb::GlobalFiberData|).

   This class has two applications: the main one is to be able to manipulate a
   given |GlobalTitsElement|, but for testing purposes an initial application
   was to build a |GlobalTitsGroup| first (from an inner class), and then to
   generate "all" valid elements for it (command 'X'). For the latter purpose
   the |square_class_gen| and associated method are included, which permits
   listing a set of initial elements from which others can be deduced. For the
   main application however this field can be left an empty list, so as not to
   waste any time computing it during construction. This dichotomy should of
   course be reimplemented through a derived class.

   The main use involves manipulating |y| values for the construction of
   possibly non-integral blocks, in which case the perspective is from that
   side; for instance the matrix |delta_tr|, called transposed, will actually
   be one that acts on character lattice $X^*$. However, in the initial use
   this class manipulates |x| values, for comparison with the results of the
   KGB construction; the perspective is then "direct". Terminology such as
   "dual side", "imaginary roots" is relative to the current perspective.
 */
class GlobalTitsGroup : public TwistedWeylGroup
{
  const PreRootDatum simple; // from DUAL side, allows W action
  WeightInvolution delta_tr; // transposed distinguished involution
  std::vector<TorusPart> alpha_v; // |simple.roots()| reduced modulo 2
  const RatWeight half_rho_v;

  // a list of gradings of the imaginary simple roots generating square classes
  std::vector<Grading> square_class_gen;

// forbid copy and assignment
  GlobalTitsGroup(const GlobalTitsGroup&);
  GlobalTitsGroup& operator= (const GlobalTitsGroup&);

 public:
  //!\brief constructor for inner class (use after latter is fully constructed)
  GlobalTitsGroup(const ComplexReductiveGroup& G);

  //!\brief constructor for dual inner class (the side of following contructor)
  GlobalTitsGroup(const ComplexReductiveGroup& G,tags::DualTag);

  //!\brief constructor from subdatum (dual side) + involution information
  GlobalTitsGroup(const subdatum::SubSystem& sub,
  		  const WeightInvolution& theta,
		  WeylWord& ww); // output: expesses |-theta^t| for |sub|

  // accessors
  size_t semisimple_rank() const { return alpha_v.size(); }
  size_t rank() const { return simple.rank(); }
  const RatWeight& torus_part_offset () const
  { return half_rho_v; }
  //  { return RatWeight(root_datum.dual_twoRho(),4); }

  //!\brief Element m_\alpha of H(2) for simple coroot \#j.
  TorusPart m_alpha(size_t j) const { return alpha_v[j]; }

  Weight parent_simple_root(weyl::Generator s) const
  { return simple.roots()[s]; }
  Coweight parent_simple_coroot(weyl::Generator s) const
  { return simple.coroots()[s]; }

  //!\brief Reflection of |TorusElement|s defined by a twisted involution
  WeightInvolution involution_matrix(const WeylElt& tw) const;

  using TwistedWeylGroup::twisted; // overloaded in this class
  using TwistedWeylGroup::twist;   // idem

  TorusElement twisted(const TorusElement& x) const;
//void twist(TorusElement& x) const { x = twisted(x); }


  const std::vector<Grading>& square_class_generators() const
  { return square_class_gen; }


// methods that only access some |GlobalTitsElement|

  bool is_valid(const GlobalTitsElement& a) const; // whether strong invovlution
  bool has_central_square(GlobalTitsElement a) const; // idem (but by-value)
  bool is_valid(const GlobalTitsElement& a, // weaker condition: square being
		const subdatum::SubSystem& sub) const; // central in subgroup

  // determine status of simple root, if assumed imaginary
  bool compact(weyl::Generator s, const GlobalTitsElement& a) const
  { return a.torus_part().negative_at(simple.coroots()[s]); }
  // compute Cayley transform
  GlobalTitsElement Cayley(weyl::Generator s, GlobalTitsElement a) const
  { leftMult(a.w,s); return a; }
  // flag length-decreasing complex cross actions and inverse Cayley transforms
  bitset::RankFlags descents(const GlobalTitsElement& a) const;

  TorusElement theta_tr_times_torus(const GlobalTitsElement& a) const;

  GlobalTitsElement prod(const GlobalTitsElement& a,
			 const GlobalTitsElement& b) const;

  bool compact(const RootSystem& rs,
	       RootNbr alpha,
	       GlobalTitsElement a) const; // whether alpha compact

// methods that manipulate a |GlobalTitsElement|

  /*!\brief
  Twisted conjugates |a| by |sigma_alpha| where |alpha| is simple root \#s,
  returns length change in $\{-1,0,+1\}$.

  This is implemented only modulo conjugation by torus elements, so there is
  no effective difference with conjugation by the inverse of |sigma_alpha|
  */
  int cross_act(weyl::Generator s, GlobalTitsElement& a) const;
  int cross_act(const WeylWord& w, GlobalTitsElement& a) const;
  int cross_act(GlobalTitsElement& a,const WeylWord& w) const;
  GlobalTitsElement cross(const WeylWord& w, GlobalTitsElement a) const
  { cross_act(w,a); return a; }

  void add(TorusPart tp,GlobalTitsElement& a) const // |tp| by value: small
  { a.t += tp; }

  // add |rw| to |a.t|
  void add(const RatWeight& rw,GlobalTitsElement& a) const;

  // modify |a| to an inverse Cayley image by (real simple root) $\alpha_s$
  void do_inverse_Cayley(weyl::Generator s,GlobalTitsElement& a) const;

 private: // this exists for pragmatic reasons only; no reason to export it
  // multiply left by either (t,sigma_ww) or delta_1(t,sigma_ww)delta_1
  void left_mult(const TorusElement& t,
		 const WeylWord& ww,
		 bool do_twist, // whether $(t,ww)$ is conjugated by $\delta_1$
		 GlobalTitsElement& b) const;
}; // |class GlobalTitsGroup|




//			     |class SubTutsGroup|


// a class that also stores a torus part induced by the parent base involution
// used only to implement the 'embedding' command, whence it has few methods
class SubTitsGroup : public GlobalTitsGroup
{
  GlobalTitsGroup parent; // a second one!
  const subdatum::SubSystem& subsys;
  TorusElement t;
 public:
  //!\brief constructor from subdatum (dual side) + involution information
  SubTitsGroup(const ComplexReductiveGroup& G,
	       const subdatum::SubSystem& sub,
	       const WeightInvolution& theta,
	       WeylWord& ww); // output: expesses |-theta^t| for |sub|
  TorusElement base_point_offset(const TwistedInvolution& tw) const;
};



//			       |class TitsElt|


/* We define two main classes, |TitsElt| and |TitsGroup|, as for Weyl groups.
   A |TitsElt| value stores both a |WeylElt| value and a |TorusPart|,
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
   declare the |TitsGroup| class a |friend| of |TitsElt|, but that's fine.
*/

/* We define |TitsElt| first, so |TitsGroup| inlines can call its methods; as
   the opposite occurs also, we postpone those inlines of |titsElt|
*/

  /*!
\brief Represents an element of a Tits group.

An element is always written $t.\sigma_w$, with $\sigma_w$ the canonical lift
to the Tits group of a Weyl group element $w$, and \f$t \in H(2)\f$.
  */
class TitsElt
{

  friend class TitsGroup; // necessary, as |TitsGroup| manipulates |TitsElt|s


  /*! \brief Factor in H(2) of the Tits group element.  */
  TorusPart d_t; // in fact thought of as to the left of |d_w|

  /*!
\brief Canonical Weyl part of the Tits group element.

If the Weyl group element $w$ has a reduced decomposition $s_1,...,s_r$, then
the canonical representative is the product of the corresponding Tits group
elements \f$\sigma_1,...,\sigma_r\f$. Because the \f$\sigma_i\f$ satisfy the
braid relations, this canonical representative is independent of the choice of
reduced decomposition.
  */
  weyl::TI_Entry d_w; // use |TI_Entry| rather than |WeylElt|, for |TE_Entry|


 public:

// constructors and destructors

// All constructors take |TitsGroup| argument, though some do not need it

/*!\brief Constructs the identity element in the group */
  explicit TitsElt(const TitsGroup& Tits);

  /*!\brief The canonical representative $\sigma_w$ of |w| in |Tits|. */
  TitsElt(const TitsGroup& Tits,const WeylElt& w); // group defines rank

  TitsElt(const TitsGroup&, TorusPart t) // pure torus part
  :  d_t(t),d_w(WeylElt())
  {}

  TitsElt(const TitsGroup&, const WeylElt& w, TorusPart t);

  TitsElt(const TitsGroup&, TorusPart t, const WeylElt& w)
  : d_t(t),d_w(w)
  {}

// copy and assignment can be left to their defaults

// accessors

// only the Weyl group component is exposed as constant reference.

/*!\brief Canonical Weyl part of the Tits group element. */
  const WeylElt& w() const { return d_w; }

/* the same componenent under another name (to make it smell sweeter); note
   however that this returns a value, not a reference (to not make the
   assumtion that |TwistedInvolution| is identical to |WeylElt|)
   */
/*!\brief twisted involution represented by canonical Weyl part */
  TwistedInvolution tw() const { return TwistedInvolution(d_w); }


// for the rest, only equality tests can bypass any use of the |TitsGroup|
  bool operator== (const TitsElt& a) const
  { return d_w == a.d_w and d_t == a.d_t; }

  bool operator!= (const TitsElt& a) const
  { return d_w != a.d_w or d_t != a.d_t; }

  bool operator< (const TitsElt& a) const // comparison, for STL use
  { // test torus part first, which is easier
    return d_t != a.d_t ? d_t < a.d_t : d_w < a.d_w;
  }

/* exceptionally we expose the raw torus part to derived classes; this should
   only be used for non-mathematical purposes like computing a hash code */
 protected:
  TorusPart t() const { return d_t; }

// to use |hashCode| method, we also give access to |d_w| as a |TI_Entry|
  const weyl::TI_Entry& ti() const { return d_w; }

/* no public manipulators: any operation defined without using the Tits group
   object would expose the implementation in an inacceptable way */
}; // |class TitsElt|




//			      |class TitsGroup|

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

  Note on independence of choices: given a root \f$\alpha\f$ for $H$ in $G$,
  the corresponding homomorphism \f$\phi_\alpha\f$ from $SL(2)$ to $G$ is
  defined only up to conjugation by $H$. This means that the generator
  \f$\sigma_\alpha\f$ of the Tits group appears to be defined only up to
  multiplication by the image of the coroot \f$\alpha^\vee\f$. But we are
  fixing a pinning, which means in particular that the maps \f$\phi_\alpha\f$
  for \f$\alpha\f$ simple are canonically defined (by the requirement that the
  standard pinning of $SL(2)$ be carried to the pinning for $G$). This means
  that the generator \f$\sigma_\alpha\f$ (still for \f$\alpha\f$ simple) is
  canonically defined.
*/
class TitsGroup : public TwistedWeylGroup
{
  /*!
\brief Dimension of the Cartan H. This is the size of torus parts of elements.
  */
  size_t d_rank;


  /*!
\brief List of the images in character lattice mod 2 of the simple roots.

Regarded as elements of order two in the dual torus \f$T^\vee\f$, these are
the elements \f$m_\alpha^\vee\f$.
  */
  std::vector<TorusPart> d_simpleRoot;

  /*!
\brief List of the elements \f$m_\alpha\f$ (for \f$\alpha\f$ simple) in $H(2)$.
  */
  std::vector<TorusPart> d_simpleCoroot;

  /*!
\brief Transpose of the reduction mod 2 of the matrix of the defining
involution of the inner class.

Gives the action of the (transposed fundamental) involution \f$\delta\f$ on
$H(2)$, for twisting TorusPart values when commuting with \f$\delta\f$.
  */
  BinaryMap d_involution;


// copy and assignment
// reserve and implement when necessary
  TitsGroup(const TitsGroup&);
  TitsGroup& operator= (const TitsGroup&);

 public:

// constructors and destructors

  //!\brief Ordinary constructor for inner class
  TitsGroup(const RootDatum&,
	    const WeylGroup& W,
	    const WeightInvolution& d);

  //!\brief Constructor for semisimple adjoint group
  TitsGroup(const int_Matrix& Cartan_matrix,
	    const WeylGroup& W,
	    const weyl::Twist& twist);

  //!\brief Constructor for Tits group relative to a subsystem
  TitsGroup(const subdatum::SubSystem& sub,
	    const WeightInvolution& theta,
	    WeylWord& ww);

  //\brief Like a copy contructor, but reference |W| rather than share or copy
  TitsGroup(const TitsGroup& Tg, const WeylGroup& W)
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
  const size_t semisimple_rank() const
    { return TwistedWeylGroup::rank(); }

  //!\brief Element m_\alpha of H(2) for simple coroot \#j.
  TorusPart m_alpha(size_t j) const { return d_simpleCoroot[j]; }

  //!\brief Image in the character lattice mod 2 of simple root \#j.
  TorusPart dual_m_alpha(size_t j) const { return d_simpleRoot[j]; }

// methods only involving a |TorusPart|

 /*!\brief Applies to the element x of H(2) simple reflection s.

    This is the same thing as conjugating |x|, viewed as embedded in the Tits
    group, by \f$\sigma_s\f$, or equivalently by its inverse (as
    \f$\sigma_s^2\f$ commutes with $x$). This is used internally to compute
    the proper commutation relations between Weyl group and torus parts of a
    Tits element.
  */
  void reflect(TorusPart& x, weyl::Generator s) const
  { if (d_simpleRoot[s].dot(x))
      x += d_simpleCoroot[s];
  }

// convert between torus parts $x$ and $y$ for which $x.w=w.y$ in Tits group
  TorusPart push_across(TorusPart x, const WeylElt& w) const;
  TorusPart pull_across(const WeylElt& w, TorusPart y) const;

  using TwistedWeylGroup::twisted; // overloaded in this class
  using TwistedWeylGroup::twist;   // idem

  //!\brief Binary matrix*vector product to compute twist on torus part
  TorusPart twisted(const TorusPart x) const { return d_involution*x; }

  //!\brief In-place imperative version of |twisted(TorusPart x)|
  void twist(TorusPart x) const { x=twisted(x); }

  //!\brief Reflection of |TorusPart|s defined by a twisted involution
  BinaryMap involutionMatrix(const WeylWord& tw) const;


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

  void left_torus_reduce
    (TitsElt& a, const subquotient::SmallSubspace& V) const
  { a.d_t=V.mod_image(a.d_t); }

  void right_torus_reduce
    (TitsElt& a, const subquotient::SmallSubspace& V) const
  { a=TitsElt(*this,a.w(),V.mod_image(right_torus_part(a))); }

// manipulators (none)


}; // |class TitsGroup|

// postponed inlines of |TitsElt|
inline TitsElt::TitsElt(const TitsGroup& Tits)
: d_t(Tits.rank()), d_w(WeylElt())
{}

inline TitsElt::TitsElt(const TitsGroup& Tits,const WeylElt& w)
: d_t(Tits.rank()), d_w(w)
{}

inline TitsElt::TitsElt
  (const TitsGroup& Tits, const WeylElt& w, TorusPart t)
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
}; // |class TE_Entry|


/*!
The class |TitsCoset| augments the |TitsGroup| class with the choice of a
basic strong involution $\delta$, which is needed to define the additional
methods |simple_grading|, |basedTwistedConjugate| and |Cayley_transform|
(although the latter happens to have an implementation independent of \delta).

All |TitsElement| values $t.\sigma_w$ get reinterpreted as $t.\sigma_w.\delta$
We are no longer in a group, but in a coset of the group. Also the images of
the operations |basedTwistedConjugate| and |Cayley_transform| are only defined
modulo the equivalence relation generated by conjugation by torus elements. If
conjugation by $w.\delta$ acts as $\theta$ on the Lie algebra $h$ of $H$, then
the equivalence class of $a=t.\sigma_w.\delta$ is obtained by modifying $t$ by
any value in the direction of the image of $\theta-1$ (the $-1$ eigenspace of
$\theta$; in the concrete representation |TorusPart|, by any reduction modulo
2 of a $-\theta$ fixed vector). In practice this means most of all that any
torus element $x$ the elements $x.a$ and $\theta(x).a=a.x$ and are equivalent.
*/
class TitsCoset
{
  TitsGroup* my_Tits_group; // pointer indicating ownership
  const TitsGroup& Tg;

  /*! \brief
    Defines conjugation action of $\delta$ on root vectors $X_\alpha$ which
    must be an element of the root space for $\beta=twist(\alpha)$. In fact we
    insist that $X_\alpha$ is set to plus or minus $X_\beta$, the sign being
    given by |grading_offset[s]| where $\alpha=\alpha_s$. The fact that
    $\delta$ is a strong involution implies that the grading offset must always
    the same for $\alpha$ as for $twist(\alpha)$.

    In term of |TitsElement|s, this relation implies the "commutation relation"
    $\sigma_\alpha.\delta=m_\alpha^grading_offset[s].\delta.\sigma_\beta$
  */
  Grading grading_offset;

  const RootSystem& rs; // needed (only) for the |grading| method

 public:
  TitsCoset(const ComplexReductiveGroup& G,
	    Grading base_grading);

  TitsCoset(const ComplexReductiveGroup& G);// adjoint case

  TitsCoset(const ComplexReductiveGroup& G,
	    tags::DualTag);// dual adjoint case

  // build Tits coset for subdatum, incorporating adapted parent base grading
  TitsCoset(const subdatum::SubDatum& sub,
	    Grading parent_base_grading);

  // this constructor computes the inner class for |sub| defined by |theta|
  TitsCoset(const subdatum::SubSystem& sub,
	    const WeightInvolution& theta,
	    Grading parent_base_grading, // by value: small
	    WeylWord& ww);

  ~TitsCoset() { delete(my_Tits_group); }

  /* accessors */

  const TitsGroup& titsGroup() const { return Tg; }
  const WeylGroup& weylGroup() const { return Tg.weylGroup(); }

  bool hasTwistedCommutation
    (weyl::Generator s, const TwistedInvolution& tw)  const
    { return Tg.hasTwistedCommutation(s,tw); }

  Grading base_grading() const { return grading_offset; }

  bool is_valid(TitsElt a) const; // whether |a| may occur at all; by-value

  // whether simple root |s| noncompact at KGB element represented by |a|
  inline bool simple_grading(const TitsElt& a, size_t s) const;

  // operation defining cross action of simple roots
  inline void basedTwistedConjugate(TitsElt& a, size_t s) const;

  void basedTwistedConjugate(TitsElt& a, const WeylWord& w) const;
  void basedTwistedConjugate(const WeylWord& w, TitsElt& a) const;

  // operation defining Cayley transform in non-compact imaginary simple roots
  void Cayley_transform(TitsElt& a, size_t s) const
    { Tg.sigma_mult(s,a); } // set |a| to $\sigma_s.a$

  // inverse Cayley transform for real simple roots
  // this requires knowing the subspace by which torus part has been reduced
  void inverse_Cayley_transform(TitsElt& a, size_t s,
				const subquotient::SmallSubspace& mod_space)
    const;

  void left_torus_reduce
    (TitsElt& a, const subquotient::SmallSubspace& V) const
  { titsGroup().left_torus_reduce(a,V); }

  void right_torus_reduce
    (TitsElt& a, const subquotient::SmallSubspace& V) const
  { titsGroup().right_torus_reduce(a,V); }

  // conjugate Tits group element by $\delta_1$
  TitsElt twisted(const TitsElt& a) const;
  // also keep a version for torus part only
  tits::TorusPart twisted(const tits::TorusPart& t) const
  { return Tg.twisted(t);}

  // general grading of a KGB element at an imaginary root: whether noncompact
  bool grading(TitsElt a, RootNbr n) const; // |a| by-value

  // various methods that provide a starting KGB element for any Cartan class
  TitsElt naive_seed (ComplexReductiveGroup& G,
		      RealFormNbr rf, size_t cn) const;
  TitsElt grading_seed (ComplexReductiveGroup& G,
			RealFormNbr rf, size_t cn) const;

}; // |class TitsCoset|

// A richer version of |TitsCoset|, with more methods
class EnrichedTitsGroup : public TitsCoset
{
  // strong real form representative
  const cartanclass::StrongRealFormRep srf;

 public:
  EnrichedTitsGroup(const RealReductiveGroup&);

  cartanclass::fiber_orbit f_orbit() const { return srf.first; }
  cartanclass::square_class square() const { return srf.second; }

  // whether arbitrary root |n| is compact for |x| in fundamental fiber
  bool is_compact(const tits::TorusPart& x, RootNbr n) const
    { return not grading(TitsElt(titsGroup(),x),n); }

  TitsElt backtrack_seed (const ComplexReductiveGroup& G,
				RealFormNbr rf, size_t cn) const;

}; // |class EnrichedTitsGroup|

// inline definitions (inlined here because of the long comments)

/* The amazingly simple way to compute the grading at simple roots, for any
   twisted involution for which they are imaginary. In short, the grading of
   an element $\sigma_w.\delta$ with trivial torus part at any imaginary
   simple root $\alpha=alpha_s$ is given by $grading_offset[s]$, regardless of
   $w$ (provided it makes $\alpha$ imaginary); the case of a nontrivial torus
   part follows easily from this basic case. Result |true| means noncompact.

   Let $w$ be the twisted involution associated to the Tits element |a| (i.e.,
   its Weyl group part) and let $\alpha=\alpha_s$ be a simple root, with image
   $\beta=\alpha_t$ under the twist associated to $\delta$ (so $t=twisted(s)$,
   and $s=t$ if $\alpha$ is imaginary for $\delta$). The precondition is that
   $s$ is imaginary for $\sigma_w.\delta$ which means that $s.w.t=w$ in $W$
   and $l(s.w)>l(w)$, so that both members of $s.w=w.t$ are reduced. Therefore
   $\sigma_s\sigma_w=\sigma_w\sigma_t$ in the Tits group, and
   $\sigma_s\sigma_w\delta=m_s^{grading_offset[s]}.\sigma_w\delta.\sigma_s$.
   Then $X_\alpha$ has eigenvalue $(-1)^grading_offset[s]$ for conjugation by
   $\sigma_w.\delta$, in other words $grading_offset[s]$ is also the grading
   of $\sigma_w.\delta$ at $\alpha$. Left multiplication of $\sigma_w.\delta$
   by a torus element $x$ will modify this grading to the opposite parity if
   and only if $\alpha(t)=-1$; since this torus part $x$ is represented as a
   bitvector, we can take the scalar product of the reduction mod 2 of
   $\alpha$ with $x$ to compute this correction. The |!=|-operation, which on
   |bool| values means exclusive or, combines the values.

   Note that by using the a left torus factor |x| rather than a right factor,
   we need not refer to $\beta$, i.e., |(twisted(s))|, in the dot product.
*/
inline
bool TitsCoset::simple_grading(const TitsElt& a, size_t s) const
{
  return grading_offset[s] != Tg.dual_m_alpha(s).dot(Tg.left_torus_part(a));
}

/* Apart from the grading methods, this is another place where we use
   |grading_offset|. The task is to compute the conjugate of the strong
   involution $a.x_0.\delta$ by $\sigma_s$ or by $\sigma_s^{-1}=\sigma_s^3$.
   These differ by conjugation by the torus element $\sigma_s^2=m_s$, which is
   not a trivial operation at the Tits group level, but as mentioned stays
   withing the equivalence relation of our interpretation of Tits elements.
   Therefore only one operation is implemented, and it defines an involution
   at the level of the KGB structure, generating for varying $s$ a $W$ action.

   For the twisted involution $w$ this operation is just twisted conjugation,
   giving $w'=s.w.t$ where $t=twisted(s)$. For the entire Tits group part
   $a=x.\sigma_w$ twisted conjugation by $\sigma_s$ in the Tits group |Tg|
   gives an element $a'=x'.\sigma_{w'}$ such that $\sigma_s.a=a'.\sigma_t$.
   Finally $\sigma_s.x.\sigma_w.\delta=x'.\sigma_{w'}.\sigma_t.\delta
   =x'.\sigma_{w'}.\delta.\sigma_s.m_s^{grading_offset[s]}$. Because of the
   equivalence in interpreting |TitsElement| values, we may replace
   multiplication on the right by the torus element $m_t^{grading_offset[s]}$
   by multiplication on the left leading to the formula implemented below.
*/
inline
void TitsCoset::basedTwistedConjugate(TitsElt& a, size_t s) const
{
  Tg.twistedConjugate(a,s);
  if (grading_offset[s])
    Tg.left_add(Tg.m_alpha(s),a);
}


} // |namespace tits|

} // |namespace atlas|

#endif
