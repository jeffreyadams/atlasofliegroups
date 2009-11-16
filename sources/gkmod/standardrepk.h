/*!
\file
\brief Class definition and function declarations for the classes
StandardRepK and KhatContext.
*/
/*
  This is standardrepk.h

  Copyright (C) 2004, 2005 Fokko du Cloux
  Copyright (C) 2008, 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups version 0.2.4

  See file main.cpp for full copyright notice
*/

#ifndef STANDARDREPK_H  /* guard against multiple inclusions */
#define STANDARDREPK_H

#include <map>
#include <set>
#include <vector>
#include <iostream>

#include "bitset.h"
#include "bitvector_fwd.h"
#include "realredgp_fwd.h"
#include "latticetypes.h"
#include "realform.h"
#include "tits.h"
#include "kgb.h"
#include "hashtable.h"
#include "free_abelian.h"

namespace atlas {



/******** type declarations **************************************************/

namespace standardrepk {

class StandardRepK;
class HechtSchmid;
class PSalgebra;

// An image of a weight of the $\rho$-cover mod $(1-\theta)X^*$
typedef std::pair
  <latticetypes::LatticeElt, // free part, after subtraction of rho
   bitset::RankFlags         // torsion part, compact representation
  > HCParam;

// a linear combination of |StandardRepK|s
typedef free_abelian::Free_Abelian<StandardRepK> Char;
typedef Char::coef_t CharCoeff;

// a $K$-type formula; first component stands for its lowest $K$-type
typedef std::pair<StandardRepK,Char> CharForm;

typedef std::pair<latticetypes::LatticeElt,tits::TitsElt> RawRep;
typedef free_abelian::Free_Abelian<RawRep> RawChar;


typedef unsigned int seq_no;
typedef unsigned int level; // unsigned latticetypes::LatticeCoeff

// a comparison object that first takes into account a grading (height)
class graded_compare
{
  const std::vector<level>* h; // pointer, so |graded_compare| can be assignable

public:
  graded_compare (const std::vector<level>& a) : h(&a) {}

  bool operator()(seq_no x, seq_no y) const
  {
    level hx=(*h)[x], hy=(*h)[y];
    return hx!=hy ? hx<hy : x<y;
  }
}; // |class graded_compare|

typedef free_abelian::Free_Abelian<seq_no,graded_compare> combination;
typedef std::pair<seq_no,combination> equation;


// class handling rewriting of StandardRepK values by directed equations
// all left hand sides are assumed distinct, so just repeated substitutions
// in fact due to ordering, it now functions as a trivial lookup table
class SR_rewrites
{
public:

  typedef std::vector<combination> sys_t;

private:

  sys_t system;

public:
  SR_rewrites() : system() {}

  const combination& lookup(seq_no n) const; // an equation must exist!

  // manipulators
  void equate(seq_no n, const combination& rhs);

}; // |SR_rewrites|



/******** function declarations *********************************************/


  matrix::Matrix<CharCoeff> triangularize
    (const std::vector<equation>& system,
     std::vector<seq_no>& new_order);

  matrix::Matrix<CharCoeff> inverse_lower_triangular
    (const matrix::Matrix<CharCoeff>& U);


/******** type definitions **************************************************/


/*!
  \brief Represents the restriction to $K$ of a (coherently) continued
  standard Harish-Chandra module.

  This is a parameter type like Tits elements; the important operations are
  modifying and comparing values, not storing additional data that facilitate
  methods. For that, auxiliary classes |SRK_context| and |KhatContext|, which
  has a role similar to |WeylGroup| with respect to |WeylElt|, will be used

  For us a standard Harish-Chandra module is attached to
  1) a real Cartan subgroup $H(R)=H(R)_c H(R)_s$, with $H(R)_c = H(R)\cap K$
     the compact part, and $H(R)_s$ (a vector group) the noncompact part;
  2) a system $\Psi_{im]$ of positive imaginary roots for $H(R)$ in $G$;
  3) a system $\Psi_{re}$ of positive real roots for $H(R)$ in $G$; and
  4) a genuine character $\gamma$ of the rho-cover $\tilde H(R)$ (here
     genuine means that $\gamma$ is in $X^*+\rho$, so unless $\rho\in X^*$,
     the character does not descend to a character of $H(R)$).

  Action of the real Weyl group preserves the meaning of these data.

  (We said "attached to" rather than parametrized, as there are subtle
  identifications and relations, which are associated to the notions of
  being "standard" (rather than continued), "final", and "normalized").

  In the atlas picture, the Cartan and complete positive root system are
  always fixed, so one does not specify 2) and 3); instead the situation
  will be conjugated to one where the positive roots are the perpetual ones,
  and what changes is the strong involution $x$ representing the real form,
  and the position of $\gamma$ with respect to the simple coroots (it need
  not be dominant for all of them). As in the KGB module, strong involutions
  are represented by Tits group elements, the precise correspondence
  depending on a "base grading" stored elsewhere. In this class we will
  always assume that the involution $\theta_x$ on $H$ is distinguished
  within its class, using $W$-conjugacy where necessary to make it so. This
  means that except for intermediate computational values, the Weyl group
  part of the Tits element can be replaced by an indication of the number
  |d_cartan| of the real Cartan we are considering (the Weyl group part is
  then the twisted involution for the distinguished involution in that class).
  The (left) torus part must is explicitly represented as |d_fiberElt|; it
  implicitly determines the real form as in the Cartan class module,
  provided the central square class (which is stored elsewhere) is known.

  Since we are interested only in HC modules restricted to $K$, we are
  interested only in the character gamma restricted to $\tilde H(R)_c$.
  Characters of the compact group $H(R)_c$ are the same as algebraic
  characters of the complexification; that is

  $\widehat{H(R)_c}$ is identified with $X^* /(1-\theta)X^*$.

  At the level of the $\rho$-cover, we get

  ($\gamma$ restricted to $K$) is identified with an element of the
  coset-quotient $(X^* + \rho)/(1-\theta) X^*$.

  This is the information held in |d_lambda|. (The name lambda refers to the
  restriction of $\gamma$ to $\tilde H(R)_c$.) This is of type |HCParam|,
  consisting of an integer vector and an integer respesenting a bit vector;
  the integer vector gives the non-torsion part of $(X^*+\rho)/(1-\theta)X^*$
  on a basis of $(1/2)X^*$ held in |KhatContext|, and the bitvector gives the
  torsion part, via a basis also given there.

*/

class StandardRepK
{

  friend class SRK_context; // which is like |WeylGroup| is for |WeylElt|
  friend class KhatContext;

/*!
  \brief Number of the Cartan to which the HC module is associated.
*/
  size_t d_cartan;

// no real form or base grading recorded in elements; they're in |KhatContext|
/*!
  \brief Element of the fiber group; left torus part of the strong involution
*/
  tits::TorusPart d_fiberElt; // a SmallBitVector

/*!
  \brief Character of the rho cover of H^theta, on basis of $(1/2)X^*$
*/
  HCParam d_lambda;

// constructors, destructors, and swap

  // main constructor is private, used by |KhatContext| methods
  // fundamental bare-bones constructor; no status is set here
  StandardRepK(size_t cn, const tits::TorusPart& x, const HCParam& lambda)
    : d_cartan(cn), d_fiberElt(x), d_lambda(lambda) {}

public:

  StandardRepK() {} // so empty slots can be created

  void swap(const StandardRepK&);

  // accessors

  bool operator< (const StandardRepK&) const;
  bool operator== (const StandardRepK&) const;

  size_t Cartan() const { return d_cartan; }

// manipulators: none (done by friend class KhatContext)

// special members required by hashtable::HashTable

  typedef std::vector<StandardRepK> Pooltype;
  bool operator!=(const StandardRepK& another) const
    { return not operator==(another); }
  size_t hashCode(size_t modulus) const;

}; // class StandardRepK


//! \brief per Cartan class information for handling |StandardRepK| values
struct Cartan_info
{
  // projection matrix to torsion free part
  latticetypes::LatticeMatrix freeProjector;
  // projection matrix to torsion part, after rho-shift and reduction mod 2
  latticetypes::BinaryMap torsionProjector;

  // matrix used to lift free part of |HCParam| back to a weight
  latticetypes::LatticeMatrix freeLift;

  // list of vectors used to lift torsion part of |HCParam| to a weight
  latticetypes::WeightList torsionLift;

  // space that fiber parts are reduced modulo
  latticetypes::SmallSubspace fiber_modulus;

}; // |struct Cartan_info|

// a data type used for storing projections to facets of fundamental chamber
struct proj_info
{
  latticetypes::LatticeMatrix projection;
  latticetypes::LatticeCoeff denom;
}; // |struct proj_info|


// a wrapper around |bitset::RankFlags| to allow a hash table indexed by them
struct bitset_entry : public bitset::RankFlags
{
  typedef bitset::RankFlags base;
  typedef std::vector<bitset_entry> Pooltype;
  bitset_entry(base b) : base(b) {}
  size_t hashCode(size_t modulus) const { return to_ulong()&(modulus-1); }
}; // |struct bitset_entry|


// This class stores the information necessary to interpret a |StandardRepK|
class SRK_context
{
  complexredgp::ComplexReductiveGroup& G;
  const tits::BasedTitsGroup& Tg;  // for getting around in KGB (unused here)
  bitmap::BitMap Cartan_set;       // marks recorded Cartan class numbers
  std::vector<Cartan_info> C_info; // indexed by number of Cartan for |GR|

// this member is precomputed to increase efficiency of certain operations
  std::vector<latticetypes::BinaryMap> simple_reflection_mod_2; // dual side

 public:
  SRK_context(realredgp::RealReductiveGroup &G);

  // accessors
  complexredgp::ComplexReductiveGroup& complexGroup() const { return G; }
  const rootdata::RootDatum& rootDatum() const { return G.rootDatum(); }
  const weyl::WeylGroup& weylGroup() const { return G.weylGroup(); }
  const weyl::TwistedWeylGroup& twistedWeylGroup() const
    { return G.twistedWeylGroup(); }
  const tits::TitsGroup& titsGroup() const { return Tg.titsGroup(); }
  const tits::BasedTitsGroup& basedTitsGroup() const { return Tg; }

  const weyl::TwistedInvolution twistedInvolution(size_t cn) const
    { return G.twistedInvolution(cn); }
  const cartanclass::Fiber& fiber(const StandardRepK& sr) const
    { return G.cartan(sr.Cartan()).fiber(); }

  const Cartan_info& info(size_t cn) const
    { return C_info[Cartan_set.position(cn)]; }
  const latticetypes::BinaryMap& dual_reflection(weyl::Generator i) const
  { return simple_reflection_mod_2[i]; }

  //!\brief Projection |Weight| (in doubled coordinates) to |HCParam|
  HCParam project(size_t cn, latticetypes::Weight lambda) const; // by value

  //!\brief A section of |project|
  latticetypes::Weight lift(size_t cn, HCParam p) const;

  //!\brief (1+theta)* lifted value; this is independent of choice of lift
  latticetypes::Weight theta_lift(size_t cn, HCParam p) const
  {
    latticetypes::Weight result=lift(cn,p);
    result += complexGroup().cartan(cn).involution().apply(result);
    return result;
  }

  latticetypes::Weight lift(const StandardRepK& s) const
  { return lift(s.d_cartan,s.d_lambda); }

  latticetypes::Weight theta_lift(const StandardRepK& s) const
  { return theta_lift(s.d_cartan,s.d_lambda); }

  StandardRepK std_rep
    (const latticetypes::Weight& two_lambda, tits::TitsElt a) const;

  StandardRepK std_rep_rho_plus
    (latticetypes::Weight lambda, tits::TitsElt a) const
    {
      (lambda *= 2) += rootDatum().twoRho();
      return std_rep(lambda,a);
    }

  RawRep Levi_rep
    (latticetypes::Weight lambda, tits::TitsElt a, bitset::RankFlags gens)
    const;

/*
  The conditions below (and Normal which is not used in tests) are defined by
   Standard: $\<\lambda,\alpha\vee>\geq0$ when $\alpha$ positive imaginary
   Normal: $\<\lambda,\alpha\vee+\theta\alpha\vee>\geq0$ when $\alpha$ simple,
     complex, and orthogonal to sums of positive imaginary resp. real roots.
   Zero: $\<\lambda,\alpha\vee>=0$ for some simple-imaginary compact $\alpha$
   Final: $\<\lambda,\alpha\vee>$ odd for all simple-real roots $\alpha$

  The condition Zero is a sufficient, but possibly not necessary, condition
  for the parameter to determine a zero representation.

  The |witness| parameter is set to an index of a root that witnesses
  the failure to be Standard, non-Zero, or Final in case of such verdicts.
  This index is into |f.simpleImaginary| for |isStandard| and |isZero|, and
  it is into |f.simpleReal| for |isFinal|, where |f| is the |Fiber| at the
  canonical twisted involution for the Cartan class of |sr|.
*/
  bool isStandard(const StandardRepK& sr, size_t& witness) const;
  bool isZero(const StandardRepK& sr, size_t& witness) const;
  bool isFinal(const StandardRepK& sr, size_t& witness) const;

  latticetypes::Weight normalize(latticetypes::Weight lambda, size_t cn) const;
  void normalize(StandardRepK& sr) const
    { sr.d_lambda=project(sr.d_cartan,normalize(lift(sr),sr.d_cartan)); }

  tits::TitsElt titsElt(const StandardRepK& s) const
  {
    return tits::TitsElt(titsGroup(),
			 twistedInvolution(s.d_cartan),
			 s.d_fiberElt);
  }

  std::ostream& print(std::ostream& strm, const StandardRepK& sr) const;
  std::ostream& print(std::ostream& strm, const Char& ch) const;

  /*!
    Returns the sum of absolute values of the scalar products of lambda
    expressed in the full basis and the positive coroots. This gives a Weyl
    group invariant limit on the size of the weights that will be needed.
  */
  level height(const StandardRepK& s) const;

}; // |SRK_context|

// This class serves to store tables of previously computed mappings from
// "bad" standard representations to good ones. Also the information
// necessary to interpret the d_lambda field in StandardRepK are stored here
class KhatContext : public SRK_context
{
  const kgb::KGB& d_KGB;

  typedef hashtable::HashTable<StandardRepK,seq_no> Hash;

  StandardRepK::Pooltype nonfinal_pool,final_pool;
  Hash nonfinals,finals;

  std::vector<level> height_of; // alongside |final_pool|
  graded_compare height_graded; // ordering object that will use |height_of|

  // a set of equations rewriting to Standard, Normal, Final, NonZero elements
  SR_rewrites d_rules; // maps from |seq_no| of |nonfinals| to |combination|

  // we cache a number of |proj_info| values, indexed by sets of generators
  bitset_entry::Pooltype proj_pool;
  hashtable::HashTable<bitset_entry,unsigned long> proj_sets;

  std::vector<proj_info> proj_data;

 public:

// constructors, destructors, and swap

  KhatContext(realredgp::RealReductiveGroup &G,const kgb::KGB& kgb);

// accessors and manipulators (manipulation only as side effect for efficiency)

  const kgb::KGB& kgb() const { return d_KGB; }

  // RepK from KGB number only, with |lambda=rho|; method is currently unused
  StandardRepK KGB_elt_rep(kgb::KGBElt z) const
    {
      return std_rep(rootDatum().twoRho(),d_KGB.titsElt(z));
    }

  seq_no nr_reps() const { return final_pool.size(); }

  StandardRepK rep_no(seq_no i) const { return final_pool[i]; }

  using SRK_context::print;
  std::ostream& print(std::ostream& strm, const combination& ch,
		      bool brief=false) const;

  using SRK_context::height;
  level height(seq_no i) const
  {
    assert(i<height_of.size());
    return height_of[i]; // which will equal |height(rep_no(i))|
  }

  const graded_compare& height_order() const { return height_graded; }

  //! Lower bound for height of representation after adding positive roots
  level height_bound(const latticetypes::Weight& lambda); // non |const|


  combination standardize(const StandardRepK& sr); // non |const|: |d_rules++|
  combination standardize(StandardRepK& sr) // non |const|, normalizes |sr|
  { normalize(sr); return standardize(static_cast<const StandardRepK&>(sr)); }
  combination standardize(const Char& chi); // non |const|

  combination truncate(const combination& c, level bound) const;


  kgb::KGBEltList sub_KGB(const PSalgebra& q) const;

  // Hecht-Schmid identity for simple-imaginary root $\alpha$
  HechtSchmid HS_id(const StandardRepK& s, rootdata::RootNbr alpha) const;

  // Hecht-Schmid identity for simple-real root $\alpha$
  HechtSchmid back_HS_id(const StandardRepK& s, rootdata::RootNbr alpha) const;

  CharForm K_type_formula
    (const StandardRepK& sr, level bound=~0u);
  equation mu_equation(seq_no, level bound=~0u); // adds equations

  std::vector<equation> saturate
    (const std::set<equation>& system, level bound);

  PSalgebra theta_stable_parabolic
    (const StandardRepK& sr, weyl::WeylWord& conjugator) const;

  // saturate and invert |system| up to |bound|, writing list into |reps|
  matrix::Matrix<CharCoeff>
    K_type_matrix(std::set<equation>& system,
		  level bound,
		  std::vector<seq_no>& reps); // non |const|

  combination branch(seq_no s, level bound); // non |const|

  void go(const StandardRepK& sr);  // used by "test" command

// manipulators

// private methods

private:
  RawChar KGB_sum(const PSalgebra& q, const latticetypes::Weight& lambda)
    const;

  const proj_info& get_projection(bitset::RankFlags gens); // non |const|

}; // class KhatContext





/* HechtSchmid identities are represented as equation with a main left hand
   member |lh|, and optional second left member |lh2|, and two optional
   right hand members |rh1| and |rh2|. Standardreps are all normalized.
*/
class HechtSchmid
{
  StandardRepK lh;
  StandardRepK* lh2; // owned pointers
  StandardRepK* rh1;
  StandardRepK* rh2;

 public:
 HechtSchmid(const StandardRepK& s, const KhatContext& khc)
    : lh(s), lh2(NULL), rh1(NULL), rh2(NULL) {}
  ~HechtSchmid() { delete lh2; delete rh1; delete rh2; }

  void add_lh(const StandardRepK& s, const KhatContext& khc)
    { lh2=new StandardRepK(s); }
  void add_rh(const StandardRepK& s, const KhatContext& khc)
    { *(rh1==NULL ? &rh1 : &rh2) = new StandardRepK(s); }

  int n_lhs() const { return lh2==NULL ? 1 : 2; }
  int n_rhs() const { return rh1==NULL ? 0 : rh2==NULL ? 1 : 2; }
  Char rhs () const; // stored right hand side: |*rh1 + *rh2|
  Char equivalent () const; // expression for initial term |lh|: |rhs() - *lh1|

}; // |class HechtSchmid|

class PSalgebra // Parabolic subalgebra
{
  tits::TitsElt strong_inv; // corresponding strong involution
  size_t cn; // number of the Cartan class
  bitset::RankFlags sub_diagram; // simple roots forming basis of Levi factor
  rootdata::RootSet nilpotents; // (positive) roots in nilpotent radical
 public:
  PSalgebra (tits::TitsElt base,
	     const kgb::KGB& kgb,
	     const complexredgp::ComplexReductiveGroup& G);

  const tits::TitsElt& strong_involution() const { return strong_inv; }
  weyl::TwistedInvolution involution() const { return strong_inv.tw(); }
  size_t Cartan_no() const { return cn; }
  bitset::RankFlags Levi_gens() const { return sub_diagram; }
  const rootdata::RootSet& radical() const { return nilpotents; }
}; // |class PSalgebra|


} // |namespace standardrepk|

} // |namespace atlas|

#endif
