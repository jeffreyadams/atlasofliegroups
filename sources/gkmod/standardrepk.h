/*!
\file
\brief Class definition and function declarations for the classes
StandardRepK and KHatComputations.
*/
/*
  This is standardrepk.h

  Copyright (C) 2004,2005 Fokko du Cloux
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
    class KHatComputations;
    class HechtSchmid;

// An image of a weight of the $\rho$-cover mod $(1-\theta)X^*$
    typedef std::pair
      <latticetypes::LatticeElt, // free part, after subtraction of rho
       bitset::RankFlags         // torsion part, compact representation
      > HCParam;

    // a linear combination of |StandardRepK|s
    typedef free_abelian::Free_Abelian<StandardRepK> Char;
    typedef long int CharCoeff; // this is what |Free_Abelian| uses

    // a character formula; first component stands for its lowest $K$-type
    typedef std::pair< StandardRepK,Char> CharForm;

    struct PSalgebra // Parabolic subalgebra
    {
      size_t cartan; // number of the Cartan class
      std::vector<atlas::rootdata::RootNbr> levi; // simple basis subsystem
      std::vector<std::pair<atlas::rootdata::RootNbr,atlas::rootdata::RootNbr>
		  > nilp; // basis of nilpotent radical, with their theta-image
    };
} // namespace standardrepk

/******** function declarations *********************************************/

namespace standardrepk {

  atlas::matrix::Matrix<CharCoeff> triangularize
    (const std::vector<CharForm>& system,
     std::vector<StandardRepK>& new_order);

  atlas::matrix::Matrix<CharCoeff> inverse_upper_triangular
    (const atlas::matrix::Matrix<CharCoeff>& U);

} // namespace standardrepk

/******** type definitions **************************************************/

namespace standardrepk {

  /*!
\brief Setting for studying the restrictions to K of standard
Harish-Chandra modules.
*/

class StandardRepK {

/*!
  \brief Represents the restriction to $K$ of a (coherently) continued
  standard Harish-Chandra module.

  This is a parameter type like Tits elements; the important operations are
  modifying and comparing values, not storing additional data that
  facilitate methods. For that auxiliary classes similar to |WeylGroup| with
  respect to |WeylElt| will be used; here mainly |KHatComputations|.

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
  the twisted involution for the distinguished involution in that class).
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
  the integer vector gives the non-torsion part of
  $(X^*+\rho)/(1-\theta)X^*$ on a basis of $(1/2)X^*$ held in
  |KHatComputations|, and the bitvector gives the torsion part, via a basis
  also given there.

*/


  friend class KHatComputations; // which is like |WeylGroup| is for |WeylElt|

 private:

  enum StatusFlagNames { IsStandard, IsNonZero, IsFinal, numFlags };
/*
   standard: $\<\lambda,\alpha\vee>\geq0$ when $\alpha$ imaginary
   normal: $\<\lambda,\alpha\vee+\theta\alpha\vee>\geq0$ when $\alpha$ simple,
     complex, and orthogonal to sums of positive imaginary resp. real roots.
   final: $\<\lambda,\alpha\vee>$ odd for all simple real roots $\alpha$
   zero: $\<\lambda,\alpha\vee>=0$ for some simple-imaginary compact $\alpha$

   The remaining three predicates can only be applied if standard is ensured
*/

  typedef bitset::BitSet<numFlags> Status;


/*!
  \brief Number of the Cartan to which the HC module is associated.
*/
  size_t d_cartan;

// no real form or base grading recorded in elements; in |KHatComputations|
/*!
  \brief Element of the fiber group; left torus part of the strong involution
*/
  tits::TorusPart d_fiberElt; // a SmallBitVector

/*!
  \brief Character of the rho cover of H^theta, on basis of $(1/2)X^*$
*/
  HCParam d_lambda;

/*!\brief
  Records whether this continued standard module is actually standard, and
  whether it is known to be in normal position, final or zero. Apart from
  isStandard, an unset bit does not mean the property necessarily lacks.

  This field is not used in comparisons.
*/
  Status d_status;

// constructors, destructors, and swap

  // main constructor is private, used by |KHatComputations| methods
  // fundamental bare-bones constructor; no status is set here
  StandardRepK(size_t cn, const tits::TorusPart& x, const HCParam& lambda)
    : d_cartan(cn), d_fiberElt(x), d_lambda(lambda), d_status() {}

public:

  StandardRepK() {} // so empty slots can be created

  void swap(const StandardRepK&);

  // accessors

  bool operator< (const StandardRepK&) const;
  bool operator== (const StandardRepK&) const;

  bool isFinal() const {
    return d_status[IsFinal];
  }

  bool isNonZero() const {
    return d_status[IsNonZero];
  }

  bool isStandard() const {
    return d_status[IsStandard];
  }


// manipulators: none (done by friend class KHatComputations)

// special members required by hashtable::HashTable

  typedef std::vector<StandardRepK> Pooltype;
  bool operator!=(const StandardRepK& another) const
    { return not operator==(another); }
  size_t hashCode(size_t modulus) const;

}; // class StandardRepK


// class handling rewriting of StandardRepK values by directed equations
// all left hand sides are assumed distinct, so just repeated substitutions
class SR_rewrites
{
public:
  typedef unsigned int seq_no;
  typedef free_abelian::Free_Abelian<seq_no> combination;
  typedef std::vector<combination> sys_t;

private:

  sys_t system;

public:
  SR_rewrites() : system() {}

  const combination& lookup(seq_no n) const; // an equation must exist!

  // manipulators
  void equate(seq_no n, const combination& rhs);

}; // SR_equations

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

}; // Cartan_info

// This class serves to store tables of previously computed mappings from
// "bad" standard representations to good ones. Also the information
// necessary to interpret the d_lambda field in StandardRepK are stored here
// (in d_realForm .. d_minusQuotient)

class KHatComputations
{
  const complexredgp::ComplexReductiveGroup* d_G;
  const kgb::KGB& d_KGB;

  typedef hashtable::HashTable<StandardRepK,SR_rewrites::seq_no> Hash;

  StandardRepK::Pooltype nonfinal_pool,final_pool;
  Hash nonfinals,finals;

  // a set of equations rewriting to Standard, Normal, Final, NonZero elements
  SR_rewrites d_rules;

  std::vector<latticetypes::BinaryMap> simple_reflection_mod_2;

// the remaining data members allow interpretation of |StandardRepK| objects
  atlas::realform::RealForm d_realForm;
  const kgb::EnrichedTitsGroup d_Tg;

  std::vector<Cartan_info> d_data;

 public:

// constructors, destructors, and swap

  KHatComputations(const realredgp::RealReductiveGroup &G,
		   const kgb::KGB& kgb);

  ~KHatComputations() {}

  void swap(KHatComputations&);

// accessors

  const complexredgp::ComplexReductiveGroup&
    complexReductiveGroup() const { return *d_G; }

  const rootdata::RootDatum& rootDatum() const { return d_G->rootDatum(); }

  StandardRepK std_rep
    (const latticetypes::Weight& two_lambda, tits::TitsElt a) const;

  StandardRepK std_rep_rho_plus
    (latticetypes::Weight lambda, tits::TitsElt a) const
    {
      (lambda *= 2) += rootDatum().twoRho();
      return std_rep(lambda,a);
    }

  StandardRepK KGB_elt_rep(kgb::KGBElt z) const
    {
      return std_rep(rootDatum().twoRho(),d_KGB.titsElt(z));
    }


  tits::TitsElt titsElt(StandardRepK s) const
  {
    return tits::TitsElt(d_Tg.titsGroup(),
			 d_G->twistedInvolution(s.d_cartan),
			 s.d_fiberElt);
  }

  std::ostream& print(std::ostream& strm,StandardRepK sr);

  //!\brief A section of |project|
  latticetypes::Weight lift(size_t cn, HCParam p) const;

  //!\brief (1+theta)* lifted value; this is independent of choice of lift
  latticetypes::Weight theta_lift(size_t cn, HCParam p) const
    {
      latticetypes::Weight result=lift(cn,p);
      result += d_G->cartan(cn).involution().apply(result);
      return result;
    }

  latticetypes::Weight lift(StandardRepK s) const
  { return lift(s.d_cartan,s.d_lambda); }

  // apply symmetry for root $\alpha_s$ to |a|
  void reflect(tits::TitsElt a, size_t s);

  //!\brief Projection |Weight| to |HCParam|
  HCParam project(size_t cn, latticetypes::Weight lambda) const; // by value


  void normalize(StandardRepK&) const;

  // Hecht-Schmid identity for simple-imaginary root $\alpha$
  HechtSchmid HS_id(StandardRepK s, rootdata::RootNbr alpha) const;

  // Hecht-Schmid identity for simple-real root $\alpha$
  HechtSchmid back_HS_id(StandardRepK s, rootdata::RootNbr alpha) const;


  CharForm character_formula(StandardRepK) const; // call by value

  std::vector<CharForm> saturate
    (std::set<CharForm> system, // call by value
     atlas::latticetypes::LatticeCoeff bound) const;

  PSalgebra theta_stable_parabolic
    (weyl::WeylWord& conjugator,
     const cartanset::CartanClassSet& cs,
     const size_t Cartan_nr) const;
  /*!
    Returns the sum of absolute values of the scalar products of lambda
    expressed in the full basis and the positive coroots. This gives a Weyl
    group invariant limit on the size of the weights that will be needed.
  */
  atlas::latticetypes::LatticeCoeff
    product_sumposroots (const StandardRepK& s) const;

  /*!
    When a standardrepK lambda is normalized this will return the product of
    (1+theta)*(the lifted lambda) with coroot |k|. A constituent in the above.
  */
  atlas::latticetypes::LatticeCoeff
    product_simpleroot(const StandardRepK& s, rootdata::RootNbr k) const;


  atlas::matrix::Matrix<CharCoeff>
    makeMULTmatrix(std::set<CharForm>& system,
		   const atlas::latticetypes::LatticeCoeff bound) const;

// manipulators

  SR_rewrites::combination standardize(StandardRepK sr); // call by value

  SR_rewrites::combination standardize(const Char& chi);

  void go(kgb::KGBElt, const latticetypes::Weight&);

}; // class KHatComputatons

/* HechtSchmid identities are represented as equation with a main left hand
   member |lh|, and optional second left member |lh2|, and two optional
   right hand members |rh1| and |rh2|
*/
struct HechtSchmid {

  StandardRepK lh;
  StandardRepK* lh2; // owned pointers
  StandardRepK* rh1;
  StandardRepK* rh2;

  HechtSchmid(StandardRepK s) : lh(s), lh2(NULL), rh1(NULL), rh2(NULL) {}
  ~HechtSchmid() { delete lh2; delete rh1; delete rh2; }

  void add_lh(const StandardRepK& s) { lh2=new StandardRepK(s); }
  void add_rh(const StandardRepK& s)
  {
    if (rh1==NULL) rh1=new StandardRepK(s); else rh2=new StandardRepK(s);
  }
}; // class HechtSchmid



} // namespace standardrepk

} // namespace atlas

#endif
