/*
  This is standardrepk.h

  Copyright (C) 2004, 2005 Fokko du Cloux
  Copyright (C) 2008, 2009, 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations version 0.2.4

  For license information see the LICENSE file
*/
/*
  Class definition and function declarations for the classes
  StandardRepK and KhatContext.
*/

#ifndef STANDARDREPK_H  /* guard against multiple inclusions */
#define STANDARDREPK_H

#include <set>
#include <vector>
#include <iostream>

#include "../Atlas.h"

#include "innerclass.h"// inlines
#include "realredgp.h"	// numerous inline methods
#include "hashtable.h"	// containment
#include "free_abelian.h"// containment and use via |Char|


namespace atlas {

/******** type declarations  and typedefs ************************************/

namespace standardrepk {


// type |HCParam|: image of a weight of the $\rho$-cover mod $(1-\theta)X^*$
// |typedef std::pair<Weight,RankFlags> HCParam;|: free part wrt rho, torsion

// |typedef Free_Abelian<StandardRepK> Char;| $\Z$-lin. comb. of |StandardRepK|
// |typedef Free_Abelian<StandardRepK,Polynomial<int> > q_Char;|

// the following cannot be in ../Atlas.h: they need free_abelian.h
typedef Char::coef_t CharCoeff;
typedef q_Char::coef_t q_CharCoeff; // i.e., |Polynomial<int>|

// remaining definitions could be in ../Atlas.h, but seem module-specific

typedef Free_Abelian<seq_no,long int,graded_compare> combination;
typedef std::pair<seq_no,combination> equation;

typedef Free_Abelian<seq_no,q_CharCoeff,graded_compare> q_combin;
typedef std::pair<seq_no,q_combin> q_equation;


/******** function declarations *********************************************/

template <typename C>
  matrix::Matrix_base<C> triangularize
    (const std::vector<
       std::pair<seq_no,
                 Free_Abelian<seq_no,C,graded_compare>
                > >& system,
     std::vector<seq_no>& new_order);

template <typename C>
  matrix::Matrix_base<C> inverse_lower_triangular
    (const matrix::Matrix_base<C>& U);

  q_combin to_q(const combination& c);
  combination q_is_1(const q_combin& c);

  q_Char to_q(const Char& chi);
  Char q_is_1(const q_Char& chi);


/******** type definitions **************************************************/


/*
  |StandardRepK| represents the restriction to $K$ of a (coherently) continued
  standard Harish-Chandra module.

  This is a parameter type like Tits elements; the important operations are
  modifying and comparing values, not storing additional data that facilitate
  methods. For that, auxiliary classes |SRK_context| and |KhatContext|, which
  has a role similar to |WeylGroup| with respect to |WeylElt|, will be used.

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

  In the Atlas picture, the Cartan and complete positive root system are
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
  The (left) torus part is explicitly represented as |d_fiberElt|; it
  implicitly determines the real form as in the Cartan class module,
  provided the central square class (which is stored elsewhere) is known.

  Since we are interested only in HC modules restricted to $K$, we are
  interested only in the character $\gamma$ restricted to $\tilde H(R)_c$.
  Characters of the compact group $H(R)_c$ are the same as algebraic
  characters of the complexification; that is

  $\widehat{H(R)_c}$ is identified with $X^* /(1-\theta)X^*$.

  At the level of the $\rho$-cover, we get

  ($\gamma$ restricted to $K$) is identified with an element of the
  coset-quotient $(X^* + \rho)/(1-\theta) X^*$.

  This is the information held in |d_lambda|. (The name $\lambda$ refers to
  the restriction of $\gamma$ to $\tilde H(R)_c$.) This is of type |HCParam|,
  consisting of an integer vector and an integer respesenting a bit vector;
  the integer vector gives the non-torsion part of $(X^*+\rho)/(1-\theta)X^*$
  on a basis of $(1/2)X^*$ held in |KhatContext|, and the bitvector gives the
  torsion part, via a basis also given there.
*/

class StandardRepK
{

  friend class SRK_context; // which is like |WeylGroup| is for |WeylElt|

  // Number of the Cartan to which the HC module is associated.
  CartanNbr d_cartan;

  // no real form or base grading recorded in elements; they're in |KhatContext|
  // Element of the fiber group; left torus part of the strong involution
  TorusPart d_fiberElt; // a SmallBitVector

  // Character of the rho cover of H^theta, on basis of $(1/2)X^*$
  HCParam d_lambda;

// constructors, destructors, and swap

  // main constructor is private, used by |SRK_context| methods
  // fundamental bare-bones constructor
  StandardRepK(CartanNbr cn, const TorusPart& x, const HCParam& lambda)
    : d_cartan(cn), d_fiberElt(x), d_lambda(lambda) {}

public:

 StandardRepK() : d_fiberElt(0) {} // so empty slots can be created

  void swap(const StandardRepK&);

  // accessors

  bool operator< (const StandardRepK&) const;
  bool operator== (const StandardRepK&) const;

  size_t Cartan() const { return d_cartan; }

// manipulators: none (done by friend class KhatContext)

// special members required by HashTable

  typedef std::vector<StandardRepK> Pooltype;
  bool operator!=(const StandardRepK& another) const
    { return not operator==(another); }
  size_t hashCode(size_t modulus) const;

}; // |class StandardRepK|


// Per Cartan class information for handling |StandardRepK| values
struct Cartan_info
{
  // projection matrix to torsion free part
  int_Matrix freeProjector;
  // projection matrix to torsion part, after rho-shift and reduction mod 2
  BinaryMap torsionProjector;

  // matrix used to lift free part of |HCParam| back to a weight
  int_Matrix freeLift;

  // list of vectors used to lift torsion part of |HCParam| to a weight
  WeightList torsionLift;

  // space that fiber parts are reduced modulo
  SmallSubspace fiber_modulus;

  // simple roots orthogonal to sums of positive imaginary and real roots
  // a chosen one from each pair of $\theta$-conjugate such simple roots
  RankFlags bi_ortho; // simple roots, and necessarily complex ones
  WeightList sum_coroots; // associated sums of 2 coroots

  Cartan_info() : torsionProjector(0) {}

  const Weight& coroot_sum(unsigned int i) const
  { assert(i<sum_coroots.size()); return sum_coroots[i]; }

}; // |struct Cartan_info|

// a data type used for storing projections to facets of fundamental chamber
struct proj_info
{
  int_Matrix projection;
  LatticeCoeff denom;
}; // |struct proj_info|


// a wrapper around |RankFlags| to allow a hash table indexed by them
struct bitset_entry : public RankFlags
{
  typedef RankFlags base;
  typedef std::vector<bitset_entry> Pooltype;
  bitset_entry(base b) : base(b) {}
  size_t hashCode(size_t modulus) const { return to_ulong()&(modulus-1); }
}; // |struct bitset_entry|


/* This class stores the information necessary to interpret a |StandardRepK|,
   but it does not store extensive tables concerning them, which is relegated
   to the derived class |KhatContext| defined below.

   Just one dynamic table is held, for projection matrices correponding to
   different subsets of simple roots; they serve to speed up the height
   computation. That computation is not a necessary part for the other
   functionality of this class, but it allows height-truncation to be built
   into for instance |K_type_formula|, which speeds up simple cases a lot.
 */
class SRK_context
{
  RealReductiveGroup& G;
  BitMap Cartan_set;       // marks recorded Cartan class numbers
  std::vector<Cartan_info> C_info; // indexed by number of Cartan for |GR|

// this member is precomputed to increase efficiency of certain operations
  std::vector<BinaryMap> simple_reflection_mod_2; // dual side

// we cache a number of |proj_info| values, indexed by sets of generators
  bitset_entry::Pooltype proj_pool;
  HashTable<bitset_entry,unsigned int> proj_sets;
  std::vector<proj_info> proj_data;

  static RootNbr offender; // |static| so |const| methods can tamper with it

 public:
  SRK_context(RealReductiveGroup &G);

  // accessors
  const InnerClass& innerClass() const { return G.innerClass(); }
  const RootDatum& rootDatum() const { return G.rootDatum(); }
  const WeylGroup& weylGroup() const { return G.weylGroup(); }
  const TwistedWeylGroup& twistedWeylGroup() const
    { return G.twistedWeylGroup(); }
  const TitsGroup& titsGroup() const { return G.titsGroup(); }
  const TitsCoset& basedTitsGroup() const { return G.basedTitsGroup(); }

  const TwistedInvolution involution_of_Cartan(CartanNbr cn) const
    { return innerClass().involution_of_Cartan(cn); }
  const Fiber& fiber(const StandardRepK& sr) const
    { return G.cartan(sr.Cartan()).fiber(); }

  const KGB& kgb() const { return G.kgb(); }

  const Cartan_info& info(CartanNbr cn) const
    { return C_info[Cartan_set.position(cn)]; }
  const BinaryMap& dual_reflection(weyl::Generator i) const
  { return simple_reflection_mod_2[i]; }

  // projection (keeping free & torsion parts): doubled |Weight| to |HCParam|
  HCParam project(CartanNbr cn, Weight lambda) const; // by value

  // a section of |project|
  Weight lift(CartanNbr cn, HCParam p) const;

  // $(1+\theta)$ times |lift(cn,p)|; this is independent of choice of lift
  Weight theta_lift(CartanNbr cn, HCParam p) const
  {
    Weight result=lift(cn,p);
    result += G.cartan(cn).involution()*result;
    return result;
  }

  Weight lift(const StandardRepK& s) const
  { return lift(s.d_cartan,s.d_lambda); }

  Weight theta_lift(const StandardRepK& s) const
  { return theta_lift(s.d_cartan,s.d_lambda); }

  StandardRepK std_rep (const Weight& two_lambda, TitsElt a) const;

  StandardRepK std_rep_rho_plus (Weight lambda_rho, TitsElt a) const
  {
    (lambda_rho *= 2) += rootDatum().twoRho();
    return std_rep(lambda_rho,a); // in fact first argument now is $2\lambda$
  }

  // RepK from KGB number only, with |lambda=rho|; method is currently unused
  StandardRepK KGB_elt_rep(KGBElt z) const;

  // RepK' for |(lambda,a)| relative to Levi subgroup given by |gens|; unused
  RawRep Levi_rep (Weight lambda, TitsElt a, RankFlags gens) const;

/*
  The conditions below are defined by
   Standard: $\<\lambda,\alpha\vee>\geq0$ when $\alpha$ positive imaginary
   Normal: $\<\alpha^\vee(1+\theta),\lambda>\geq0$ when $\alpha$ simple,
     complex, and orthogonal to sums of positive imaginary resp. real roots.
   Zero: $\<\lambda,\alpha\vee>=0$ for some simple-imaginary compact $\alpha$
   Final: $\<\lambda,\alpha\vee>$ odd for all simple-real roots $\alpha$

  The condition Zero is a sufficient, but possibly not necessary, condition
  for the parameter to determine a zero representation.

  On a negative result, these functions set the value of |offender| to the
  number of a root that witnesses the failure to be Standard, non-Zero, final.
  This number can be retrieved (and is then cleared) by calling |witness|.
  For |isNormal| the stored number indexes |sum_coroots|, so is not a root
  number (use |bi_ortho.n_th_bit| for that) but otherwise use is the same.
*/
  bool isStandard(const StandardRepK& sr) const;
  bool isNormal(Weight lambda, CartanNbr cn) const;
  bool isNormal(const StandardRepK& sr) const
    { return isNormal(lift(sr),sr.Cartan()); }
  bool isZero(const StandardRepK& sr) const;
  bool isFinal(const StandardRepK& sr) const;

  // after one of the above fail, |witness| tells which root caused it to
  RootNbr witness() const { auto result=offender; offender=-1; return result; }

  void normalize(StandardRepK& sr) const;

  q_Char q_normalize_eq (const StandardRepK& sr,size_t witness) const;
  q_Char q_reflect_eq (const StandardRepK& sr,size_t i,
		       Weight lambda,
		       const Weight& cowt) const;

  TitsElt titsElt(const StandardRepK& s) const
  {
    return TitsElt(titsGroup(), involution_of_Cartan(s.d_cartan), s.d_fiberElt);
  }

  // list of elements below element representing |q|, by Levi (real) generators
  KGBEltList sub_KGB(const PSalgebra& q) const;

  PSalgebra theta_stable_parabolic
    (const StandardRepK& sr, WeylWord& conjugator) const;

  // Hecht-Schmid identity for simple-imaginary root $\alpha$
  HechtSchmid HS_id(const StandardRepK& s, RootNbr alpha) const;

  // Hecht-Schmid identity for simple-real root $\alpha$
  HechtSchmid back_HS_id(const StandardRepK& s, RootNbr alpha) const;

  // equivalent by $q$-Hecht-Schmid identity for simple-imaginary root $\alpha$
  q_Char q_HS_id_eq(const StandardRepK& s, RootNbr alpha) const;

  // no need for |q_back_HS_id_eq|, it would not involve $q$; use |back_HS_id|

/*
  Half the sum of absolute values on the positive coroots of $(1+theta)\lambda$
  or equivalently that weight made dominant paired with $\rho^\vee$. This gives
  a Weyl group invariant limit on the size of the weights that will be needed.
*/
  level height(const StandardRepK& s) const;

  std::ostream& print(std::ostream& strm, const StandardRepK& sr) const;
  std::ostream& print(std::ostream& strm, const Char& ch) const;
  std::ostream& print(std::ostream& strm, const q_Char& ch) const;

  // manipulators

  // lower bound for height of representation after adding positive roots
  level height_bound(const Weight& lambda); // non |const|, see |get_projection|

  CharForm K_type_formula
    (const StandardRepK& sr, level bound=~0u); // non-|const| (|height_bound|)
  q_CharForm q_K_type_formula
    (const StandardRepK& sr, level bound=~0u); // non-|const| (|height_bound|)

// private methods
 private:
  // alternating sum over |sub_KGB(q)| of |RawRep| values deduced from |lambda|
  RawChar KGB_sum(const PSalgebra& q, const Weight& lambda) const;

  Raw_q_Char q_KGB_sum(const PSalgebra& q, const Weight& lambda) const;

  const proj_info& get_projection(RankFlags gens); // non |const|

}; // |SRK_context|

// a comparison object that first takes into account a grading (height)
class graded_compare
{
  const std::vector<level>* h; // pointer, so |graded_compare| can be assignable

public:
  graded_compare (const std::vector<level>& a) : h(&a) {}

  bool operator()(seq_no x, seq_no y) const
  {
    level hx=(*h)[x], hy=(*h)[y];
    return hx!=hy ? hx<hy : x<y; // compare levels; sequence nbrs break ties
  }
}; // |class graded_compare|

/*
  This class serves to store tables of previously computed mappings from "bad"
  standard representations to good ones. Also the information necessary to
  interpret the |d_lambda| field in |StandardRepK| are stored here
*/
class KhatContext : public SRK_context
{
  typedef HashTable<StandardRepK,seq_no> Hash;

  StandardRepK::Pooltype nonfinal_pool,final_pool;
  Hash nonfinals,finals;

  std::vector<level> height_of; // alongside |final_pool|
  graded_compare height_graded; // ordering object that will use |height_of|

  // a set of equations rewriting to Standard, Normal, Final, NonZero elements
  std::vector<combination> expanded; // equivalents for |nonfinal| elements

 public:

// constructors, destructors, and swap

  KhatContext(RealReductiveGroup &G);

// accessors and manipulators (manipulation only as side effect for efficiency)

  seq_no nr_reps() const { return final_pool.size(); }

  StandardRepK rep_no(seq_no i) const { return final_pool[i]; }


  using SRK_context::height;
  level height(seq_no i) const
  {
    assert(i<height_of.size());
    return height_of[i]; // which will equal |height(rep_no(i))|
  }

  const graded_compare& height_order() const { return height_graded; }

  seq_no match_final(const StandardRepK& sr); // look up |sr|, assumed final
  combination standardize(const StandardRepK& sr); // non |const|: |expanded++|
  combination standardize(const Char& chi); // non |const|

  combination truncate(const combination& c, level bound) const;


  equation mu_equation(seq_no, level bound=~0u); // adds equations

  std::vector<equation> saturate
    (const std::set<equation>& system, level bound);

 // saturate and invert |system| up to |bound|, writing list into |reps|
  matrix::Matrix_base<CharCoeff>
    K_type_matrix(std::set<equation>& system,
		  level bound,
		  std::vector<seq_no>& reps,
		  matrix::Matrix_base<CharCoeff>* direct_p); // non |const|

  combination branch(seq_no s, level bound); // non |const|

  void go(const StandardRepK& sr);  // used by "test" command

  using SRK_context::print;
  std::ostream& print(std::ostream& strm, const combination& ch,
		      bool brief=false) const;
// manipulators
 private:
  const combination& equate(seq_no n, const combination& rhs); // nonfinal #n

}; // |class KhatContext|

// This class serves to store tables of previously computed mappings from
// "bad" standard representations to good ones. Also the information
// necessary to interpret the d_lambda field in StandardRepK are stored here
class qKhatContext : public SRK_context
{
  typedef HashTable<StandardRepK,seq_no> Hash;

  StandardRepK::Pooltype nonfinal_pool,final_pool;
  Hash nonfinals,finals;

  std::vector<level> height_of; // alongside |final_pool|
  graded_compare height_graded; // ordering object that will use |height_of|

  // a set of equations rewriting to Standard, Normal, Final, NonZero elements
  std::vector<q_combin> expanded; // equivalents for |nonfinal| elements

 public:

// constructors, destructors, and swap

  qKhatContext(RealReductiveGroup &G);

// accessors and manipulators (manipulation only as side effect for efficiency)

  seq_no nr_reps() const { return final_pool.size(); }

  StandardRepK rep_no(seq_no i) const { return final_pool[i]; }


  using SRK_context::height;
  level height(seq_no i) const
  {
    assert(i<height_of.size());
    return height_of[i]; // which will equal |height(rep_no(i))|
  }

  const graded_compare& height_order() const { return height_graded; }

  q_combin standardize(const StandardRepK& sr); // non |const|: |expanded++|
  q_combin standardize(const q_Char& chi); // non |const|

  q_combin truncate(const q_combin& c, level bound) const;

  q_equation mu_equation(seq_no, level bound=~0u); // adds equations

  std::vector<q_equation> saturate
    (const std::set<q_equation>& system, level bound);

  matrix::Matrix_base<q_CharCoeff>
    K_type_matrix(std::set<q_equation>& system,
		  level bound,
		  std::vector<seq_no>& reps,
		  matrix::Matrix_base<q_CharCoeff>* direct_p); // non |const|

  q_combin branch(seq_no s, level bound); // non |const|

  using SRK_context::print;
  std::ostream& print(std::ostream& strm, const q_combin& ch, bool brief=false)
    const;

// manipulators
 private:
  const q_combin& equate(seq_no n, const q_combin& rhs); // nonfinal #n

}; // |class qKhatContext|


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
 HechtSchmid(const StandardRepK& s)
    : lh(s), lh2(NULL), rh1(NULL), rh2(NULL) {}
  ~HechtSchmid() { delete lh2; delete rh1; delete rh2; }

  void add_lh(const StandardRepK& s) { lh2=new StandardRepK(s); }
  void add_rh(const StandardRepK& s)
    { *(rh1==NULL ? &rh1 : &rh2) = new StandardRepK(s); }

  int n_lhs() const { return lh2==NULL ? 1 : 2; }
  int n_rhs() const { return rh1==NULL ? 0 : rh2==NULL ? 1 : 2; }
  Char rhs () const; // stored right hand side: |*rh1 + *rh2|
  Char equivalent () const; // expression for initial term |lh|: |rhs() - *lh1|

}; // |class HechtSchmid|

class PSalgebra // Parabolic subalgebra
{
  TitsElt strong_inv; // corresponding strong involution
  CartanNbr cn; // number of the Cartan class
  RankFlags sub_diagram; // simple roots forming basis of Levi factor
  RootNbrSet nilpotents; // (positive) roots in nilpotent radical
 public:
  PSalgebra (TitsElt base,
	     const InnerClass& G);

  const TitsElt& strong_involution() const { return strong_inv; }
  TwistedInvolution involution() const { return strong_inv.tw(); }
  CartanNbr Cartan_no() const { return cn; }
  RankFlags Levi_gens() const { return sub_diagram; }
  const RootNbrSet& radical() const { return nilpotents; }
}; // |class PSalgebra|


} // |namespace standardrepk|

} // |namespace atlas|

#endif
