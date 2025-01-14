/*       Implementation of the classes TitsGroup and TitsElt.

  This module contains an implementation of a slight variant of the
  Tits group (also called extended Weyl group) as defined in J. Tits,
  J. of Algebra 4 (1966), pp. 96-116.

  The slight variant is that we include all elements of order two in the
  torus, instead of just the subgroup generated by the $m_\alpha$ (denoted
  $h_\alpha$ in Tits' paper.) Tits' original group may be defined by
  generators $\sigma_\alpha$ for $\alpha$ simple, subject to the braid
  relations and to $\sigma_\alpha^2= m_\alpha$; to get our group we just
  add a basis of elements of $H(2)$ as additional generators, and express the
  $m_\alpha$ in this basis. This makes for a simpler implementation, where
  torus parts are just elements of the $Z/2Z$-vector space $H(2)$.

  On a practical level, because the $\sigma_\alpha$ satisfy the braid
  relations, any element of the Weyl group has a canonical lifting in the Tits
  group; so we may uniquely represent any element of the Tits group as a pair
  $(t,w)$, with $t$ in $T(2)$ and $w$ in $W$ (the latter representing the
  canonical lift $\sigma_w$. The multiplication rules have to be
  thoroughly changed, though, to take into account the new relations.

  We have not tried to be optimally efficient here, as it is not
  expected that Tits computations will be significant computationally.
*/
/*
  This is tits.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include <string> // used implicitly in throwing |std::runtime_error|

#include "tits.h"

#include "lattice.h"
#include "rootdata.h"
#include "permutations.h"
#include "arithmetic.h"

#include "innerclass.h"
#include "realredgp.h"
#include "subquotient.h"
#include "subsystem.h"

#include <cassert>

namespace atlas {

namespace tits {

namespace {

std::vector<Grading> compute_square_classes
  (const InnerClass& G);

} // |namespace|

/****************************************************************************

	  Chapter I -- Classes |GlobalTitsElement|, |GlobalTitsGroup|

*****************************************************************************/



//              |GlobalTitsElement|

GlobalTitsElement GlobalTitsElement::simple_imaginary_cross
  (const RootDatum& dual_rd, // dual for pragmatic reasons
   RootNbr alpha) const // any simple-imaginary root
{
  GlobalTitsElement a(*this); // make a copy in all cases
  if (not t.negative_at(dual_rd.coroot(alpha))) // |alpha| is a noncompact root
    a.torus_part() += TorusPart(dual_rd.root(alpha)); // so add $m_\alpha$
  return a;
}


//              |GlobalTitsGroup|


GlobalTitsGroup::GlobalTitsGroup(const InnerClass& G)
  : TwistedWeylGroup(G.twistedWeylGroup())
  , prd(G.root_datum(),tags::DualTag()) // viewed from the dual side
  , delta_tr(G.distinguished().transposed())
  , alpha_v()
  , half_rho_v(G.root_datum().dual_twoRho(),4)
  , square_class_gen(compute_square_classes(G))
{
  alpha_v.reserve(G.semisimple_rank());
  for (size_t i=0; i<G.semisimple_rank(); ++i) // reduce vectors mod 2
    alpha_v.push_back(TorusPart(prd.simple_root(i)));
}

WeightInvolution
  GlobalTitsGroup::involution_matrix(const WeylElt& tw) const
{
  WeightInvolution M = delta_tr;
  M.negate(); // we need the involution |-^delta| corresponding to |delta|
  Weyl_group().act(prd,tw,M); // twisted involution, so left action is OK
  return M;
}

TorusElement GlobalTitsGroup::twisted(const TorusElement& x) const
{
  RatWeight rw = x.log_pi(false);
  return y_values::exp_pi(RatWeight(delta_tr*rw.numerator(),
				    rw.denominator()));
}

 // act by $-w_0\delta$; the three parts are mutually commuting involutions
TorusElement GlobalTitsGroup::dual_twisted(const TorusElement& x) const
{
  RatWeight rw = x.log_pi(false);
  const WeylGroup& W = Weyl_group();
  W.act(prd,W.longest(),rw);
  return y_values::exp_pi(RatWeight(delta_tr*-rw.numerator(),
				    rw.denominator()));
}

TorusElement GlobalTitsGroup::theta_tr_times_torus(const GlobalTitsElement& a)
  const
{ RatWeight rw = a.torus_part().log_pi(false);
  RatWeight delta_rw(delta_tr*rw.numerator(),rw.denominator());
  Weyl_group().act(prd,a.tw(),delta_rw);
  return y_values::exp_pi(delta_rw);
}

// assuming |a.tw()| is a twisted involution, shifted square is easy to compute
TorusElement GlobalTitsGroup::square_shifted(const GlobalTitsElement& a) const
{ return a.torus_part()+theta_tr_times_torus(a); }

bool GlobalTitsGroup::is_valid(const GlobalTitsElement& a,bool check_tw) const
{
  if (check_tw and TwistedWeylGroup::prod(twisted(a.tw()),a.tw())!=WeylElt())
    return false;
  return is_central(prd.simple_coroots_mat(),square_shifted(a));
}

 // weaker condition: square being central in subgroup
bool GlobalTitsGroup::is_valid(const GlobalTitsElement& a,
			       const SubSystem& sub) const
{
  const auto sr=sub.rank();
  LatticeMatrix alpha(sub.parent_datum().rank(),sr);
  for (weyl::Generator s=0; s<sr; ++s)
    alpha.set_column(s,sub.parent_datum().coroot(sub.parent_nr_simple(s)));
  return is_central(alpha,square_shifted(a));
}


// find simple roots giving length-decreasing links (complex descent and real)
// this could have been a method of |WeylGroup| at |TwistedInvolution|s.
RankFlags GlobalTitsGroup::descents(const GlobalTitsElement& a) const
{
  RankFlags result;
  for (weyl::Generator s=0; s<semisimple_rank(); ++s)
    result.set(s,has_descent(s,a.tw())); // this covers all cases precisely!
  return result;
}


void
GlobalTitsGroup::imaginary_cross_act(weyl::Generator s,TorusElement& t) const
{
  RatNum r = t.evaluate_at(prd.simple_coroot(s)) - RatNum(1); // $\rho_{im}$ shift
  if (r.numerator()!=0) // compact imaginary case needs no action
    add(RatWeight // now reflect for |s|, shifted to fix $\rho_{im}$
	(prd.simple_root(s)*-r.numerator(),2*r.denominator()),t);
}

/* The code below uses Tits element representation as $t*\sigma_w*\delta_1$
   and $\delta_1*\sigma_\alpha^{-1} = \sigma_\beta*\delta_1$, simple $\alpha$
   so conjugation by $\sigma_\alpha$ gives $\sigma_\alpha*(t,w)*\sigma_\beta$
*/
int // length change in $\{-1,0,+1\}$, half that of change for |x.tw()|
GlobalTitsGroup::cross_act(weyl::Generator s, GlobalTitsElement& x) const
{
  int d=twistedConjugate(x.w,s); // Tits group must add $m_\alpha$ iff $d=0$
  if (d!=0) // complex cross action: reflect torus part by |s|
  {
    complex_cross_act(s,x.t);
    return d/2; // report half of length change in Weyl group
  }

  if (not has_descent(s,x.w)) // imaginary cross action
    imaginary_cross_act(s,x.t);
  return 0; // no length change this case
}

int GlobalTitsGroup::cross_act(const WeylWord& w, GlobalTitsElement& a)
  const
{
  int d=0;
  for (size_t i=w.size(); i-->0; )
    d += cross_act(w[i],a);
  return d; // total length change, in case the caller is interested
}

int GlobalTitsGroup::cross_act(GlobalTitsElement& a,const  WeylWord& w)
   const
{
  int d=0;
  for (size_t i=0; i<w.size(); ++i )
    d += cross_act(w[i],a);
  return d; // total length change, in case the caller is interested
}

/*
  When trying to invert a Cayley transform, we may need to change the torus
  element so that its evaluation on $\alpha$ becomes integer (half-integer
  means we have a potentially valid Tits element, but the requirement that the
  simple root |alpha| become a noncompact root means the value should in fact
  be integer). We may modify by a rational multiple of $\alpha^\vee$, since
  being $-\theta$-fixed after the Cayley transform such a coweight has zero
  evaluation on $\theta$-fixed weights (so the evaluation on imaginary roots
  stays half-integral), and modulo coweights that are $-\theta$-fixed before
  the Cayley transform (which do not affect $\alpha$, and which we continue to
  not care about) multiples of $\alpha^\vee$ are the only freedom we have to
  modify the evaluation at $\alpha$. If not already integral, we simply make
  the evaluation at $\alpha$ zero, which is as good as any integer. Note that
  in the reduction modulo 2 that will be used in |TitsGroup| below, this kind
  of correction is no longer possible, and we must perform a more tedious
  search for a valid correction.

  This is Lemma 14.6 in the Algorithms paper.
 */
void
GlobalTitsGroup::do_inverse_Cayley(weyl::Generator s,TorusElement& t) const
{
  const Coweight& eval_pt=prd.simple_coroot(s);
  RatWeight w=t.log_2pi();
  int num = eval_pt.dot(w.numerator()); // should become multiple of denominator

  if (num % w.denominator()!=0) // correction needed if alpha compact
  {
    const Weight& shift_vec= prd.simple_root(s);
    // |eval_pt.dot(shift_vec)==2|, so correct numerator by |(num/2d)*shift_vec|
    w -= RatWeight(shift_vec*num,2*w.denominator());
    assert(eval_pt.dot(w.numerator())==0);
    t=y_values::exp_2pi(w);
  }
}

void GlobalTitsGroup::do_inverse_Cayley(weyl::Generator s,GlobalTitsElement& a)
  const
{
  do_inverse_Cayley(s,a.t);
  left_multiply(a.w,s);
}

// Sometimes we need to compute the grading at non-simple imaginary roots.
// This could be computed using expression in simple-imaginary roots |alpha|,
// for which grading is \emph{compact} iff |torus_part().negative_at(alpha)|
// however that is not easy to implement; conjugating to simple is easier.
// Root system |rs| necessary to interpret |alpha|, dual makes no difference
bool GlobalTitsGroup::compact(const RootSystem& rs,
			      RootNbr alpha,
			      GlobalTitsElement a)
  const
{
  make_positive(rs,alpha);
  weyl::Generator s; // declare outside loop, allow inspection of final value
  while (alpha!=rs.simpleRootNbr(s=rs.find_descent(alpha)))
  {
    cross_act(s,a);
    rs.simple_reflect_root(s,alpha);
  }

  return compact(s,a); // since now |alpha==rs.simpleRootNbr(s)|
}

namespace {

/*
 The |square_class_gen| field of a |GlobalTitsGroup| contains a list of
 generators of gradings of the simple roots that generate the square classes.
 Each generator is a grading that flags a single $\delta$-fixed simple root.
 A $\Z/2\Z$-linear combination of these generators will later be translated
 into a central |TorusElement| value by taking the corresponding combination
 $c$ of fundamental coweights, and applying the map $c\mapsto\exp(2\pi i c)$.
 As representative torus element for the class one can take $t=\exp(\pi i c)$.
 The function |compute_square_classes| computes this list of generators.

 Actually the valid gradings do not form a $\Z/2\Z$ vector space, but an
 affine space, for which one can take the grading of $\delta_1$ as base point
 (it grades all simple-imaginary roots as non-compact). These gradings apply
 to $\delta$-fixed (imaginary) roots, and by extension to the set
 $(\<\Phi>_\Z)^\delta$ of $\delta$-fixed elements of the root lattice. The
 grading for the strong involution representative $t.\delta$, with torus part
 $t=\exp(\pi i c)$, differs at $\delta$-fixed element $v$ from the grading for
 $\delta$ by a factor $(-1)^{\<v,c>}$ (the exponent will be integer).
 Projecting $c$ to the $\delta$-fixed subspace of $X_*\tensor\Q$ gives a
 different representative of the same strong involution, for which the strong
 involution condition $t=\exp(\pi i(1+delta)c)\in Z(G)$ says that $\<.c>$
 takes integer values on all roots. Then $\<(1+\delta)x,c>\in2\Z$ for all
 $x\in\<\Phi>_\Z$; this means all obtained gradings are invariant under
 translation by any element of the $1+\delta$-image of the root lattice.
 Because all $v\in(\<\Phi>_\Z)^\delta$ are sums of $\delta$-fixed simple roots
 on one hand and of sums of $\delta$-exchanged simple root pairs on the other,
 it suffices to know the grading on $\delta$-fixed \emph{simple} roots.

 The list is formed of generators of the quotient of all possible gradings of
 those simple roots by the subgroup of gradings coming from $\delta$-fixed
 elements of $X_*$ (since gradings in the same coset for that subgroup define
 the same square class). We get a basis of the subgroup in |Vplus|, and then
 get gradings by taking scalar products with $\delta$-fixed simple roots.

 The generators can all be taken to be canonical basis elements (gradings with
 exactly one bit set) at $\delta$-fixed simple roots. Once the subgroup to
 quotient by is constructed, we can simply take those canonical basis vectors
 not in the "support" of the subgroup (so that reduction would be trivial).
 */
std::vector<Grading> compute_square_classes
  (const InnerClass& G)
{
  const RootDatum& rd = G.root_datum();
  const WeightInvolution& delta = G.distinguished();
  const weyl::Twist& twist = G.twistedWeylGroup().twist();

  size_t r= rd.rank();
  assert(delta.n_rows()==r);
  assert(delta.n_columns()==r);

  int_Matrix roots(0,r); // rows: coordinates of $\delta$-fixed simple roots
  RankFlags fixed; // the set of $\delta$-fixed simple roots
  for (size_t i=0; i<rd.semisimple_rank(); ++i)
    if (twist[i]==i)
    {
      fixed.set(i);
      roots.add_row(rd.simpleRoot(i));
    }

  BinaryMap A(lattice::eigen_lattice(delta.transposed(),1)); // $(X_*)^\delta$
  SmallSubspace Vplus(A); // mod 2: modding subgroup, in $X_*$ coordinates
  BinaryMap to_grading(roots); // evaluation at $\delta$-fixed simple roots
  Vplus.apply(to_grading); // converts |Vplus| to grading coordinates

  RankFlags supp = Vplus.support(); // pivot positions
  supp.complement(roots.n_rows()); // non-pivot positions among fixed ones
  supp.unslice(fixed);  // bring bits back to the simple-root positions

  std::vector<Grading> result; result.reserve(supp.count());

  for (RankFlags::iterator it=supp.begin(); it(); ++it)
  {
    result.push_back(Grading());
    result.back().set(*it);
  }

  return result;
}

} // |namespace|



/****************************************************************************

        Chapter II -- The TitsGroup class

*****************************************************************************/

/*!
  Constructs the Tits group corresponding to the root datum |rd|, and
  the fundamental involution |d| (which also defines the twist).
*/
TitsGroup::TitsGroup(const RootDatum& rd,
		     const WeylGroup& W,
		     const WeightInvolution& delta)
  : TwistedWeylGroup(W,weyl::make_twist(rd,delta))
  , d_rank(rd.rank())
  , d_simpleRoot()   // set number of vectors, but not yet
  , d_simpleCoroot() // their size (which will be |d_rank|)
  , d_involution(delta.transposed())
  , dual_involution(rd.rank()) // set below
{
  d_simpleRoot.reserve(rd.semisimple_rank());
  d_simpleCoroot.reserve(rd.semisimple_rank());
  for (size_t i = 0; i<rd.semisimple_rank(); ++i) // reduce vectors mod 2
  {
    d_simpleRoot.push_back(TorusPart(rd.simpleRoot(i)));
    d_simpleCoroot.push_back(TorusPart(rd.simpleCoroot(i)));
  }

  // we could compute $w_0*d_invulution$ to give |dual_involution|
  // but computing the integer matrix and then reducing modulo 2 is easier
  int_Matrix M = delta;
  W.act(rd,W.longest(),M);
  // M.negate(); // morally it is $-w_0$, but we will reduce modulo 2
  dual_involution = BinaryMap(M.transposed()); // (-w_0*delta).transposed()
}


/*!
  Constructs the Tits group corresponding to the adjoint group of type given
  by |Cartan_matrix| and distinguished involution given by |twist|
*/
TitsGroup::TitsGroup(const int_Matrix& Cartan_matrix,
		     const WeylGroup& W,
		     const weyl::Twist& twist)
  : TwistedWeylGroup(W,twist)
  , d_rank(Cartan_matrix.n_rows())
  , d_simpleRoot()
  , d_simpleCoroot()
  , d_involution(d_rank)   // square matrix, initially zero
  , dual_involution(d_rank)   // square matrix, initially zero
{
  d_simpleRoot.reserve(d_rank); d_simpleCoroot.reserve(d_rank);
  for (size_t i=0; i<d_rank; ++i)
  {
    d_simpleRoot.push_back(TorusPart(d_rank,i)); // $e_i$
    d_simpleCoroot.push_back(TorusPart
			     (Cartan_matrix.column(i))); // adjoint coroot
    d_involution.set(i,twist[i]); // (transpose of) |twist| permutation matrix
    dual_involution.set(i,W.Chevalley_dual(twist[i]));
  }
}

// Switching between left and right torus parts is a fundamental tool.
/*
  Find torus part $x'$ so that $x.w=w.x'$

  Algorithm: for simple reflections this is |reflect|; repeat left-to-right

  Note: a more sophisticated algorithm would involve precomputing the
  conjugation matrices for the various pieces of w. Don't do it unless it
  turns out to really matter.
*/
TorusPart TitsGroup::push_across(TorusPart x, const WeylElt& w) const
{
  WeylWord ww=Weyl_group().word(w);

  for (size_t i = 0; i < ww.size(); ++i)
    reflect(x,ww[i]);

  return x;
}

// find torus part $x'$ so that $w.x=x'.w$; inverse to |push_across|
TorusPart TitsGroup::pull_across(const WeylElt& w, TorusPart x) const
{
  WeylWord ww=Weyl_group().word(w);
  for (size_t i=ww.size(); i-->0; )
    reflect(x,ww[i]);
  return x;
}

/*
  Left multiply |a| by the canonical generator $\sigma_s$.

  This is the basic case defining the group structure in the Tits group (since
  left-multiplication by an element of $T(2)$ just adds to the torus part).

  Algorithm: This is surprisingly simple: multiplying by $\sigma_s$ just
  amounts to reflecting the torus part through |s|, then left-multiplying the
  Weyl group part $w$ by |s| in the Weyl group, and if the length goes down in
  the latter step, add a factor of $(\sigma_s)^2=m_s$ to the torus part.
  (To see this, write in the length-decreasing case $w=s.w'$ with $w'$
  reduced; then $\sigma_s\sigma_w=m_s\sigma_{w'}$, so we need to add $m_s$
  to the reflected left torus part in that case.

  The upshot is a multiplication algorithm almost as fast as in the Weyl group!
*/
void TitsGroup::sigma_mult(weyl::Generator s,TitsElt& a) const
{
  reflect(a.d_t,s); // commute (left) torus part with $\sigma_s$
  if (Weyl_group().left_multiply(a.d_w,s)<0) // on length decrease
    left_add(d_simpleCoroot[s],a);     // adjust torus part
}

void TitsGroup::sigma_inv_mult(weyl::Generator s,TitsElt& a) const
{
  reflect(a.d_t,s); // commute (left) torus part with $\sigma_s$
  if (Weyl_group().left_multiply(a.d_w,s)>0) // on length increase
    left_add(d_simpleCoroot[s],a);     // adjust torus part
}


/*
  Right multiply |a| by the canonical generator $sigma_s$.

  Algorithm: similar to above, but omitting the |reflect| (since the torus
  part is at the left), and using |weyl::mult| instead of |weyl::left_multiply|,
  and |right_add| to add the possible contribution $m_s$ instead of
  |left_add|; the latter implicitly involves a call to |pull_across|.
*/
void TitsGroup::mult_sigma(TitsElt& a, weyl::Generator s) const
{
// |WeylGroup::mult| multiplies $w$ by $s$, returns sign of the length change
  if (Weyl_group().mult(a.d_w,s)<0)  // on length decrease adjust torus part
    right_add(a,d_simpleCoroot[s]);
}

void TitsGroup::mult_sigma_inv(TitsElt& a, weyl::Generator s) const
{
  if (Weyl_group().mult(a.d_w,s)>0) // on length increase adjust torus part
    right_add(a,d_simpleCoroot[s]);
}


/*
  Product of general Tits group elements

  Algorithm: The algorithm is to multiply the second factors successively by
  the generators in a reduced expression for |a.w()|, then left-add |a.t()|.

  Since we have torus parts on the left, we prefer left-multiplication here.
*/
TitsElt TitsGroup::prod(const TitsElt& a, TitsElt b) const
{
  WeylWord ww=Weyl_group().word(a.w());

  // first incorporate the Weyl group part
  for (size_t i = ww.size(); i-->0; )
    sigma_mult(ww[i],b);

  // and finally the torus part
  left_add(a.d_t,b);
  return b;
}

// here the we basically do Weyl group action at torus side, including twist
// the (binary) matrix produced represents $delta.ww^{-1}$ which, being an
// involution (and delta too), could also have been computed as $ww.delta$
BinaryMap
TitsGroup::involutionMatrix(const WeylWord& ww) const
{
  BinaryMap result(rank(),0);
  for (size_t i=0; i<rank(); ++i)
  {
    SmallBitVector v(rank(),i); // basis vector $e_i$
    // act by letters of |ww| from left to right (because transposed action)
    for (size_t j=0; j<ww.size(); ++j)
      reflect(v,ww[j]);
    result.addColumn(twisted(v)); // finally (since we are on dual side) twist
  }
  return result;
}






/*
 *
 *				TitsCoset
 *
 */

TitsCoset::TitsCoset(const InnerClass& G, Grading base_grading)
  : my_Tits_group(nullptr) // no ownership in this case
  , Tg(G.Tits_group())
  , grading_offset(base_grading)
  , rs(G.root_datum())
{
}

// Based Tits group for the adjoint group
TitsCoset::TitsCoset(const InnerClass& G)
  : my_Tits_group(new TitsGroup(G.root_datum().Cartan_matrix(),
				      G.Weyl_group(),
				      G.twistedWeylGroup().twist()))
  , Tg(*my_Tits_group)
  , grading_offset()
  , rs(G.root_datum())
{ // make all imaginary simple roots noncompact
  for (unsigned i=0; i<G.semisimple_rank(); ++i)
    grading_offset.set(i,Tg.twisted(i)==i);
}

// Based Tits group for the adjoint dual group
TitsCoset::TitsCoset(const InnerClass& G,tags::DualTag)
  : my_Tits_group(new TitsGroup(G.root_datum().Cartan_matrix().transposed(),
				      G.Weyl_group(),
				      G.twistedWeylGroup().dual_twist()))
  , Tg(*my_Tits_group)
  , grading_offset()
  , rs(G.dual_root_datum())
{ // make all imaginary simple roots noncompact
  for (unsigned i=0; i<G.semisimple_rank(); ++i)
    grading_offset.set(i,Tg.twisted(i)==i);
}

void TitsCoset::basedTwistedConjugate
  (TitsElt& a, const WeylWord& w) const
{
  for (size_t i=0; i<w.size(); ++i)
    basedTwistedConjugate(a,w[i]);
}

void TitsCoset::basedTwistedConjugate
  (const WeylWord& w, TitsElt& a) const
{
  for (size_t i=w.size(); i-->0; )
    basedTwistedConjugate(a,w[i]);
}


/* Inverse Cayley transform is more delicate than Cayley transform, since the
   torus part has been modularly reduced when passing from an the involution
   $\theta$ to $\theta_\alpha$. While before the Cayley transform the value
   ($\pm1$) of each $\theta$-stable weight on the torus part was well defined,
   the modular reduction retains only the values for $\theta_\alpha$-stable
   weights. When going back, the values of all $\theta$-stable weights must be
   defined again, but some might have changed by the reduction. It turns out
   that our reconstruction of the torus part is valid if and only if the
   evaluation of the simple root $\alpha=\alpha_i$ through which we transform
   (and which is imaginary before Cayley transform) is $-1$, so that $\alpha$
   becomes noncompact. Indeed, if the transformation was type I, then the
   value at this root is already determined by the values at weights fixed by
   the new involution $\theta_\alpha$: the modular reduction was by
   $m_\alpha$, for which $\alpha$ evaluates to $+1$, and which forms the
   difference between the torus parts of the two elements which the same value
   of the Cayley transform, from which we have to make a choice anyway. If the
   transformation was type II, then the sublattice of $\theta$-stable weights
   is the direct sum of the $\theta_\alpha$-stable weights and the multiples
   of $\alpha$, and lifting the torus part means determining the evaluation of
   $\alpha$ at it; since there are always two possible lifts, one of them
   makes $\alpha$ noncompact.
 */
void TitsCoset::inverse_Cayley_transform
  (TitsElt& a, size_t i,
   const SmallSubspace& mod_space) const
{
  Tg.sigma_inv_mult(i,a); // set |a| to $\sigma_i^{-1}.a$
  if (not simple_grading(a,i)) // $\alpha_i$ should not become compact!
  {
    // we must find a vector in |mod_space| affecting grading at $\alpha_i$
    for (size_t j=0; j<mod_space.dimension(); ++j) // a basis vector will do
      if (Tg.dual_m_alpha(i).dot(mod_space.basis(j)))
      { // scalar product true means grading is affected: we found it
	Tg.left_add(mod_space.basis(j),a); // adjust |a|'s torus part
	break;
      }

    assert(simple_grading(a,i)); // the problem must now be corrected
  }
}


TitsElt TitsCoset::twisted(const TitsElt& a) const
{
  TitsElt result(Tg,Tg.twisted(Tg.left_torus_part(a)));
  WeylWord ww=Tg.word(a.w());
  for (size_t i=0; i<ww.size(); ++i)
  {
    weyl::Generator s=Tg.twisted(ww[i]);
    Tg.mult_sigma(result,s);
    if (grading_offset[s]) // $s$ noncompact imaginary: add $m_\alpha$
      Tg.right_add(result,Tg.m_alpha(s));
  }
  return result;
}

bool TitsCoset::grading(TitsElt a, RootNbr alpha) const
{
  make_positive(rs,alpha);
  weyl::Generator s; // declare outside loop, allow inspection of final value
  while (alpha!=rs.simpleRootNbr(s=rs.find_descent(alpha)))
  {
    basedTwistedConjugate(a,s);
    rs.simple_reflect_root(s,alpha);
  }

  return simple_grading(a,s);
}

/*
   A somewhat easier special case of previous: |alpha| simple-imaginary.
   However this is mostly a proof-of-possibility implementation: in this
   special case we can do with just the TorusPart |t| of a TitsElement |a|.

   If the base grading is 1 on all simple roots, then the result can be
   obtained by just evaluating the root |alpha| at |t|, and adding 1: this is
   the property that at the base point all simple-imaginary roots are
   noncompact. It can be explained by the fact that when repeatedly reflecting
   |alpha| to become a simple root, and transforming |a| correspongingly as in
   |grading|, each |basedTwistedConjugate| in the latter operation will be by
   a complex root (if it were real it would be orthogonal to the imaginary
   root |alpha|, and if it were imaginary this would contradict the
   simple-imaginary status of |alpha|), and this just reflects |t|; in the end
   the evaluation of the simple root obtained on the transformed |t| is the
   same as the evaluation of |alpha| had at |t|. The base grading at the
   simple root that was reached gives the remaining term 1, independently of
   which simple root that is. If the base grading is different from all-one,
   we need to compensate with a contribution for its complement. One can think
   of the complement as a sum of fundamental coweights for its values 1, which
   can be added to the torus part (since left- and right-multiplication by
   torus parts are equivalent), and we can evaluate |alpha| on this sum to get
   the compensating term. This evaluation can be obtained by expressing
   |alpha| as a sum of simple roots, reducing the coefficients modulo 2, and
   taking the scalar product with the complement of the base grading. For the
   evaluation of |alpha| at |t|, we can again use the bits of the mentioned
   reduction modulo 2, but now we must multiply by the dual $m_\alpha$'s
   (simple roots modulo 2) before taking the scalar product with |t|.
*/

bool TitsCoset::simple_imaginary_grading(TorusPart t, RootNbr alpha) const
{
  assert(rs.is_posroot(alpha));

  RankFlags re_mod2(rs.root_expr(alpha));

  bool evaluation = not // incorporate the term 1 for base grading $\delta_1$
    re_mod2.dot(base_grading().complement(constants::RANK_MAX)); // and rest

  for (RankFlags::iterator it=re_mod2.begin(); it(); ++it)
    evaluation ^= Tg.dual_m_alpha(*it).dot(t);

  return evaluation;
}

// the  difference with |basedTwistedConjugate| is that we |right_add| here
void TitsCoset::strict_based_twisted_conjugate(TitsElt& a, size_t s) const
{
  Tg.twistedConjugate(a,s);
  if (grading_offset[s])
    Tg.right_add(a,Tg.m_alpha(Tg.twisted(s)));
}


bool TitsCoset::is_valid(TitsElt a) const
{
  static TitsElt e(Tg);  // identity
  Tg.mult(a,twisted(a));
  return a==e;
}


EnrichedTitsGroup::EnrichedTitsGroup(const RealReductiveGroup& GR)
  : TitsCoset(GR.innerClass(),
	      grading_of_simples(GR.innerClass(),GR.g_rho_check()))
  , srf(GR.innerClass().sample_strong_form(GR.realForm()))
{}


/*
  The purpose of the method |backtrack_seed| is to find some |TitsElt| in the
  strong real form associated with our |EnrichedTitsGroup|, at (the canonical
  involution of) the Cartan class |cn|. We simulate the KGB construction back
  from the fundamental fiber to the one for which we try to find a seed,
  trying all the representatives in the fundamental fiber of the strong real
  form, until finding one that, along the chosen path of cross actions and
  Cayley transforms, proves to be suited for every necessary Cayley transform
  (meaning that it makes all the simple roots involved noncompact). Rather
  than set out for the journey with the full fundamental fiber and have the
  unfit ones die along the way, we select a fit element from the fundamental
  fiber and transform that one; but this does require some preliminary work.
*/
TitsElt EnrichedTitsGroup::backtrack_seed
  (const InnerClass& G, RealFormNbr rf, size_t cn) const
{
  const TitsGroup& Tgr= Tits_group();
  // a name chosen to avoid warnings about shadowing (the inaccessible) |Tg|

  const TwistedInvolution& tw=G.involution_of_Cartan(cn);

  RootNbrSet rset;
  WeylWord cross;
  innerclass::Cayley_and_cross_part(rset,cross,tw,G.root_datum(),Tgr);

  /* at this point we can get from the fundamental fiber to |tw| by first
     applying cross actions according to |cross|, and then applying Cayley
     transforms in the strongly orthogonal set |rset|.
  */

  // transform strong orthogonal set |Cayley| back to distinguished involution
  RootNbrList Cayley(rset.begin(),rset.end()); // convert to |RootNbrList|

  for (auto it=Cayley.begin(); it!=Cayley.end(); ++it)
    *it = G.root_datum().permuted_root(cross,*it);

  /* at this point we can get from the fundamental fiber to |tw| by first
     applying Cayley transforms in the strongly orthogonal set |Cayley|, and
     then applying cross actions according to |cross|
  */

/* Now find an element in the stored strong real form |srf| at the fundamental
   fiber that has noncompact grading on all the roots of |Cayley| (which are
   all imaginary for $\delta$). Methods |square| and |f_orbit| use this |srf|.
*/
  TitsElt result(Tgr);

  const Partition& srp = G.fundamental_fiber_partition(square());
  for (unsigned long x=0; x<srp.size(); ++x)
    if (srp.class_of(x)==f_orbit()) // test membership of strong real form
    {
      TorusPart t = G.lift_from_fundamental_fiber(x);
      for (size_t i=0; i<Cayley.size(); ++i)
	if (is_compact(t,Cayley[i]))
	  goto again; // none of the |Cayley[i]| should be compact

      // if we get here, |t| is OK as torus part
      result = TitsElt(Tits_group(),t); // pure torus part
      goto found;
    again: {}
    }
  assert(false); // getting here means none of the orbit elements is appropriate

found:

  /* Now we must apply the Cayley transforms and cross actions to |result|.
     However, Cayley transforms by non-simple roots are not implemented, and
     so we reorder the operations as in |Tgr.involution_expr(tw)|, which gives
     the same cross actions, but interspersed with simple Cayley transforms.
   */

  // transform |result| via Cayley transforms and cross actions
  std::vector<signed char> dec=Tgr.involution_expr(tw);
  for (size_t j=dec.size(); j-->0; )
    if (dec[j]>=0)
    {
      assert(simple_grading(result,dec[j])); // must be noncompact
      Cayley_transform(result,dec[j]);
    }
    else
      basedTwistedConjugate(result,~dec[j]);

  assert(result.tw()==tw);

  return result;  // result should be reduced immediately by caller
} // |EnrichedTitsGroup::backtrack_seed|


//				Functions

// torus parts: modulo the mod-2 reduction of the $-\theta$-fixed sublattice
SmallSubspace fiber_denom(const WeightInvolution& theta)
{
  BinaryMap A(lattice::eigen_lattice(theta.transposed(),-1));
  return SmallSubspace(A);
}

Grading compact_simples
  (const RootDatum& rd, const TorusElement& t, RankFlags imag)
{
  Grading result;
  for (auto it=imag.begin(); it(); ++it)
    result.set(*it,t.negative_at(rd.simpleRoot(*it)));
  return result;
}

Grading compact_simples(const TitsCoset& Tc, const TitsElt& a, RankFlags imag)
{
  Grading result;
  for (auto it=imag.begin(); it(); ++it)
    result.set(*it,not Tc.simple_grading(a,*it));
  return result;
}

} // |namespace tits|


} // |namespace atlas|
