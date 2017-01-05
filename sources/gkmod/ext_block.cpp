/*
  This is ext_block.cpp

  Copyright (C) 2013-2016 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_block.h"

#include <cassert>
#include <vector>

#include "bitmap.h"
#include "polynomials.h"
#include "matreduc.h"
#include "sl_list.h"

#include "innerclass.h"
#include "weyl.h"
#include "kgb.h"
#include "blocks.h"
#include "repr.h"
#include "prettyprint.h" 

/*
  For an extended group, the block structure is more complicated than an
  ordinary block, because each link in fact represents a local part of the
  parent block structure that links a twist-fixed block element to another
  twist-fixed element, via intermediate elements that are non twist-fixed.
 */

namespace atlas {

namespace ext_block {

bool is_complex(DescValue v)
{
  static unsigned long mask =
    1ul << one_complex_ascent   | 1ul << one_complex_descent |
    1ul << two_complex_ascent   | 1ul << two_complex_descent |
    1ul << three_complex_ascent | 1ul << three_complex_descent;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool is_unique_image (DescValue v)
{
  static unsigned long mask =
      1ul << one_real_pair_fixed    | 1ul << one_imaginary_pair_fixed
    | 1ul << two_semi_imaginary     | 1ul << two_semi_real
    | 1ul << two_real_double_double | 1ul << two_imaginary_double_double
    | 1ul << three_semi_imaginary   | 1ul << three_real_semi
    | 1ul << three_imaginary_semi   | 1ul << three_semi_real;

  return (1ul << v & mask) != 0 // whether |v| is one of the above
    or is_complex(v);  // these are also unique images
}

bool has_double_image(DescValue v)
{
  static unsigned long mask =
      1ul << one_real_pair_fixed         | 1ul << one_imaginary_pair_fixed
    | 1ul << two_real_double_double      | 1ul << two_imaginary_double_double
    | 1ul << two_imaginary_single_double_fixed
    | 1ul << two_real_single_double_fixed;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}


bool is_like_nonparity(DescValue v) // ascents with 0 neighbours
{
  static unsigned long mask =
      1ul << one_real_nonparity | 1ul << one_imaginary_pair_switched
    | 1ul << two_real_nonparity | 1ul << two_imaginary_single_double_switched
    | 1ul << three_real_nonparity;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool is_like_compact(DescValue v) // descents with 0 neighbours
{
  static unsigned long mask =
      1ul << one_imaginary_compact | 1ul << one_real_pair_switched
    | 1ul << two_imaginary_compact | 1ul << two_real_single_double_switched
    | 1ul << three_imaginary_compact;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool is_like_type_1(DescValue v)
{
  static unsigned long mask =
      1ul << one_imaginary_single | 1ul << one_real_pair_fixed
    | 1ul << two_imaginary_single_single  | 1ul << two_real_double_double;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool is_like_type_2(DescValue v)
{
  static unsigned long mask =
      1ul << one_imaginary_pair_fixed | 1ul << one_real_single
    | 1ul << two_imaginary_double_double  | 1ul << two_real_single_single;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool has_defect(DescValue v)
{
  static unsigned long mask =
      1ul << two_semi_imaginary   | 1ul << two_semi_real
    | 1ul << three_semi_imaginary | 1ul << three_real_semi
    | 1ul << three_imaginary_semi | 1ul << three_semi_real;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool has_quadruple(DescValue v)
{
  static const unsigned long mask =
      1ul << two_imaginary_single_double_fixed
    | 1ul << two_real_single_double_fixed;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

// It came to light in October 2016 that links giving an even length difference
// had to be singled out in certain circumstances; this had not been expected
bool has_october_surprise(DescValue v)
{
  return generator_length(v) == (has_defect(v) ? 3 : 2);
}

bool is_proper_ascent(DescValue v)
{
  return not(is_descent(v) or is_like_nonparity(v));
}

unsigned int generator_length(DescValue v)
{ return v<two_complex_ascent ? 1 : v<three_complex_ascent ? 2 : 3; }

unsigned int link_count(DescValue v)
{ switch(v)
  {
    // zero valued Cayleys: nothing recorded (cross action is trivial)
  case one_real_nonparity: case one_imaginary_compact:
  case one_imaginary_pair_switched: case one_real_pair_switched:
  case two_real_nonparity: case two_imaginary_compact:
  case two_imaginary_single_double_switched:
  case two_real_single_double_switched:
  case three_real_nonparity: case three_imaginary_compact:
    return 0;

    // complex cases (record only the cross action)
  case one_complex_ascent:
  case one_complex_descent:
  case two_complex_ascent:
  case two_complex_descent:
  case three_complex_ascent:
  case three_complex_descent:
    return 1;

    // semi cases do not record thei (trivial) cross action
  case two_semi_imaginary: case two_semi_real:
  case three_semi_imaginary: case three_real_semi:
  case three_imaginary_semi: case three_semi_real:
    return 1;

    // some single valued extended Cayleys use second link for cross action
  case one_imaginary_single:
  case one_real_single:
  case two_imaginary_single_single:
  case two_real_single_single:
    return 2;

    // double valued Cayleys also have unrecorded trivial cross actions
  case one_real_pair_fixed: case one_imaginary_pair_fixed:
  case two_real_double_double: case two_imaginary_double_double:
  case two_imaginary_single_double_fixed: case two_real_single_double_fixed:
    return 2;
  }
  assert(false); return -1; // keep compiler happy
}

// number of links that are ascents or descents (not real/imaginary cross)
unsigned int scent_count(DescValue v)
{ return has_double_image(v) ? 2 : link_count(v)==0 ? 0 : 1; }

// find element |n| such that |z(n)>=zz|
BlockElt ext_block::element(BlockElt zz) const
{
  BlockElt min=0, max=size();
  // loop invariant: |(min==0 or z(min-1)<zz) and (max<size() or z(max)>=zz)|
  while (max>min)
  {
    BlockElt x=(min+max)/2;
    if (z(x)<zz)
      min=x+1;
    else
      max=x;
  }
  assert(min==max);
  return min;
}

// compute |bgv-(bgv+t_bits)*(1+theta)/2 == (bgv-t_bits-(bgv+t_bits)*theta)/2|
Coweight ell (const KGB& kgb, KGBElt x)
{ auto diff= (kgb.base_grading_vector()-kgb.torus_factor(x)).normalize();
  assert(diff.denominator()==1);
  return Coweight(diff.numerator().begin(),diff.numerator().end());
}

void validate(const param& E)
{
  const auto& i_tab = E.rc().innerClass().involution_table();
  const auto& rd = E.rc().innerClass().rootDatum();
  const auto& theta = i_tab.matrix(E.tw);
  const auto& delta = E.ctxt.delta();
  assert(delta*theta==theta*delta);
  assert((delta-1)*E.lambda_rho==(1-theta)*E.tau);
  if (not ((delta-1).right_prod(E.l)==(theta+1).right_prod(E.t)))
    {
      auto ec = E.ctxt;
      Weight gamma_numer(ec.gamma().numerator().begin(),
	     ec.gamma().numerator().end());
      unsigned int gamma_denom = ec.gamma().denominator();
      prettyprint::printVector(std::cout << "gamma_denom = "
			       << gamma_denom << ", gamma_numer = "
			       ,gamma_numer);
      std::cout << std::endl;
    }
  //  assert((delta-1).right_prod(E.l)==(theta+1).right_prod(E.t));
  if (not ((E.ctxt.g_rho_check()-E.l)*(1-theta)).numerator().isZero())
    {
      auto ec = E.ctxt;
      prettyprint::printVector(std::cout << "l = ",E.l);
      std::cout << std::endl;
      Coweight g_numer(ec.g_rho_check().numerator().begin(),
		       ec.g_rho_check().numerator().end());
      unsigned int g_denom = ec.g_rho_check().denominator();
      prettyprint::printVector(std::cout << "g_rho_check_denom = "
			       << g_denom << ", g_rho_check_numer = "
			       ,g_numer);
      std::cout << std::endl;
    }
  // assert(((E.ctxt.g_rho_check()-E.l)*(1-theta)).numerator().isZero());
  assert(((theta+1)*(E.ctxt.gamma()-E.lambda_rho-rho(rd)))
	 .numerator().isZero());
  ndebug_use(delta); ndebug_use(theta); ndebug_use(rd);
}

param::param (const context& ec, const StandardRepr& sr, bool flipped)
  : ctxt(ec)
  , tw(ec.rc().kgb().involution(sr.x()))
  , l(ell(ec.realGroup().kgb(),sr.x()))
  , lambda_rho(ec.rc().lambda_rho(sr))
  , tau(matreduc::find_solution(1-theta(),(delta()-1)*lambda_rho))
  , t(matreduc::find_solution
	(theta().transposed()+1,(delta()-1).right_prod(l)))
  , flipped(flipped)
{
  validate(*this);
}

param::param (const context& ec,
	      KGBElt x, const Weight& lambda_rho, bool flipped)
  : ctxt(ec)
  , tw(ec.realGroup().kgb().involution(x))
  , l(ell(ec.realGroup().kgb(),x))
  , lambda_rho(lambda_rho)
  , tau(matreduc::find_solution(1-theta(),(delta()-1)*lambda_rho))
  , t(matreduc::find_solution
	(theta().transposed()+1,(delta()-1).right_prod(l)))
  , flipped(flipped)
{
  validate(*this);
}

param::param (const context& ec, const TwistedInvolution& tw,
	      Weight lambda_rho, Weight tau, Coweight l, Coweight t,
	      bool flipped)
  : ctxt(ec), tw(tw)
  , l(std::move(l))
  , lambda_rho(std::move(lambda_rho))
  , tau(std::move(tau))
  , t(std::move(t))
  , flipped(flipped)
{
  validate(*this);
}

bool in_L_image(Weight beta,WeightInvolution&& A)
{ int_Matrix L,R;
  auto inv_fact = matreduc::diagonalise(std::move(A),L,R);
  int_Vector image = L*beta;

  unsigned i;
  for (i=0; i<inv_fact.size(); ++i)
    if (image[i]%inv_fact[i]!=0)
      return false;
  for (/* continue with |i| */ ; i<image.size(); ++i)
    if (image[i]!=0)
      return false;
  return true;
}

bool in_R_image(WeightInvolution&& A,Coweight b)
{ int_Matrix L,R;
  auto inv_fact = matreduc::diagonalise(std::move(A),L,R);
  int_Vector image = R.right_prod(b);

  unsigned i;
  for (i=0; i<inv_fact.size(); ++i)
    if (image[i]%inv_fact[i]!=0)
      return false;
  for (/* continue with |i| */ ; i<image.size(); ++i)
    if (image[i]!=0)
      return false;
  return true;
}

// whether |E| and |F| lie over equivalent |StandardRepr| values
bool same_standard_reps (const param& E, const param& F)
{
  if (&E.ctxt!=&F.ctxt)
  { if (&E.ctxt.innerClass()!=&F.ctxt.innerClass())
      throw std::runtime_error
	("Comparing extended parameters from different inner classes");
    if (E.delta()!=F.delta()
	or E.ctxt.g_rho_check()!=F.ctxt.g_rho_check()
	or E.ctxt.gamma()!=F.ctxt.gamma())
      return false;
  } // otherwise there might still be a match, so fall through
  return E.theta()==F.theta()
    and in_R_image(E.theta()+1,E.l-F.l)
    and in_L_image(E.lambda_rho-F.lambda_rho,E.theta()-1);
}

KGBElt param::x() const
{ TitsElt a(ctxt.innerClass().titsGroup(),TorusPart(l),tw);
  return rc().kgb().lookup(a);
}

/*
  This function serves to replace and circumvent |Rep_context::make_dominant|,
  which maps any ordinary parameter to one with a dominant |gamma| component,
  and moreover descends through simgular complex descents in the block to the
  lowest parameter equivalent to the inital parameter. The difference with
  that method is that here we keep track of all extended parameter components,
  transforming them from the default choices at the initial elemnt, and at the
  end comparing with the default choices at the final paremeter, recording the
  sign in |flipped|.

  This is intended for use in deformation, and the initial extended parameter
  components are those inherited from |sr| before scaling its |nu| part by
  |factor|. The user should make sure |sr| itself has dominant |gamma|, which
  moreover is assumed to be fixed by |delta| (if not, don't use this function).
 */
StandardRepr scaled_extended_dominant // result will have its |gamma()| dominant
(const Rep_context rc,
 const StandardRepr& sr, const WeightInvolution& delta,
 Rational factor, // |z.nu()| is scaled by |factor| first
 bool& flipped // records whether an extended flip was recorded
 )
{ const RootDatum& rd=rc.rootDatum(); const KGB& kgb = rc.kgb();
  const ext_gens orbits = rootdata::fold_orbits(rd,delta);
  assert(is_dominant_ratweight(rd,sr.gamma())); // dominant
  assert(((delta-1)*sr.gamma().numerator()).isZero()); // $\delta$-fixed
  // First approximation to result is scaled input; will later be overwritten
  StandardRepr result = rc.sr(sr.x(),rc.lambda_rho(sr),sr.gamma()*factor);

  // it will be convenent to have a working copy of the numerator of |gamma|
  Weight gamma_numer(result.gamma().numerator().begin(),
		     result.gamma().numerator().end());
  context old_ctxt(rc,delta,result.gamma());
  // class |param| cannot change its |gamma|, so work on separate components
  Weight lr, tau; Coweight l,t;
  { // context ctxt(rc,delta,result.gamma()); // scaffolding for construction
    param E(old_ctxt,result); // default extend |result| to extended parameter
    lr=E.lambda_rho; tau=E.tau; l=E.l; t=E.t; // and copy fields to variables
   }
   WeightInvolution theta = kgb.involution_matrix(sr.x());
   const RatCoweight& g_r = rc.realGroup().g_rho_check();

  assert((delta-1).right_prod(l)==(theta+1).right_prod(t));
  assert(((g_r-l)*(1-theta)).numerator().isZero());
  KGBElt x = result.x(); // another variable, for convenience

  const auto grc = old_ctxt.g_rho_check();
  int denom=grc.denominator();
  l*=denom; // scale to make shift applied to |l| below an integer vector
  Coweight l_offset(grc.numerator().begin(),grc.numerator().end()); // convert
  l-=l_offset; // shift so that reflections can apply directly

  int_Vector r_g_eval (rd.semisimpleRank()); // evaluations at |-gr|
  { // const RatCoweight& g_r=rc.realGroup().g_rho_check();
    for (unsigned i=0; i<r_g_eval.size(); ++i)
      r_g_eval[i] = -g_r.dot(rd.simpleRoot(i));
  }
  // since |gamma| reflects along, our action with be affine about $-\rho$
  const int_Vector ones(rd.semisimpleRank(),1);

  { unsigned i; // index into |orbits|
    do
      for (i=0; i<orbits.size(); ++i)
	if (kgb.status(x).isComplex(orbits[i].s0))
	{ const auto& s=orbits[i];
	  const auto& alpha_v = rd.simpleCoroot(s.s0);
	  int v=alpha_v.dot(gamma_numer);
	  if (v<0 or (v==0 and kgb.isDescent(s.s0,x)))
	  { if (v<0)
	      rd.act(s.w_kappa,gamma_numer); // change inf.char representative
            else // 3Cr excluded, as it would make |s.s0+s.s1|
	         // real singular parity
              assert(s.length()!=3 or kgb.cross(s.w_kappa,x)!=x);
	    if (v<0 or s.length()!=2 or kgb.cross(s.s0,x)!=kgb.cross(s.s1,x))
	    {
	      rd.shifted_act(s.w_kappa,lr,ones);
	      rd.act(s.w_kappa,tau);
	      // rd.shifted_dual_act(l,s.w_kappa,r_g_eval);
	      rd.dual_act(l,s.w_kappa); // this was in old master,
	      // where l was shifted first.
	      rd.dual_act(t,s.w_kappa);
	      x = kgb.cross(s.w_kappa,x);
	      theta = kgb.involution_matrix(x);
	      assert((delta-1).right_prod(l)==(theta+1).right_prod(t));
	      assert((1-theta).right_prod(l).isZero());
	    }
	    else // we have a singular 2Cr descent; do just one reflection
	    { // corrections to |tau| and |t| are as in 2Cr case of |star| below
	      const int f = alpha_v.dot(lr)+1;
	      const auto& alpha = rd.simpleRoot(s.s0);
	      lr -= alpha*f; // equivalently |rd.simple_reflect(s.s0,lr,1)|
	      rd.simple_reflect(s.s0,tau,-f); // shift |-f| adds extra |alpha*f|
	      const int df = l.dot(alpha)+r_g_eval[s.s0];
	      //	      l -= alpha_v*df;
	      // |rd.simple_coreflect(l,s.s0,r_g_eval[s.s0]);|
	      rd.simple_coreflect(l,s.s0);
	      rd.simple_coreflect(t,s.s0,df); // |df| subs extra |alpha_v*df|
	      x = kgb.cross(s.s0,x);
	      theta = kgb.involution_matrix(x);
	      assert((delta-1).right_prod(l)==(theta+1).right_prod(t));
	      assert((1-theta).right_prod(l).isZero());
	    }
	    break;
	  }
	} // |for(s)|, if |isComplex|
    while(i<orbits.size()); // continue until above |for| runs to completion
  } // end of transformation of extended parameter components
  l+=l_offset;
  l/=denom; //shift and scale back to original

  // since |gamma| may have changed, we only now build our |context|
  context ctxt(rc,delta, RatWeight(gamma_numer,result.gamma().denominator()));
  // now ensure that |E| gets matching |gamma| and |theta| (for flipped test)
  param E(ctxt,kgb.involution(x),lr,tau,l,t);

  // finally extract |StandardRepr| from |E|, overwriting |result|
  result = rc.sr_gamma(x,E.lambda_rho,ctxt.gamma());

  // but the whole point of this function is to record the relative flip too!
  flipped = not same_sign(E,param(ctxt,result)); // compare |E| to default ext.
  return result;

} // |scaled_extended_dominant|

/*
 This function is destined to be used after |scaled_extended_dominant|, to
 express the standard representation as linear combination of ones without
 singular descents, while keeping track (unlike |Rep_context::expand_final|)
 of flips that might occur during the process.
 */
containers::sl_list<std::pair<StandardRepr,bool> > extended_finalise
  (const repr::Rep_context& rc,
   const StandardRepr& sr, const WeightInvolution& delta)
{ // in order that |singular_generators| generate the whole singular system:
  assert(is_dominant_ratweight(rc.rootDatum(),sr.gamma()));
  // must assume gamma dominant, DON'T call make_dominant here
  context ctxt(rc,delta,sr.gamma());
  const ext_gens orbits = rootdata::fold_orbits(ctxt.id(),delta);
  const RankFlags singular_orbits =
    reduce_to(orbits,singular_generators(ctxt.id(),sr.gamma()));

  containers::queue<param> to_do { param(ctxt,sr) };
  containers::sl_list<std::pair<StandardRepr,bool> > result;

  do
  { const param E= to_do.front();
    to_do.pop(); // we are done with |head|
    auto s = first_descent_among(singular_orbits,orbits,E);
    if (s>=orbits.size()) // no singular descents, so append to result
      result.emplace_back(std::make_pair(E.restrict(),not is_default(E)));
    else // |s| is a singular descent orbit
    { containers::sl_list<param> links;
      auto type = star(E,orbits[s],links);
      if (not is_like_compact(type)) // some descent, push to front of |to_do|
      { bool flip = has_october_surprise(type); // to undo extra flip |star|
	auto l_it=links.begin();
	l_it->flip(flip);
	to_do.push(*l_it);
	if (has_double_image(type)) // then append a second node after |head|
	{ ++l_it;
	  l_it->flip(flip);
	  to_do.push(*l_it);
	}
      }
    }
  }
  while(not to_do.empty());

  return result;
} // |extended_finalise|

#if 0 // unused code, but the formula is referred to in the comment below
int z (const param& E) // value modulo 4, exponent of imaginary unit $i$
{ return
    (E.l.dot((E.delta()-1)*E.tau) + 2*E.t.dot(E.lambda_rho)) % 4;
}
#endif

/*
  The quotient of |z| values across a Cayley transform will only be used when
  value of |t| is the same upstairs as downstairs. Then in the quotient of |z|
  values the second term contributes |t.dot(E.lambda_rho-F.lambda_rho)|.
  That difference is often zero for the Cayley transform by a \emph{simple}
  root |alpha| or (length 1) changes by |alpha| with |t.cdot(alpha)==0|; in
  those cases the second term in |z| contributes nothing. For Cayley by a
  non-simple root, although the quotient of |z| values might pick up a
  contribution from the second term due to a |Cayley_shift| added to
  |lambda_rho|, it turns out that the Right Thing is not to use that quotient,

  hmmm?

  but the quotient evaluated at the simple situation. So for these cases we
  should just compute the contribution to the difference of values |z| that
  would come \emph{from its first term} above only.

  This function does that, and possibly flips the second parameter accordingly
*/
void z_align (const param& E, param& F)
{ // assert(E.t==F.t); // we require preparing |t| upstairs to get this
  // assert(E.t.dot(E.lambda_rho-F.lambda_rho) == 0);
  int d = E.l.dot((E.delta()-1)*E.tau) - F.l.dot((F.delta()-1)*F.tau)
    // add correction so always works?
    + 2*E.t.dot(E.lambda_rho) - 2*F.t.dot(F.lambda_rho);
  assert(d%2==0);
  F.flipped &= (d%4!=0);
}

/*
  In some cases, notably 2i12, one or both of the Cayley transforms may
  involve (in the simple root case) a change |mu| in |lambda_rho| that does
  not necessarily satisfy |t.dot(mu)==0|. In such cases the previous version
  of |z_align| is insufficient, and we should include a contribution from the
  second term of the formula for |z|. But retrieving |mu| from the parameters
  |E| and |F| themselves is complicated by the posssible contribution from
  |Cayley_shift|, which contribution should be ignored;

  hmmm???

  however at the place
  of call the value of |mu| is explicitly available, so we ask here to pass
  |t.dot(mu)| as third argument |t_mu|.
 */
void z_align (const param& E, param& F, int t_mu)
{
  assert(E.t.dot(E.lambda_rho-F.lambda_rho) == t_mu);
  z_align(E,F);
  F.flip(t_mu%2!=0);
}



// this implements (comparison using) the formula from Proposition 16 in
// "Parameters for twisted repressentations" (with $\delta-1 = -(1-\delta)$
// the relation is symmetric in |E|, |F|, although not obviously so
bool same_sign (const param& E, const param& F)
{
  assert(same_standard_reps(E,F));
  const WeightInvolution& delta = E.delta();
  Weight kappa1=E.tau, kappa2=F.tau;
  kappa1 -= delta*kappa1;
  kappa2 -= delta*kappa2;
  int i_exp = E.l.dot(kappa1) - F.l.dot(kappa2);
  assert(i_exp%2==0);
  int n1_exp =
    (F.l-E.l).dot(E.tau) + F.t.dot(F.lambda_rho-E.lambda_rho);
  return ((i_exp/2+n1_exp)%2==0)!=(E.is_flipped()!=F.is_flipped());
}

bool same_sign_with_one_of (const param& E, const param& F1, const param& F2)
{ return
    same_standard_reps(E,F1) ? same_sign(E,F1)
    : same_standard_reps(E,F2) ? same_sign(E,F2)
    : throw std::runtime_error("Neither candidate has same standard repn");
}

void ext_block::add_neighbours
  (BlockEltList& dst, weyl::Generator s, BlockElt n) const
{
  const BlockEltPair& links = data[s][n].links;
  if (links.first==UndefBlock)
    return;
  dst.push_back(links.first);
  if (links.second==UndefBlock)
    return;
  dst.push_back(links.second);
}

void ext_block::flip_edge(weyl::Generator s, BlockElt x, BlockElt y)
{
  BlockEltPair p= data[s][x].links;
  int i= p.first==y ? 0 : p.second==y ? 1 : -1;
  assert(i>=0);
  info[x].flips[i].flip(s);
}

// whether link for |s| from |x| to |y| has a sign flip attached
int ext_block::epsilon(weyl::Generator s, BlockElt x, BlockElt y ) const
{
  BlockEltPair p= data[s][x].links;
  int i= p.first==y ? 0 : 1;
  assert(i==0 or p.second==y);

  return info[x].flips[i][s] ? -1 : 1;
}

BlockEltList ext_block::down_set(BlockElt n) const
{
  BlockEltList result; result.reserve(rank());
  for (weyl::Generator s=0; s<rank(); ++s)
  {
    const DescValue type = descent_type(s,n);
    if (is_descent(type) and not is_like_compact(type))
    {
      result.push_back(data[s][n].links.first);
      if (has_double_image(type))
	result.push_back(data[s][n].links.second);
    }
  }
  return result;
}

context::context
  (const repr::Rep_context& rc,
   const WeightInvolution& delta,
   const RatWeight& gamma)
    : d_rc(rc)
    , d_delta(delta)
    , d_gamma(gamma)
    , integr_datum(integrality_datum(rc.rootDatum(),gamma))
    , sub(SubSystem::integral(rc.rootDatum(),gamma))
    , pi_delta(rc.rootDatum().rootPermutation(d_delta))
    , twist()
    , lambda_shifts (integr_datum.semisimpleRank())
    , l_shifts (integr_datum.semisimpleRank())
{
  const RootDatum& rd = rc.rootDatum();
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    twist[s] = rd.simpleRootIndex(delta_of(rd.simpleRootNbr(s)));

  // the reflections for |E.lambda_rho| pivot around $\gamma-\rho$
  const RatWeight gamma_rho = gamma - rho(rd);
  for (unsigned i=0; i<lambda_shifts.size(); ++i)
    lambda_shifts[i] = -gamma_rho.dot(integr_datum.simpleCoroot(i));
  // the reflections for |E.l| pivot around |g_rho_check()|
  const RatCoweight& g_rho_check = this->g_rho_check();
  for (unsigned i=0; i<l_shifts.size(); ++i)
    l_shifts[i] = -g_rho_check.dot(integr_datum.simpleRoot(i));
}


// old version of |extended_type| below, this one uses |Hermitian_dual| method
DescValue extended_type(const Block_base& block, BlockElt z, const ext_gen& p,
			BlockElt& link)
{
  switch (p.type)
  {
  case ext_gen::one:
    switch (block.descentValue(p.s0,z))
    {
    case DescentStatus::ComplexAscent:
      link=block.cross(p.s0,z); return one_complex_ascent;
    case DescentStatus::ComplexDescent:
      link=block.cross(p.s0,z); return one_complex_descent;
    case DescentStatus::RealNonparity:
      link=UndefBlock; return one_real_nonparity;
    case DescentStatus::ImaginaryCompact:
      link=UndefBlock; return one_imaginary_compact;
    case DescentStatus::ImaginaryTypeI:
      link=block.cayley(p.s0,z).first; return one_imaginary_single;
    case DescentStatus::RealTypeII:
      link=block.inverseCayley(p.s0,z).first; return one_real_single;
    case DescentStatus::ImaginaryTypeII:
      link=block.cayley(p.s0,z).first;
      if (link!=UndefBlock and block.Hermitian_dual(link)==link)
	return one_imaginary_pair_fixed;
      link=UndefBlock; return one_imaginary_pair_switched;
    case DescentStatus::RealTypeI:
      link=block.inverseCayley(p.s0,z).first;
      if (link!=UndefBlock and block.Hermitian_dual(link)==link)
	return one_real_pair_fixed;
      link=UndefBlock; return one_real_pair_switched;
    }
  case ext_gen::two:
    switch (block.descentValue(p.s0,z))
    {
    case DescentStatus::ComplexAscent:
      link=block.cross(p.s0,z);
      if (link==block.cross(p.s1,z))
	return two_semi_imaginary; // just a guess if |link==UndefBlock|
      if (link!=UndefBlock)
	link=block.cross(p.s1,link);
      return two_complex_ascent;
    case DescentStatus::ComplexDescent:
      link=block.cross(p.s0,z);
      if (link==block.cross(p.s1,z))
	return two_semi_real; // just a guess if |link==UndefBlock|
      if (link!=UndefBlock)
	link=block.cross(p.s1,link);
      return two_complex_descent;
    case DescentStatus::RealNonparity:
      link=UndefBlock; return two_real_nonparity;
    case DescentStatus::ImaginaryCompact:
      link=UndefBlock; return two_imaginary_compact;
    case DescentStatus::ImaginaryTypeI:
      link=block.cayley(p.s0,z).first;
      if (link==UndefBlock)
	return two_imaginary_single_single; // really just a guess
      link=block.cayley(p.s1,link).first;
      if (link==UndefBlock)
	return two_imaginary_single_single; // really just a guess
      assert(block.Hermitian_dual(link)==link);
      return block.descentValue(p.s0,link)==DescentStatus::RealTypeI
	? two_imaginary_single_single : two_imaginary_single_double_fixed;
    case DescentStatus::RealTypeII:
      link=block.inverseCayley(p.s0,z).first;
      if (link==UndefBlock)
	return two_real_single_single; // really just a guess
      link=block.inverseCayley(p.s1,link).first;
      if (link==UndefBlock)
	return two_real_single_single; // really just a guess
      assert(block.Hermitian_dual(link)==link);
      return block.descentValue(p.s0,link)==DescentStatus::ImaginaryTypeII
	? two_real_single_single : two_real_single_double_fixed;
    case DescentStatus::ImaginaryTypeII:
      link=block.cayley(p.s0,z).first;
      if (link==UndefBlock)
	return two_imaginary_double_double;
      link=block.cayley(p.s1,link).first;
      if (link==UndefBlock)
	return two_imaginary_double_double;
      if (block.Hermitian_dual(link)!=link)
      {
	link=block.cross(p.s0,link);
	assert(link==UndefBlock or block.Hermitian_dual(link)==link);
      }
      return two_imaginary_double_double;
    case DescentStatus::RealTypeI:
      link=block.inverseCayley(p.s0,z).first;
      if (link==UndefBlock)
	return two_real_double_double;
      link=block.inverseCayley(p.s1,link).first;
      if (link==UndefBlock)
	return two_real_double_double;
      if (block.Hermitian_dual(link)!=link)
      {
	link=block.cross(p.s0,link);
	assert(link==UndefBlock or block.Hermitian_dual(link)==link);
      }
      return two_real_double_double;
    }
  case ext_gen::three:
    switch (block.descentValue(p.s0,z))
    {
    case DescentStatus::RealNonparity:
      link=UndefBlock; return three_real_nonparity;
    case DescentStatus::ImaginaryCompact:
      link=UndefBlock; return three_imaginary_compact;
    case DescentStatus::ComplexAscent:
      link=block.cross(p.s0,z);
      if (link==UndefBlock)
	return three_complex_ascent; // just a guess
      if (link==block.cross(p.s1,link))
      {
	assert(block.descentValue(p.s1,link)==
	       DescentStatus::ImaginaryTypeII);
	link=block.cayley(p.s1,link).first;
	if (link!=UndefBlock and block.Hermitian_dual(link)!=link)
	{
	  link=block.cross(p.s1,link);
	  assert(link==UndefBlock or block.Hermitian_dual(link)==link);
	}
	return three_semi_imaginary;
      }
      link=block.cross(p.s1,link);
      if (link!=UndefBlock)
	link=block.cross(p.s0,link);
      if (link!=UndefBlock)
	assert(block.Hermitian_dual(link)==link);
      return three_complex_ascent;
    case DescentStatus::ComplexDescent:
      link=block.cross(p.s0,z);
      if (link==UndefBlock)
	return three_complex_descent; // just a guess
      if (link==block.cross(p.s1,link))
      {
	assert(block.descentValue(p.s1,link)==DescentStatus::RealTypeI);
	link=block.inverseCayley(p.s1,link).first;
	if (link!=UndefBlock and block.Hermitian_dual(link)!=link)
	{
	  link=block.cross(p.s1,link);
	  assert(link==UndefBlock or block.Hermitian_dual(link)==link);
	}
	return three_semi_real;
      }
      link=block.cross(p.s1,link);
      if (link!=UndefBlock)
	link=block.cross(p.s0,link);
      if (link!=UndefBlock)
	assert(block.Hermitian_dual(link)==link);
      return three_complex_descent;
    case DescentStatus::ImaginaryTypeI:
      link=block.cayley(p.s0,z).first;
      if (link!=UndefBlock)
      {
	link=block.cross(p.s1,link);
	if (block.cayley(p.s1,z).first!=UndefBlock)
	  assert(link==block.cross(p.s0,block.cayley(p.s1,z).first));
      }
      else if ((link=block.cayley(p.s1,z).first)!=UndefBlock)
	link=block.cross(p.s0,link);
      if (link!=UndefBlock)
	assert(block.Hermitian_dual(link)==link);
      return three_imaginary_semi;
    case DescentStatus::RealTypeII:
      link=block.inverseCayley(p.s0,z).first;
      if (link!=UndefBlock)
      {
	link=block.cross(p.s1,link);
	if (block.inverseCayley(p.s1,z).first!=UndefBlock)
	  assert(link==block.cross(p.s0,block.inverseCayley(p.s1,z).first));
      }
      else if ((link=block.inverseCayley(p.s1,z).first)!=UndefBlock)
	link=block.cross(p.s0,link);
      if (link!=UndefBlock)
	assert(block.Hermitian_dual(link)==link);
      return three_real_semi;
    case DescentStatus::ImaginaryTypeII: case DescentStatus::RealTypeI:
      assert(false); // these cases should never occur
    }
  } // |switch (p.type)|
  assert(false); return one_complex_ascent; // keep compiler happy
} // |extended_type|

// auxiliary function to recognise local situation in |ext_block| construction
// the following function assumes a full block, and precomputed |fixed_points|
DescValue extended_type(const Block_base& block, BlockElt z, const ext_gen& p,
			BlockElt& link, const BitMap& fixed_points)
{
  switch (p.type)
  {
  case ext_gen::one:
    switch (block.descentValue(p.s0,z))
    {
    case DescentStatus::ComplexAscent:
      return link=block.cross(p.s0,z), one_complex_ascent;
    case DescentStatus::ComplexDescent:
      return link=block.cross(p.s0,z), one_complex_descent;
    case DescentStatus::RealNonparity:
      return link=UndefBlock, one_real_nonparity;
    case DescentStatus::ImaginaryCompact:
      return link=UndefBlock, one_imaginary_compact;
    case DescentStatus::ImaginaryTypeI:
      return link=block.cayley(p.s0,z).first, one_imaginary_single;
    case DescentStatus::RealTypeII:
      return link=block.inverseCayley(p.s0,z).first, one_real_single;
    case DescentStatus::ImaginaryTypeII:
      { const BlockElt t=block.cayley(p.s0,z).first;
	if (fixed_points.isMember(t))
	  return link=t, one_imaginary_pair_fixed;
	return link=UndefBlock, one_imaginary_pair_switched;
      }
    case DescentStatus::RealTypeI:
      { const BlockElt t=block.inverseCayley(p.s0,z).first;
	if (fixed_points.isMember(t))
	  return link=t, one_real_pair_fixed;
	return link=UndefBlock, one_real_pair_switched;
      }
    }
  case ext_gen::two:
    switch (block.descentValue(p.s0,z))
    {
    case DescentStatus::ComplexAscent:
      { const BlockElt t=block.cross(p.s0,z);
	if (t==block.cross(p.s1,z))
	  return link=t, two_semi_imaginary;
	return link=block.cross(p.s1,t),  two_complex_ascent;
      }
    case DescentStatus::ComplexDescent:
     { const BlockElt t=block.cross(p.s0,z);
	if (t==block.cross(p.s1,z))
	  return link=t, two_semi_real;
	return link=block.cross(p.s1,t), two_complex_descent;
     }
    case DescentStatus::RealNonparity:
      return link=UndefBlock, two_real_nonparity;
    case DescentStatus::ImaginaryCompact:
      return link=UndefBlock, two_imaginary_compact;
    case DescentStatus::ImaginaryTypeI:
      { const BlockElt t=block.cayley(p.s0,z).first;
        if (block.descentValue(p.s1,t)==DescentStatus::ImaginaryTypeI)
	  return link=block.cayley(p.s1,t).first, two_imaginary_single_single;
        return fixed_points.isMember(link=block.cayley(p.s1,t).first)
	  ? two_imaginary_single_double_fixed
	  : (link=UndefBlock, two_imaginary_single_double_switched);
      }
    case DescentStatus::RealTypeII:
      { const BlockElt t=block.inverseCayley(p.s0,z).first;
	if (block.descentValue(p.s1,t)==DescentStatus::RealTypeII)
	  return link=block.inverseCayley(p.s1,t).first,two_real_single_single;
	return fixed_points.isMember(link=block.inverseCayley(p.s1,t).first)
	  ? two_real_single_double_fixed
	  : (link=UndefBlock, two_real_single_double_switched);
      }
    case DescentStatus::ImaginaryTypeII:
      link=block.cayley(p.s1,block.cayley(p.s0,z).first).first;
      if (not fixed_points.isMember(link))
	link=block.cross(p.s0,link), assert(fixed_points.isMember(link));
      return two_imaginary_double_double;
    case DescentStatus::RealTypeI:
      link=block.inverseCayley(p.s1,block.inverseCayley(p.s0,z).first).first;
      if (not fixed_points.isMember(link))
	link=block.cross(p.s0,link), assert(fixed_points.isMember(link));
      return two_real_double_double;
    }
  case ext_gen::three:
    switch (block.descentValue(p.s0,z))
    {
    case DescentStatus::RealNonparity:
      return link=UndefBlock, three_real_nonparity;
    case DescentStatus::ImaginaryCompact:
      return link=UndefBlock, three_imaginary_compact;
    case DescentStatus::ComplexAscent:
      { const BlockElt t=block.cross(p.s0,z);
	if (t==block.cross(p.s1,t))
	{
	  assert(block.descentValue(p.s1,t)==DescentStatus::ImaginaryTypeII);
	  link=block.cayley(p.s1,t).first;
	  if (not fixed_points.isMember(link))
	    link=block.cross(p.s1,link), assert(fixed_points.isMember(link));
	  return three_semi_imaginary;
	}
	link=block.cross(p.s0,block.cross(p.s1,t));
	assert(fixed_points.isMember(link));
	return three_complex_ascent;
      }
    case DescentStatus::ComplexDescent:
      link=block.cross(p.s0,z);
      { const BlockElt t=block.cross(p.s0,z);
	if (t==block.cross(p.s1,t))
	{
	  assert(block.descentValue(p.s1,t)==DescentStatus::RealTypeI);
	  link=block.inverseCayley(p.s1,link).first;
	  if (not fixed_points.isMember(link))
	    link=block.cross(p.s1,link), assert(fixed_points.isMember(link));
	  return three_semi_real;
	}
	link=block.cross(p.s0,block.cross(p.s1,link));
	assert(fixed_points.isMember(link));
	return three_complex_descent;
      }
    case DescentStatus::ImaginaryTypeI:
      link=block.cross(p.s1,block.cayley(p.s0,z).first);
      assert(link==block.cross(p.s0,block.cayley(p.s1,z).first));
      assert(fixed_points.isMember(link));
      return three_imaginary_semi;
    case DescentStatus::RealTypeII:
      link=block.cross(p.s1,block.inverseCayley(p.s0,z).first);
      assert(link==block.cross(p.s0,block.inverseCayley(p.s1,z).first));
      assert(fixed_points.isMember(link));
      return three_real_semi;
    case DescentStatus::ImaginaryTypeII: case DescentStatus::RealTypeI:
      assert(false); // these cases should never occur
    }
  } // |switch (p.type)|
  assert(false); return one_complex_ascent; // keep compiler happy
} // |extended_type|

// involution for dual KGB is just $\delta$ transposed, no factor $-w_0$ here
// both |kgb| and |dual_kgb| must be known to be twist-stable
BlockElt twisted (const Block& block,
		  const KGB& kgb, const KGB& dual_kgb, // all are needed
		  BlockElt z,
		  const WeightInvolution& delta)
{ return block.element
    (kgb.twisted(block.x(z),delta),
     dual_kgb.twisted(block.y(z),delta.transposed()));
}

BlockElt twisted
  (const param_block& block, BlockElt z, const WeightInvolution& delta)
{
  return block.lookup(block.context().twisted(block.sr(z),delta));
}


/*
  An auxiliary routine to compute extended parameters across complex links.
  The situation is complicated by the fact that the cross action is by a
  generator of the folded integral system, so we need to exapand it first into
  a product of |length<=3| integral generators, and then have those generators
  act on the components of |E|. For the purpose of changing |E.tw| we further
  develop those generators into reflection words for the full root datum, but
  the reflection action on the othe components can be done more directly.

  However, though the integral generators are complex, they action of those
  reflection words need not be purely complex, which implies that the effect
  on the |lambda_rho| and |l| components are not purely reflection. The
  difference with respect to pure reflection action can be computed comparing
  |rho_r| values (half sums of positive real roots in the full system) with
  the reflected image of that value taken at the starting point, respectively
  (for |l|) the same thing with |rho_check_imaginary|. This is done by the
  "correction" terms below.
 */
param complex_cross(const ext_gen& p, param E) // by-value for |E|, modified
{ const RootDatum& rd = E.rc().rootDatum();
  const RootDatum& id = E.ctxt.id();
  const InvolutionTable& i_tab = E.rc().innerClass().involution_table();
  auto &tW = E.rc().twistedWeylGroup(); // caution: |p| refers to integr. datum

  InvolutionNbr theta = i_tab.nr(E.tw);
  Weight rho_r_shift = rd.twoRho(i_tab.real_roots(theta));
  Coweight dual_rho_im_shift = rd.dual_twoRho(i_tab.imaginary_roots(theta));

  for (unsigned i=p.w_kappa.size(); i-->0; ) // at most 3 letters, right-to-left
  { weyl::Generator s=p.w_kappa[i]; // generator for integrality datum
    tW.twistedConjugate(E.ctxt.subsys().reflection(s),E.tw);
    id.simple_reflect(s,E.lambda_rho,E.ctxt.lambda_shift(s));
    id.simple_reflect(s,rho_r_shift);
    id.simple_reflect(s,E.tau);
    id.simple_coreflect(E.l,s,E.ctxt.l_shift(s));
    id.simple_coreflect(dual_rho_im_shift,s);
    id.simple_coreflect(E.t,s);
  }

  rho_r_shift -= rd.twoRho(i_tab.real_roots(i_tab.nr(E.tw)));
  //  assert(rho_r_shift==(rho_r_shift/2) + (rho_r_shift/2));
  // rho_r_shift/=2; // now it is just a sum of (real) roots
  E.lambda_rho -= rho_r_shift; //add real on new side, subtract on original

  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // diff of $\delta$-fixed

  dual_rho_im_shift -= rd.dual_twoRho(i_tab.imaginary_roots(i_tab.nr(E.tw)));
  // dual_rho_im_shift/=2; // now it is just a sum of (imaginary) coroots
  assert(dual_rho_im_shift.isZero());
  E.l -= dual_rho_im_shift;

  assert(E.ctxt.delta().right_prod(dual_rho_im_shift)==dual_rho_im_shift);

  validate(E);
  return E;
} // |complex_cross|


WeylWord fixed_conjugate_simple (const context& ctxt, RootNbr& alpha)
{ const RootDatum& rd = ctxt.innerClass().rootDatum();

  WeylWord result;
  while (not rd.is_simple_root(alpha)) // also |break| halfway is possible
  {
    weyl::Generator s = rd.descent_set(alpha)
      .andnot(rd.ascent_set(ctxt.delta_of(alpha))).firstBit();
    assert(s<rd.semisimpleRank()); // exists for positive non-simple roots
    weyl::Generator t = ctxt.twisted(s);
    if (rd.simple_reflected_root(s,alpha)==rd.simpleRootNbr(t))
      break; // |alpha| is sum of (non-commuting) simple roots |s|,|twisted(s)|
    result.push_back(s);
    rd.simple_reflect_root(s,alpha);
    if (s!=t) // second generator for cases of length 2,3
    { result.push_back(t);
      rd.simple_reflect_root(t,alpha);
      if (rd.diagram_linked(s,t)) // then simple reflections |s,t| don't commute
      { // we have $sts=tst$, re-apply |s| to symmetrise the reflection word
	result.push_back(s);
	rd.simple_reflect_root(s,alpha);
      }
    }
  }
  std::reverse(result.begin(),result.end());
  return result;
} // |fixed_conjugate_simple|

/*
  for real Cayley transforms, one will subtract $\rho_r$ from |lambda_rho|
  before projecting it parallel to |alpha| so as to make |alpha_v| vanish on
  |gamma-lambda_rho-rho|. Here we compute from |E.lambda_rho|, corrected by
  that |shift|, the multiple of $\alpha/2$ that such a projection would add
  to |lambda_rho| (or equivalently, subtract from |gamma-lambda_rho-rho|).
*/
int level_a (const param& E, const Weight& shift, RootNbr alpha)
{
  const RootDatum& rd = E.rc().rootDatum();
  return (E.ctxt.gamma() - E.lambda_rho + shift).dot(rd.coroot(alpha))
    - rd.colevel(alpha); // final term $<\alpha^\vee,\rho>$
}

/*
  For the upstairs and downstairs conjugation scenario where we call
  |repr::Cayley_shift| to make a shift to |lambda_rho|, we also need the sign
  by which |delta| acts on wedge product for |repr::Cayley_shift| roots, more
  precisely those positive roots that |to_simple| maps to negative, and which
  are complex downstairs while they were real upstairs. In fact this is just
  full set of those positive-to-negative roots that are real upstairs.
 */
bool Cayley_shift_flip
  (const context& ec,
   InvolutionNbr theta_upstairs, // top of the link (more split Cartan)
   InvolutionNbr theta_downstairs, // at the bottom of the link
   const WeylWord& to_simple)
{ const RootDatum& rd = ec.rc().rootDatum();
  const InvolutionTable& i_tab = ec.innerClass().involution_table();
  RootNbrSet S = pos_to_neg(rd,to_simple) & i_tab.real_roots(theta_upstairs);
  RootNbrSet T = pos_to_neg(rd,to_simple) & i_tab.real_roots(theta_downstairs);
  unsigned countup=0; // will count 2-element |delta|-orbits upstairs
  for (auto it=S.begin(); it(); ++it)
    if (*it!=ec.delta_of(*it) and not rd.sumIsRoot(*it,ec.delta_of(*it)))
      ++countup;
  assert(countup%2==0); // since |S| is supposed to be $\delta$-stable
  unsigned countdown=0; // will count 2-element |delta|-orbits downstairs
  for (auto it=T.begin(); it(); ++it)
    if (*it!=ec.delta_of(*it) and not rd.sumIsRoot(*it,ec.delta_of(*it)))
      ++countdown;
  //  Weight gamma_numer(ec.gamma().numerator().begin(),
  //		     ec.gamma().numerator().end());
  //  unsigned int gamma_denom = ec.gamma().denominator();
  if (not (countdown==0)) std::cout << "countup = " << countup
				    << ", countdown = " << countdown
				    << ", thetaup = " << theta_upstairs
				    << ", thetadown = " << theta_downstairs
				    << std::endl;
  // if (not (countdown==0)) prettyprint::printVector(std::cout
  // << "gamma_denom = "
  //						   << gamma_denom
  //						   << ", gamma_numer = "
  //						   ,gamma_numer);
  // if (not (countdown==0)) std::cout << std::endl;
    //std::cout << "gamma_numer = " << gamma_numer
    //				    <<std::endl;
  // assert(countdown==0); // since that's what we see
  return (countup-countdown)%4!=0;
} // Cayley_shift_flip

  /* In the case of Cayley_shift_flip, the tau parameter of the target must also be modified by the sum of one real root from each delta orbit.
   */
  /* Weight tau_shift
  (const context& ec,
   InvolutionNbr theta_upstairs, // top of the link (more split Cartan)
   InvolutionNbr theta_downstairs, // at the bottom of the link
   const WeylWord& to_simple)
{ const RootDatum& rd = ec.rc().rootDatum();
  const InvolutionTable& i_tab = ec.innerClass().involution_table();
  RootNbrSet S = pos_to_neg(rd,to_simple) & i_tab.real_roots(theta_upstairs);
  RootNbrSet T = pos_to_neg(rd,to_simple) & i_tab.real_roots(theta_downstairs);
  Weight result(rd.rank(),0); // this will be one root from each delta pair
  for (auto it=S.begin(); it(); ++it)
    if (*it < ec.delta_of(*it))
      // (*it!=ec.delta_of(*it) and
	// not rd.sumIsRoot(*it,ec.delta_of(*it)))
      result += rd.root(*it); // add smaller root in delta orbit
  for (auto it=T.begin(); it(); ++it)
    if (*it < ec.delta_of(*it))
      result -= rd.root(*it); // subtract downstairs
  return result;
} // tau_shift
*/

// version of |type| that will also export signs for every element of |links|
DescValue star (const param& E,	const ext_gen& p,
		containers::sl_list<param>& links)
{
  param E0=E; // a copy of |E| that might be modified below to "normalise"
  DescValue result;

  const TwistedWeylGroup& tW = E.rc().twistedWeylGroup();
  const InnerClass& ic = E.rc().innerClass();
  const InvolutionTable& i_tab = ic.involution_table();
  const RootDatum& rd = E.rc().rootDatum();
  const RootDatum& integr_datum = E.ctxt.id();
  const SubSystem& subs = E.ctxt.subsys();
  const InvolutionNbr theta = i_tab.nr(E.tw);
  switch (p.type)
  {
  case ext_gen::one:
    { const Weight& alpha = integr_datum.simpleRoot(p.s0);
      const Coweight& alpha_v = integr_datum.simpleCoroot(p.s0);
      const RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);

      if (theta_alpha==n_alpha) // length 1 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	if (tf_alpha%2!=0) // then $\alpha$ is compact
	  return one_imaginary_compact; // quit here, do not collect \$200

	// noncompact case
	bool monoflip = false;
	const TwistedInvolution new_tw= tW.prod(subs.reflection(p.s0),E.tw);
	unsigned int d = i_tab.length(new_tw) - i_tab.length(E.tw);
	if (d%2==0) monoflip = not monoflip;
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1; // upstairs

	int tau_coef = alpha_v.dot(E.tau); // take $\tau_\alpha$ of table 2

	// try to make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs
	const Weight rho_r_shift = repr::Cayley_shift(ic,theta_p,theta,ww);
	//	const Weight tau_shift = tau_shift(E.ctxt,theta_p,theta,ww);
	const bool flipped = Cayley_shift_flip(E.ctxt,theta_p,theta,ww);
        if(flipped) std::cout << "1i flip" << std::endl;

	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(E.t.dot(alpha)==0); // follows from $\delta*\alpha=\alpha$

	Weight first; // maybe a root with |(1-delta)*first==alpha|
	RootNbr alpha_0 = 0;
	bool has_first=false;
	if (rd.is_simple_root(alpha_simple))
	  first = Weight(rd.rank(),0); // effectively not used in this case
	else
	{
	  has_first=true;
	  std::cout << "can't make simple alpha = " << n_alpha << std::endl;
	  --tau_coef; // the parity change and decrease are both relevant
	  weyl::Generator s = // first switched root index
	    rd.find_descent(alpha_simple);
	  alpha_0 = rd.permuted_root(rd.simpleRootNbr(s),ww);
	  first = // corresponding root summand, conjugated back
	      rd.root(rd.permuted_root(rd.simpleRootNbr(s),ww));
	  assert(alpha == first + E.ctxt.delta()*first );
	} // with this set-up, |alpha_simple| needs no more inspection

	// now separate cases; based on type 1 or 2 first
	if (matreduc::has_solution(th_1,alpha))
	{ // type 1, so extended type is 1i1
	  result = one_imaginary_single;

	  Weight diff = // called $-\sigma$ in table 2 of [Ptr] (NOTE MINUS)
	      matreduc::find_solution(th_1,alpha); // solutions are equivalent

	  param F(E.ctxt,new_tw,
		  E.lambda_rho + first + rho_r_shift,
		  E0.tau+diff*tau_coef,
		  E.l+alpha_v*(tf_alpha/2), E.t,
		  flipped);

 	  E0.l = tf_alpha%4==0 ? F.l+alpha_v : F.l; // for cross
	  assert(not same_standard_reps(E,E0));
	  F.lambda_rho-=first + rho_r_shift;
	  z_align(E,F);
	  //	  z_align(F,E0);
	  z_align(E,E0);
	  F.lambda_rho+=first + rho_r_shift;
	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E0)); // cross link
	} // end of 1i1 case
	else
	{ // imaginary type 2; now we need to distinguish 1i2f and 1i2s

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
	    assert(not matreduc::has_solution
		   (th_1, (E.ctxt.delta()-1)*(E.lambda_rho+rho_r_shift)));
	    return one_imaginary_pair_switched; // case 1i2s
	  }
	  result = one_imaginary_pair_fixed;  // what remains is case 1i2f
	  param F0(E.ctxt,new_tw,
		   E.lambda_rho + first + rho_r_shift,
		   E.tau - alpha*(tau_coef/2) - first, // + tau_shift,
		   E.l + alpha_v*(tf_alpha/2), E.t,
		   flipped);
	  param F1(E.ctxt,new_tw,
		   E.lambda_rho + first + rho_r_shift
		   + alpha, F0.tau, F0.l, E.t,
		   flipped);
	  if(has_first and (((-F0.lambda_rho - first - rho_r_shift)
			      //    rho_r_shift_old)
			      .dot(rd.coroot(alpha_0)) -
			     rd.colevel(alpha_0))%2 != 0 ))
	    F0.flipped = not F0.flipped; //test alpha0 nonparity at nu=0
	    if(has_first) std::cout << "parity test0 = " <<
			    (-F0.lambda_rho - rho_r_shift - first)
			    .dot(rd.coroot(alpha_0)) -
			    rd.colevel(alpha_0) << std::endl
				    << "alpha_0 = " << alpha_0
				    << ", alpha = " << n_alpha << std::endl;

	    if(has_first and (((-F1.lambda_rho - first - rho_r_shift)
			      // - rho_r_shift_old)
				.dot(rd.coroot(alpha_0)) -
			     rd.colevel(alpha_0))%2 != 0 ))
	    F1.flipped = not F1.flipped; //test alpha0 nonparity at nu=0
	    if(has_first) std::cout << "parity test1 = " <<
			    (-F1.lambda_rho + rho_r_shift)
			    .dot(rd.coroot(alpha_0)) -
			    rd.colevel(alpha_0) << std::endl;
	  F0.lambda_rho-=first+rho_r_shift;
	  F1.lambda_rho-=first+rho_r_shift;
	  z_align(E,F0);
	  z_align(E,F1);
	  F0.lambda_rho+=first+rho_r_shift;
	  F1.lambda_rho+=first+rho_r_shift;
	  links.push_back(std::move(F0)); // Cayley link
	  links.push_back(std::move(F1)); // Cayley link
	} // end of type 2 case
      } // end of length 1 imaginary case

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 1 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s0),E.tw);
	bool monoflip=false;
	unsigned int d = i_tab.length(E.tw) - i_tab.length(new_tw);
	if (d%2==0) monoflip = not monoflip;
	const auto theta_p = i_tab.nr(new_tw); // downstairs
	Weight rho_r_shift = repr::Cayley_shift(ic,theta,theta_p,ww);
	const bool flipped = Cayley_shift_flip(E.ctxt,theta,theta_p,ww);
	//	const Weight tau_shift = tau_shift(E.ctxt,theta,theta_p,ww);
        if(flipped) std::cout << "1r flip" << std::endl;
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	RootNbr alpha_0 = // maybe one of |alpha==alpha_0+alpha_1|
	  rd.is_simple_root(alpha_simple) ? 0 // unused
	  : rd.permuted_root(rd.simpleRootNbr(rd.find_descent(alpha_simple)),
			     ww);

	// test parity, taking into account modifications that will be applied
	bool shift_correct = // whether |alpha_0| is defined and real at |theta|
	  not rd.is_simple_root(alpha_simple) and
	  i_tab.root_involution(theta,alpha_0)==rd.rootMinus(alpha_0);
	bool flipped_correct = // shift_correct, and alpha_0 nonparity
			       // at nu=0
	  shift_correct and ( ((- E.lambda_rho + rho_r_shift)
				// - rho_r_shift_new)
			       .dot(rd.coroot(alpha_0))
			       - rd.colevel(alpha_0))%2!=0);
			      bool flipped1 = flipped_correct ?
			      not flipped : flipped;

	const int level = level_a(E,rho_r_shift,n_alpha) +
	   (shift_correct ? 1 : 0 ); // add 1 if |alpha_0| is defined and real

	if (level%2!=0) // nonparity
	   return one_real_nonparity; // no link added here

	const WeightInvolution& th_1 = i_tab.matrix(E.tw)-1; // upstairs
	bool type1 = matreduc::has_solution(th_1,alpha);

	Weight tau_correction; // adapt to integrality based change of lambda
	if (rd.is_simple_root(alpha_simple))
	  tau_correction = Weight(rd.rank(),0); // no correction needed here
	else
	{
	  const Weight a0 = rd.root(alpha_0);
	  assert(alpha == a0 + E.ctxt.delta()*a0);
	  if (shift_correct)
	  {
	    rho_r_shift += a0; // non delta-fixed contribution

	    // now we must add $d$ to $\tau$ with $(1-\theta')d=(1-\delta)*a0$
	    // since $\theta'*a0 = a1 = \delta*a_0$, we can take $d=a0$
	    tau_correction = a0;
	    assert((i_tab.matrix(new_tw)-1)*tau_correction
		   ==(E.ctxt.delta()-1)*a0);
	  }
	}

	const Weight new_lambda_rho =
	  E.lambda_rho - rho_r_shift + alpha*(level/2);
	assert((E.ctxt.gamma()-new_lambda_rho).dot(alpha_v)
	       ==rd.colevel(n_alpha)); // check that |level_a| did its work

	const int t_alpha = E.t.dot(alpha);
	if (type1)
	{ // now distinguish 1r1f and 1r1s
	  if (t_alpha%2!=0) // no effect of |alpha_simple|, unlike 1i2 cases
	    return one_real_pair_switched;
	  result = one_real_pair_fixed; // what remains is case 1r1f

	  E0.t -= alpha_v*(t_alpha/2);
	  assert(same_sign(E,E0)); // since only |t| changes

	  param F0(E.ctxt,new_tw,
		   new_lambda_rho, E.tau + tau_correction,
		   E.l, E0.t, flipped1);
	  param F1(E.ctxt,new_tw,
		   new_lambda_rho, F0.tau, E.l + alpha_v, E0.t,
		   flipped1);

	  F0.lambda_rho+=rho_r_shift;
	  F1.lambda_rho+=rho_r_shift;
	  z_align(E0,F0);
	  z_align(E0,F1);
	  F0.lambda_rho-=rho_r_shift;
	  F1.lambda_rho-=rho_r_shift;
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley

	} // end of 1r1 case
	else
	{ // type 1r2
	  result = one_real_single;
	  Coweight diff = // called $s$ in table 2 of [Ptr]
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v);

	  E0.t -= diff*t_alpha;
	  assert(same_sign(E,E0)); // since only |t| changes

	  param E1 = E0; // for cross neighbour; share updated value of |t|
	  E1.lambda_rho += alpha;
	  assert(not same_standard_reps(E0,E1));

	  param F(E.ctxt,new_tw,
		  new_lambda_rho, E.tau + tau_correction,
		  E.l, E0.t, flipped1);
	  F.lambda_rho+=rho_r_shift;
	  z_align(E0,F);
	  //	  z_align(F,E1);
	  z_align(E0,E1);
	  F.lambda_rho-=rho_r_shift;
	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E1)); // cross link
	} // end of 1r2 case
      }
      else // length 1 complex case
      { result = rd.is_posroot(theta_alpha)
	  ? one_complex_ascent : one_complex_descent ;

	TwistedInvolution new_tw = E.tw;
	bool monoflip=false;
	unsigned int d = rd.is_posroot(theta_alpha) ?
	  i_tab.length(new_tw) - i_tab.length(E.tw) :
	  i_tab.length(E.tw) - i_tab.length(new_tw);
	if (d%2==0) monoflip = not monoflip;
	tW.twistedConjugate(subs.reflection(p.s0),new_tw);
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	//	assert(rd.is_simple_root(alpha_simple));
	// Haven't thought whether there's an issue
	const auto theta_q = i_tab.nr(new_tw);
	//	if (rd.is_posroot(theta_alpha)) assert (theta_q > theta);
	// Weight rho_r_shift = rd.is_posroot(theta_alpha) ?
	//   repr::Cayley_shift(ic,theta_q,ww) :
	//   repr::Cayley_shift(ic,theta,ww);

	const bool flipped = rd.is_posroot(theta_alpha) ?
	  Cayley_shift_flip(E.ctxt,theta_q,theta,ww) :
	  Cayley_shift_flip(E.ctxt,theta,theta_q,ww);
	if(flipped) std::cout << "1C flip" << std::endl;
	auto E1=complex_cross(p,E0);
	E1.flipped = flipped;
	//  E1.lambda_rho -= rho_r_shift; already done in complex_cross
	links.push_back(E1);
      }
    }
    break;

  case ext_gen::two:
    { const Weight& alpha = integr_datum.simpleRoot(p.s0);
      const Coweight& alpha_v = integr_datum.simpleCoroot(p.s0);
      RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
      const Weight& beta = integr_datum.simpleRoot(p.s1);
      const Coweight& beta_v = integr_datum.simpleCoroot(p.s1);
      RootNbr n_beta = subs.parent_nr_simple(p.s1);
      // RootNbr theta_beta = i_tab.root_involution(theta,n_beta);

      if (theta_alpha==n_alpha) // length 2 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	int tf_beta = (E.ctxt.g_rho_check() - E.l).dot(beta);
	assert((tf_alpha-tf_beta)%2==0); // same compactness
	if (tf_alpha%2!=0) // then $\alpha$ and $\beta$ are compact
	  return two_imaginary_compact;

	// noncompact case
	const TwistedInvolution new_tw =
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));
	bool monoflip=false;
	unsigned int d = i_tab.length(new_tw) - i_tab.length(E.tw);
	if (d%2!=0) monoflip = not monoflip;
	// make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const Weight rho_r_shift = repr::Cayley_shift(ic,theta_p,theta,ww);
	const bool flipped = Cayley_shift_flip(E.ctxt,theta_p,theta,ww);
	//	const Weight tau_shift = tau_shift(E.ctxt,theta,theta_p,ww);
	if(flipped) std::cout << "2i flip" << std::endl;
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 2

	int at = alpha_v.dot(E.tau); int bt = beta_v.dot(E.tau);
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1;

	if (matreduc::has_solution(th_1,alpha)) // then type 2i11
	{ result = two_imaginary_single_single;
	  const Weight sigma = matreduc::find_solution(th_1,alpha*at+beta*bt);

	  param F (E.ctxt, new_tw,
		   E.lambda_rho + rho_r_shift,
		   E.tau + sigma, E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2),
		   E.t, flipped);

	  E0.l += alpha_v+beta_v;
	  F.lambda_rho-=rho_r_shift;
	  z_align(E,F); // no 3rd arg, since |E.lambda_rho| unchanged
	  //	  z_align(F,E0);
	  F.lambda_rho-=rho_r_shift;
	  z_align(E,E0);
	  F.lambda_rho+=rho_r_shift;
	  links.push_back(std::move(F));  // Cayley link
	  links.push_back(std::move(E0)); // cross link
	}
	else if (matreduc::has_solution(th_1,alpha+beta)) // case 2i12
	{
	  if ((at+bt)%2!=0)
	    return two_imaginary_single_double_switched; // 2i12s
	  result = two_imaginary_single_double_fixed; // 2i12f
	  const int m =  unsigned(at)%2; // safe modular reduction
          const int mm=1-m;

	  // one of the $\tau$ requires upstairs solution for an odd-odd pair:
	  const Weight sigma =
	    matreduc::find_solution(th_1,alpha*(at+mm)+beta*(bt-mm));

	  const Weight new_tau0 = E.tau - alpha*((at+m)/2) - beta*((bt-m)/2);
	  //  + tau_shift;
          const Coweight new_l = E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2);

	  // first Cayley link |F0| will be the one that does not need |sigma|
	  param F0(E.ctxt, new_tw,
		   E.lambda_rho + rho_r_shift
		   + alpha*m, new_tau0, new_l, E.t, flipped);
	  param F1(E.ctxt, new_tw,
		   E.lambda_rho + rho_r_shift
		   + alpha*mm, E.tau + sigma, new_l, E.t, flipped);

	  //	  int t_alpha=E.t.dot(alpha);
	  // z_align(E,F0,m*t_alpha);
	  F0.lambda_rho-=rho_r_shift;
	  F1.lambda_rho-=rho_r_shift;
	  z_align(E,F0);
	  // z_align(E,F1,mm*t_alpha);
	  z_align(E,F1);
	  F0.lambda_rho+=rho_r_shift;
	  F1.lambda_rho+=rho_r_shift;
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley
	} // end of case 2i12f
	else
	{ // type 2i22
	  result = two_imaginary_double_double;
	  // $\alpha^\vee$ and $\beta^\vee$ are even on $(X^*)^\theta$ and
	  // $(1-\delta)\tau\in(X^*)^\theta+2X^*$ so $<av-bv,\tau>$ is even
	  assert((at-bt)%2==0);
	  int m = static_cast<unsigned int>(at)%2; // safe modular reduction
	  assert (m==0 or m==1);
	  param F0(E.ctxt, new_tw,
		   E.lambda_rho + rho_r_shift
		   + alpha*m,
		   E.tau - alpha*((at+m)/2) - beta*((bt-m)/2),
		   E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t,
		   flipped);
	  param F1(E.ctxt, new_tw,
		   E.lambda_rho + rho_r_shift
		   + alpha*(1-m) + beta,
		   E.tau - alpha*((at-m)/2) - beta*((bt+m)/2),
		   F0.l,E.t,
		   flipped);

	  //	  int ta = E.t.dot(alpha), tb=E.t.dot(beta);
	  // z_align(E,F0,ta*m);
	  F0.lambda_rho -= rho_r_shift;
	  F1.lambda_rho -= rho_r_shift;
	  z_align(E,F0);
	  // z_align(E,F1,ta*(1-m)+tb);
	  z_align(E,F1);
	  F0.lambda_rho += rho_r_shift;
	  F1.lambda_rho += rho_r_shift;
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley
	} // end type 2i22 case
      }

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 2 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));
	bool monoflip=false;
	unsigned int d = i_tab.length(E.tw) - i_tab.length(new_tw);
	if (d%2!=0) monoflip = not monoflip;
	const auto theta_p = i_tab.nr(new_tw); // downstairs
	const Weight rho_r_shift = repr::Cayley_shift(ic,theta,theta_p,ww);
	const bool flipped = Cayley_shift_flip(E.ctxt,theta,theta_p,ww);
	if(flipped) std::cout << "2r flip" << std::endl;
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return two_real_nonparity; // no link added here

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	WeightInvolution theta_1 = i_tab.matrix(theta)-1; // upstairs

	const Weight new_lambda_rho = E.lambda_rho - rho_r_shift
	  + alpha*(a_level/2) + beta*(b_level/2);

	int ta = E.t.dot(alpha); int tb = E.t.dot(beta);
	param E1=E; // another modifiable copy, like |E0|

	if (matreduc::has_solution(theta_1,alpha))
	{ // type 2r11
	  result = two_real_double_double;
	  // $\alpha$ is even on $(X_*)^{-\theta'}$ (so is $\beta$), and
	  // $t(1-\delta)\in(X_*)^{-\theta'}+2X_*$ so $<t,alpha-beta>$ is even
	  assert((ta-tb)%2==0);
	  int m =  static_cast<unsigned int>(ta)%2;

	  // set two values for |t|; actually the same value in case |m==0|
	  E0.t -= alpha_v*((ta+m)/2) + beta_v*((tb-m)/2);
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E0.t.dot(alpha)==-m and E0.t.dot(beta)==m);

	  E1.t -= alpha_v*((ta-m)/2) + beta_v*((tb+m)/2);
	  assert(same_sign(E,E1)); // since only |t| changes
	  assert(E1.t.dot(alpha)==m and E1.t.dot(beta)==-m);

	  param F0(E.ctxt, new_tw,
		   new_lambda_rho, E.tau, E.l+alpha_v*m, E0.t,
		   flipped);
	  param F1(E.ctxt, new_tw,
		   new_lambda_rho, E.tau,
		   E.l+alpha_v*(1-m)+beta_v,E1.t, flipped);
	  // z_align(E0,F0,m*((b_level-a_level)/2));
	  F0.lambda_rho -= rho_r_shift;
	  F1.lambda_rho -= rho_r_shift;
	  z_align(E0,F0);
	  //	  z_align(E1,F1,m*((a_level-b_level)/2));
	  // z_align(E0,F1,m*((a_level-b_level)/2));
	  z_align(E0,F1);
	  F0.lambda_rho -= rho_r_shift;
	  F1.lambda_rho -= rho_r_shift;
	  // Cayley links
	  links.push_back(std::move(F0));
	  links.push_back(std::move(F1));
	} // end 2r11 case
	else if (matreduc::has_solution(theta_1,alpha+beta))
	{ // type 2r21
	  if ((ta+tb)%2!=0)
	    return two_real_single_double_switched; // 2r21s
	  result = two_real_single_double_fixed; // 2r21f
	  const int m =  static_cast<unsigned int>(ta)%2;
	  const int mm=1-m;

	  // one of the $t$ requires downstairs solution for an odd-odd pair:
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v*(ta+mm)+beta_v*(tb-mm));

	  // E0 is parameter adapted to Cayley transform that does not need |s|
	  E0.t -= alpha_v*((ta+m)/2) + beta_v*((tb-m)/2);
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E0.t.dot(alpha)==-m and E0.t.dot(beta)==m);

	  E1.t -= s;
	  assert(same_sign(E,E1)); // since only |t| changes
	  assert(E1.t.dot(alpha)==-mm and E1.t.dot(beta)==mm);

	  // Cayley links
	  param F0(E.ctxt, new_tw,
		   new_lambda_rho, E.tau, E.l+alpha_v*m, E0.t,
		   flipped);
	  param F1(E.ctxt, new_tw,
		   new_lambda_rho, E.tau, E.l+alpha_v*mm, E1.t,
		   flipped);

	  // z_align(E0,F0,m *((b_level-a_level)/2));
	  F0.lambda_rho += rho_r_shift;
	  F1.lambda_rho += rho_r_shift;
	  z_align(E0,F0);
	  //	  z_align(E1,F1,mm*((b_level-a_level)/2));
	  // z_align(E0,F1,mm*((b_level-a_level)/2));
	  z_align(E0,F1);
	  F0.lambda_rho -= rho_r_shift;
	  F1.lambda_rho -= rho_r_shift;
	  links.push_back(std::move(F0));
	  links.push_back(std::move(F1));
	} // end of case 2r21f
	else // case 2r22
	{ result = two_real_single_single;
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v*ta+beta_v*tb);

	  E0.t -= s; // parameter adapted to Cayley transform |F|
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E.t.dot(alpha)==0 and E.t.dot(beta)==0);

	  E1.lambda_rho += alpha+beta;
	  E1.t = E0.t; // cross action, keeps adaptation of |t| to |F| below
	  assert(not same_standard_reps(E0,E1));

	  param F(E.ctxt, new_tw, new_lambda_rho, E.tau,
		  E.l, E0.t, flipped);
	  F.lambda_rho += rho_r_shift;
	  z_align(E0,F); // no 3rd arg, as |E.t.dot(alpha)==0| etc.
	  //	  z_align(F,E1);
	  F.lambda_rho -= rho_r_shift;
	  z_align(E0,E1);
	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E1)); // cross link
	} // end of case 2r22
      }
      else // length 2 complex case
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha != (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // twisted non-commutation with |s0.s1|
	  result = ascent ? two_complex_ascent : two_complex_descent;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw);
	  tW.twistedConjugate(subs.reflection(p.s1),new_tw);
	  bool monoflip=false;
	  unsigned int d = ascent ?
	  i_tab.length(new_tw) - i_tab.length(E.tw) :
	  i_tab.length(E.tw) - i_tab.length(new_tw);
	  if (d%2!=0) monoflip = not monoflip;

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here
	  const auto theta_q = i_tab.nr(new_tw);
	  Weight rho_r_shift = ascent ?
	    repr::Cayley_shift(ic,theta_q,theta,ww) :
	    repr::Cayley_shift(ic,theta,theta_q,ww);
	  const bool flipped = ascent ?
	    Cayley_shift_flip(E.ctxt,theta_q,theta,ww) :
	    Cayley_shift_flip(E.ctxt,theta,theta_q,ww);
	  if(flipped) std::cout << "2C flip" << std::endl;
	auto E1=complex_cross(p,E0);
	E1.flipped = flipped;
	E1.lambda_rho += rho_r_shift;
	z_align(E0,E1);
	E1.lambda_rho -= rho_r_shift;
	links.push_back(E1);
	}
	else if (ascent)
	{ // twisted commutation with |s0.s1|: 2Ci
	  result = two_semi_imaginary;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|
	  bool monoflip=false;
	  unsigned int d = i_tab.length(new_tw) - i_tab.length(E.tw);
	  if (d%2==0) monoflip = not monoflip;
	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p = i_tab.nr(new_tw); // upstairs
	  const Weight rho_r_shift = repr::Cayley_shift(ic,theta_p,theta,ww);
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  const bool flipped = Cayley_shift_flip(E.ctxt,theta_p,theta,ww);
	  if(flipped) std::cout << "2Ci flip" << std::endl;

	  // downstairs cross by |ww| only has imaginary and complex steps, so
	  // $\alpha_v.(\gamma-\lambda_\rho)$ is unchanged across |ww|
	  const int f = // number of times $\alpha$ is added to $\lambda_\rho$
	    (E.ctxt.gamma() - E.lambda_rho).dot(alpha_v)- rd.colevel(n_alpha);

	  const Weight new_lambda_rho = E.lambda_rho + alpha*f + rho_r_shift;
	  // both $\gamma-\lambda$ and $\tau$ get $f*alpha$ subtracted by
	  // $\alpha$-reflection; adapt $\tau$ for vanishing $1-\delta$ image
	  const Weight new_tau = rd.reflection(n_alpha,E.tau) + alpha*f;

	  // but |dual_v| needs correction by |ell_shift|
	  const int dual_f = (E.ctxt.g_rho_check() - E.l).dot(alpha);

	  const Coweight new_l = E.l + alpha_v*dual_f;
          const Coweight new_t =
	    rd.coreflection(E.t,n_alpha) - alpha_v*dual_f;
	  param F (E.ctxt, new_tw, new_lambda_rho,
		   new_tau, new_l, new_t, flipped);
	  F.lambda_rho-=rho_r_shift;
	  int ab_tau = (alpha_v+beta_v).dot(E.tau);
	  assert (ab_tau%2==0);
	  // F.flip((F.flip)&((ab_tau*dual_f)%4!=0));
	  F.flip((ab_tau*dual_f)%4!=0);
	  F.lambda_rho+=rho_r_shift;
	  links.push_back(std::move(F));  // "Cayley" link
	}
	else // twisted commutation with |s0.s1|, and not |ascent|: 2Cr
	{ result = two_semi_real;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|
	  const auto theta_p = i_tab.nr(new_tw); // downstairs
	  bool monoflip=false;
	  unsigned int d = i_tab.length(E.tw) - i_tab.length(new_tw);
	  if (d%2==0) monoflip = not monoflip;

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here
	  const Weight rho_r_shift = repr::Cayley_shift(ic,theta,theta_p,ww);
	  const bool flipped = Cayley_shift_flip(E.ctxt,theta,theta_p,ww);
	  if(flipped) std::cout << "2Cr flip" << std::endl;
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  const int f = level_a(E,rho_r_shift,n_alpha);

	  const Weight new_lambda_rho = // \emph{reflect} parallel to alpha
	    E.lambda_rho - rho_r_shift + alpha*f;
	  const Weight new_tau = rd.reflection(n_alpha,E.tau) - alpha*f;

	  const int dual_f = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	  const Coweight new_l = E.l + alpha_v*dual_f;
          const Coweight new_t =
	    rd.coreflection(E.t,n_alpha) + alpha_v*dual_f;

	  param F (E.ctxt, new_tw, new_lambda_rho, new_tau, new_l,
		   new_t,  flipped);

	  F.lambda_rho+=rho_r_shift;
	  int t_ab = E.t.dot(beta-alpha);
	  assert(t_ab%2==0);
	  // F.flip((F.flip)&((t_ab * (f+alpha_v.dot(E.tau)))%4!=0));
	  F.flip((t_ab * (f+alpha_v.dot(E.tau)))%4!=0);
	  F.lambda_rho-=rho_r_shift;
	  links.push_back(std::move(F));  // "Cayley" link
	}
      }
    }
    break;
  case ext_gen::three:
    { const Weight& alpha = integr_datum.simpleRoot(p.s0);
      const Coweight& alpha_v = integr_datum.simpleCoroot(p.s0);
      RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
      const Weight& beta = integr_datum.simpleRoot(p.s1);
      RootNbr n_beta = subs.parent_nr_simple(p.s1);

      RootNbr n_kappa =integr_datum.simple_reflected_root
	 (p.s1, integr_datum.simpleRootNbr(p.s0));
      WeylWord s_kappa = subs.reflection(integr_datum.posRootIndex(n_kappa));

      const Weight& kappa = integr_datum.root(n_kappa);
      const Coweight& kappa_v = integr_datum.coroot(n_kappa);

      const TwistedInvolution new_tw = tW.prod(s_kappa,E.tw); // when applicable
      bool monoflip=false;
      int e = i_tab.length(new_tw) - i_tab.length(E.tw);
      unsigned int d = (e>0) ? e : -e;
      if (d%2==0) monoflip = not monoflip;

      if (theta_alpha==n_alpha) // length 3 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	int tf_beta = (E.ctxt.g_rho_check() - E.l).dot(beta);
	assert((tf_alpha-tf_beta)%2==0); // same compactness
	if (tf_alpha%2!=0) // then $\alpha$ and $\beta$ are compact
	  return three_imaginary_compact;

	// noncompact case
	result = three_imaginary_semi; //3i

	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const Weight rho_r_shift = repr::Cayley_shift(ic,theta_p,theta,ww);
	const bool flipped = Cayley_shift_flip(E.ctxt,theta_p,theta,ww);
	if(flipped) std::cout << "3Ci flip" << std::endl;
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 3

	param F(E.ctxt, new_tw,
		E.lambda_rho, // + rho_r_shift,
		E.tau - alpha*kappa_v.dot(E.tau),
		E.l + kappa_v*((tf_alpha+tf_beta)/2), E.t,
		flipped);
	z_align(E,F); // |lambda_rho| unchanged at simple Cayley
	F.lambda_rho+=rho_r_shift;
	links.push_back(std::move(F)); // Cayley link
      }
      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 3 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here
	const auto theta_p = i_tab.nr(new_tw); // downstairs

	const Weight rho_r_shift = repr::Cayley_shift(ic,theta,theta_p,ww);
	const bool flipped = Cayley_shift_flip(E.ctxt,theta,theta_p,ww);
	if(flipped) std::cout << "3Cr flip" << std::endl;
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return three_real_nonparity; // no link added here

	// parity case
	result = three_real_semi; // 3r

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	const Weight new_lambda_rho = // make level for |kappa| zero
	  E.lambda_rho-rho_r_shift + kappa*((a_level+b_level)/2);

	E0.t -= alpha_v*kappa.dot(E.t); // makes |E.t.dot(kappa)==0|
	assert(same_sign(E,E0)); // since only |t| changes

	param F(E.ctxt, new_tw,	new_lambda_rho + rho_r_shift,
		E.tau,E.l,E0.t, flipped);

	z_align(E0,F); // no 3rd arg since |E.t.dot(kappa)==0|
	F.lambda_rho-=rho_r_shift;
	links.push_back(std::move(F)); // Cayley link
      }
      else // length 3 complex case 3Ci or 3Cr or 3C
      { const bool ascent = rd.is_posroot(theta_alpha);
	const RootNbr n_beta = subs.parent_nr_simple(p.s1);
	if (theta_alpha == (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // reflection by |alpha+beta| twisted commutes with |E.tw|: 3Ci or 3Cr
	  result = ascent ? three_semi_imaginary : three_semi_real;

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_upstairs = ascent ? i_tab.nr(new_tw) : theta;
	  const auto theta_downstairs = ascent ? theta : i_tab.nr(new_tw);
	  Weight rho_r_shift =
	    repr::Cayley_shift(ic,theta_upstairs,theta_downstairs,ww);
	  rho_r_shift = ascent ? rho_r_shift : -rho_r_shift;
	  const bool flipped =
	    Cayley_shift_flip(E.ctxt,theta_upstairs,theta_downstairs,ww);
	  if(flipped) std::cout << "3Ci flip" << std::endl;
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	  int dtf_alpha = (E.ctxt.gamma() - E.lambda_rho).dot(alpha_v)
	    - rd.colevel(n_alpha);
	  Weight new_lambda_rho = E.lambda_rho; // + rho_r_shift; // for now

	  if (ascent) // 3Ci
	  { param F(E.ctxt,new_tw,
		    dtf_alpha%2==0 ? new_lambda_rho : new_lambda_rho + kappa,
		    E.tau - kappa*(kappa_v.dot(E.tau)/2),
		    E.l + kappa_v*tf_alpha, E.t,
		    flipped);

	    assert(E.t.dot(kappa)==0);
	    // since it is half of |t*(1+theta)*kappa=l*(delta-1)*kappa==0|
	    z_align(E,F); // may ignore possible shift by |kappa|
	    F.lambda_rho+=rho_r_shift;
	    links.push_back(std::move(F)); // Cayley link
	  }
	  else // descent, so 3Cr
	  { // make |E.t.dot(kappa)==0| using |kappa_v|
	    E0.t -= kappa_v*(kappa.dot(E0.t)/2);
	    assert(same_sign(E,E0)); // since only |t| changes

	    param F(E.ctxt, new_tw,
		    new_lambda_rho + kappa*dtf_alpha, E.tau,
		    tf_alpha%2==0 ? E.l : E.l+kappa_v, E0.t,
		    flipped);

	    z_align(E0,F); // no 3rd arg since |E.t.dot(kappa)==0|
	    F.lambda_rho+=rho_r_shift;
	    links.push_back(std::move(F)); // Cayley link
	  }

	}
	else // twisted non-commutation: 3C+ or 3C-
	{
	  result = ascent ? three_complex_ascent : three_complex_descent;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw);
	  tW.twistedConjugate(subs.reflection(p.s1),new_tw);
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw);
	  bool monoflip=false;
	  unsigned int d = rd.is_posroot(theta_alpha) ?
	  i_tab.length(new_tw) - i_tab.length(E.tw) :
	  i_tab.length(E.tw) - i_tab.length(new_tw);
	  if (d%2==0) monoflip = not monoflip;

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here
	  const auto theta_q = i_tab.nr(new_tw);

	  //	  Weight rho_r_shift = ascent ?
	  //	    repr::Cayley_shift(ic,theta_q,ww) :
	  //	    repr::Cayley_shift(ic,theta,ww);
	  const bool flipped = ascent ?
	    Cayley_shift_flip(E.ctxt,theta_q,theta,ww) :
	    Cayley_shift_flip(E.ctxt,theta,theta_q,ww);
	  auto E1=complex_cross(p,E0);
	  E1.flipped = flipped;
	  //	  E1.lambda_rho -= rho_r_shift;
	  links.push_back(E1);
	}
      }
    }
    break;
  }

  // October surprise: add a flip to links with a length difference of 2
  if (p.length()-(has_defect(result)?1:0)==2)
  { auto it=links.begin(); auto c=scent_count(result);
    for (unsigned i=0; i<c; ++i,++it) // only affect ascent/descent links
      it->flip(); // do the flip
  }

  return result;
} // |star|

bool is_descent (const ext_gen& kappa, const param& E)
{ // like in |star|, generators are not simple for |E.rc().twistedWeylGroup()|
  const TwistedWeylGroup& tW = E.rc().twistedWeylGroup();
  const InnerClass& ic = E.rc().innerClass();
  const InvolutionTable& i_tab = ic.involution_table();
  const InvolutionNbr theta = i_tab.nr(E.tw); // so use root action of |E.tw|
  const RootDatum& rd = E.rc().rootDatum();
  const SubSystem& subs = E.ctxt.subsys();
  const RootNbr n_alpha = E.ctxt.subsys().parent_nr_simple(kappa.s0);
  const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
  const Weight& alpha = E.ctxt.id().simpleRoot(kappa.s0);

  // we don't need to inspect |kappa.type|, it does not affect descent status
  if (theta_alpha==n_alpha) // imaginary case, return whether compact
    return (E.ctxt.g_rho_check()-E.l).dot(alpha) %2!=0;
  TwistedInvolution new_tw = tW.prod(subs.reflection(kappa.s0),E.tw);
  if (theta_alpha==rd.rootMinus(n_alpha)) // real, return whether parity
  {
    if (not (kappa.type == ext_gen::one))
      new_tw = tW.prod(subs.reflection(kappa.s1),new_tw);
    if (kappa.type == ext_gen::three)
      new_tw = tW.prod(subs.reflection(kappa.s0),new_tw);
    InvolutionNbr theta_downstairs = i_tab.nr(new_tw);
    RootNbr alpha_simple = n_alpha; // copy to be made simple
    const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
    const auto level = level_a(E,repr::Cayley_shift(ic,theta,theta_downstairs,
						    ww),n_alpha);
    if (rd.is_simple_root(alpha_simple)) // |fixed_conjugate_simple| succeeded
      return level%2==0; // parity if |level| is even
    else
    {
      RootNbr alpha_0 = // one of |alpha==alpha_0+alpha_1|
	rd.permuted_root(rd.simpleRootNbr(rd.find_descent(alpha_simple)),ww);

      if (i_tab.root_involution(theta,alpha_0)==rd.rootMinus(alpha_0))
	// when |alpha_0| is real for theta, the parity condition is flipped:
	return level%2!=0; // parity if |level| is odd
      return level%2==0; // parity if |level| is even
    }
  }
  else // complex
    return rd.is_negroot(theta_alpha);
}

weyl::Generator first_descent_among
  (RankFlags singular_orbits, const ext_gens& orbits, const param& E)
{ for (auto it=singular_orbits.begin(); it(); ++it)
    if (is_descent(orbits[*it],E))
      return *it;
  return orbits.size(); // no singular descents found
}

ext_block::ext_block // for external twist; old style blocks
  (const InnerClass& G,
   const Block& block,
   const KGB& kgb, const KGB& dual_kgb, // all are needed
   const WeightInvolution& delta)
  : parent(block)
  , orbits(fold_orbits(G.rootDatum(),delta))
  , folded(orbits,block.Dynkin())
  , d_delta(delta)
  , info()
  , data(orbits.size()) // create that many empty vectors
  , l_start(parent.length(parent.size()-1)+2,0)
{
  BitMap fixed_points(block.size());

 // compute |child_nr| and |parent_nr| tables
  weyl::Twist twist(orbits);

  if (kgb.twisted(0,delta)==UndefKGB or
      dual_kgb.twisted(0,delta.transposed())==UndefKGB)
    return; // if one or other not delta-stable, leave |size==0| and quit

  for (BlockElt z=0; z<block.size(); ++z)
    if (twisted(block,kgb,dual_kgb,z,delta)==z)
      fixed_points.insert(z);

  complete_construction(fixed_points);
  // FIXME cannot call |check| here, although setting sign flips depends on it

} // |ext_block::ext_block|

ext_block::ext_block // for an external twist
  (const InnerClass& G,
   const param_block& block, const WeightInvolution& delta,
   bool verbose)
  : parent(block)
  , orbits(block.fold_orbits(delta))
  , folded(orbits,block.Dynkin())
  , d_delta(delta)
  , info()
  , data(orbits.size()) // create that many empty vectors
  , l_start(parent.length(parent.size()-1)+2,0)
{
  BitMap fixed_points(block.size());

  // compute the delta-fixed points of the block

  // test if twisting some block element lands in the same block
  if (twisted(block,0,delta)==UndefBlock)
    return; // if block not delta-stable, leave |size==0| and quit

  for (BlockElt z=0; z<block.size(); ++z)
    if (twisted(block,z,delta)==z)
      fixed_points.insert(z);

  complete_construction(fixed_points);
  if (not check(block,verbose)) // this sets the edge signs, not just a check!
    throw std::runtime_error("Failure detected in extended block construction");

} // |ext_block::ext_block|, from a |param_block|

void ext_block::complete_construction(const BitMap& fixed_points)
{
  unsigned int folded_rank = orbits.size();
  std::vector<BlockElt> child_nr(parent.size(),UndefBlock);
  std::vector<BlockElt> parent_nr(fixed_points.size());
  { BlockElt x=0; unsigned cur_len=0;
    for (auto it=fixed_points.begin(); it(); ++it,++x)
    {
      parent_nr[x]=*it;
      child_nr[*it]=x;
      while (cur_len<parent.length(*it)) // for new length level(s) reached
	l_start[++cur_len]=x; // mark |x| as first of length at least |cur_len|
    }
    while (++cur_len<l_start.size())
      l_start[cur_len]=fixed_points.size(); // allow |l_start[length(...)+1]|
  }

  info.reserve(parent_nr.size());  // reserve size of (smaller) extended block
  for (weyl::Generator s=0; s<folded_rank; ++s)
    data[s].reserve(parent_nr.size()); // same for each |data[s]|.

  for (BlockElt n=0; n<parent_nr.size(); ++n)
  {
    BlockElt z=parent_nr[n];
    info.push_back(elt_info(z));
    for (weyl::Generator oi=0; oi<orbits.size(); ++oi) // |oi|: orbit index
    {
      const weyl::Generator s = orbits[oi].s0, t=orbits[oi].s1;
      BlockElt link, second = UndefBlock; // these index parent block elements
      DescValue type = extended_type(parent,z,orbits[oi],link,fixed_points);
      data[oi].push_back(block_fields(type)); // create entry

      if (link==UndefBlock)
	continue; // done with |s| for imaginary compact, real nonparity cases

      // now maybe set |second|, depending on case
      switch (type)
      {
      default: break;

      case one_imaginary_single:
      case one_real_single: // in these cases: parent cross neighbour for |s|
	second = parent.cross(s,z);
	break;

      case one_real_pair_fixed:
      case one_imaginary_pair_fixed: // in these cases get second Cayley image
	second = parent.cross(s,link);
	break;

      case two_imaginary_single_single:
      case two_real_single_single: // here: double cross neighbour for |s|
	second = parent.cross(t,parent.cross(s,z));
	assert(second==parent.cross(s,parent.cross(t,z)));
	break;

      case two_imaginary_single_double_fixed:
      case two_real_single_double_fixed: // find second Cayley image, which is
	second = parent.cross(s,link); // parent cross link
	assert(second==parent.cross(t,link)); // (for either generator)
	if (link>second) // to make sure ordering is same for a twin pair
	  std::swap(link,second); // we order both by block number (for now)
	break;

      case two_imaginary_double_double:
      case two_real_double_double: // find second Cayley image
	second = parent.cross(t,parent.cross(s,link));
	assert(second==parent.cross(s,parent.cross(t,link)));
	break;
      } // |switch(type)|

      // enter translations of |link| and |second| to child block numbering
      BlockEltPair& dest = data[oi].back().links;
      dest.first=child_nr[link];
      if (second!=UndefBlock)
	dest.second = child_nr[second];
    }
  } // |for(n)|
} // |ext_block::complete_construction|

// we compute $\max\{l\mid l_start[l]\leq n\}$, i.e. |upper_bound(,,n)-1|
unsigned ext_block::length(BlockElt n) const
{
  unsigned min=0, max=l_start.size()-1; // the answer will lie in $[min,max)$
  while (max>min+1) // body strictly reduces |max-min| in all cases
  {
    unsigned l=(min+max)/2;
    if (l_start[l]>n)
      max=l; // preserves invariant |l_start[max]>n|
    else
      min=l; // preserves invariant |l_start[min]<=n|
  }
  assert(min+1==max);
  return min;
}

BlockElt ext_block::cross(weyl::Generator s, BlockElt n) const
{
  switch (descent_type(s,n))
  {
  case one_complex_ascent:
  case one_complex_descent:
  case two_complex_ascent:
  case two_complex_descent:
  case three_complex_ascent:
  case three_complex_descent:
    return data[s][n].links.first;

    // zero valued Cayleys have trivial cross actions
  case one_real_nonparity: case one_imaginary_compact:
  case one_imaginary_pair_switched: case one_real_pair_switched:
  case two_real_nonparity: case two_imaginary_compact:
  case two_imaginary_single_double_switched:
  case two_real_single_double_switched:
  case three_real_nonparity: case three_imaginary_compact:

    // double valued Cayleys also have trivial cross actions
  case one_real_pair_fixed: case one_imaginary_pair_fixed:
  case two_real_double_double: case two_imaginary_double_double:

    // cases with back-and-forth cross actions
  case two_semi_imaginary: case two_semi_real:
  case two_imaginary_single_double_fixed: case two_real_single_double_fixed:
  case three_semi_imaginary: case three_real_semi:
  case three_imaginary_semi: case three_semi_real:
    return n;

    // some single valued extended Cayleys use second link for cross action
  case one_imaginary_single:
  case one_real_single:
  case two_imaginary_single_single:
  case two_real_single_single:
    return data[s][n].links.second;

  }
  assert(false); return UndefBlock; // keep compiler happy
} // |ext_block::cross|

BlockElt ext_block::some_scent(weyl::Generator s, BlockElt n) const
{
  const BlockElt c = data[s][n].links.first;
  assert(c!=UndefBlock);
  return c;
}

BlockElt ext_block::Cayley(weyl::Generator s, BlockElt n) const
{
  return  is_complex(descent_type(s,n)) ? UndefBlock : data[s][n].links.first;
}

BlockEltPair ext_block::Cayleys(weyl::Generator s, BlockElt n) const
{
  assert(has_double_image(descent_type(s,n)));
  return data[s][n].links;
}


// check validity, by comparing with results found using extended parameters
// the signs are recorded in |eb|, and printed to |cout| if |verbose| holds.
bool ext_block::check(const param_block& block, bool verbose)
{
  context ctxt (block.context(),delta(),block.gamma());
  containers::sl_list<param> links;
  for (BlockElt n=0; n<size(); ++n)
  { auto z=this->z(n);
    for (weyl::Generator s=0; s<rank(); ++s)
    { param E(ctxt,block.x(z),block.lambda_rho(z)); // re-init each iteration
      ext_gen p=orbit(s); links.clear(); // output arguments for |star|
      auto tp = star(E,p,links);
      if (tp!=descent_type(s,n))
	return false;

      auto it = links.begin();

      switch (tp)
      {
      case one_imaginary_pair_switched: case one_real_pair_switched:
      case one_real_nonparity: case one_imaginary_compact:
      case two_imaginary_single_double_switched:
      case two_real_single_double_switched:
      case two_real_nonparity: case two_imaginary_compact:
      case three_real_nonparity: case three_imaginary_compact:
	assert(links.empty()); break;
      case one_complex_ascent: case one_complex_descent:
      case two_complex_ascent: case two_complex_descent:
      case three_complex_ascent: case three_complex_descent:
	{ assert(links.size()==1);
	  BlockElt m=cross(s,n); // cross neighbour as bare element of |*this|
	  BlockElt cz = this->z(m); // corresponding element of (parent) |block|
	  param F(ctxt,block.x(cz),block.lambda_rho(cz)); // default extension
	  assert(same_standard_reps(*it,F)); // must lie over same
	  if (not same_sign(*it,F))
	  {
	    flip_edge(s,n,m);
	    if (verbose)
	      std::cout << "Flip at cross link " << unsigned{s}
                        << " from " << z << " to " << cz << '.' << std::endl;
	  }
	} break;
      case one_imaginary_single: case one_real_single:
      case two_imaginary_single_single: case two_real_single_single:
	{ assert(links.size()==2);
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  BlockElt Cz = this->z(m); // corresponding element of block
	  param F(ctxt,block.x(Cz),block.lambda_rho(Cz));
	  assert(same_standard_reps(*it,F));
	  if (not same_sign(*it,F))
	  {
	    flip_edge(s,n,m);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
	                << " from " << z << " to " << Cz << '.' << std::endl;
	  }
	  ++it;
	  m=cross(s,n); BlockElt cz = this->z(m);
	  param Fc(ctxt,block.x(cz),block.lambda_rho(cz));
	  assert(same_standard_reps(*it,Fc));
	  if (not same_sign(*it,Fc))
	  {
	    flip_edge(s,n,m);
	    if (verbose)
	      std::cout << "Flip at cross link " << unsigned{s}
	                << " from " << z << " to " << cz << '.' << std::endl;
	  }
	} break;
      case two_semi_imaginary: case two_semi_real:
      case three_semi_imaginary: case three_real_semi:
      case three_imaginary_semi: case three_semi_real:
	{ assert(links.size()==1);
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  BlockElt Cz = this->z(m); // corresponding element of block
	  param F(ctxt,block.x(Cz),block.lambda_rho(Cz));
	  assert(same_standard_reps(*it,F));
	  if (not same_sign(*it,F))
	  {
	    flip_edge(s,n,m);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
		      << " from " << z << " to " << Cz << '.' << std::endl;
	  }
	} break;
      case one_imaginary_pair_fixed: case one_real_pair_fixed:
      case two_imaginary_double_double: case two_real_double_double:
	{ assert(links.size()==2);
	  BlockEltPair m=Cayleys(s,n);
	  BlockElt Cz0 = this->z(m.first); BlockElt Cz1= this->z(m.second);
	  param F0(ctxt,block.x(Cz0),block.lambda_rho(Cz0));
	  param F1(ctxt,block.x(Cz1),block.lambda_rho(Cz1));
	  bool straight=same_standard_reps(*it,F0);
          const auto& node0 = straight ? *it : *std::next(it);
          const auto& node1 = straight ? *std::next(it) : *it;
	  if (not straight)
	    assert(same_standard_reps(node0,F0));
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node0,F0))
	  {
	    flip_edge(s,n,m.first);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
			<< " from " << z << " to " << Cz0 << '.' << std::endl;
	  }
	  if (not same_sign(node1,F1))
	  {
	    flip_edge(s,n,m.second);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
			<< " from " << z << " to " << Cz1 << '.' << std::endl;
	  }
	} break;
      case two_imaginary_single_double_fixed: case two_real_single_double_fixed:
	{ assert(links.size()==2);
	  BlockEltPair m=Cayleys(s,n);
	  BlockElt Cz0 = this->z(m.first); BlockElt Cz1= this->z(m.second);
	  param F0(ctxt,block.x(Cz0),block.lambda_rho(Cz0));
	  param F1(ctxt,block.x(Cz1),block.lambda_rho(Cz1));
	  bool straight=same_standard_reps(*it,F0);
          const auto& node0 = straight ? *it : *std::next(it);
          const auto& node1 = straight ? *std::next(it) : *it;
	  if (not straight)
	    assert(same_standard_reps(node0,F0));
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node0,F0))
	  {
	    flip_edge(s,n,m.first);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
			<< " from " << z << " to " << Cz0 << '.' << std::endl;
	  }
	  if (not same_sign(node1,F1))
	  {
	    flip_edge(s,n,m.second);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
			<< " from " << z << " to " << Cz1 << '.' << std::endl;
	  }
	} break;
      } // |switch(tp)|
    } // |for(s)|
  } // |for(n)|
  return true; // report success if we get here
} // |check|

RankFlags reduce_to(const ext_gens orbits, RankFlags gen_set)
{ RankFlags result;
  for (weyl::Generator s=0; s<orbits.size(); ++s)
    result.set(s,gen_set[orbits[s].s0]);
  return result;
}


RankFlags ext_block::singular_orbits (const param_block& parent) const
{ return reduce_to(orbits,parent.singular_simple_roots()); }

weyl::Generator
ext_block::first_descent_among(RankFlags singular_orbits, BlockElt y) const
{ auto it=singular_orbits.begin();
  for (; it(); ++it)
    if (is_descent(descent_type(*it,y)))
      return *it;

  return rank();
}

// reduce matrix to rows for extended block elements without singular descents
// the other rows are not removed, but the result lists the rows to retain
template<typename C> // matrix coefficient type (signed)
containers::simple_list<BlockElt> // returns list of elements selected
  ext_block::condense(matrix::Matrix<C>& M, const param_block& parent) const
{
  RankFlags sing_orbs = singular_orbits(parent);
  containers::simple_list<BlockElt> result;

  for (BlockElt y=M.numRows(); y-->0; ) // reverse loop is essential here
  { auto s = first_descent_among(sing_orbs,y);
    if (s==rank())
      result.push_front(y); // no singular descents, so a survivor
    else // a singular descent found, so not a survivor
    { // we contribute row |y| to all its descents by |s| with sign |-1|
      // then conceptually we clear row |y|, but don't bother: it gets ignored
      auto type=descent_type(s,y);
      if (is_like_compact(type))
	continue; // no descents, |y| represents zero; nothing to do for |y|

      // length difference and October surprise combine to always give -1
      const C c(-1); // surprise!
      if (has_double_image(type)) // 1r1f, 2r11
      { auto pair = Cayleys(s,y);
	M.rowOperation(pair.first,y,c*epsilon(s,pair.first,y));
	M.rowOperation(pair.second,y,c*epsilon(s,pair.second,y));
      }
      else
      { auto x = some_scent(s,y);
	M.rowOperation(x,y,c*epsilon(s,x,y));
      }
    }
  }
  return result;
} // |ext_block::condense|

// coefficient of neighbour |sx| of |x| in the action $(T_s+1)*a_x$
Pol ext_block::T_coef(weyl::Generator s, BlockElt sx, BlockElt x) const
{
  DescValue v = descent_type(s,x);
  if (not is_descent(v))
  {
    if (x==sx) // diagonal coefficient
      if (has_defect(v))
	return Pol(1,1)+Pol(1); // $q+1$
      else  // $0$, $1$, or $2$
	return Pol(is_like_nonparity(v) ? 0 : has_double_image(v) ? 2 : 1);
    else if (is_like_type_1(v) and sx==cross(s,x)) // type 1 imaginary cross
    {
      BlockElt y = Cayley(s,x); // pass via this element for signs
      int sign = epsilon(s,x,y)*epsilon(s,sx,y); // combine two Cayley signs
      return Pol(sign);
    }
    else // below diagonal coefficient
    {
      assert (data[s][x].links.first==sx or data[s][x].links.second==sx);
      int sign = epsilon(s,x,sx);
      if (has_defect(v)) // $\pm(q+1)$
	return Pol(1,sign)+Pol(sign);
      else
	return Pol(sign); // $\pm1$
    }
  } // |if (not is_descent(v))|

  int k = orbit(s).length();
  Pol result(k,1); // start with $q^k$
  if (x==sx) // diagonal coefficient
  {
    if (has_double_image(v))  // diagonal coefficient
      result[0] = -1; // $q^k-1$
    else if (is_like_compact(v))
      result[0] = 1; // $q^k+1$
    else if (has_defect(v))
      result[1] = -1; // $q^k-q$
    // |else| leave $q^k$
  }
  else if (is_like_type_2(v) and sx==cross(s,x)) // type 2 real cross
  {
    BlockElt y = Cayley(s,x); // pass via Cayley descent for signs
    int sign = epsilon(s,y,x)*epsilon(s,y,sx); // combine two Cayley signs
    return Pol(-sign); // forget term $q^k$, return $\mp 1$ instead
  }
  else // remaining cases involve descending edge (above-diagonal coefficient)
  {
    assert (data[s][x].links.first==sx or data[s][x].links.second==sx);
    if (not is_complex(v))
      result[has_defect(v) ? 1 : 0] = -1; // change into $q^k-1$ or $q^k-q$
    result *= epsilon(s,x,sx); // flip sign according to chosen edge
  }

  return result;
} // |ext_block::T_coef|

void set(matrix::Matrix<Pol>& dst, unsigned i, unsigned j, Pol P)
{
  // P.print(std::cerr << "M[" << i << ',' << j << "] =","q") << std::endl;
  dst(i,j)=P;
}

void show_mat(std::ostream& strm,const matrix::Matrix<Pol> M,unsigned inx)
{
  strm << "T_" << inx+1 << " [";
  for (unsigned i=0; i<M.numRows(); ++i)
    {        strm << " " << std::endl;
    for (unsigned j=0; j<M.numColumns(); ++j)
      //      if (not M(i,j).isZero())
      //	M(i,j).print(strm << i << ',' << j << ':',"q")  << ';';
      M(i,j).print(strm  << ' ',"q");
    //  strm << " " << std::endl;
    }
}

bool check_braid
  (const ext_block& b, weyl::Generator s, weyl::Generator t, BlockElt x,
   BitMap& cluster)
{
  if (s==t)
    return true;
  static const unsigned int cox_entry[] = {2, 3, 4, 6};
  unsigned int len = cox_entry[b.Dynkin().edge_multiplicity(s,t)];

  BitMap to_do(b.size()),used(b.size());
  to_do.insert(x);
  for (unsigned int i=0; i<len; ++i) // repeat |len| times, |i| is not used
    for (BitMap::iterator it=to_do.begin(); it(); ++it)
    {
      used.insert(*it);
      to_do.remove(*it);
      BlockEltList l; l.reserve(4); // for neighbours of |*it| by |s| and |t|
      b.add_neighbours(l,s,*it);
      b.add_neighbours(l,t,*it);
      for (unsigned j=0; j<l.size(); ++j)
	if (not used.isMember(l[j]))
	  to_do.insert(l[j]);
    }

  unsigned int n=used.size();
  matrix::Matrix<Pol> Ts(n,n,Pol()), Tt(n,n,Pol());

  unsigned int j=0;
  for (BitMap::iterator jt=used.begin(); jt(); ++jt,++j)
  {
    BlockElt y = *jt;
    set(Ts,j,j, b.T_coef(s,y,y)-Pol(1)); set(Tt,j,j, b.T_coef(t,y,y)-Pol(1));
    BlockEltList l; l.reserve(2);
    b.add_neighbours(l,s,*jt);
    for (unsigned int i=0; i<l.size(); ++i)
      if (used.isMember(l[i]))
	set(Ts,used.position(l[i]),j, b.T_coef(s,l[i],y));
    l.clear();
    b.add_neighbours(l,t,*jt);
    for (unsigned int i=0; i<l.size(); ++i)
      if (used.isMember(l[i]))
	set(Tt,used.position(l[i]),j, b.T_coef(t,l[i],y));
  }
  matrix::Vector<Pol> v(n,Pol()), w;
  v[used.position(x)]=Pol(1); w=v;

  // finally compute braid relation
  for (unsigned i=0; i<len; ++i)
    if (i%2==0)
      Ts.apply_to(v), Tt.apply_to(w);
    else
      Tt.apply_to(v), Ts.apply_to(w);

  cluster |= used;

  static bool verbose = false;
  bool success = v==w;
  if (verbose and (not success or b.z(x)==59))
  {
    //    std::cout << "success: " << success << std::endl;
    show_mat(std::cout,Ts,s);
    std::cout << std::endl;
    show_mat(std::cout,Tt,t);
  }
  return success;
} // |check_braid|

template containers::simple_list<BlockElt> ext_block::condense
  (matrix::Matrix<int>& M, const param_block& parent) const;
template containers::simple_list<BlockElt> ext_block::condense
  (matrix::Matrix<Split_integer>& M, const param_block& parent) const;

} // |namespace ext_block|

} // |namespace atlas|
