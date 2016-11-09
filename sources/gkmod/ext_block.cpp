/*
  This is ext_block.cpp

  Copyright (C) 2013-2016 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_block.h"

#include <cassert>
#include <vector>
#include <iostream> // |std::cout| used in methods like |ext_block::list_edges|

#include "innerclass.h"
#include "weyl.h"
#include "kgb.h"
#include "blocks.h"
#include "repr.h"

#include "bitmap.h"
#include "polynomials.h"
#include "matreduc.h"
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
  while (max>min) // invar: |(m==0 or z(m-1)<zz) and (max<size() or z(max)>=zz)|
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

unsigned int ext_block::list_edges()
{
  std::set<BlockEltPair>::iterator it;
  std::cout << "flipped edges:" << std::endl;
  unsigned int count=0;
  for (BlockElt x=0; x<info.size(); ++x)
    for (unsigned i=0; i<2; ++i) // two groups of links
      for (auto it=info[x].flips[i].begin(); it(); ++it)
      {
	auto &l = data[*it][x].links;
	std::cout << z(x) << "->" << z(i==0 ? l.first : l.second)
		  << ", s=" << *it+1 << std::endl;
	++count;
      }
  std::cout << std::endl;
  return count;
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
  assert((delta-1)*E.lambda_rho()==(1-theta)*E.tau());
  assert((delta-1).right_prod(E.l())==(theta+1).right_prod(E.t()));
  assert(((E.ctxt.g()-E.l()-rho_check(rd))*(1-theta)).numerator().isZero());
  assert(((theta+1)*(E.ctxt.gamma()-E.lambda_rho()-rho(rd)))
	 .numerator().isZero());
  ndebug_use(delta); ndebug_use(theta); ndebug_use(rd);
}

param::param (const context& ec, const StandardRepr& sr)
  : ctxt(ec)
  , tw(ec.rc().kgb().involution(sr.x()))
  , d_l(ell(ec.realGroup().kgb(),sr.x()))
  , d_lambda_rho(ec.rc().lambda_rho(sr))
  , d_tau(matreduc::find_solution(1-theta(),(delta()-1)*lambda_rho()))
  , d_t(matreduc::find_solution
	(theta().transposed()+1,(delta()-1).right_prod(l())))
{
  validate(*this);
}

param::param (const context& ec, KGBElt x, const Weight& lambda_rho)
  : ctxt(ec)
  , tw(ec.realGroup().kgb().involution(x))
  , d_l(ell(ec.realGroup().kgb(),x))
  , d_lambda_rho(lambda_rho)
  , d_tau(matreduc::find_solution(1-theta(),(delta()-1)*lambda_rho))
  , d_t(matreduc::find_solution
	(theta().transposed()+1,(delta()-1).right_prod(l())))
{
  validate(*this);
}

param::param (const context& ec, const TwistedInvolution& tw,
	      Weight lambda_rho, Weight tau, Coweight l, Coweight t)
  : ctxt(ec), tw(tw)
  , d_l(std::move(l))
  , d_lambda_rho(std::move(lambda_rho))
  , d_tau(std::move(tau))
  , d_t(std::move(t))
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
	or E.ctxt.g()!=F.ctxt.g()
	or E.ctxt.gamma()!=F.ctxt.gamma())
      return false;
  } // otherwise there might still be a match, so fall through
  return E.theta()==F.theta()
    and in_R_image(E.theta()+1,E.l()-F.l())
    and in_L_image(E.lambda_rho()-F.lambda_rho(),E.theta()-1);
}

KGBElt param::x() const
{ TitsElt a(ctxt.innerClass().titsGroup(),TorusPart(l()),tw);
  return rc().kgb().lookup(a);
}

StandardRepr scaled_extended_dominant // result will have its |gamma()| dominant
(const Rep_context rc,
 const StandardRepr& sr, const WeightInvolution& delta,
 Rational factor, // |z.nu()| is scaled by |factor| first
 bool& flipped // records whether and extended flip was recorded
 )
{ const RootDatum& rd=rc.rootDatum();
  assert(is_dominant_ratweight(rd,sr.gamma())); // dominant
  assert(((delta-1)*sr.gamma().numerator()).isZero()); // $\delta$-fixed
  StandardRepr result = rc.sr(sr.x(),rc.lambda_rho(sr),sr.gamma()*factor);
  Weight gamma_numer(result.gamma().numerator().begin(),
		     result.gamma().numerator().end());
  context ctxt(rc,delta,result.gamma());
  const ext_gens orbits = rootdata::fold_orbits(rd,delta);

  Weight lr, tau; Coweight l,t;
  { param E(ctxt,result); // comput fields as for extended parameter
    lr=E.lambda_rho(); tau=E.tau(); l=E.l(); t=E.t();
  }
  KGBElt x = result.x(); // another variable, for convenience

  const auto grc = ctxt.g_rho_check(); int denom=grc.denominator();
  l*=denom; // scale to make shift applied to |l| below an integer vector
  Coweight l_offset(grc.numerator().begin(),grc.numerator().end()); // convert
  l-=l_offset; // shift so that reflections can apply directly
  { unsigned i; // index into |orbits|
    do
      for (i=0; i<orbits.size(); ++i)
	if (rc.kgb().status(x).isComplex(orbits[i].s0))
	{ const auto& s=orbits[i];
	  int v=rd.simpleCoroot(s.s0).dot(gamma_numer);
	  if (v<0)
	  { rd.act(s.w_kappa,gamma_numer); // change inf.char representative
	    lr = rd.image_by(s.w_kappa,lr) - rho_minus_w_rho(rd,s.w_kappa);
	    rd.act(s.w_kappa,tau);
	    rd.dual_act(l,s.w_kappa);
	    rd.dual_act(t,s.w_kappa);
	    x = rc.kgb().cross(s.w_kappa,x);
	    break;
	  }
	} // |for(s)|, if |isComplex|
    while(i<orbits.size()); // continue until above |for| runs to completion
  }
  l+=l_offset;
  l/=denom; // shift and scale back to original size

  // since |gamma| may have changed, we need to buid a new |context|
  context new_ctxt(rc,delta,
		   RatWeight(gamma_numer,result.gamma().denominator()));
  // now ensure that |E| gets matching |gamma| and |theta| for flipped test
  param E(new_ctxt,rc.kgb().involution(x),lr,tau,l,t);
  result = rc.sr_gamma(x,E.lambda_rho(),new_ctxt.gamma());
  flipped = not same_sign(E,param(new_ctxt,result));
  return result;
}

containers::sl_list<std::pair<StandardRepr,bool> > finalise
  (const repr::Rep_context& rc,
   StandardRepr sr, const WeightInvolution& delta)
{ // in order that |singular_generators| generate the whole singular system:
  rc.make_dominant(sr); // ensure that |sr.gamma()| is dominant
  context ctxt(rc,delta,sr.gamma());
  const ext_gens orbits = rootdata::fold_orbits(ctxt.id(),delta);
  const RankFlags singular_orbits =
    reduce_to(orbits,singular_generators(ctxt.id(),sr.gamma()));

  containers::sl_list<std::pair<param,bool> >
    to_do(1,std::make_pair(param(ctxt,sr),false));
  containers::sl_list<std::pair<StandardRepr,bool> > result;

  do
  { auto& head=to_do.front(); // but don't pop just yet
    const param E= std::move(head.first);
    bool flipped = head.second;
    to_do.pop_front(); // we are done with |head|
    auto s = first_descent_among(singular_orbits,orbits,E);
    if (s>=orbits.size()) // no singular descents, so append to result
      result.emplace_back(std::make_pair(E.restrict(),flipped==is_default(E)));
    else // |s| is a singular descent orbit
    { containers::sl_list<std::pair<int,param> > links;
      auto type = star(E,orbits[s],links);
      if (not is_like_compact(type)) // some descent, push to front of |to_do|
      { if(has_october_surprise(type))  flipped = not flipped;
	auto it = to_do.begin(); auto l_it=links.begin();
	to_do.insert(it,std::make_pair
		     (std::move(l_it->second),flipped==(l_it->first>0)));
	if (has_double_image(type)) // then append a second node after |head|
	{ ++it, ++l_it;
	  to_do.insert(it,std::make_pair
		       (std::move(l_it->second),flipped==(l_it->first>0)));
	}
      }
    }
  }
  while(not to_do.empty());

  return result;
}

#if 0 // unused code, but the formula is referred to in the comment below
int z (const param& E) // value modulo 4, exponent of imaginary unit $i$
{ return
    (E.l().dot((E.delta()-1)*E.tau()) + 2*E.t().dot(E.lambda_rho())) % 4;
}
#endif

/*
  The quotient of |z| values across a Cayley transform will only be used when
  value of |t| is the same upstairs as downstairs. Then in the quotient of |z|
  values the second term contributes |t.dot(E.lambda_rho()-F.lambda_rho())|.
  That difference is often zero for the Cayley transform by a \emph{simple}
  root |alpha| or (length 1) changes by |alpha| with |t.cdot(alpha)==0|; in
  those cases the second term in |z| contributes nothing. For Cayley by a
  non-simple root, although the quotient of |z| values might pick up a
  contribution from the second term due to a |Cayley_shift| added to
  |lambda_rho|, it turns out that the Right Thing is not to use that quotient,
  but the quotient evaluated at the simple situation. So for these cases we
  should just compute the contribution to the difference of values |z| that
  would come \emph{from its first term} above only. This function does that:
 */
int z_quot (const param& E, const param& F)
{ assert(E.t()==F.t()); // we require preparing |t| upstairs to get this
  int d = E.l().dot((E.delta()-1)*E.tau()) - F.l().dot((F.delta()-1)*F.tau());
  return arithmetic::exp_i(d); // asserts |d| is even, and returns $(-1)^(d/2)$
}

/*
  In some cases, notably 2i12, one or both of the Cayley transforms may
  involve (in the simple root case) a change |mu| in |lambda_rho| that does
  not necessarily satisfy |t.dot(mu)==0|. In such cases the previous version
  of |z_quot| is insufficient, and we should include a contribution from the
  second term. But retrieving |mu| from the parameters |E| and |F| themselves
  is complicated by the posssible contribution from |Cayley_shift| that should
  be ignored; however at the place of call the value of |mu| is explicitely
  available, so we ask here to pass |t.dot(mu)| as third argument |t_mu|.
 */

int z_quot (const param& E, const param& F, int t_mu)
{ return z_quot(E,F)*arithmetic::exp_minus_1(t_mu); }

// this implements (comparison using) the formula from Proposition 16 in
// "Parameters for twisted repressentations" (with $\delta-1 = -(1-\delta)$
// the relation is symmetric in |E|, |F|, although not obviously so
bool same_sign (const param& E, const param& F)
{
  assert(same_standard_reps(E,F));
  const WeightInvolution& delta = E.delta();
  Weight kappa1=E.tau(), kappa2=F.tau();
  kappa1 -= delta*kappa1;
  kappa2 -= delta*kappa2;
  int i_exp = E.l().dot(kappa1) - F.l().dot(kappa2);
  assert (i_exp%2==0);
  int n1_exp =
    (F.l()-E.l()).dot(E.tau()) + F.t().dot(F.lambda_rho()-E.lambda_rho());
  return (i_exp/2+n1_exp)%2==0;
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
  (const repr::Rep_context& rc, WeightInvolution delta, const RatWeight& gamma)
    : d_rc(rc)
    , d_delta(std::move(delta)), d_gamma(gamma)
    , d_g(rc.kgb().base_grading_vector()+rho_check(rc.rootDatum()))
    , integr_datum(integrality_datum(rc.rootDatum(),gamma))
    , sub(SubSystem::integral(rc.rootDatum(),gamma))
{}


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

// act by external twist |delta| on KGB element |x|
KGBElt twisted (const KGB& kgb, KGBElt x,
		const WeightInvolution& delta, const weyl::Twist& twist)
{
  RatCoweight g_rho_l = kgb.torus_factor(x);
  delta.right_mult(g_rho_l.numerator());
  const RatCoweight diff = (g_rho_l - kgb.base_grading_vector()).normalize();
  if (diff.denominator()!=1)
    return UndefKGB;

  TorusPart delta_t(diff.numerator());

  const WeylGroup& W = kgb.innerClass().weylGroup();
  TitsElt te = kgb.titsElt(x);
  WeylElt delta_w = W.translation(te.w(),twist); // act on twisted involution
  TitsElt delta_te (kgb.innerClass().titsGroup(),delta_t,delta_w);
  return kgb.lookup(delta_te);
}

// involution for dual KGB is just $\delta$ transposed, no factor $-w_0$ here
// both |kgb| and |dual_kgb| must be known to be twist-stable
BlockElt twisted (const Block& block,
		  const KGB& kgb, const KGB& dual_kgb, // all are needed
		  BlockElt z,
		  const WeightInvolution& delta,
		  const weyl::Twist& twist)
{ return block.element
    (twisted(kgb,block.x(z),delta,twist),
     twisted(dual_kgb,block.y(z),delta.transposed(),twist));
}

BlockElt twisted (const param_block& block, const KGB& kgb,
		  BlockElt z,
		  const WeightInvolution& delta,
		  const weyl::Twist& twist)
{
  KGBElt x = block.x(z);
  TorusElement y = block.y_rep(block.y(z));
  KGBElt xx= twisted(kgb,x,delta,twist);
  if (xx==UndefKGB)
    return UndefBlock;
  y.act_by(delta);
  BlockElt result=block.lookup(xx,y);
  return result;
}


/*
  An auxiliary routine to compute extended parameters across complex links.
  The situation is complicated by the fact that the cross action is by a
  generator of the folded integral system, so already to compute the cross
  action at the level of involutions involves "unfolding" the generator to a
  |length<=3| word in the integral generators, and then translating those
  integrally-simple reflections into a word on the simple generators. As a
  consequence many simple reflections for the full root datum can be involved,
  and it is not certain that all of them will be complex.
 */
param complex_cross(ext_gen p, const param& E)
{ const RootDatum& rd = E.rc().rootDatum();
  const InvolutionTable& i_tab = E.rc().innerClass().involution_table();
  auto &tW = E.rc().twistedWeylGroup(); // caution: |p| refers to integr. datum

  TwistedInvolution tw=E.tw;
  InvolutionNbr theta = i_tab.nr(tw);
  const RatWeight gamma_rho = E.ctxt.gamma() - rho(rd);
  RatWeight gamma_lambda =  gamma_rho - E.lambda_rho();
  auto& ga_la_num = gamma_lambda.numerator();
  Weight rho_r_shift = rd.twoRho(i_tab.real_roots(theta));

  const RatCoweight g_rho_check = E.ctxt.g() - rho_check(rd);
  RatCoweight torus_factor =  g_rho_check - E.l();
  auto& tf_num = torus_factor.numerator();
  Coweight dual_rho_im_shift = rd.dual_twoRho(i_tab.imaginary_roots(theta));

  Weight tau=E.tau();
  Coweight t=E.t();
  const RootDatum& id = E.ctxt.id();
  for (unsigned i=p.w_kappa.size(); i-->0; )
  { weyl::Generator s=p.w_kappa[i]; // generator for integrality datum
    tW.twistedConjugate(E.ctxt.subsys().reflection(s),tw);
    id.simple_reflect(s,ga_la_num);
    id.simple_reflect(s,rho_r_shift);
    id.simple_reflect(s,tau);
    id.simple_coreflect(tf_num,s);
    id.simple_coreflect(t,s);
    id.simple_coreflect(dual_rho_im_shift,s);
  }
  const RatWeight lr_ratvec = (gamma_rho - gamma_lambda).normalize();
  assert(lr_ratvec.denominator()==1);
  Weight lambda_rho(lr_ratvec.numerator().begin(),
		    lr_ratvec.numerator().end()); // convert to |Weight|
  rho_r_shift -= rd.twoRho(i_tab.real_roots(i_tab.nr(tw)));
  rho_r_shift/=2; // now it is just a sum of (real) roots
  Weight tau_corr = ((E.ctxt.delta()-1)*rho_r_shift)/2; // hope it divides

  const RatWeight l_ratvec = (g_rho_check - torus_factor).normalize();
  assert(l_ratvec.denominator()==1);
  Coweight l(l_ratvec.numerator().begin(), l_ratvec.numerator().end());
  dual_rho_im_shift -= rd.dual_twoRho(i_tab.imaginary_roots(i_tab.nr(tw)));
  dual_rho_im_shift/=2; // now it is just a sum of (imaginary) coroots
  Coweight t_corr = ((E.ctxt.delta()-1).right_prod(dual_rho_im_shift))/2;

  return param(E.ctxt, tw,
	       lambda_rho-rho_r_shift, tau+tau_corr,
	       l-dual_rho_im_shift, t+t_corr);
} // |complex_cross|


WeylWord fixed_conjugate_simple (const context& ctxt, RootNbr& alpha)
{ const RootDatum& rd = ctxt.innerClass().rootDatum();
  std::vector<weyl::Generator> delta (rd.semisimpleRank());
  std::vector<bool> is_length_3 (delta.size());

  // tabulate action of |delta| on simple roots
  // this could and should be precomputed in |ctxt|
  for (weyl::Generator s=0; s<delta.size(); ++s)
  { weyl::Generator t =
      rd.simpleRootIndex(rd.root_index(ctxt.delta()*rd.simpleRoot(s)));
    delta[s]=t;
    if (s==t)
      is_length_3[s]=false;
    else
    { delta[t]=s;
      is_length_3[s] = rd.cartan(s,t)<0;
    }
  }

  RootNbr delta_alpha = rd.root_index(ctxt.delta()*rd.root(alpha));
  WeylWord result;
  while (not rd.is_simple_root(alpha)) // also |break| halfway is possible
  {
    weyl::Generator s =
      rd.descent_set(alpha).andnot(rd.ascent_set(delta_alpha)).firstBit();
    assert(s<rd.semisimpleRank()); // exists for positive non-simple roots
    if (is_length_3[s] and //"sum of swapped non-commuting roots" case:
	rd.simple_reflected_root(s,alpha)==rd.simpleRootNbr(delta[s]))
      break;
    result.push_back(s);
    rd.simple_reflect_root(s,alpha);
    rd.simple_reflect_root(delta[s],delta_alpha);
    if (delta[s]!=s) // second generator for cases of length 2,3
    { result.push_back(delta[s]);
      rd.simple_reflect_root(delta[s],alpha);
      rd.simple_reflect_root(s,delta_alpha);
      if (is_length_3[s])
      { // final generator for cases of length 3
	result.push_back(s);
	rd.simple_reflect_root(s,alpha);
	rd.simple_reflect_root(delta[s],delta_alpha);
      }
    }
  }
  std::reverse(result.begin(),result.end());
  return result;
} // |fixed_conjugate_simple|

/*
  for real Cayley transforms, one will subtract $\rho_r$ from |lambda_rho|
  before projecting it parallel to |alpha| so as to make |alpha_v| vanish on
  |gamma-lambda_rho-rho|. Here we compute from |E.lambda_rho()|, corrected by
  that |shift|, the multiple of $\alpha/2$ that such a projection would add
  to |lambda_rho| (or equivalently, subtract from |gamma-lambda_rho-rho|).
*/
int level_a (const param& E, const Weight& shift, RootNbr alpha)
{
  const RootDatum& rd = E.rc().rootDatum();
  return (E.ctxt.gamma() - E.lambda_rho() + shift).dot(rd.coroot(alpha))
    - rd.colevel(alpha); // final term $<\alpha^\vee,\rho>$
}

// version of |type| that will also export signs for every element of |links|
DescValue star (const param& E,
		const ext_gen& p,
		containers::sl_list<std::pair<int,param> >& links)
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
  const WeightInvolution delta_1 = E.ctxt.delta()-1;
  switch (p.type)
  {
  case ext_gen::one:
    { const Weight& alpha = integr_datum.simpleRoot(p.s0);
      const Coweight& alpha_v = integr_datum.simpleCoroot(p.s0);
      const RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);

      if (theta_alpha==n_alpha) // length 1 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g() - E.l()).dot(alpha)-rd.level(n_alpha);
	if (tf_alpha%2!=0) // then $\alpha$ is compact
	  return one_imaginary_compact; // quit here, do not collect \$200

	// noncompact case
	const TwistedInvolution new_tw= tW.prod(subs.reflection(p.s0),E.tw);
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1; // upstairs

	int tau_coef = alpha_v.dot(E.tau()); // take $\tau_\alpha$ of table 2

	// try to make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const Weight rho_r_shift = repr::Cayley_shift(ic,i_tab.nr(new_tw),ww);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(E.t().dot(alpha)==0); // follows from $\delta*\alpha=\alpha$

	Weight first; // maybe a root with |(1-delta)*first==alpha|
	if (rd.is_simple_root(alpha_simple))
	  first = Weight(rd.rank(),0); // effectively not used in this case
	else
	{
	  --tau_coef; // the parity change and decrease are both relevant
	  weyl::Generator s = // first switched root index
	    rd.find_descent(alpha_simple);
	  first = // corresponding root summand, conjugated back
	      rd.root(rd.permuted_root(rd.simpleRootNbr(s),ww));
	} // with this set-up, |alpha_simple| needs no more inspection

	// now separate cases; based on type 1 or 2 first
	if (matreduc::has_solution(th_1,alpha))
	{ // type 1, so extended type is 1i1
	  result = one_imaginary_single;

	  Weight diff = // called $-\sigma$ in table 2 of [Ptr] (NOTE MINUS)
	      matreduc::find_solution(th_1,alpha); // solutions are equivalent

	  param F(E.ctxt,new_tw,
		  E.lambda_rho() + first + rho_r_shift, E0.tau()+diff*tau_coef,
		  E.l()+alpha_v*(tf_alpha/2), E.t());

 	  E0.set_l(tf_alpha%4==0 ? F.l()+alpha_v : F.l()); // for cross
	  assert(not same_standard_reps(E,E0));
	  int sign  = z_quot(E,F);
	  int sign0 = sign*z_quot(E0,F);

	  links.push_back(std::make_pair(sign,std::move(F))); // Cayley link
	  links.push_back(std::make_pair(sign0,std::move(E0))); // cross link
	} // end of 1i1 case
	else
	{ // imaginary type 2; now we need to distinguish 1i2f and 1i2s

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
	    assert(not matreduc::has_solution
		   (th_1, delta_1*(E.lambda_rho()+rho_r_shift)));
	    return one_imaginary_pair_switched; // case 1i2s
	  }
	  result = one_imaginary_pair_fixed;  // what remains is case 1i2f

	  param F0(E.ctxt,new_tw,
		   E.lambda_rho() + first + rho_r_shift,
		   E.tau() - alpha*(tau_coef/2) - first,
		   E.l() + alpha_v*(tf_alpha/2), E.t());
	  param F1(E.ctxt,new_tw,
		   F0.lambda_rho() + alpha, F0.tau(), F0.l(), E.t());

	  int sign0 = z_quot(E,F0);
	  int sign1 = z_quot(E,F1);

	  links.push_back(std::make_pair(sign0,std::move(F0))); // Cayley link
	  links.push_back(std::make_pair(sign1,std::move(F1))); // Cayley link
	} // end of type 2 case
      } // end of length 1 imaginary case

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 1 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s0),E.tw);

	Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	RootNbr alpha_0 = // maybe one of |alpha==alpha_0+alpha_1|
	  rd.is_simple_root(alpha_simple) ? 0 // unused
	  : rd.permuted_root(rd.simpleRootNbr(rd.find_descent(alpha_simple)),
			     ww);

	// test parity, taking into account modifications that will be applied
	bool shift_correct = // whether |alpha_0| is defined and real at |theta|
	  not rd.is_simple_root(alpha_simple) and
	  i_tab.root_involution(theta,alpha_0)==rd.rootMinus(alpha_0);
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
	  if (shift_correct)
	  {
	    rho_r_shift += a0; // non delta-fixed contribution

	    // now we must add $d$ to $\tau$ with $(1-\theta')d=(1-\delta)*a0$
	    // since $\theta'*a0 = a1 = \delta*a_0$, we can take $d=a0$
	    tau_correction = a0;
	    assert((i_tab.matrix(new_tw)-1)*tau_correction==delta_1*a0);
	  }
	}

	const Weight new_lambda_rho =
	  E.lambda_rho() - rho_r_shift + alpha*(level/2);
	assert((E.ctxt.gamma()-new_lambda_rho).dot(alpha_v)
	       ==rd.colevel(n_alpha)); // check that |level_a| did its work

	const int t_alpha = E.t().dot(alpha);
	if (type1)
	{ // now distinguish 1r1f and 1r1s
	  if (t_alpha%2!=0) // no effect of |alpha_simple|, unlike 1i2 cases
	    return one_real_pair_switched;
	  result = one_real_pair_fixed; // what remains is case 1r1f

	  E0.set_t(E.t() - alpha_v*(t_alpha/2));
	  assert(same_sign(E,E0)); // since only |t| changes

	  param F0(E.ctxt,new_tw,
		   new_lambda_rho, E.tau() + tau_correction, E.l(), E0.t());
	  param F1(E.ctxt,new_tw,
		   new_lambda_rho, F0.tau(), E.l() + alpha_v, E0.t());

	  int sign0 = z_quot(E0,F0), sign1 = z_quot(E0,F0);

	  links.push_back(std::make_pair(sign0,std::move(F0))); // first Cayley
	  links.push_back(std::make_pair(sign1,std::move(F1))); // second Cayley

	} // end of 1r1 case
	else // real type 2
	{
	  result = one_real_single;
	  Coweight diff = // called $s$ in table 2 of [Ptr]
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v);

	  E0.set_t(E.t() - diff*t_alpha);
	  assert(same_sign(E,E0)); // since only |t| changes

	  param E1 = E0; // for cross neighbour; share updated value of |t|
	  E1.set_lambda_rho(E.lambda_rho()+alpha);
	  assert(not same_standard_reps(E0,E1));

	  param F(E.ctxt,new_tw,
		  new_lambda_rho, E.tau() + tau_correction, E.l(), E0.t());
	  int sign0 = z_quot(E0,F);
	  int sign1 = sign0*z_quot(E1,F);

	  links.push_back(std::make_pair(sign0,std::move(F ))); // Cayley link
	  links.push_back(std::make_pair(sign1,std::move(E1))); // cross link
	}
      }
      else // length 1 complex case
      { result = rd.is_posroot(theta_alpha)
	  ? one_complex_ascent : one_complex_descent ;
	links.push_back(std::make_pair(1,complex_cross(p,E)));
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
	int tf_alpha = (E.ctxt.g() - E.l()).dot(alpha)-rd.level(n_alpha);
	int tf_beta = (E.ctxt.g() - E.l()).dot(beta)-rd.level(n_alpha);
	assert((tf_alpha-tf_beta)%2==0); // same compactness
	if (tf_alpha%2!=0) // then $\alpha$ and $\beta$ are compact
	  return two_imaginary_compact;

	// noncompact case
	const TwistedInvolution new_tw =
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));
	// make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const Weight rho_r_shift = repr::Cayley_shift(ic,i_tab.nr(new_tw),ww);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 2

	int at = alpha_v.dot(E.tau()); int bt = beta_v.dot(E.tau());
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1;

	if (matreduc::has_solution(th_1,alpha)) // then type 2i11
	{ result = two_imaginary_single_single;
	  const Weight sigma = matreduc::find_solution(th_1,alpha*at+beta*bt);

	  param F (E.ctxt, new_tw,
		   E.lambda_rho() + rho_r_shift,  E.tau() + sigma,
		   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t());

	  E0.set_l(E.l()+alpha_v+beta_v);
	  int sign = z_quot(E,F); // no 3rd arg, since |E.lambda_rho| unchanged
	  int sign0 = sign * z_quot(E0,F);
	  links.push_back(std::make_pair(sign,std::move(F)));	// Cayley link
	  links.push_back(std::make_pair(sign0,std::move(E0))); // cross link
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

	  const Weight new_tau0 = E.tau() - alpha*((at+m)/2) - beta*((bt-m)/2);
          const Coweight new_l = E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2);

	  param F0(E.ctxt, new_tw,
		   E.lambda_rho() + rho_r_shift + alpha*m, new_tau0,
		   new_l, E.t());
	  param F1(E.ctxt, new_tw,
		   E.lambda_rho() + rho_r_shift + alpha*mm, E.tau() + sigma,
		   new_l, E.t());

	  // compute signs before invoking |std::move|
	  int t_alpha=E.t().dot(alpha);
	  int sign0=z_quot(E,F0,m*t_alpha), sign1=z_quot(E,F1,mm*t_alpha);

	  // first Cayley link will be the one that does not need |sigma|
	  links.push_back(std::make_pair(sign0,std::move(F0))); // first Cayley
	  links.push_back(std::make_pair(sign1,std::move(F1))); // second Cayley
	} // end of case 2i12f
	else
	{ // type 2i22
	  result = two_imaginary_double_double;
	  // $\alpha^\vee$ and $\beta^\vee$ are even on $(X^*)^\theta$ and
	  // $(1-\delta)\tau\in(X^*)^\theta+2X^*$ so $<av-bv,\tau>$ is even
	  assert((at-bt)%2==0);
	  int m = static_cast<unsigned int>(at)%2; // safe modular reduction

	  param F0(E.ctxt, new_tw,
		   E.lambda_rho() + rho_r_shift + alpha*m,
		   E.tau() - alpha*((at+m)/2) - beta*((bt-m)/2),
		   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t());
	  param F1(E.ctxt, new_tw,
		   E.lambda_rho() + rho_r_shift + alpha*(1-m) + beta,
		   E.tau() - alpha*((at-m)/2) - beta*((bt+m)/2),
		   F0.l(),E.t());
	  // get signs before invoking the |std::move| from |F0|, |F1|
	  int ta = E.t().dot(alpha), tb=E.t().dot(beta);
	  int sign0=z_quot(E,F0,ta*m)
	    , sign1=z_quot(E,F1,ta*(1-m)+tb);

	  links.push_back(std::make_pair(sign0,std::move(F0))); // first Cayley
	  links.push_back(std::make_pair(sign1,std::move(F1))); // second Cayley
	} // end type 2i22 case
      }

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 2 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here

	const Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return two_real_nonparity; // no link added here

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	WeightInvolution theta_1 = i_tab.matrix(theta)-1; // upstairs
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));

	const Weight new_lambda_rho = E.lambda_rho() - rho_r_shift
	  + alpha*(a_level/2) + beta*(b_level/2);

	int ta = E.t().dot(alpha); int tb = E.t().dot(beta);
	param E1=E; // another modifiable copy, like |E0|

	if (matreduc::has_solution(theta_1,alpha))
	{ // type 2r11
	  result = two_real_double_double;
	  // $\alpha$ is even on $(X_*)^{-\theta'}$ (so is $\beta$), and
	  // $t(1-\delta)\in(X_*)^{-\theta'}+2X_*$ so $<t,alpha-beta>$ is even
	  assert((ta-tb)%2==0);
	  int m =  static_cast<unsigned int>(ta)%2;

	  // set two values for |t|; actually the same value in case |m==0|
	  E0.set_t(E.t() - alpha_v*((ta+m)/2) - beta_v*((tb-m)/2));
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E0.t().dot(alpha)==-m and E0.t().dot(beta)==m);

	  E1.set_t(E.t() - alpha_v*((ta-m)/2) - beta_v*((tb+m)/2));
	  assert(same_sign(E,E1)); // since only |t| changes
	  assert(E1.t().dot(alpha)==m and E1.t().dot(beta)==-m);

	  param F0(E.ctxt, new_tw,
		   new_lambda_rho,E.tau(), E.l()+alpha_v*m, E0.t());
	  param F1(E.ctxt, new_tw,
		   new_lambda_rho,E.tau(), E.l()+alpha_v*(1-m)+beta_v,E1.t());

	  int sign0=z_quot(E0,F0,m*((b_level-a_level)/2));
	  int sign1=z_quot(E1,F1,m*((a_level-b_level)/2));

	  // Cayley links
	  links.push_back(std::make_pair(sign0,std::move(F0)));
	  links.push_back(std::make_pair(sign1,std::move(F1)));
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
	  E0.set_t(E.t() - alpha_v*((ta+m)/2) - beta_v*((tb-m)/2));
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E0.t().dot(alpha)==-m and E0.t().dot(beta)==m);

	  E1.set_t(E.t() - s);
	  assert(same_sign(E,E1)); // since only |t| changes
	  assert(E1.t().dot(alpha)==-mm and E1.t().dot(beta)==mm);

	  param F0(E.ctxt, new_tw,
		   new_lambda_rho, E.tau(), E.l()+alpha_v*m, E0.t());
	  param F1(E.ctxt, new_tw,
		   new_lambda_rho, E.tau(), E.l()+alpha_v*mm, E1.t());

	  int sign0=z_quot(E0,F0,m *((b_level-a_level)/2));
	  int sign1=z_quot(E1,F1,mm*((b_level-a_level)/2));

	  // Cayley links
	  links.push_back(std::make_pair(sign0,std::move(F0)));
	  links.push_back(std::make_pair(sign1,std::move(F1)));

	} // end of case 2r21f
	else // case 2r22
	{ result = two_real_single_single;
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v*ta+beta_v*tb);

	  E0.set_t(E.t() - s); // parameter adapted to Cayley transform |F|
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E.t().dot(alpha)==0 and E.t().dot(beta)==0);

	  E1.set_lambda_rho(E.lambda_rho()+alpha+beta);
	  E1.set_t(E0.t()); // cross action, keeps adaption of |t| to |F| below
	  assert(not same_standard_reps(E0,E1));

	  param F(E.ctxt, new_tw, new_lambda_rho, E.tau(), E.l(), E0.t());

	  int sign0=z_quot(E0,F); // no 3rd arg, as |E.t().dot(alpha)==0| etc.
	  int sign1=sign0*z_quot(E1,F); // total sign from |E| to its cross |E1|

	  links.push_back(std::make_pair(sign0,std::move(F ))); // Cayley link
	  links.push_back(std::make_pair(sign1,std::move(E1))); // cross link
	} // end of case 2r22
      }
      else // length 2 complex case
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha != (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // twisted non-commutation with |s0.s1|
	  result = ascent ? two_complex_ascent : two_complex_descent;
	  links.push_back(std::make_pair(1,complex_cross(p,E)));
	}
	else if (ascent)
	{ // twisted commutation with |s0.s1|: 2Ci
	  result = two_semi_imaginary;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p = i_tab.nr(new_tw); // upstairs
	  const Weight rho_r_shift = repr::Cayley_shift(ic,theta_p,ww);
	  assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	  // downstairs cross by |ww| only has imaginary and complex steps, so
	  // $\alpha_v.(\gamma-\lambda_\rho)$ is unchanged across |ww|
	  const int f = // number of times $\alpha$ is added to $\lambda_\rho$
	    (E.ctxt.gamma() - E.lambda_rho()).dot(alpha_v)- rd.colevel(n_alpha);

	  const Weight new_lambda_rho = E.lambda_rho() + alpha*f + rho_r_shift;
	  // both $\gama-\lambda$ and $\tau$ get $f*alpha$ subtracted by
	  // $\alpha$-reflection; adapt $\tau$ for vanishing $1-\delta$ image
	  const Weight new_tau = rd.reflection(n_alpha,E.tau()) + alpha*f;

	  // but |dual_v| needs correction by |ell_shift|
	  const int dual_f =
	    (E.ctxt.g() - E.l()).dot(alpha) - rd.level(n_alpha);

	  const Coweight new_l = E.l() + alpha_v*dual_f;
          const Coweight new_t =
	    rd.coreflection(E.t(),n_alpha) - alpha_v*dual_f;
	  param F (E.ctxt, new_tw, new_lambda_rho, new_tau, new_l, new_t);

	  int ab_tau = (alpha_v+beta_v).dot(E.tau());
	  assert (ab_tau%2==0);
	  int sign = arithmetic::exp_i(ab_tau * dual_f);
	  links.push_back(std::make_pair(sign,std::move(F)));  // "Cayley" link
	}
	else // twisted commutation with |s0.s1|, and not |ascent|: 2Cr
	{ result = two_semi_real;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	  assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	  const int f = level_a(E,rho_r_shift,n_alpha);

	  const Weight new_lambda_rho = E.lambda_rho() - rho_r_shift + alpha*f;
	  const Weight new_tau = rd.reflection(n_alpha,E.tau()) - alpha*f;

	  const int dual_f =
	    (E.ctxt.g() - E.l()).dot(alpha) - rd.level(n_alpha);
	  const Coweight new_l = E.l() + alpha_v*dual_f;
          const Coweight new_t =
	    rd.coreflection(E.t(),n_alpha) + alpha_v*dual_f;

	  param F (E.ctxt, new_tw, new_lambda_rho, new_tau, new_l, new_t);

	  int t_ab = E.t().dot(beta-alpha);
	  assert (t_ab%2==0);
	  int sign = arithmetic::exp_i(t_ab * (f+alpha_v.dot(E.tau())));
	  links.push_back(std::make_pair(sign,std::move(F)));  // "Cayley" link
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

      if (theta_alpha==n_alpha) // length 3 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g() - E.l()).dot(alpha)-rd.level(n_alpha);
	int tf_beta = (E.ctxt.g() - E.l()).dot(beta)-rd.level(n_alpha);
	assert((tf_alpha-tf_beta)%2==0); // same compactness
	if (tf_alpha%2!=0) // then $\alpha$ and $\beta$ are compact
	  return three_imaginary_compact;

	// noncompact case
	result = three_imaginary_semi;

	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const Weight rho_r_shift = repr::Cayley_shift(ic,i_tab.nr(new_tw),ww);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 3

	param F(E.ctxt, new_tw,
		E.lambda_rho() + rho_r_shift,
		E.tau() - alpha*kappa_v.dot(E.tau()),
		E.l() + kappa_v*((tf_alpha+tf_beta)/2), E.t());

	int sign = z_quot(E,F); // |lambda_rho| unchanged at simple Cayley

	links.push_back(std::make_pair(sign,std::move(F))); // Cayley link
      }
      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 3 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here

	const Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return three_real_nonparity; // no link added here

	// parity case
	result = three_real_semi;

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	const Weight new_lambda_rho = // make level for |kappa| zero
	  E.lambda_rho()-rho_r_shift + kappa*((a_level+b_level)/2);

	E0.set_t(E.t()-alpha_v*kappa.dot(E.t())); // makes |E.t().dot(kappa)==0|
	assert(same_sign(E,E0)); // since only |t| changes

	param F(E.ctxt, new_tw,	new_lambda_rho,E.tau(), E.l(), E0.t() );

	int sign = z_quot(E0,F); // no 3rd arg since |E.t().dot(kappa)==0|
	links.push_back(std::make_pair(sign,std::move(F))); // Cayley link
      }
      else // length 3 complex case
      { const bool ascent = rd.is_posroot(theta_alpha);
	const RootNbr n_beta = subs.parent_nr_simple(p.s1);
	if (theta_alpha == (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // reflection by |alpha+beta| twisted commutes with |E.tw|: 3Ci or 3Cr
	  result = ascent ? three_semi_imaginary : three_semi_real;

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const Weight rho_r_shift =
	    repr::Cayley_shift(ic,ascent ? i_tab.nr(new_tw) : theta,ww);
	  assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	  int tf_alpha = (E.ctxt.g() - E.l()).dot(alpha) - rd.level(n_alpha);
	  int dtf_alpha = (E.ctxt.gamma() - E.lambda_rho()).dot(alpha_v)
	    - rd.colevel(n_alpha);
	  Weight new_lambda_rho = E.lambda_rho() + rho_r_shift; // for now

	  if (ascent) // 3Ci
	  { param F(E.ctxt,new_tw,
		    dtf_alpha%2==0 ? new_lambda_rho : new_lambda_rho + kappa,
		    E.tau() - kappa*(kappa_v.dot(E.tau())/2),
		    E.l() + kappa_v*tf_alpha, E.t());

	    assert(E.t().dot(kappa)==0);
	    // since it is half of |t*(1+theta)*kappa=l*(delta-1)*kappa==0|
	    int sign  = z_quot(E,F); // may ignore possible shift by |kappa|
	    links.push_back(std::make_pair(sign,std::move(F))); // Cayley link
	  }
	  else // 3Cr
	  {
	    E0.set_t // make |E.t().dot(kappa)==0| using |kappa_v|
	      (E.t() - kappa_v*(kappa.dot(E.t())/2));
	    assert(same_sign(E,E0)); // since only |t| changes

	    param F(E.ctxt, new_tw,
		    new_lambda_rho + kappa*dtf_alpha, E.tau(),
		    tf_alpha%2==0 ? E.l() : E.l()+kappa_v, E0.t());

	    int sign = z_quot(E0,F); // no 3rd arg since |E.t().dot(kappa)==0|
	    links.push_back(std::make_pair(sign,std::move(F))); // Cayley link
	  }

	}
	else // twisted non-commutation: 3C+ or 3C-
	{
	  result = ascent ? three_complex_ascent : three_complex_descent;
	  links.push_back(std::make_pair(1,complex_cross(p,E)));
	}
      }
    }
    break;
  }

  // add a flip to links with a length difference of 2
  if (p.length()-(has_defect(result)?1:0)==2)
  { auto it=links.begin(); auto c=scent_count(result);
    for (unsigned i=0; i<c; ++i,++it) // only affect ascent/descent links
      it->first  = -it->first; // do the flip
  }

  return result;
} // |star|

bool is_descent (const ext_gen& kappa, const param& E)
{ // like in |star|, generators are not simple for |E.rc().twistedWeylGroup()|
  const InnerClass& ic = E.rc().innerClass();
  const InvolutionTable& i_tab = ic.involution_table();
  const InvolutionNbr theta = i_tab.nr(E.tw); // so use root action of |E.tw|
  const RootDatum& rd = E.rc().rootDatum();

  const RootNbr n_alpha = E.ctxt.subsys().parent_nr_simple(kappa.s0);
  const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
  const Weight& alpha = E.ctxt.id().simpleRoot(kappa.s0);

  // we don't need to inspect |kappa.type|, it does not affect descent status
  if (theta_alpha==n_alpha) // imaginary case, return whether compact
    return ( ((E.ctxt.g()-E.l()).dot(alpha)-rd.level(n_alpha)) %2!=0 );
  if (theta_alpha==rd.rootMinus(n_alpha)) // real, return whether parity
  {
    RootNbr alpha_simple = n_alpha; // copy to be made simple
    const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
    const auto level = level_a(E,repr::Cayley_shift(ic,theta,ww),n_alpha);
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

  if (twisted(kgb,0,delta,twist)==UndefKGB or
      twisted(dual_kgb,0,delta.transposed(),twist)==UndefKGB)
    return; // if one or other not delta-stable, leave |size==0| and quit

  for (BlockElt z=0; z<block.size(); ++z)
    if (twisted(block,kgb,dual_kgb,z,delta,twist)==z)
      fixed_points.insert(z);

  complete_construction(fixed_points);
  // FIXME cannot call |check| here, although setting sign flips depends on it

} // |ext_block::ext_block|

// we use these prediates to flip edges

bool is_2ir (DescValue v)
{ return generator_length(v)==2 and not is_complex(v) and not has_defect(v); }

bool is_2C (DescValue v)
{ return generator_length(v)==2 and is_complex(v); }

bool is_3Cir (DescValue v)
{ return generator_length(v)==3 and has_defect(v); }

ext_block::ext_block // for an external twist
  (const InnerClass& G,
   const param_block& block, const KGB& kgb,
   const WeightInvolution& delta,
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

  // the following is NOT |twist(orbits)|, which would be for subsystem
  weyl::Twist twist(fold_orbits(block.rootDatum(),delta));

  // test if twisting some block element lands in the same block
  if (twisted(block,kgb,0,delta,twist)==UndefBlock)
    return; // if block not delta-stable, leave |size==0| and quit

  for (BlockElt z=0; z<block.size(); ++z)
    if (twisted(block,kgb,z,delta,twist)==z)
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
// the signs are recorded in |eb|, and printed to |cout| is |verbose| holds.
bool ext_block::check(const param_block& block, bool verbose)
{
  context ctxt (block.context(),delta(),block.gamma());
  containers::sl_list<std::pair<int,param> > links;
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
	  assert(same_standard_reps(it->second,F)); // must lie over same
	  if (it->first!=sign_between(it->second,F)) // here != means XOR
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
	  assert(same_standard_reps(it->second,F));
	  if (it->first!=sign_between(it->second,F))
	  {
	    flip_edge(s,n,m);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
	                << " from " << z << " to " << Cz << '.' << std::endl;
	  }
	  ++it;
	  m=cross(s,n); BlockElt cz = this->z(m);
	  param Fc(ctxt,block.x(cz),block.lambda_rho(cz));
	  assert(same_standard_reps(it->second,Fc));
	  if (it->first!=sign_between(it->second,Fc))
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
	  assert(same_standard_reps(it->second,F));
	  if (it->first!=sign_between(it->second,F))
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
	  bool straight=same_standard_reps(it->second,F0);
          const auto& node0 = straight ? *it : *std::next(it);
          const auto& node1 = straight ? *std::next(it) : *it;
	  if (not straight)
	    assert(same_standard_reps(node0.second,F0));
	  assert(same_standard_reps(node1.second,F1));
	  if (node0.first!=sign_between(node0.second,F0))
	  {
	    flip_edge(s,n,m.first);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
			<< " from " << z << " to " << Cz0 << '.' << std::endl;
	  }
	  if (node1.first!=sign_between(node1.second,F1))
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
	  bool straight=same_standard_reps(it->second,F0);
          const auto& node0 = straight ? *it : *std::next(it);
          const auto& node1 = straight ? *std::next(it) : *it;
	  if (not straight)
	    assert(same_standard_reps(node0.second,F0));
	  assert(same_standard_reps(node1.second,F1));
	  if (node0.first!=sign_between(node0.second,F0))
	  {
	    flip_edge(s,n,m.first);
	    if (verbose)
	      std::cout << "Flip at Cayley link " << unsigned{s}
			<< " from " << z << " to " << Cz0 << '.' << std::endl;
	  }
	  if (node1.first!=sign_between(node1.second,F1))
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
  return true; // report sucess if we get here
} // |check|

void ext_block::flip_edges(extended_predicate match)
{
  DescValue t;
  for (weyl::Generator s=0; s<rank(); ++s)
    for (BlockElt n=0; n<size(); ++n)
      if (match(t=descent_type(s,n)))
      { if (is_like_compact(t) or is_like_nonparity(t))
	  continue;
	if (has_double_image(t))
	{ BlockEltPair Cs = Cayleys(s,n);
	  flip_edge(s,n,Cs.first);
	  flip_edge(s,n,Cs.second);
	}
	else
	  flip_edge(s,n,some_scent(s,n));
      }
}

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

// reduce a matrix to elements without descents among singular generators
template<typename C> // matrix coefficient type (signed)
containers::simple_list<BlockElt> // returns list of elements selected
  ext_block::condense(matrix::Matrix<C>& M, const param_block& parent) const
{
  RankFlags sing_orbs = singular_orbits(parent);
  containers::simple_list<BlockElt> result;

  for (BlockElt y=M.numColumns(); y-->0; ) // reverse loop is essential here
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

  BitMap todo(b.size()),used(b.size());
  todo.insert(x);
  for (unsigned int i=0; i<len; ++i)
    for (BitMap::iterator it=todo.begin(); it(); ++it)
    {
      used.insert(*it);
      todo.remove(*it);
      BlockEltList l; l.reserve(4);
      b.add_neighbours(l,s,*it);
      b.add_neighbours(l,t,*it);
      for (unsigned j=0; j<l.size(); ++j)
	if (not used.isMember(l[j]))
	  todo.insert(l[j]);
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
