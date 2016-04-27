/*
  This is ext_block.cpp

  Copyright (C) 2013-2016 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_block.h"

#include <cassert>
#include <vector>

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

bool has_double_image(DescValue v)
{
  static unsigned long mask =
      1ul << one_real_pair_fixed         | 1ul << one_imaginary_pair_fixed
    | 1ul << two_real_double_double      | 1ul << two_imaginary_double_double
    | 1ul << two_imaginary_single_double_fixed
    | 1ul << two_real_single_double_fixed;

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

bool is_like_nonparity(DescValue v)
{
  static unsigned long mask =
      1ul << one_imaginary_pair_switched | 1ul << one_real_nonparity
    | 1ul << two_real_nonparity          | 1ul << three_real_nonparity;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool is_like_compact(DescValue v)
{
  static unsigned long mask =
      1ul << one_real_pair_switched | 1ul << one_imaginary_compact
    | 1ul << two_imaginary_compact  | 1ul << three_imaginary_compact;

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

bool has_quadruple(DescValue v)
{
  static const unsigned long mask =
      1ul << two_imaginary_single_double_fixed
    | 1ul << two_real_single_double_fixed;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

bool is_proper_ascent(DescValue v)
{
  return not(is_descent(v) or is_like_nonparity(v));
}

bool has_defect(DescValue v)
{
  static unsigned long mask =
      1ul << two_semi_imaginary   | 1ul << two_semi_real
    | 1ul << three_semi_imaginary | 1ul << three_real_semi
    | 1ul << three_imaginary_semi | 1ul << three_semi_real;

  return (1ul << v & mask) != 0; // whether |v| is one of the above
}

int generator_length(DescValue v)
{ return v<two_complex_ascent ? 1 : v<three_complex_ascent ? 2 : 3; }

context::context
  (const repr::Rep_context& rc, WeightInvolution delta, const RatWeight& gamma)
    : d_rc(rc)
    , d_delta(std::move(delta)), d_gamma(gamma)
    , d_g(rc.kgb().base_grading_vector()+rho_check(rc.rootDatum()))
    , integr_datum(integrality_datum(rc.rootDatum(),gamma))
    , sub(SubSystem::integral(rc.rootDatum(),gamma))
{}

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
  ndebug_use(delta); ndebug_use(theta); ndebug_use(rd);
  assert(((theta+1)*(E.ctxt.gamma()-E.lambda_rho()-rho(rd)))
	 .numerator().isZero());
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

// whether |E| and |F| lie over equivalent |StandrdRepr| values
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

KGBElt x(const param& E)
{ TorusPart tp(E.l());
  TitsElt a(E.ctxt.innerClass().titsGroup(),tp,E.tw);
  return E.rc().kgb().lookup(a);
}

// this implements (comparison using) the formula from Propodition 16 in
// "Parameters for twisted repressentations" (with $\delta-1=-(1-\delta)$
// the relation is symmetric in |E|, |F|, although not obviously so
bool signs_differ (const param& E, const param& F)
{
  const WeightInvolution& delta = E.delta();
  Weight kappa1=E.tau(), kappa2=F.tau();
  kappa1 -= delta*kappa1;
  kappa2 -= delta*kappa2;
  int i_exp = E.l().dot(kappa1) - F.l().dot(kappa2);
  assert (i_exp%2==0);
  int n1_exp =
    (F.l()-E.l()).dot(E.tau()) + F.t().dot(F.lambda_rho()-E.lambda_rho());
  return (i_exp/2+n1_exp)%2!=0;
}

bool sign_differs_with_one_of (const param& E, const param& F1, const param& F2)
{ return
    same_standard_reps(E,F1) ? signs_differ(E,F1)
    : same_standard_reps(E,F2) ? signs_differ(E,F2)
    : throw std::runtime_error("Neither candidate has same standard repn");
}

bool is_default (const param& E)
{ return not signs_differ(E,param(E.ctxt,x(E),E.lambda_rho())); }

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
	    link=block.cross(p.s0,link), assert(fixed_points.isMember(link));
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
	    link=block.cross(p.s0,link), assert(fixed_points.isMember(link));
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
  for (unsigned i=p.w_tau.size(); i-->0; )
  { weyl::Generator s=p.w_tau[i]; // generator for integrality datum
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

// for real Cayley transforms, $\lambda$ needs $\rho_r$ shift, then projection;
// so compute scalar product of |alpha_check| with dual torus factor plus shift
int level_a (const param& E, const Weight& shift, RootNbr alpha)
{
  const RootDatum& rd = E.rc().rootDatum();
  return (E.ctxt.gamma() - E.lambda_rho() + shift).dot(rd.coroot(alpha))
    - rd.colevel(alpha); // final term $<\alpha^\vee,\rho>$
}

DescValue type (const param& E, const ext_gen& p,
		containers::sl_list<param>& links)
{
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

	// now separate cases; based on type 1/2 first
	if (matreduc::has_solution(th_1,alpha))
	{ // type 1, so extended type is 1i1
	  result = one_imaginary_single;

	  Weight diff = // called $-\sigma$ in table 2 of [Ptr] (NOTE MINUS)
	      matreduc::find_solution(th_1,alpha); // solutions are equivalent
	  links.push_back(param // Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + first + rho_r_shift,
			   E.tau() + diff*tau_coef - first,
			   E.l()+alpha_v*(tf_alpha/2),
			   E.t()
			   ));
	  links.push_back(param // cross link
			  (E.ctxt,E.tw,E.lambda_rho(),
			   E.tau(), E.l()+alpha_v, E.t()));
	} // end of type 1 case
	else
	{ // type 2; now we need to distinguish 1i2f and 1i2s

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
	    assert(not matreduc::has_solution
		   (th_1, delta_1*(E.lambda_rho()+rho_r_shift)));
	    return one_imaginary_pair_switched; // case 1i2s
	  }
	  result = one_imaginary_pair_fixed;  // what remains is case 1i2f

	  links.push_back(param // first Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + first + rho_r_shift,
			   E.tau() - alpha*(tau_coef/2) - first,
			   E.l()+alpha_v*(tf_alpha/2),
			   E.t()
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + first + rho_r_shift + alpha,
			   E.tau() - alpha*(tau_coef/2) - first,
			   E.l()+alpha_v*(tf_alpha/2),
			   E.t()
			   ));
	} // end of type 2 case
      } // end of length 1 imaginary case

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 1 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);

	const Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	const int level = level_a(E,rho_r_shift,n_alpha);

	if (level%2!=0) // nonparity
	   return one_real_nonparity; // no link added here

	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s0),E.tw);

	Weight new_lambda_rho = E.lambda_rho()-rho_r_shift + alpha*(level/2);
	assert((E.ctxt.gamma()-new_lambda_rho).dot(alpha_v)
	       ==rd.colevel(n_alpha)); // check that |level_a| did its work

	const WeightInvolution& th_1 = i_tab.matrix(E.tw)-1; // upstairs
	bool type1 = matreduc::has_solution(th_1,alpha);

	Weight tau_correction;
	if (rd.is_simple_root(alpha_simple))
	  tau_correction = Weight(rd.rank(),0); // no correction needed here
	else
	{
	  weyl::Generator s = // first switched root index
	    rd.find_descent(alpha_simple);
	  Weight first = // corresponding root summand, conjugated back
	    rd.root(rd.permuted_root(rd.simpleRootNbr(s),ww));
	  new_lambda_rho -= first; // final, non-delta-fixed contribution

	  if (type1)
	    // $(1-\theta)(s+f)=(\delta-1)*(-f)$ when $(\theta-1)s=\alpha$
	    tau_correction = matreduc::find_solution(th_1,alpha)+first;
	  else // type 2; now |matrix(new_tw)*first==delta*first| (second root)
	    tau_correction = first; // since $(1-\theta)f = (\delta-1)*(-f)$
	  assert((i_tab.matrix(new_tw)-1)*tau_correction==delta_1*-first);
	}

	const int t_alpha = E.t().dot(alpha);
	if (type1)
	  { // now distinguish 1r1f and 1r1s
	    if (t_alpha%2!=0) // no effect |alpha_simple|, unlike 1i2 cases
	      return one_real_pair_switched;
	    result = one_real_pair_fixed; // what remains is case 1r1f
	    const Coweight new_t = E.t() - alpha_v*(t_alpha/2);
	    links.push_back(param // first Cayley link
			    (E.ctxt,new_tw,
			     new_lambda_rho, E.tau() + tau_correction,
			     E.l(), new_t
			     ));
	    links.push_back(param // second Cayley link
			    (E.ctxt,new_tw,
			     new_lambda_rho, E.tau() + tau_correction,
			     E.l() + alpha_v, new_t
			     ));
	  } // end of real type 1 case
	  else // real type 2
	  {
	    result = one_real_single;
	    Coweight diff = // called $s$ in table 2 of [Ptr]
	      matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				      alpha_v);
	    links.push_back(param // Cayley link
			    (E.ctxt,new_tw,
			     new_lambda_rho,
			     E.tau() + tau_correction,
			     E.l(),
			     E.t() - diff*t_alpha
			     ));
	    links.push_back(param // cross link
	        (E.ctxt,E.tw, E.lambda_rho()+alpha ,E.tau(),E.l(),E.t()));
	  }
      }
      else // length 1 complex case
      { result = rd.is_posroot(theta_alpha)
	  ? one_complex_ascent : one_complex_descent ;
	links.push_back(complex_cross(p,E));
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
	  links.push_back(param // Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + rho_r_shift,
			   E.tau() + sigma,
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2),
			   E.t()
			   ));
	  links.push_back(param // cross link
			  (E.ctxt,E.tw,
			   E.lambda_rho(),E.tau(),
			   E.l()+alpha_v+beta_v, E.t()));
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

	  // first Cayley link will be the one that does not need |sigma|
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + rho_r_shift + alpha*m,
			   E.tau() - alpha*((at+m)/2) - beta*((bt-m)/2),
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t()
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + rho_r_shift + alpha*mm,
			   E.tau() + sigma,
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t()
			   ));
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(E.tw).transposed()+1,
				    beta_v-alpha_v);
	  links.push_back(param // false cross action link
			  (E.ctxt,E.tw,
			   E.lambda_rho(),E.tau(), E.l()+alpha_v, E.t()+s));
	} // end of case 2i12f
	else
	{ // type 2i22
	  result = two_imaginary_double_double;
	  // $\alpha^\vee$ and $\beta^\vee$ are even on $(X^*)^\theta$ and
	  // $(1-\delta)\tau\in(X^*)^\theta+2X^*$ so $<av-bv,\tau>$ is even
	  assert((at-bt)%2==0);
	  int m =  static_cast<unsigned int>(at)%2; // safe modular reduction
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + alpha*m + rho_r_shift,
			   E.tau() - alpha*((at+m)/2) - beta*((bt-m)/2),
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2),
			   E.t()
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + alpha*(1-m) + beta + rho_r_shift,
			   E.tau() - alpha*((at-m)/2) - beta*((bt+m)/2),
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2),
			   E.t()
			   ));
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
	const TwistedInvolution new_tw =
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));

	const Weight new_lambda_rho = E.lambda_rho()-rho_r_shift
	  + alpha*(a_level/2) + beta*(b_level/2);

	int ta = E.t().dot(alpha); int tb = E.t().dot(beta);

	if (matreduc::has_solution(theta_1,alpha))
	{ // type 2r11
	  result = two_real_double_double;
	  // $\alpha$ and $\beta$ are even on $(X_*)^{-\theta'}$ and
	  // $t(1-\delta)\in(X_*)^{-\theta'}+2X_*$ so $<t,alpha-beta>$ is even
	  assert((ta-tb)%2==0);
	  int m =  static_cast<unsigned int>(ta)%2;
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*m,
			   E.t() - alpha_v*((ta+m)/2)+beta_v*((tb-m)/2)
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*(1-m)+beta_v,
			   E.t() - alpha_v*((ta-m)/2)+beta_v*((tb+m)/2)
			   ));
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

	  // first Cayley link will be the one that does not need |sigma|
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*m,
			   E.t() - alpha_v*((ta+m)/2) - beta_v*((tb-m)/2)
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*mm,
			   E.t() - s
			   ));
	  const Weight sigma = matreduc::find_solution(theta_1,alpha-beta);
	  links.push_back(param // false cross action link
			  (E.ctxt,E.tw,
			   E.lambda_rho()+alpha,E.tau()+sigma, E.l(), E.t()));
	} // end of case 2r21f
	else // case 2r22
	{ result = two_real_single_single;
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v*ta+beta_v*tb);
	  links.push_back(param // Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l(), E.t() - s
			   ));
	  links.push_back(param // cross link
			  (E.ctxt,E.tw,
			   E.lambda_rho()+alpha+beta ,E.tau(),E.l(),E.t()));
	}
      }
      else // length 2 complex case
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha == (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // twisted commutation with |s0.s1|: 2Ci or 2Cr
	  result = ascent ? two_semi_imaginary : two_semi_real;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  const int f = (E.ctxt.gamma() - E.lambda_rho()).dot(alpha_v)
	    - rd.colevel(n_alpha);
	  const int dual_f = (E.ctxt.g() - E.l()).dot(alpha)
	    - rd.level(n_alpha);
	  const Weight new_lambda_rho = E.lambda_rho() + alpha*f;
	  const Weight new_tau =
	    rd.reflection(n_alpha,E.tau()) + alpha*(ascent ? f : -f);
	  const Coweight new_l = E.l() + alpha_v*dual_f;
          const Coweight new_t = rd.coreflection(E.t(),n_alpha)
	    + alpha_v*(ascent ? -dual_f : dual_f);
	  links.push_back(param (E.ctxt, new_tw, // Cayley link
				 new_lambda_rho, new_tau, new_l, new_t));
	}
	else // twisted non-commutation with |s0.s1|
	{
	  result = ascent ? two_complex_ascent : two_complex_descent;
	  links.push_back(complex_cross(p,E));
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
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 2

	const Weight new_tau = E.tau() - alpha*kappa_v.dot(E.tau());
 	const Coweight new_l = E.l() + kappa_v*((tf_alpha+tf_beta)/2);
	links.push_back(param // Cayley link
			(E.ctxt, new_tw,
			 E.lambda_rho() + rho_r_shift, new_tau, new_l, E.t()));
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

	const Weight new_lambda_rho =
	  E.lambda_rho()-rho_r_shift + kappa*((a_level+b_level)/2);
	links.push_back(param // Cayley link
			(E.ctxt, new_tw,
			 new_lambda_rho,E.tau(),
			 E.l(),E.t()-alpha_v*kappa.dot(E.t())));
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

	  const Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	  assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	  int tf_alpha = (E.ctxt.g() - E.l()).dot(alpha) - rd.level(n_alpha);
	  int dtf_alpha = (E.ctxt.gamma() - E.lambda_rho()).dot(alpha_v)
	    - rd.colevel(n_alpha);
	  Weight new_lambda_rho = E.lambda_rho() + rho_r_shift; // for now

	  if (ascent) // 3Ci
	    links.push_back(param // Cayley link
			    (E.ctxt, new_tw,
			     dtf_alpha%2==0 ? new_lambda_rho
			      : new_lambda_rho + kappa,
			     E.tau() - kappa*(kappa_v.dot(E.tau())/2),
			     E.l() + kappa_v*tf_alpha,
			     E.t()));
	  else // 3Cr
	    links.push_back(param // Cayley link
			    (E.ctxt, new_tw,
			     new_lambda_rho + kappa*dtf_alpha, E.tau(),
			     tf_alpha%2==0 ? E.l() : E.l()+kappa_v,
			     E.t() - kappa_v*(kappa.dot(E.t())/2)));

	}
	else // twisted non-commutation: 3C+ or 3C-
	{
	  result = ascent ? three_complex_ascent : three_complex_descent;
	  links.push_back(complex_cross(p,E));
	}
      }
    }
    break;
  }
  return result;
} // |type|

// version of |type| that will also export a sign, from normalisation of |E|
DescValue star (param& E, // parameter will be normalised relative to |p|
		const ext_gen& p,
		containers::sl_list<param>& links,
		int& sign)
{
  const param E0=E; // retain initial value for comparison
  DescValue result; sign=1; // default value, may be changed

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

	  // normalise
	  E.set_tau(E.tau()+diff*tau_coef);
	  E.set_l(E.l()+alpha_v*(tf_alpha/2));

	  links.push_back(param // Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + first + rho_r_shift, E.tau(),
			   E.l(), E.t()
			   ));
	  links.push_back(param // cross link
			  (E.ctxt,E.tw,
			   E.lambda_rho(),E.tau(), E.l()+alpha_v,E.t()));
	} // end of type 1 case
	else
	{ // type 2; now we need to distinguish 1i2f and 1i2s

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
	    assert(not matreduc::has_solution
		   (th_1, delta_1*(E.lambda_rho()+rho_r_shift)));
	    return one_imaginary_pair_switched; // case 1i2s
	  }
	  result = one_imaginary_pair_fixed;  // what remains is case 1i2f

	  // normalise
	  E.set_tau(E.tau()+alpha*(tau_coef/2));
	  E.set_l(E.l()+alpha_v*(tf_alpha/2));

	  links.push_back(param // first Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + first + rho_r_shift,
			   E.tau() - first,
			   E.l(), E.t()
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + first + rho_r_shift + alpha,
			   E.tau() - first,
			   E.l(), E.t()
			   ));
	} // end of type 2 case
      } // end of length 1 imaginary case

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 1 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);

	const Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	const int level = level_a(E,rho_r_shift,n_alpha);

	if (level%2!=0) // nonparity
	   return one_real_nonparity; // no link added here

	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s0),E.tw);

	E.set_lambda_rho(E.lambda_rho() + alpha*(level/2));

	Weight new_lambda_rho = E.lambda_rho()-rho_r_shift;
	assert((E.ctxt.gamma()-new_lambda_rho).dot(alpha_v)
	       ==rd.colevel(n_alpha)); // check that |level_a| did its work

	const WeightInvolution& th_1 = i_tab.matrix(E.tw)-1; // upstairs
	bool type1 = matreduc::has_solution(th_1,alpha);

	Weight tau_correction; // adapt to integrality based change of lambda
	if (rd.is_simple_root(alpha_simple))
	  tau_correction = Weight(rd.rank(),0); // no correction needed here
	else
	{
	  weyl::Generator s = // first switched root index
	    rd.find_descent(alpha_simple);
	  Weight first = // corresponding root summand, conjugated back
	    rd.root(rd.permuted_root(rd.simpleRootNbr(s),ww));
	  new_lambda_rho -= first; // final, non delta-fixed contribution

	  // we must add $d$ to $\tau$ with $(\theta'-1)d=(1-\delta)*first$
	  if (type1)
	    // $(1-\theta)(s+f)=(1-\delta)*f$ when $(\theta-1)s=\alpha$
	    tau_correction = matreduc::find_solution(th_1,alpha)+first;
	  else // type 2; now |matrix(new_tw)*first==delta*first| (second root)
	    tau_correction = first; // since $(1-\theta)f = (1-\delta)*f$
	  assert((i_tab.matrix(new_tw)-1)*tau_correction==delta_1*-first);
	}

	const int t_alpha = E.t().dot(alpha);
	if (type1)
	{ // now distinguish 1r1f and 1r1s
	  if (t_alpha%2!=0) // no effect of |alpha_simple|, unlike 1i2 cases
	    return one_real_pair_switched;
	  result = one_real_pair_fixed; // what remains is case 1r1f

	  E.set_t(E.t() - alpha_v*(t_alpha/2));

	  links.push_back(param // first Cayley link
			  (E.ctxt,new_tw,
			   new_lambda_rho, E.tau() + tau_correction,
			   E.l(), E.t()
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   new_lambda_rho, E.tau() + tau_correction,
			   E.l() + alpha_v, E.t()
			   ));
	} // end of real type 1 case
	else // real type 2
	{
	  result = one_real_single;
	  Coweight diff = // called $s$ in table 2 of [Ptr]
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v);

	  E.set_t(E.t() - diff*t_alpha);

	  links.push_back(param // Cayley link
			  (E.ctxt,new_tw,
			   new_lambda_rho,
			   E.tau() + tau_correction,
			   E.l(), E.t()
			   ));
	  links.push_back(param // cross link
	        (E.ctxt,E.tw, E.lambda_rho()+alpha ,E.tau(),E.l(),E.t()));
	}
      }
      else // length 1 complex case
      { result = rd.is_posroot(theta_alpha)
	  ? one_complex_ascent : one_complex_descent ;
	links.push_back(complex_cross(p,E));
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
	  links.push_back(param // Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + rho_r_shift,
			   E.tau() + sigma,
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2),
			   E.t()
			   ));
	  links.push_back(param // cross link
			  (E.ctxt,E.tw,
			   E.lambda_rho(),E.tau(),
			   E.l()+alpha_v+beta_v, E.t()));
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

	  // first Cayley link will be the one that does not need |sigma|
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + rho_r_shift + alpha*m,
			   E.tau() - alpha*((at+m)/2) - beta*((bt-m)/2),
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t()
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + rho_r_shift + alpha*mm,
			   E.tau() + sigma,
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t()
			   ));
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(E.tw).transposed()+1,
				    beta_v-alpha_v);
	  links.push_back(param // false cross action link
			  (E.ctxt,E.tw,
			   E.lambda_rho(),E.tau(), E.l()+alpha_v, E.t()+s));
	} // end of case 2i12f
	else
	{ // type 2i22
	  result = two_imaginary_double_double;
	  // $\alpha^\vee$ and $\beta^\vee$ are even on $(X^*)^\theta$ and
	  // $(1-\delta)\tau\in(X^*)^\theta+2X^*$ so $<av-bv,\tau>$ is even
	  assert((at-bt)%2==0);
	  int m =  static_cast<unsigned int>(at)%2; // safe modular reduction
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   E.lambda_rho() + alpha*m + rho_r_shift,
			   E.tau() - alpha*((at+m)/2) - beta*((bt-m)/2),
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2),
			   E.t()
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho() + alpha*(1-m) + beta + rho_r_shift,
			   E.tau() - alpha*((at-m)/2) - beta*((bt+m)/2),
			   E.l()+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2),
			   E.t()
			   ));
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
	const TwistedInvolution new_tw =
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));

	const Weight new_lambda_rho = E.lambda_rho()-rho_r_shift
	  + alpha*(a_level/2) + beta*(b_level/2);

	int ta = E.t().dot(alpha); int tb = E.t().dot(beta);

	if (matreduc::has_solution(theta_1,alpha))
	{ // type 2r11
	  result = two_real_double_double;
	  // $\alpha$ and $\beta$ are even on $(X_*)^{-\theta'}$ and
	  // $t(1-\delta)\in(X_*)^{-\theta'}+2X_*$ so $<t,alpha-beta>$ is even
	  assert((ta-tb)%2==0);
	  int m =  static_cast<unsigned int>(ta)%2;
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*m,
			   E.t() - alpha_v*((ta+m)/2)+beta_v*((tb-m)/2)
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*(1-m)+beta_v,
			   E.t() - alpha_v*((ta-m)/2)+beta_v*((tb+m)/2)
			   ));
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

	  // first Cayley link will be the one that does not need |sigma|
	  links.push_back(param // first Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*m,
			   E.t() - alpha_v*((ta+m)/2) - beta_v*((tb-m)/2)
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l()+alpha_v*mm,
			   E.t() - s
			   ));
	  const Weight sigma = matreduc::find_solution(theta_1,alpha-beta);
	  links.push_back(param // false cross action link
			  (E.ctxt,E.tw,
			   E.lambda_rho()+alpha,E.tau()+sigma, E.l(), E.t()));
	} // end of case 2r21f
	else // case 2r22
	{ result = two_real_single_single;
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v*ta+beta_v*tb);
	  links.push_back(param // Cayley link
			  (E.ctxt, new_tw,
			   new_lambda_rho, E.tau(),
			   E.l(), E.t() - s
			   ));
	  links.push_back(param // cross link
			  (E.ctxt,E.tw,
			   E.lambda_rho()+alpha+beta ,E.tau(),E.l(),E.t()));
	}
      }
      else // length 2 complex case
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha == (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // twisted commutation with |s0.s1|: 2Ci or 2Cr
	  result = ascent ? two_semi_imaginary : two_semi_real;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  const int f = (E.ctxt.gamma() - E.lambda_rho()).dot(alpha_v)
	    - rd.colevel(n_alpha);
	  const int dual_f = (E.ctxt.g() - E.l()).dot(alpha)
	    - rd.level(n_alpha);
	  const Weight new_lambda_rho = E.lambda_rho() + alpha*f;
	  const Weight new_tau =
	    rd.reflection(n_alpha,E.tau()) + alpha*(ascent ? f : -f);
	  const Coweight new_l = E.l() + alpha_v*dual_f;
          const Coweight new_t = rd.coreflection(E.t(),n_alpha)
	    + alpha_v*(ascent ? -dual_f : dual_f);
	  links.push_back(param (E.ctxt, new_tw, // Cayley link
				 new_lambda_rho, new_tau, new_l, new_t));
	}
	else // twisted non-commutation with |s0.s1|
	{
	  result = ascent ? two_complex_ascent : two_complex_descent;
	  links.push_back(complex_cross(p,E));
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
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 2

	const Weight new_tau = E.tau() - alpha*kappa_v.dot(E.tau());
 	const Coweight new_l = E.l() + kappa_v*((tf_alpha+tf_beta)/2);
	links.push_back(param // Cayley link
			(E.ctxt, new_tw,
			 E.lambda_rho() + rho_r_shift, new_tau, new_l, E.t()));
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

	const Weight new_lambda_rho =
	  E.lambda_rho()-rho_r_shift + kappa*((a_level+b_level)/2);
	links.push_back(param // Cayley link
			(E.ctxt, new_tw,
			 new_lambda_rho,E.tau(),
			 E.l(),E.t()-alpha_v*kappa.dot(E.t())));
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

	  const Weight rho_r_shift = repr::Cayley_shift(ic,theta,ww);
	  assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	  int tf_alpha = (E.ctxt.g() - E.l()).dot(alpha) - rd.level(n_alpha);
	  int dtf_alpha = (E.ctxt.gamma() - E.lambda_rho()).dot(alpha_v)
	    - rd.colevel(n_alpha);
	  Weight new_lambda_rho = E.lambda_rho() + rho_r_shift; // for now

	  if (ascent) // 3Ci
	    links.push_back(param // Cayley link
			    (E.ctxt, new_tw,
			     dtf_alpha%2==0 ? new_lambda_rho
			      : new_lambda_rho + kappa,
			     E.tau() - kappa*(kappa_v.dot(E.tau())/2),
			     E.l() + kappa_v*tf_alpha,
			     E.t()));
	  else // 3Cr
	    links.push_back(param // Cayley link
			    (E.ctxt, new_tw,
			     new_lambda_rho + kappa*dtf_alpha, E.tau(),
			     tf_alpha%2==0 ? E.l() : E.l()+kappa_v,
			     E.t() - kappa_v*(kappa.dot(E.t())/2)));

	}
	else // twisted non-commutation: 3C+ or 3C-
	{
	  result = ascent ? three_complex_ascent : three_complex_descent;
	  links.push_back(complex_cross(p,E));
	}
      }
    }
    break;
  }
  sign *= sign_between(E0,E);
  return result;
} // |star|


// move to standard form, return sign of change; root is integrally-simple
int normalize_1i1 (param& E, RootNbr n_alpha)
{ const param E0=E;
  const InnerClass& ic = E.rc().innerClass();
  const InvolutionTable& i_tab = ic.involution_table();
  const RootDatum& rd = E.rc().rootDatum();
  const Weight alpha=rd.root(n_alpha);
  const Weight alpha_v=rd.coroot(n_alpha);
  const int tf_alpha = (E.ctxt.g() - E.l()).dot(alpha)-rd.level(n_alpha);
  assert(E.t().dot(alpha)==0); // follows from $\delta*\alpha=\alpha$
  // too hard to do: |assert(level_a(E,Cayley_shift(theta,ww),n_alpha)==0);|

  const WeightInvolution th_1 = // theta upstairs minus identity
    rd.root_reflection(n_alpha) * i_tab.matrix(E.tw) - 1;

  Weight diff = matreduc::find_solution(th_1,alpha); // MINUS $\sigma$

  E.set_l(E.l()+alpha_v*(tf_alpha/2));
  E.set_tau(E.tau()+diff*alpha_v.dot(E.tau()));
  return signs_differ(E0,E) ? -1 : 1;
}

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
  , flipped_edges()
{
  BitMap fixed_points(block.size());

  { // compute |child_nr| and |parent_nr| tables
    weyl::Twist twist(orbits);

    if (twisted(kgb,0,delta,twist)==UndefKGB or
	twisted(dual_kgb,0,delta.transposed(),twist)==UndefKGB)
      return; // if one or other not delta-stable, leave |size==0| and quit

    for (BlockElt z=0; z<block.size(); ++z)
    {
      BlockElt tz = twisted(block,kgb,dual_kgb,z,delta,twist);
      if (tz==z)
	fixed_points.insert(z);
    }
  }

  complete_construction(fixed_points);
}

ext_block::ext_block // for an external twist
  (const InnerClass& G,
   const param_block& block, const KGB& kgb,
   const WeightInvolution& delta)
  : parent(block)
  , orbits(block.fold_orbits(delta))
  , folded(orbits,block.Dynkin())
  , d_delta(delta)
  , info()
  , data(orbits.size()) // create that many empty vectors
  , flipped_edges()
{
  BitMap fixed_points(block.size());

  { // compute the delta-fixed points of the block
    // the following is NOT |twist(orbits)|, which would be for subsystem
    weyl::Twist twist(fold_orbits(block.rootDatum(),delta));

    // test if twisting some block element lands in the same block
    if (twisted(block,kgb,0,delta,twist)==UndefBlock)
      return; // if block not delta-stable, leave |size==0| and quit

    for (BlockElt z=0; z<block.size(); ++z)
    {
      BlockElt tz = twisted(block,kgb,z,delta,twist);
      if (tz==z)
	fixed_points.insert(z);
    }
  }

  complete_construction(fixed_points);
}

void ext_block::complete_construction(const BitMap& fixed_points)
{
  unsigned int folded_rank = orbits.size();
  std::vector<BlockElt> child_nr(parent.size(),UndefBlock);
  std::vector<BlockElt> parent_nr(fixed_points.size());
  { KGBElt x=0;
    for (auto it=fixed_points.begin(); it(); ++it,++x)
    {
      parent_nr[x]=*it;
      child_nr[*it]=x;
    }
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

      if (is_descent(type)) // then set or check length from predecessor
      {
	if (info.back().length==0)
	  info.back().length=info[child_nr[link]].length+1;
	else
	  assert(info.back().length==info[child_nr[link]].length+1);
      }

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
} // |ext_block::ext_block|

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

BlockEltPair ext_block::Cayleys(weyl::Generator s, BlockElt n) const
{
  const DescValue type = descent_type(s,n);
  if (has_double_image(type))
    return data[s][n].links;
  if (is_complex(type) or is_like_compact(type) or is_like_nonparity(type))
    return {UndefBlock,UndefBlock};
  return {data[s][n].links.first,UndefBlock};
}



// check validity, by comparing with results found using extended parameters
bool check(const ext_block eb, const param_block& block)
{
  context ctxt (block.context(),eb.delta(),block.gamma());
  containers::sl_list<param> links;
  for (BlockElt n=0; n<eb.size(); ++n)
  { auto z=eb.z(n);
    for (weyl::Generator s=0; s<eb.rank(); ++s)
    { ext_gen p=eb.orbit(s);
      param E(ctxt,block.x(z),block.lambda_rho(z)); // re-init each iteration
      int sign; links.clear(); // output arguments for |star|
      auto tp = star(E,p,links,sign);
      if (tp!=eb.descent_type(s,n))
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
	  BlockElt m=eb.cross(s,n);
	  BlockElt cz = eb.z(m); // corresponding element of block
	  param F(ctxt,block.x(cz),block.lambda_rho(cz));
	  assert(same_standard_reps(*it,F));
	  if (sign!=sign_between(*it,F))
	    std::cout << "Flip at cross link " << unsigned{s} << " from " << z
		      << " to " << cz << '.' << std::endl;
	} break;
      case one_imaginary_single: case one_real_single:
      case two_imaginary_single_single: case two_real_single_single:
	{ assert(links.size()==2);
	  BlockElt m=eb.some_scent(s,n); // the unique (inverse) Cayley
	  BlockElt Cz = eb.z(m); // corresponding element of block
	  param F(ctxt,block.x(Cz),block.lambda_rho(Cz));
	  assert(same_standard_reps(*it,F));
	  if (sign!=sign_between(*it,F))
	    std::cout << "Flip at Cayley link " << unsigned{s} << " from " << z
		      << " to " << Cz << '.' << std::endl;
	  ++it;
	  m=eb.cross(s,n); BlockElt cz = eb.z(m);
	  param Fc(ctxt,block.x(cz),block.lambda_rho(cz));
	  assert(same_standard_reps(*it,Fc));
	  if (sign!=sign_between(*it,Fc))
	    std::cout << "Flip at cross link " << unsigned{s} << " from " << z
		      << " to " << cz << '.' << std::endl;
	} break;
      case two_semi_imaginary: case two_semi_real:
      case three_semi_imaginary: case three_real_semi:
      case three_imaginary_semi: case three_semi_real:
	{ assert(links.size()==1);
	  BlockElt m=eb.some_scent(s,n); // the unique (inverse) Cayley
	  BlockElt Cz = eb.z(m); // corresponding element of block
	  param F(ctxt,block.x(Cz),block.lambda_rho(Cz));
	  assert(same_standard_reps(*it,F));
	  if (sign!=sign_between(*it,F))
	    std::cout << "Flip at Cayley link " << unsigned{s} << " from " << z
		      << " to " << Cz << '.' << std::endl;
	} break;
      case one_imaginary_pair_fixed: case one_real_pair_fixed:
      case two_imaginary_double_double: case two_real_double_double:
	{ assert(links.size()==2);
	  BlockEltPair m=eb.Cayleys(s,n);
	  BlockElt Cz0 = eb.z(m.first); BlockElt Cz1= eb.z(m.second);
	  param F0(ctxt,block.x(Cz0),block.lambda_rho(Cz0));
	  param F1(ctxt,block.x(Cz1),block.lambda_rho(Cz1));
	  bool straight=same_standard_reps(*it,F0);
	  if (straight)
	    assert(same_standard_reps(*std::next(it),F1));
	  else // it must be crossed
	    assert(same_standard_reps(*it,F1) and
		   same_standard_reps(*std::next(it),F0));
	  unsigned iF0 = straight ? 0 : 1; // |links| index that pairs with |F0|
	  if (sign!=sign_between(*std::next(it,iF0),F0))
	    std::cout << "Flip at Cayley link " << unsigned{s} << " from " << z
		      << " to " << Cz0 << '.' << std::endl;
	  if (sign!=sign_between(*std::next(it,1-iF0),F1))
	    std::cout << "Flip at Cayley link " << unsigned{s} << " from " << z
		      << " to " << Cz1 << '.' << std::endl;
	} break;
      case two_imaginary_single_double_fixed: case two_real_single_double_fixed:
	{ assert(links.size()==3);
	  BlockEltPair m=eb.Cayleys(s,n);
	  BlockElt Cz0 = eb.z(m.first); BlockElt Cz1= eb.z(m.second);
	  param F0(ctxt,block.x(Cz0),block.lambda_rho(Cz0));
	  param F1(ctxt,block.x(Cz1),block.lambda_rho(Cz1));
	  bool straight=same_standard_reps(*it,F0);
	  if (straight)
	    assert(same_standard_reps(*std::next(it),F1));
	  else // it must be crossed
	    assert(same_standard_reps(*it,F1) and
		   same_standard_reps(*std::next(it),F0));
	  unsigned iF0 = straight ? 0 : 1; // |links| index that pairs with |F0|
	  if (sign!=sign_between(*std::next(it,iF0),F0))
	    std::cout << "Flip at Cayley link " << unsigned{s} << " from " << z
		      << " to " << Cz0 << '.' << std::endl;
	  if (sign!=sign_between(*std::next(it,1-iF0),F1))
	    std::cout << "Flip at Cayley link " << unsigned{s} << " from " << z
		      << " to " << Cz1 << '.' << std::endl;
	  // check the false cross link
	  std::advance(it,2);
	  m = eb.Cayleys(s,m.first);
	  BlockElt fcz = eb.z(m.first==n ? m.second : m.first);
	  assert(fcz==block.cross(p.s0,z) and fcz==block.cross(p.s1,z));
	  param F(ctxt,block.x(fcz),block.lambda_rho(fcz));
	  assert(same_standard_reps(*it,F));
	  if (sign!=sign_between(*it,F))
	    std::cout << "Flip at false cross link " << unsigned{s}
		      << " from " << z << " to " << fcz << '.' << std::endl;
	} break;
      } // |switch(tp)|
    } // |for(s)|
  } // |for(n)|
  return true;
}

} // |namespace ext_block|

} // |namespace atlas|
