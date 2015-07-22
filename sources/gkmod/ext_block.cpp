/*
  This is ext_block.cpp

  Copyright (C) 2013-2014 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_block.h"

#include <cassert>
#include <vector>

#include "complexredgp.h"
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
    |  1ul << two_real_double_double      | 1ul << two_imaginary_double_double
    |  1ul << two_imaginary_single_double | 1ul << two_real_single_double;

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
    1ul << two_imaginary_single_double | 1ul << two_real_single_double;

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

// find element |n| such that |z(n)>=zz|
BlockElt extended_block::element(BlockElt zz) const
{
  BlockElt n=0;
  while (n<size() and z(n)<zz)
    ++n;
  return n;
}

context::context
  (repr::Rep_context& rc, WeightInvolution delta, const RatWeight& gamma)
    : d_rc(rc)
    , d_delta(std::move(delta)), d_gamma(gamma)
    , d_g(rc.kgb().base_grading_vector()-rho_check(rc.rootDatum()))
    , integr_datum(integrality_datum(rc.rootDatum(),gamma))
    , sub(SubSystem::integral(rc.rootDatum(),gamma))
{}

// compute |bgv-(bgv+t_bits)*(1+theta)/2 == (bgv-t_bits-(bgv+t_bits)*theta)/2|
Coweight ell (const KGB& kgb, KGBElt x)
{ const RatCoweight bgv=kgb.base_grading_vector();
  const Coweight t_bits (kgb.torus_part(x)); // lift from bitvector to vector
  RatCoweight result=bgv+t_bits;
  twisted_act(kgb.complexGroup(),result,kgb.involution(x));
  result.numerator().negate();
  result += bgv;
  result -= t_bits;
  result /= 2;
  result.normalize();
  assert(result.denominator()==1);
  return Coweight(result.numerator().begin(),result.numerator().end());
}

param::param (const context& ec, const StandardRepr& sr)
  : ctxt(ec)
  , tw(ec.rc().kgb().involution(sr.x()))
  , l(ell(ec.realGroup().kgb(),sr.x()))
  , lambda_rho(ec.rc().lambda_rho(sr))
  , tau(matreduc::find_solution(1-theta(),(delta()-1)*lambda_rho))
  , t(matreduc::find_solution(1-theta().transposed(),(delta()-1).right_prod(l)))
{}

param::param (const context& ec, KGBElt x, const Weight& lambda_rho)
  : ctxt(ec)
  , tw(ec.realGroup().kgb().involution(x))
  , l(ell(ec.realGroup().kgb(),x))
  , lambda_rho(lambda_rho)
  , tau(matreduc::find_solution(1-theta(),(delta()-1)*lambda_rho))
  , t(matreduc::find_solution(1-theta().transposed(),(delta()-1).right_prod(l)))
{}

bool in_L_image(Weight beta,WeightInvolution&& A)
{ int_Matrix L,R;
  auto inv_fact = matreduc::diagonalise(A,L,R);
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
  auto inv_fact = matreduc::diagonalise(A,L,R);
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
  { if (&E.ctxt.complexGroup()!=&F.ctxt.complexGroup())
      throw std::runtime_error
	("Comparing extended parameters from different inner classes");
    if (E.delta()!=F.delta()
	or E.ctxt.g()!=F.ctxt.g()
	or E.ctxt.gamma()!=F.ctxt.gamma())
      return false;
  } // otherwise the might still be a match, so fall through
  return E.theta()==F.theta()
    and in_R_image(E.theta()+1,E.l-F.l)
    and in_L_image(E.lambda_rho-F.lambda_rho,E.theta()+1);
}

KGBElt x(const param& E)
{ TorusPart tp(E.l);
  TitsElt a(E.ctxt.complexGroup().titsGroup(),tp,E.tw);
  return E.rc().kgb().lookup(a);
}

// this implements (comparison using) the formula from Propodition 16 in
// "Parameters for twisted repressentations" (with $\delta-1=-(1-\delta)$
// the relation is symmetric in |E|, |F|, although not obviously so
bool signs_differ (const param& E, const param& F)
{
  const WeightInvolution& delta = E.delta();
  Weight kappa1=E.tau, kappa2=F.tau;
  kappa1 -= delta*kappa1;
  kappa2 -= delta*kappa2;
  int i_exp = E.l.dot(kappa1) - F.l.dot(kappa2);
  assert (i_exp%2==0);
  int n1_exp = (F.l-E.l).dot(E.tau) + F.t.dot(F.lambda_rho-E.lambda_rho);
  return (i_exp/2+n1_exp)%2!=0;
}

bool sign_differs_with_one_of (const param& E, const param& F1, const param& F2)
{ return
    same_standard_reps(E,F1) ? signs_differ(E,F1)
    : same_standard_reps(E,F2) ? signs_differ(E,F2)
    : throw std::runtime_error("Neither candidate has same standard repn");
}

bool is_default (const param& E)
{ return not signs_differ(E,param(E.ctxt,x(E),E.lambda_rho)); }

bool extended_block::toggle_edge(BlockElt x,BlockElt y, bool verbose)
{
  x = element(x); y=element(y);
  assert (x!=UndefBlock and y!=UndefBlock);
  BlockEltPair p= x<y ? std::make_pair(x,y) : std::make_pair(y,x);
  std::pair<std::set<BlockEltPair>::iterator,bool>
    inserted = flipped_edges.insert(p);
  if (not inserted.second)
    flipped_edges.erase(inserted.first);

  if (verbose)
    std::cerr << (inserted.second ? "Set" : "Unset") << " edge ("
	      << z(p.first) << ',' << z(p.second) << ')' << std::endl;

  return inserted.second;
}

void extended_block::order_quad
  (BlockElt x,BlockElt y, BlockElt p, BlockElt q, int s, bool verbose)
{
  x = element(x); y=element(y); // decipher user friendly numbering
  p = element(p); q=element(q);
  assert (x!=UndefBlock and y!=UndefBlock and p!=UndefBlock and q!=UndefBlock);
  assert (descent_type(s-1,x)==two_imaginary_single_double);
  assert (descent_type(s-1,y)==two_imaginary_single_double);
  assert (descent_type(s-1,p)==two_real_single_double);
  assert (descent_type(s-1,q)==two_real_single_double);
  const BlockEltPair xy(x,y);
  const BlockEltPair pq(p,q);
  std::vector<block_fields>& data_s = data[s-1];
  data_s[x].links = data_s[y].links = pq;
  data_s[p].links = data_s[q].links = xy;
  if (verbose)
    std::cerr << "Ordering (" << z(x) << ',' << z(y) << ';'
	      << z(p) << ',' << z(q) << ") for generator " << s << std::endl;
}

BlockElt extended_block::cross(weyl::Generator s, BlockElt n) const
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
  case two_real_nonparity: case two_imaginary_compact:
  case three_real_nonparity: case three_imaginary_compact:

    // double valued Cayleys also have trivial cross actions
  case one_real_pair_fixed: case one_real_pair_switched:
  case one_imaginary_pair_fixed: case one_imaginary_pair_switched:
  case two_real_double_double: case two_imaginary_double_double:

    // cases with back-and-forth cross actions
  case two_semi_imaginary: case two_semi_real:
  case two_imaginary_single_double: case two_real_single_double:
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
} // |extended_block::cross|

BlockElt extended_block::Cayley(weyl::Generator s, BlockElt n) const
{
  const DescValue type = descent_type(s,n);
  return
    is_descent(type) or is_complex(type) ? UndefBlock : data[s][n].links.first;
}

BlockElt extended_block::inverse_Cayley(weyl::Generator s, BlockElt n) const
{
  const DescValue type = descent_type(s,n);
  return
    not is_descent(type) or is_complex(type) ? UndefBlock
    : data[s][n].links.first;
}

BlockElt extended_block::some_scent(weyl::Generator s, BlockElt n) const
{
  const BlockElt c = data[s][n].links.first;
  assert(c!=UndefBlock);
  return c;
}

void extended_block::add_neighbours
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

BlockEltPair extended_block::Cayleys(weyl::Generator s, BlockElt n) const
{
  const DescValue type = descent_type(s,n);
  BlockEltPair result(UndefBlock,UndefBlock);
  if (not is_descent(type) and not is_complex(type))
  {
    result.first = data[s][n].links.first;
    if (has_double_image(type))
      result.second = data[s][n].links.second;
  }
  return result;
}

BlockEltPair
extended_block::inverse_Cayleys(weyl::Generator s, BlockElt n) const
{
  const DescValue type = descent_type(s,n);
  BlockEltPair result(UndefBlock,UndefBlock);
  if (is_descent(type) and not is_complex(type))
  {
    result.first = data[s][n].links.first;
    if (has_double_image(type))
      result.second = data[s][n].links.second;
  }
  return result;
}

// whether link for |s| from |x| to |y| has a signe flip attached
int extended_block::epsilon(weyl::Generator s, BlockElt x, BlockElt y ) const
{
  BlockEltPair p= x<y ? std::make_pair(x,y) : std::make_pair(y,x);
  int sign = flipped_edges.count(p)==0 ? 1 : -1;

  // each 2i12/21r21 quadruple has one negative sign not using |flipped_edges|
  if (has_quadruple(descent_type(s,x)) and
      data[s][x].links.second==y and data[s][y].links.second==x)
    sign = -sign; // it is between second elements in both pairs of the quad

  return sign;
}

BlockEltList extended_block::down_set(BlockElt n) const
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

DescValue extended_type(const Block_base& block, BlockElt z, ext_gen p,
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
	? two_imaginary_single_single : two_imaginary_single_double;
    case DescentStatus::RealTypeII:
      link=block.inverseCayley(p.s0,z).first;
      if (link==UndefBlock)
	return two_real_single_single; // really just a guess
      link=block.inverseCayley(p.s1,link).first;
      if (link==UndefBlock)
	return two_real_single_single; // really just a guess
      assert(block.Hermitian_dual(link)==link);
      return block.descentValue(p.s0,link)==DescentStatus::ImaginaryTypeII
	? two_real_single_single : two_real_single_double;
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

param complex_cross(ext_gen p, const param& E)
{ const RootDatum& rd = E.rc().rootDatum();
  const InvolutionTable& i_tab = E.rc().complexGroup().involution_table();
  auto &tW = E.rc().twistedWeylGroup(); // caution: |p| refers to integr. datum

  TwistedInvolution tw=E.tw;
  InvolutionNbr theta = i_tab.nr(tw);
  const RatWeight gamma_rho = E.ctxt.gamma() - rho(rd);
  RatWeight gamma_lambda =  gamma_rho - E.lambda_rho;
  auto& ga_la_num = gamma_lambda.numerator();
  Weight rho_r_shift = rd.twoRho(i_tab.real_roots(theta));

  const RatCoweight g_rho_check = E.ctxt.g() - rho_check(rd);
  RatCoweight torus_factor =  g_rho_check - E.l;
  auto& tf_num = torus_factor.numerator();
  Coweight dual_rho_im_shift = rd.dual_twoRho(i_tab.imaginary_roots(theta));

  Weight tau=E.tau;
  Coweight t=E.t;
  const RootDatum& id = E.ctxt.id();
  for (unsigned i=p.w_tau.size(); i-->0; )
  { weyl::Generator s=p.w_tau[i]; // generator for integrality datum
    tW.twistedConjugate(E.ctxt.subsys().reflection(s),tw);
    id.reflect(s,ga_la_num);
    id.reflect(s,rho_r_shift);
    id.reflect(s,tau);
    id.coreflect(tf_num,s);
    id.coreflect(t,s);
    id.coreflect(dual_rho_im_shift,s);
  }
  RatWeight lr_ratvec = (gamma_rho - gamma_lambda).normalize();
  assert(lr_ratvec.denominator()==1);
  Weight lambda_rho(lr_ratvec.numerator().begin(),
		    lr_ratvec.numerator().end()); // convert to |Weight|
  rho_r_shift -= rd.twoRho(i_tab.real_roots(i_tab.nr(tw)));
  rho_r_shift/=2; // now it is just a sum of (real) roots
  Weight tau_corr = ((E.ctxt.delta()-1)*rho_r_shift)/2; // hope it divides

  RatWeight l_ratvec = (g_rho_check - torus_factor).normalize();
  assert(l_ratvec.denominator()==1);
  Coweight l(l_ratvec.numerator().begin(), l_ratvec.numerator().end());
  dual_rho_im_shift -= rd.dual_twoRho(i_tab.imaginary_roots(i_tab.nr(tw)));
  dual_rho_im_shift/=2; // now it is just a sum of (imaginary) coroots
  Coweight t_corr = ((E.ctxt.delta()-1).right_prod(dual_rho_im_shift))/2;

  return param(E.ctxt, tw,
	       lambda_rho-rho_r_shift, tau+tau_corr,
	       l-dual_rho_im_shift, t+t_corr);
}


WeylWord fixed_conjugate_simple (const context& ctxt, RootNbr& alpha)
{ const RootDatum& rd = ctxt.complexGroup().rootDatum();
  std::vector<weyl::Generator> delta (rd.semisimpleRank());
  std::vector<bool> is_length_3 (delta.size());
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
  while (true) // termination condition after definition of |s| (readability)
  {
    weyl::Generator s =
      rd.descent_set(alpha).andnot(rd.ascent_set(delta_alpha)).firstBit();
    if (alpha==rd.simpleRootNbr(s) or
	(is_length_3[s] and //"sum of swapped non-commuting roots" case:
	 rd.simple_reflected_root(s,alpha)==rd.simpleRootNbr(delta[s])))
      break;
    result.push_back(s);
    rd.simple_reflect_root(s,alpha);
    rd.simple_reflect_root(delta[s],delta_alpha);
    if (delta[s]!=s) // second gernerator for cases of length 2,3
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
}

DescValue type (const param& E, ext_gen p, std::vector<param>& links)
{
  DescValue result;
  const TwistedWeylGroup& tW = E.rc().twistedWeylGroup();
  const InvolutionTable& i_tab = E.rc().complexGroup().involution_table();
  const RootDatum& rd = E.rc().rootDatum();
  const RootDatum& integr_datum = E.ctxt.id();
  const SubSystem& subs = E.ctxt.subsys();
  const InvolutionNbr theta = i_tab.nr(E.tw);
  const WeightInvolution delta_1 = E.ctxt.delta()-1;
  switch (p.type)
  {
  case ext_gen::one:
    { const Weight& alpha = integr_datum.root(p.s0);
      const Coweight& alpha_v = integr_datum.coroot(p.s0);
      const RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);

      if (theta_alpha==n_alpha) // imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g() - E.l).dot(alpha)-rd.level(n_alpha);
	if (tf_alpha%2!=0) // then $\alpha$ is compact
	  return one_imaginary_compact; // quit here, do not collect \$200

	// noncompact case
	const TwistedInvolution new_tw= tW.prod(E.tw,subs.reflection(p.s0));
	const WeightInvolution& th_1 = i_tab.matrix(new_tw)-1;

	int tau_coef = alpha_v.dot(E.tau); // $\tau_\alpha$ of table 2

	// try to make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const Weight rho_r_shift = repr::Cayley_shift
	  (E.rc().complexGroup(),i_tab.nr(new_tw),ww);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	Weight first; // maybe a root with |(1-delta)*first==alpha|
	if (rd.is_simple_root(alpha_simple))
	  first = Weight(rd.rank(),0); // effectively not used in this case
	else
	{
	  --tau_coef; // parity change and decrease both relevant
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
	      matreduc::find_solution(th_1,alpha); // solutions equivalent ?
	  links.push_back(param // Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho + first + rho_r_shift,
			   E.tau + diff*tau_coef - first,
			   E.l+alpha_v*(tf_alpha/2),
			   E.t
			   ));
	  links.push_back(param // cross link
			  (E.ctxt,E.tw,E.lambda_rho,E.tau, E.l+alpha_v, E.t));
	} // end of type 1 case
	else
	{ // type 2; now we need to distinguish 1i2f and 1i2s

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
	    assert(not matreduc::has_solution
		   (th_1, delta_1*(E.lambda_rho+rho_r_shift)));
	    return one_imaginary_pair_switched; // case 1i2s
	  }
	  result = one_imaginary_pair_fixed;  // what remains is case 1i2f

	  links.push_back(param // first Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho + first + rho_r_shift,
			   E.tau - alpha*(tau_coef/2) - first,
			   E.l+alpha_v*(tf_alpha/2),
			   E.t
			   ));
	  links.push_back(param // second Cayley link
			  (E.ctxt,new_tw,
			   E.lambda_rho + first + rho_r_shift + alpha,
			   E.tau - alpha*(tau_coef/2) - first,
			   E.l+alpha_v*(tf_alpha/2),
			   E.t
			   ));
	} // end of type 2 case
      } // end of imaginary case
      else if (theta_alpha==rd.rootMinus(n_alpha)) // real case
      { // the folowing is coherehnt with using |Cayley_shift|, but easier
	const int parity_n = (E.ctxt.gamma() - E.lambda_rho).dot(alpha_v)
	  - rd.colevel(n_alpha) // <alpha_v,rho>
	  + rd.twoRho(i_tab.real_roots(theta)).dot(alpha_v)/2;
	if (parity_n%2==0) // nonparity
	   return one_real_nonparity; // no link added here

	const WeightInvolution& th_1 = i_tab.matrix(E.tw)-1; // at more split
	bool type1 = matreduc::has_solution(th_1,alpha);

	const TwistedInvolution new_tw= tW.prod(E.tw,subs.reflection(p.s0));

	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);

	const Weight rho_r_shift =
	  repr::Cayley_shift(E.rc().complexGroup(),theta,ww);
	assert((delta_1*rho_r_shift).isZero()); // since $ww\in W^\delta$

	Weight new_lambda_rho = E.lambda_rho-rho_r_shift;
	Weight tau_correction;
	if (rd.is_simple_root(alpha_simple))
	  tau_correction = Weight(rd.rank(),0); // no correction needed here
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
	{ int d = // amount by which $\gamma-\tilde\lambda$ needs correction
	    (E.ctxt.gamma()-new_lambda_rho).dot(alpha_v)- rd.colevel(n_alpha);
	  assert(d%2==0); // parity condition says this
	  new_lambda_rho -= alpha*(d/2); // project to correct $\tilde\lambda$
	}

	int t_alpha = E.t.dot(alpha);
	if (type1)
	  { // now distinguish 1r1f and 1r1s
	    if (t_alpha%2!=0)
	      return one_real_pair_switched;
	    result = one_real_pair_fixed; // what remains is case 1r1f
	    Coweight new_t = E.t - alpha_v*(t_alpha/2);
	    links.push_back(param // first Cayley link
			    (E.ctxt,new_tw,
			     new_lambda_rho,
			     E.tau + tau_correction,
			     E.l,
			     new_t
			     ));
	    links.push_back(param // second Cayley link
			    (E.ctxt,new_tw,
			     new_lambda_rho,
			     E.tau + tau_correction,
			     E.l + alpha_v,
			     new_t
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
			     E.tau + tau_correction,
			     E.l,
			     E.t - diff*t_alpha
			     ));
	    links.push_back(param // cross link
	        (E.ctxt,E.tw, E.lambda_rho+alpha ,E.tau,E.l,E.t));
	  }
      }
      else // complex case
      { result = rd.is_posroot(theta_alpha)
	  ? one_complex_ascent : one_complex_descent ;
	links.push_back(complex_cross(p,E));
      }
    }
    break;
  case ext_gen::two:
    { const Weight& alpha = integr_datum.root(p.s0);
      const Coweight& alpha_v = integr_datum.coroot(p.s0);
      RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
      const Weight& beta = integr_datum.root(p.s1);
      const Coweight& beta_v = integr_datum.coroot(p.s1);
      // RootNbr n_beta = subs.parent_nr_simple(p.s1);
      // RootNbr theta_beta = i_tab.root_involution(theta,n_beta);
      if (theta_alpha==n_alpha) //imaginary case
      { const RatCoweight grading_coweight = E.ctxt.g() - E.l -
	  RatCoweight(rd.dual_twoRho()
		      -rd.dual_twoRho(i_tab.imaginary_roots(theta)),2);
	if (grading_coweight.dot(alpha)%2==0) // compact
	  result=two_imaginary_compact;
	else // noncompact
	{ CoweightInvolution theta1_t = i_tab.matrix(theta).transposed()+1;
	  if (matreduc::has_solution(theta1_t, alpha_v)) // type 2
	  {
	    result = two_imaginary_double_double;
	  }
	  else // type 1
	    if (matreduc::has_solution(theta1_t, alpha_v+beta_v))
	    { // |cross(alpha,E) == cross(beta,E)|, so case 2i12
	      result = two_imaginary_single_double;
	    }
	    else // case 2i11
	    {
	      result = two_imaginary_single_single;
	      links.push_back(param
		(E.ctxt,E.tw,E.lambda_rho,E.tau, E.l+alpha_v+beta_v, E.t));
	    }
	}
      }
      else if (theta_alpha==rd.rootMinus(n_alpha)) // real case
      { const RatWeight parity_weight = E.ctxt.gamma() - E.lambda_rho -
	  RatWeight(rd.twoRho()-rd.twoRho(i_tab.real_roots(theta)),2);
	if (parity_weight.dot(alpha_v)%2==0) // nonparity
	  result = two_real_nonparity; // no link added here
	else // parity
	{ WeightInvolution theta_1 = i_tab.matrix(theta)-1;
	  if (matreduc::has_solution(theta_1,alpha))
	  { // type 1
	    result = two_real_double_double;
	  }
	  else // real type 2
	    if (matreduc::has_solution(theta_1,alpha+beta))
	    { // |cross(alpha,E) == cross(beta,E)|, so case 2r21
	      result = two_real_single_double;
	    }
	    else // case 2r22
	    {
	      result = two_real_single_single;
	      links.push_back(param
	        (E.ctxt,E.tw, E.lambda_rho+alpha+beta ,E.tau,E.l,E.t));
	    }
	}
      }
      else // complex case
      { TwistedInvolution tw = E.tw;
	tW.twistedConjugate(tw,p.s0);
	tW.twistedConjugate(tw,p.s1);
	if (tw==E.tw) // twisted commutation with |s0.s1|
	{
	  result = rd.is_posroot(theta_alpha)
	    ? two_semi_imaginary : two_semi_real;
	}
	else // twisted non-commutation with |s0.s1|
	{
	  result = rd.is_posroot(theta_alpha)
	    ? two_complex_ascent : two_complex_descent;
	  links.push_back(complex_cross(p,E));
	}
      }
    }
    break;
  case ext_gen::three:
    break;
  }
  return result;
}

extended_block::extended_block
  (const Block_base& block,const TwistedWeylGroup& W)
  : parent(block)
  , tW(W)
  , folded(blocks::folded(block.Dynkin(),block.fold_orbits()))
  , info()
  , data(parent.folded_rank())
  , l_start(parent.length(parent.size()-1)+2)
  , flipped_edges()
{
  unsigned int folded_rank = data.size();
  if (folded_rank==0 or parent.Hermitian_dual(0)==UndefBlock)
    return; // block not stable under twist, so leave extended block empty

  std::vector<BlockElt> child_nr(parent.size(),UndefBlock);
  std::vector<BlockElt> parent_nr;

  { // compute |child_nr| and |parent_nr| tables, and the |l_start| vector
    size_t cur_len=0; l_start[cur_len]=0;
    for (BlockElt z=0; z<parent.size(); ++z)
      if (parent.Hermitian_dual(z)==z)
      {
	child_nr[z]=parent_nr.size();
	parent_nr.push_back(z);
	while (cur_len<parent.length(z)) // for new length level(s) reached
	  l_start[++cur_len]=child_nr[z]; // mark element as first of |cur_len|
      }
    assert(cur_len+1<l_start.size());
    l_start[++cur_len]=parent.size(); // makes |l_start[length(...)+1]| legal
  }

  info.reserve(parent_nr.size());
  for (weyl::Generator s=0; s<folded_rank; ++s)
    data[s].reserve(parent_nr.size());

  for (BlockElt n=0; n<parent_nr.size(); ++n)
  {
    BlockElt z=parent_nr[n];
    info.push_back(elt_info(z));
    info.back().length = parent.length(z);
    for (weyl::Generator s=0; s<folded_rank; ++s)
    {
      BlockElt link, second = UndefBlock;
      DescValue type = extended_type(block,z,parent.orbit(s),link);
      data[s].push_back(block_fields(type)); // create entry for block element
      if (link==UndefBlock)
	continue; // |s| done for imaginary compact and real nonparity cases

      switch (type) // maybe set |second|, depending on case
      {
      default: break;

      case one_imaginary_single:
      case one_real_single: // in these cases: parent cross neighbour for |s|
	second = parent.cross(parent.orbit(s).s0,z);
	break;

      case one_real_pair_fixed:
      case one_imaginary_pair_fixed: // in these cases: second Cayley image
	second = parent.cross(parent.orbit(s).s0,link);
	break;

      case two_imaginary_single_single:
      case two_real_single_single: // here: double cross neighbour for |s|
	if (parent.cross(parent.orbit(s).s0,z)==UndefBlock)
	  if (parent.cross(parent.orbit(s).s1,z)==UndefBlock)
	    continue; // can't get there, let's hope it doesn't exist
	  else second = parent.cross(parent.orbit(s).s0,
				     parent.cross(parent.orbit(s).s1,z));
	else
	{
	  second = parent.cross(parent.orbit(s).s1,
				parent.cross(parent.orbit(s).s0,z));
	  if (parent.cross(parent.orbit(s).s1,z)!=UndefBlock)
	    assert(parent.cross(parent.orbit(s).s0,
				parent.cross(parent.orbit(s).s1,z))==second);
	}
	break;

      case two_imaginary_single_double:
      case two_real_single_double: // find second Cayley image, which is
	second = parent.cross(parent.orbit(s).s0,link); // parent cross link
	assert(parent.cross(parent.orbit(s).s1,link)==second); // (either gen)
	if (link>second) // to make sure ordering is same for a twin pair
	  std::swap(link,second); // we order both by block number (for now)
	break;

      case two_imaginary_double_double:
      case two_real_double_double: // find second Cayley image
	if (parent.cross(parent.orbit(s).s0,link)==UndefBlock)
	  if (parent.cross(parent.orbit(s).s1,link)==UndefBlock)
	    continue; // can't get there, let's hope it doesn't exist
	  else second = parent.cross(parent.orbit(s).s0,
				     parent.cross(parent.orbit(s).s1,link));
	else
	{
	  second = parent.cross(parent.orbit(s).s1,
				parent.cross(parent.orbit(s).s0,link));
	  if (parent.cross(parent.orbit(s).s1,link)!=UndefBlock)
	    assert(parent.cross(parent.orbit(s).s0,
				parent.cross(parent.orbit(s).s1,link))
		   ==second);
	}
	break;
      } // |switch(type)|

      // enter translations of |link| and |second| to child block numbering
      BlockEltPair& dest = data[s].back().links;
      dest.first=child_nr[link];
      if (second!=UndefBlock)
	dest.second = child_nr[second];
    }
  } // |for(n)|
} // |extended_block::extended_block|

// coefficient of neighbour |sx| for $s$ in action $(T_s+1)*a_x$
Pol extended_block::T_coef(weyl::Generator s, BlockElt sx, BlockElt x) const
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
    BlockElt y = inverse_Cayley(s,x); // pass via this element for signs
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
} // |extended_block::T_coef|

} // |namespace ext_block|

} // |namespace atlas|
