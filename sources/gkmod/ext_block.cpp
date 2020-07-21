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
#include "common_blocks.h"
#include "repr.h"
#include "ext_kl.h"

/*
  For an extended group, the block structure is more complicated than an
  ordinary block, because each link in fact represents a local part of the
  parent block structure that links a twist-fixed block element to another
  twist-fixed element, via intermediate elements that are non twist-fixed.
 */

namespace atlas {

namespace ext_block {

  // Declaration of a local function
WeylWord fixed_conjugate_simple (const context& c, RootNbr& alpha);

  // Function definitions

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

// for these ascent types a link may remain undefined if at edge of partial block
bool might_be_uncertain(DescValue v) // which might make type itself uncertain
{
  static const unsigned long mask =
      1ul << one_complex_ascent // type itself cannot be wrong here
    | 1ul << two_complex_ascent // maybe |two_semi_imaginary|
    | 1ul << three_complex_ascent // maybe |three_semi_imaginary|
    | 1ul << one_imaginary_pair_fixed          // maybe switched
    | 1ul << two_imaginary_single_double_fixed // maybe switched
    | 1ul << two_imaginary_double_double // type itself cannot be wrong here
    | 1ul << three_imaginary_semi        // type itself cannot be wrong here
    | 1ul << three_semi_imaginary;       // type itself cannot be wrong here

  return (1ul << v & mask) != 0; // whether |v| is one of the above
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

    // semi cases do not record their (trivial) cross action
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
} // |link_count|

// number of links that are ascents or descents (not real/imaginary cross)
unsigned int scent_count(DescValue v)
{ return has_double_image(v) ? 2 : link_count(v)==0 ? 0 : 1; }

// find |n| such that |z(n)>=zz| (for |zz==parent.size()| returns |n==size()|)
BlockElt ext_block::element(BlockElt zz) const
{
  BlockElt min=0, max=size();
  // loop invariant: |(min==0 or z(min-1)<zz) and (max>=size() or z(max)>=zz)|
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

// compute the components in a default extended parameter for parameter |sr|
void set_default_extended
(const Rep_context& rc, const StandardRepr& sr, const WeightInvolution& delta,
   Weight& lambda_rho, Weight& tau, Coweight& l, Coweight& t)
{
  const auto& kgb = rc.kgb(); const auto x=sr.x();
  WeightInvolution theta = rc.inner_class().matrix(kgb.involution(x));

  lambda_rho=rc.lambda_rho(sr);
  tau=matreduc::find_solution(1-theta,(delta-1)*lambda_rho);
  l=ell(kgb,x);
  t=matreduc::find_solution(theta.transposed()+1,(delta-1).right_prod(l));
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
    , integr_datum(integrality_datum(rc.root_datum(),gamma))
    , sub(SubSystem::integral(rc.root_datum(),gamma))
    , pi_delta(rc.root_datum().rootPermutation(d_delta))
    , delta_fixed_roots(fixed_points(pi_delta))
    , twist()
    , lambda_shifts (integr_datum.semisimpleRank())
    , l_shifts (integr_datum.semisimpleRank())
{
  const RootDatum& rd = rc.root_datum();
  assert(is_dominant_ratweight(rd,d_gamma)); // this is a class invariant

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

bool context::is_very_complex (InvolutionNbr theta, RootNbr alpha) const
{ const auto& i_tab = inner_class().involution_table();
  const auto& rd = root_datum();
  assert (rd.is_posroot(alpha)); // this is a precondition
  auto image = i_tab.root_involution(theta,alpha);
  make_positive(rd,image);
  return image!=alpha and image!=delta_of(alpha);
}

Weight context::to_simple_shift
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ const InvolutionTable& i_tab = inner_class().involution_table();
  S &= (i_tab.real_roots(theta) ^i_tab.real_roots(theta_p));
  return root_sum(root_datum(),S);
}

/*
  For the conjugation to simple scenario, we compute a set |pos_neg| of
  positive roots that become negative under an element of $W^\delta$ that
  makes the integrally-simple root(s) in question simple. The function
  |shift_flip| computes from this set, and the involutions at both ends of the
  link in the block, whether an additional flip is to be added to the link.

  This comes from an action of |delta| on a certain top wedge product of
  root spaces, and the formula below tells whether that action is by $-1$.
*/
bool context::shift_flip
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ S.andnot(delta_fixed()); // $\delta$-fixed roots won't contribute

  unsigned count=0; // will count 2-element |delta|-orbit elements
  for (auto it=S.begin(); it(); ++it)
    if (is_very_complex(theta,*it) != is_very_complex(theta_p,*it) and
	not root_datum().sumIsRoot(*it,delta_of(*it)))
      ++count;

  assert(count%2==0); // since |pos_to_neg| is supposed to be $\delta$-stable
  return count%4!=0; // whether number of 2-element orbits (themselves) is odd
}

// auxiliary function to recognise local situation in |ext_block| construction
// assumes precomputed |fixed_points|; block may be partial or complete
// for partial blocks, some boundary elements may return an uncertain type
DescValue extended_type(const Block_base& block, BlockElt z, const ext_gen& p,
			BlockElt& link, const BitMap& fixed_points)
{
  switch (p.type)
  {
  case ext_gen::one:
    switch (block.descentValue(p.s0,z))
    {
    case DescentStatus::ComplexAscent:
      return link=block.cross(p.s0,z), one_complex_ascent; // maybe |UndefBlock|
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
	// now if |t==UndefBlock| we are uncertain; tentatively return "fixed"
	if (t==UndefBlock or fixed_points.isMember(t))
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
	if (t==UndefBlock)
	  return link=t, two_complex_ascent; // uncertain
	else if(t==block.cross(p.s1,z))
	  return link=t, two_semi_imaginary;
	return link=block.cross(p.s1,t), two_complex_ascent; // maybe undefined
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
      { const BlockElt t=block.cayley(p.s0,z).first; // unique Cayley ascent
	if (t==UndefBlock)
	  return link=t, two_imaginary_single_double_fixed; // uncertain
	if (block.descentValue(p.s1,t)==DescentStatus::ImaginaryTypeI)
	  return link=block.cayley(p.s1,t).first, two_imaginary_single_single;
	link=block.cayley(p.s1,t).first; // uncertain when |link==UndefBlock|
	return link==UndefBlock or fixed_points.isMember(link)
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
      { BlockElt tmp=block.cayley(p.s0,z).first;
	if (tmp==UndefBlock) // first Cayley ascent crossed edge of partial block
	  // since both our links have it as ascent, they are beyond edge too
	  return link=tmp, two_imaginary_double_double; // certain, unset |link|
	auto pair = block.cayley(p.s1,tmp);
	if (pair.first==UndefBlock or // then both components are |UndefBlock|
	    (not fixed_points.isMember(pair.first) and pair.second==UndefBlock))
	{ // try again with other pair of Cayley ascent by |s0| of |z|
	  if ((tmp=block.cayley(p.s0,z).second)==UndefBlock)
	    return link=tmp, two_imaginary_double_double; // crt, unset |link|
	  pair = block.cayley(p.s1,tmp); // try other pair
	  if (pair.first==UndefBlock)
	    return link=UndefBlock, two_imaginary_double_double;
	}
	link = fixed_points.isMember(pair.first) ? pair.first : pair.second;
	assert(link==UndefBlock or fixed_points.isMember(link));
	return two_imaginary_double_double;
      }
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
	if (t==UndefBlock)
	  return link=t, three_complex_ascent; // uncertain
	if (t==block.cross(p.s1,t))
	{
	  assert(block.descentValue(p.s1,t)==DescentStatus::ImaginaryTypeII);
	  link=block.cayley(p.s1,t).first;
	  if (link!=UndefBlock and not fixed_points.isMember(link))
	    link=block.cayley(p.s1,t).second, // choose the door without a goat
	      assert(link==UndefBlock or fixed_points.isMember(link));
	  return three_semi_imaginary; // certain, but link may be |UndefBlock|
	}
	link=block.cross(p.s1,t);
	if (link!=UndefBlock)
	  link=block.cross(p.s0,link), // continue up the third link
	    assert(link==UndefBlock or fixed_points.isMember(link));
	return three_complex_ascent; // certain, but link may be |UndefBlock|
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
      { const BlockElt t=block.cayley(p.s0,z).first;
	if (t==UndefBlock)
	  return link=t, three_imaginary_semi; // certain, but with unset |link|
	link=block.cross(p.s1,t); // could be |UndefBlock|; then leave it
	if (link!=UndefBlock) // then |block.cayley(p.s1,z)| is defined, and
	  assert(fixed_points.isMember(link) and
		 link==block.cross(p.s0,block.cayley(p.s1,z).first));
	return three_imaginary_semi;
      }
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


/* Try to conjugate |alpha| by product of folded-generators for the (full)
   root system of |c| to a simple root, and return the left-conjugating word
   that was applied. This may fail, if after some conjugation one ends up with
   the long root of a nontrivially folded A2 subsystem (in which case there
   cannot be any solution because |alpha| is fixed by the involution but none
   of the simple roots in its component of the root system is). In this case
   |alpha| is left as that non simple root, and the result conjugates to it.
 */
WeylWord fixed_conjugate_simple (const context& ctxt, RootNbr& alpha)
{ const RootDatum& rd = ctxt.root_datum();

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





ext_block::~ext_block() = default;

ext_block::ext_block // for external twist; old style blocks
  (const InnerClass& G,
   const Block& block,
   const KGB& kgb, const KGB& dual_kgb, // all are needed
   const WeightInvolution& delta)
  : parent(block)
  , orbits(fold_orbits(G.rootDatum(),delta))
  , info()
  , data(orbits.size()) // create that many empty vectors
  , l_start(parent.length(parent.size()-1)+2,0)
  , pol_hash(nullptr)
  , KL_ptr(nullptr)
  , folded_diagram(block.Dynkin().folded(orbits))
  , delta(delta)
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

// create tables defining extended block structure
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

  for (BlockElt n=0; n<parent_nr.size(); ++n) // |n| is index in extended block
  {
    BlockElt z=parent_nr[n]; // |z| is index in parent block
    info.push_back(elt_info(z));
    for (weyl::Generator oi=0; oi<orbits.size(); ++oi) // |oi|: orbit index
    {
      const weyl::Generator s = orbits[oi].s0, t=orbits[oi].s1;
      BlockElt link, second = UndefBlock; // these index parent block elements
      DescValue type = extended_type(parent,z,orbits[oi],link,fixed_points);
      data[oi].push_back(block_fields(type)); // create entry

      if (is_like_compact(type) or is_like_nonparity(type))
	continue; // leave both link fields |UndefBlock| in those cases

      // now maybe set |second|, depending on case
      switch (type)
      {
      default: break;

	// cases where second link is cross neighbour for |s|
      case one_imaginary_single:
      case one_real_single:
	second = parent.cross(s,z);
	break;

	// cases where second link is second Cayley image, cross of |link|
      case one_real_pair_fixed:
      case one_imaginary_pair_fixed:
	if (link!=UndefBlock)
	  second = parent.cross(s,link);
	break;

	// cases where second link is double cross neighbour for |s| of |z|
      case two_imaginary_single_single:
      case two_real_single_single:
	{
	  BlockElt tmp = parent.cross(s,z);
	  if (tmp!=UndefBlock)
	  {
	    second = parent.cross(t,tmp);
	    assert((tmp=parent.cross(t,z))==UndefBlock or
		   second==parent.cross(s,tmp));
	  }
	  else if ((tmp=parent.cross(t,z))!=UndefBlock) // try alternative route
	    second=parent.cross(s,tmp);
	  else if (type==two_real_single_single and // try to pass from above
		   (tmp=parent.cross(t,parent.inverseCayley(s,z).first))
		    !=UndefBlock)
	  {
	    auto pair = parent.cayley(s,tmp);
	    second = pair.first!=UndefBlock and fixed_points.isMember(pair.first)
	      ? pair.first : pair.second;
	  }
	  // for |two_imaginary_single_single| a nasty case remains: though the
	  // double cross neighbour may be in the block, all intermediates could
	  // be absent. Then leave |second| undefined, hoping it is never needed
	}
	break;

	// pair-to-pair link cases; second link is second Cayley, and sort
      case two_imaginary_single_double_fixed:
      case two_real_single_double_fixed:
	if (link!=UndefBlock)
	{
	  second = parent.cross(s,link); // second Cayley image is cross of first
	  assert(second==parent.cross(t,link)); // (for either generator)
	  if (link>second) // to make sure ordering is same for a twin pair
	    std::swap(link,second); // we order both by block number (for now)
	}
	break;

	// cases where second link is second Cayley image, double cross of |link|
      case two_imaginary_double_double:
      case two_real_double_double:
	if (link!=UndefBlock)
	{
	  BlockElt tmp = parent.cross(s,link);
	  if (tmp!=UndefBlock)
	  {
	    second = parent.cross(t,tmp);
	    assert((tmp=parent.cross(t,link))==UndefBlock or
		   second==parent.cross(s,tmp));
	  }
	  else if ((tmp=parent.cross(t,link))!=UndefBlock)
	    second = parent.cross(s,tmp);
	  else if ((tmp=parent.cayley(s,z).second)!=UndefBlock)
	  { // in |two_imaginary_double_double| case, try again from above
	    auto pair = parent.cayley(t,tmp);
	    if (pair.first!=UndefBlock)
	      second = fixed_points.isMember(pair.first)
		? pair.first : pair.second;
	  }
	  if (link>second) // make sure single |UndefBlock| is ranked second
	      std::swap(link,second); // by ordering by block number
	}
	break;
      } // |switch(type)|

      // enter translations of |link| and |second| to child block numbering
      BlockEltPair& dest = data[oi].back().links;
      if (link!=UndefBlock)
	dest.first=child_nr[link];
      if (second!=UndefBlock)
	dest.second = child_nr[second];
    }
  } // |for(n)|
} // |ext_block::complete_construction|

// we compute $\max\{l: l_start[l]\leq n\}$, i.e. |upper_bound(,,n)-1|
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

BlockElt ext_block::Cayley(weyl::Generator s, BlockElt n) const
{
  return  is_complex(descent_type(s,n)) ? UndefBlock : data[s][n].links.first;
}

BlockEltPair ext_block::Cayleys(weyl::Generator s, BlockElt n) const
{
  assert(has_double_image(descent_type(s,n)));
  return data[s][n].links;
}


// flag those among |orbits| whose elements are flagged in |gen_set|
// here |gen_set| is supposed a union of orbits, so (any flagged => all flagged)
RankFlags reduce_to(const ext_gens& orbits, RankFlags gen_set)
{ RankFlags result;
  for (weyl::Generator s=0; s<orbits.size(); ++s)
    result.set(s,gen_set[orbits[s].s0]); // set whether |s0| element in |gen_set|
  return result;
}


weyl::Generator
ext_block::first_descent_among(RankFlags singular_orbits, BlockElt y) const
{ auto it=singular_orbits.begin();
  for (; it(); ++it)
    if (is_descent(descent_type(*it,y)))
      return *it;

  return rank();
}

const ext_kl::KL_table& ext_block::kl_table
  (BlockElt limit, ext_KL_hash_Table* pol_hash)
{
  if (KL_ptr==nullptr)
    KL_ptr.reset(new ext_kl::KL_table(*this,pol_hash));
  KL_ptr->fill_columns(limit);
  return *KL_ptr;
}

void ext_block::swallow // integrate older partial block, using |embed| mapping
  (ext_block&& sub, const BlockEltList& embed)
{
  if (sub.KL_ptr!=nullptr)
  {
    if (KL_ptr==nullptr)
      KL_ptr.reset(new ext_kl::KL_table(*this,pol_hash));

    // restrict the translation vector to its |sub| elements (|delta|-fixed ones)
    BlockEltList eblock_embed; eblock_embed.reserve(sub.size());
    for (BlockElt x=0; x<sub.size(); ++x)
    {
      eblock_embed.push_back(element(embed[sub.z(x)]));
      assert(eblock_embed.back()!=size()); // check that lookup succeeded
    }
    KL_ptr->swallow(std::move(*sub.KL_ptr),eblock_embed);
  }
}


// reduce matrix to rows for extended block elements without singular descents
// the other rows are not removed, but the result lists the rows to retain
template<typename C> // matrix coefficient type (signed)
containers::simple_list<BlockElt> // returns list of elements selected
  ext_block::condense(matrix::Matrix<C>& M, RankFlags sing_orbs) const
{
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

bool ext_block::add_neighbours
  (containers::sl_list<BlockElt>& dst, weyl::Generator s, BlockElt n) const
{
  const BlockEltPair& links = data[s][n].links;
  if (links.first==UndefBlock)
    return 0<link_count(descent_type(s,n)); // whether too short
  dst.push_back(links.first);
  if (links.second==UndefBlock)
    return 1<link_count(descent_type(s,n)); // whether too short;
  dst.push_back(links.second);
  return false; // success, |link_count| cannot exceed 2
}

bool check_quadratic (const ext_block& b, weyl::Generator s, BlockElt x0)
{ containers::sl_list<BlockElt> l;

  if (b.add_neighbours(l,s,x0))
    return true;

  if (l.empty()) // compact or nonparity cases, there is nothing to check
    return true;

  // check symmetry of link signs
  for (BlockElt y : l)
    if (b.epsilon(s,x0,y)!=b.epsilon(s,y,x0))
      return false;

  auto tp = b.descent_type(s,x0);
  if (l.size()==1) // cases without cycles; we're done
    return true;
  assert(l.size()==2);
  BlockElt y0 = l.front(); l.pop_front();
  BlockElt y1 = l.front(); l.pop_front();

  if (has_quadruple(tp))
  { if (b.add_neighbours(l,s,y0))
      return true;
    if (x0==l.front())
      l.pop_front(); // so that |l.front()| won't be |x0|
    BlockElt x1 = l.front();
    assert (x0!=x1);
    return b.epsilon(s,x0,y0)*b.epsilon(s,x0,y1) != // negative product here!
      b.epsilon(s,y0,x1)*b.epsilon(s,y1,x1);
  }
  else
    return b.epsilon(s,x0,y0)*b.epsilon(s,x0,y1)==b.epsilon(s,y0,y1);
}

bool check_braid
  (const ext_block& b, weyl::Generator s, weyl::Generator t, BlockElt x,
   BitMap& cluster)
{
  if (s==t)
    return true;
  static const unsigned int cox_entry[] = {2, 3, 4, 6};
  unsigned int len = cox_entry[b.folded_diagram.edge_multiplicity(s,t)];

  BitMap used(b.size());
  containers::queue<BlockElt> to_do { x };
  do
  {
    BlockElt z=to_do.front();
    to_do.pop();
    used.insert(z);
    containers::sl_list<BlockElt> l;
    if (b.add_neighbours(l,s,z) or b.add_neighbours(l,t,z))
      return true;
    for (BlockElt y : l)
      if (not used.isMember(y))
	to_do.push(y);
  }
  while (not to_do.empty());

  unsigned int n=used.size();
  matrix::Matrix<Pol> Ts(n,n,Pol()), Tt(n,n,Pol());

   unsigned int j=0; // track index of |y|
  for (const BlockElt y : used)
  {
    set(Ts,j,j, b.T_coef(s,y,y)-Pol(1));
    set(Tt,j,j, b.T_coef(t,y,y)-Pol(1));
    containers::sl_list<BlockElt> l;
    if (b.add_neighbours(l,s,y))
      return true;
    for (BlockElt z : l)
      if (used.isMember(z))
	set(Ts,used.position(z),j, b.T_coef(s,z,y));
    l.clear();
    if (b.add_neighbours(l,t,y))
      return true;
    for (BlockElt z : l)
      if (used.isMember(z))
	set(Tt,used.position(z),j, b.T_coef(t,z,y));
    ++j; // keep |j| in phase with |y|
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
  if (verbose and not success)
  {
    //    std::cout << "success: " << success << std::endl;
    show_mat(std::cout,Ts,s);
    std::cout << std::endl;
    show_mat(std::cout,Tt,t);
  }
  return success;
} // |check_braid|

template containers::simple_list<BlockElt> ext_block::condense
  (matrix::Matrix<int>& M, RankFlags sing_orbs) const;

} // |namespace ext_block|

} // |namespace atlas|
