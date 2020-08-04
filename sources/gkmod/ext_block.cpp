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
#include "ext_kl.h"

/*
  For an extended group, the block structure is more complicated than an
  ordinary block, because each link in fact represents a local part of the
  parent block structure that links a twist-fixed block element to another
  twist-fixed element, via intermediate elements that are non twist-fixed.
 */

namespace atlas {

namespace ext_block {

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
WeylWord fixed_conjugate_simple
  (const repr::Ext_common_context& ctxt, RootNbr& alpha)
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


//			|ext_block| methods



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

ext_block::ext_block
  (const blocks::common_block& block, const WeightInvolution& delta,
   ext_KL_hash_Table* pol_hash)
  : parent(block)
  , orbits(block.fold_orbits(delta))
  , info()
  , data(orbits.size()) // create that many empty vectors
  , l_start(parent.length(parent.size()-1)+2,0)
  , pol_hash(pol_hash)
  , KL_ptr(nullptr)
  , folded_diagram(block.Dynkin().folded(orbits))
  , delta(delta)
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
  if (not tune_signs(block))
    throw std::runtime_error("Failure detected in extended block construction");

} // |ext_block::ext_block|, from a |common_block|

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


void validate(const ext_param& E)
{
#ifndef NDEBUG // make sure this is a no-op when debugging is disabled
  const auto& i_tab = E.ctxt.inner_class().involution_table();
  const auto& theta = i_tab.matrix(E.tw);
  const auto& delta = E.ctxt.delta();
  assert(delta*theta==theta*delta);
  const Weight diff = E.gamma_lambda.integer_diff<int>(delta*E.gamma_lambda);
  assert(diff == (1-theta)*E.tau);
  assert((delta-1).right_prod(E.l)==(theta+1).right_prod(E.t));
  assert(((E.ctxt.g_rho_check()-E.l)*(1-theta)).numerator().isZero());
  assert(((theta+1)*E.gamma_lambda).numerator().isZero());
#endif
}
/*
  An auxiliary routine to compute extended parameters across complex links.
  The situation is complicated by the fact that the cross action is by a
  generator of the folded integral system, so we need to expand it first into
  a product of |length<=3| integral generators, and then have those generators
  act on the components of |E|. For the purpose of changing |E.tw| we further
  develop those generators into reflection words for the full root datum, but
  the reflection action on the other components can be done more directly.

  However, although the integral generators are complex, the action of those
  reflection words need not be purely complex, which implies that the effect
  on the |gamma_lambda| and |l| components are not purely reflections. The
  difference with respect to pure reflection action can be computed comparing
  |rho_r| values (half sums of positive real roots in the full system) with
  the reflected image of that value taken as the starting point, respectively
  (for |l|) the same thing with |rho_check_imaginary|. This is done by the
  "correction" terms below.
 */
ext_param complex_cross(const repr::Ext_common_context& ctxt,
		      const ext_gen& p, ext_param E) // by-value for |E|, modified
{ const RootDatum& rd = E.rc().root_datum();
  const RootDatum& id = ctxt.id();
  const InvolutionTable& i_tab = E.rc().inner_class().involution_table();
  auto &tW = E.rc().twisted_Weyl_group(); // caution: |p| refers to integr. datum
  const SubSystem& subs=ctxt.subsys();

  InvolutionNbr theta = i_tab.nr(E.tw);
  const RootNbrSet& theta_real_roots = i_tab.real_roots(theta);
  Weight rho_r_shift = rd.twoRho(theta_real_roots);
  Coweight dual_rho_im_shift = rd.dual_twoRho(i_tab.imaginary_roots(theta));

  for (unsigned i=p.w_kappa.size(); i-->0; ) // at most 3 letters, right-to-left
  { weyl::Generator s=p.w_kappa[i]; // generator for integrality datum
    tW.twistedConjugate(subs.reflection(s),E.tw);
    id.simple_reflect(s,E.gamma_lambda.numerator());
    id.simple_reflect(s,rho_r_shift);
    id.simple_reflect(s,E.tau);
    id.simple_coreflect(E.l,s,ctxt.l_shift(s));
    id.simple_coreflect(dual_rho_im_shift,s);
    id.simple_coreflect(E.t,s);
  }

  InvolutionNbr new_theta = i_tab.nr(E.tw);
  const RootNbrSet& new_theta_real_roots = i_tab.real_roots(new_theta);
  rho_r_shift -= rd.twoRho(new_theta_real_roots);
  rho_r_shift/=2; // now it is just a sum of (real) roots
  E.gamma_lambda += rho_r_shift;

  assert(ctxt.delta()*rho_r_shift==rho_r_shift); // diff of $\delta$-fixed

  dual_rho_im_shift -= rd.dual_twoRho(i_tab.imaginary_roots(new_theta));
  dual_rho_im_shift/=2; // now it is just a sum of (imaginary) coroots
  E.l -= dual_rho_im_shift;

  assert(ctxt.delta().right_prod(dual_rho_im_shift)==dual_rho_im_shift);
  validate(E);

  RootNbr alpha_simple = subs.parent_nr_simple(p.s0);
  const WeylWord to_simple = fixed_conjugate_simple(ctxt,alpha_simple);
  // by symmetry by $\delta$, |to_simple| conjugates $\delta(\alpha)$ to simple:
  assert(p.length()==1 or rd.is_simple_root(rd.permuted_root(to_simple,
				                subs.parent_nr_simple(p.s1))));
  // apply flip for $\delta$ acting on root set for |to_simple|, as elsewhere
  E.flip(ctxt.shift_flip(theta,new_theta,pos_to_neg(rd,to_simple)));

  E.flip(p.length()==2); // to parallel the 2i,2r flips

  return E;
} // |complex_cross|

// whether |E| and |F| lie over equivalent |StandardRepr| values
bool same_standard_reps (const ext_param& E, const ext_param& F)
{
  if (&E.ctxt!=&F.ctxt)
  { if (&E.ctxt.inner_class()!=&F.ctxt.inner_class())
      throw std::runtime_error
	("Comparing extended parameters from different inner classes");
    if (E.delta()!=F.delta() or E.ctxt.g_rho_check()!=F.ctxt.g_rho_check())
      return false;
    // otherwise (contexts differ, but agree on vitals), match could still occur
  } // so fall through
  return E.theta()==F.theta()
    and in_R_image(E.theta()+1,E.l-F.l)
    and in_L_image(E.gamma_lambda.integer_diff<int>(F.gamma_lambda),E.theta()-1);
} // |same_standard_reps|

// this implements (comparison using) the formula from Proposition 16 in
// "Parameters for twisted repressentations" (with $\delta-1 = -(1-\delta)$
// the relation is symmetric in |E|, |F|, although not obviously so
bool same_sign (const ext_param& E, const ext_param& F)
{
  assert(same_standard_reps(E,F));
  const WeightInvolution& delta = E.delta();
  Weight kappa1=E.tau, kappa2=F.tau;
  kappa1 -= delta*kappa1;
  kappa2 -= delta*kappa2;
  int i_exp = E.l.dot(kappa1) - F.l.dot(kappa2);
  assert(i_exp%2==0);
  int n1_exp =
    (F.l-E.l).dot(E.tau) + (E.gamma_lambda-F.gamma_lambda).dot(F.t);
  return ((i_exp/2+n1_exp)%2==0)!=(E.is_flipped()!=F.is_flipped());
}

inline bool is_default (const ext_param& E)
  { return same_sign(E,ext_param(E.ctxt,E.x(),E.gamma_lambda)); }


/*
  for real Cayley transforms, one will add $\rho_r$ to |E.gamma_lambda|
  before projecting it parallel to |alpha| so as to make |alpha_v| vanish on
  |E.gamma_lambda|. Here we compute from |E.gamma_lambda|, corrected by
  that |shift|, the multiple of $\alpha/2$ that such a projection would
  subtract from |E.gamma_lambda|).
*/
int level_a (const ext_param& E, const Weight& shift, RootNbr alpha)
{
  const RootDatum& rd = E.rc().root_datum();
  return (E.gamma_lambda + shift).dot(rd.coroot(alpha));
}



void z_align (const ext_param& E, ext_param& F, bool extra_flip)
{ assert(E.t==F.t); // we require preparing |t| upstairs to get this
  int d = E.l.dot((E.delta()-1)*E.tau) - F.l.dot((F.delta()-1)*F.tau);
  assert(d%2==0);
  F.flipped = E.is_flipped()^(d%4!=0)^extra_flip; // XOR 3 Booleans into |F|
}

/*
  In some cases, notably 2i12, one or both of the Cayley transforms may
  involve (in the simple root case) a change |mu| in |lambda_rho| that does
  not necessarily satisfy |t.dot(mu)==0|. In such cases the previous version
  of |z_align| is insufficient, and we should include a contribution from the
  second term of the formula for |z|. But retrieving |mu| from the parameters
  |E| and |F| themselves is complicated by the posssible contribution from
  |Cayley_shift|, which contribution should be ignored; however at the place
  of call the value of |mu| is explicitly available, so we ask here to pass
  |t.dot(mu)| as third argument |t_mu|.
 */
void z_align (const ext_param& E, ext_param& F, bool extra_flip, int t_mu)
{ z_align(E,F,extra_flip^(t_mu%2!=0)); }


// compute type of |p| for |E|, and export adjacent |ext_param| values in |links|
DescValue star (const repr::Ext_common_context& ctxt,
		const ext_param& E, const ext_gen& p,
		containers::sl_list<ext_param>& links)
{
  ext_param E0=E; // a copy of |E| that might be modified below to "normalise"
  DescValue result;

  const TwistedWeylGroup& tW = E.rc().twisted_Weyl_group();
  const InnerClass& ic = E.rc().inner_class();
  const InvolutionTable& i_tab = ic.involution_table();
  const RootDatum& rd = E.rc().root_datum();
  const RootDatum& integr_datum = ctxt.id();
  const SubSystem& subs = ctxt.subsys();
  const InvolutionNbr theta = i_tab.nr(E.tw);
  switch (p.type)
  {
  default: assert(false);
    result=one_complex_ascent; // shut up "maybe uninitialsed" warning
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
	const TwistedInvolution new_tw= tW.prod(subs.reflection(p.s0),E.tw);
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1; // upstairs

	int tau_coef = alpha_v.dot(E.tau); // take $\tau_\alpha$ of table 2

	// try to make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);

	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(E.t.dot(alpha)==0); // follows from $\delta*\alpha=\alpha$

	// now separate cases; based on type 1 or 2 first
	if (matreduc::has_solution(th_1,alpha))
	{ // type 1, so extended type is 1i1
	  result = one_imaginary_single;

	  /* if imaginary integrally simple |alpha| were sum of two simple
	     roots these cannot be imaginary (hence integral), so they must be
	     interchanged by $\theta$, but then we are in type 2
	   */
	  assert(rd.is_simple_root(alpha_simple));

	  Weight diff = // called $-\sigma$ in table 2 of [Ptr] (NOTE MINUS)
	      matreduc::find_solution(th_1,alpha); // solutions are equivalent

	  ext_param F(E.ctxt,new_tw,
		    E.gamma_lambda- rho_r_shift,
		    E0.tau+diff*tau_coef,
		    E.l+alpha_v*(tf_alpha/2), E.t);

	  E0.l = tf_alpha%4==0 ? F.l+alpha_v : F.l; // for cross
	  assert(not same_standard_reps(E,E0));
	  z_align(E,F,flipped);
	  z_align(F,E0,flipped);

	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E0)); // cross link
	} // end of 1i1 case
	else
	{ // imaginary type 2; now we need to distinguish 1i2f and 1i2s

	  auto new_gamma_lambda = E.gamma_lambda; auto new_tau = E.tau;
	  if (not rd.is_simple_root(alpha_simple))
	  {
	    --tau_coef; // the parity change and decrease are both relevant
	    weyl::Generator s = // first switched root index
	      rd.find_descent(alpha_simple);
	    RootNbr first = // corresponding root summand, conjugated back
	      rd.permuted_root(rd.simpleRootNbr(s),ww);
	    assert(alpha == (E.ctxt.delta()+1)*rd.root(first));
	    new_gamma_lambda -= rd.root(first);
	    new_tau -= rd.root(first);
	  }

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
#ifndef NDEBUG
	    auto ratv = (E.ctxt.delta()-1)*(E.gamma_lambda-rho_r_shift);
	    ratv.normalize(); // having a multiple of the weight here won't do!
	    const Weight target
	      { ratv.numerator().begin(),ratv.numerator().end() };
	    assert(not matreduc::has_solution (th_1,target));
#endif
	    return one_imaginary_pair_switched; // case 1i2s
	  }
	  result = one_imaginary_pair_fixed;  // what remains is case 1i2f


	  ext_param F0(E.ctxt,new_tw,
		     new_gamma_lambda - rho_r_shift,
		     new_tau - alpha*(tau_coef/2),
		     E.l + alpha_v*(tf_alpha/2), E.t);
	  ext_param F1(E.ctxt,new_tw,
		     F0.gamma_lambda - alpha, F0.tau, F0.l, E.t);

	  if (not rd.is_simple_root(alpha_simple))
	    flipped = not flipped;

	  z_align(E,F0,flipped);
	  z_align(E,F1,flipped);
	  links.push_back(std::move(F0)); // Cayley link
	  links.push_back(std::move(F1)); // Cayley link
	} // end of type 2 case
      } // end of length 1 imaginary case

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 1 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s0),E.tw);

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	const int t_alpha = E.t.dot(alpha);
	if (matreduc::has_solution(i_tab.matrix(E.tw)-1,alpha)) // then type 1
	{ // length 1 type 1 real case
	  // for the same reason as |alpha| must be simple in case 1i1, we have:
	  assert(rd.is_simple_root(alpha_simple));
	  const int level = level_a(E,rho_r_shift,n_alpha);
	  if (level%2!=0) // nonparity
	    return one_real_nonparity; // case 1rn, no link added here

	  // now distinguish 1r1f and 1r1s
	  if (t_alpha%2!=0)
	    return one_real_pair_switched; // case 1r1s
	  result = one_real_pair_fixed; // what remains is case 1r1f

	  const RatWeight new_gamma_lambda =
	    E.gamma_lambda + rho_r_shift - alpha*(level/2);
	  assert(new_gamma_lambda.dot(alpha_v)==0); // check effect of |level_a|

	  E0.t -= alpha_v*(t_alpha/2);
	  assert(same_sign(E,E0)); // since only |t| changes

	  ext_param F0(E.ctxt,new_tw,new_gamma_lambda, E.tau, E.l          , E0.t);
	  ext_param F1(E.ctxt,new_tw,new_gamma_lambda, E.tau, E.l + alpha_v, E0.t);

	  z_align(E0,F0,flipped);
	  z_align(E0,F1,flipped);
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley

	} // end of 1r1f case
	else // type 2
	{ // length 1 type 2 real

	  int level = level_a(E,rho_r_shift,n_alpha);

	  Weight new_tau=E.tau; // maybe modified below
	  if (not rd.is_simple_root(alpha_simple))
	  { // adapt to integrality based change of lambda
	    RootNbr first = // one of summand roots in |alpha==first+second|
	      rd.permuted_root(rd.simpleRootNbr(rd.find_descent(alpha_simple)),
			       ww);
	    assert(alpha == (E.ctxt.delta()+1)*rd.root(first));
	    assert(i_tab.real_roots(theta).isMember(first));

	    rho_r_shift += rd.root(first); // non delta-fixed contribution
	    ++level; // the change in |rho_r_shift| augments its $\alpha$-level

	    // now we must add $d$ to $\tau$ with $(1-\theta')d=(1-\delta)*a0$
	    // where |a0=rd.root(first)|
	    // since $\theta'*a0 = a1 = \delta*a_0$, we can take $d=a0$
	    assert((i_tab.matrix(new_tw))*rd.root(first) ==
		   (E.ctxt.delta()      )*rd.root(first));
	    new_tau += rd.root(first);

	    flipped = not flipped; // flip the Cayley links in this case
	  }
	  if (level%2!=0) // nonparity
	    return one_real_nonparity; // case 1rn, no link added here
	  result = one_real_single; // case 1r2

	  const RatWeight new_gamma_lambda =
	    E.gamma_lambda + rho_r_shift - alpha*(level/2);
	  assert(new_gamma_lambda.dot(alpha_v)==0); // check effect of |level_a|

	  const Coweight diff = // called $s$ in table 2 of [Ptr]
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v);
	  E0.t -= diff*t_alpha;
	  assert(same_sign(E,E0)); // since only |t| changes

	  ext_param E1 = E0; // for cross neighbour; share updated value of |t|
	  E1.gamma_lambda -= alpha;
	  assert(not same_standard_reps(E0,E1));

	  ext_param F(E.ctxt,new_tw, new_gamma_lambda, new_tau, E.l, E0.t);


	  z_align(E0,F,flipped);
	  z_align(F,E1,flipped);

	  // since |z_align| ignores |gamma_lambda|, we must have equal flips:
	  assert(E0.is_flipped()==E1.is_flipped());

	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E1)); // cross link
	} // end of 1r2 case
      }
      else // length 1 complex case
      { result = rd.is_posroot(theta_alpha)
	  ? one_complex_ascent : one_complex_descent ;
	links.push_back(complex_cross(ctxt,p,E));
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
	// make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 2

	flipped = not flipped; // because of wedge correction for 2i/2r cases

	int at = alpha_v.dot(E.tau); int bt = beta_v.dot(E.tau);
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1;

	if (matreduc::has_solution(th_1,alpha)) // then type 2i11
	{ result = two_imaginary_single_single;
	  const Weight sigma = matreduc::find_solution(th_1,alpha*at+beta*bt);
	  ext_param F (E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift,  E.tau + sigma,
		     E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t);

	  E0.l += alpha_v+beta_v;
	  z_align(E,F,flipped); // no 4th arg, since |E.gamma_lambda| unchanged
	  z_align(F,E0,flipped);
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
	  const Coweight new_l = E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2);

	  // first Cayley link |F0| will be the one that does not need |sigma|
	  ext_param F0(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift - alpha*m, new_tau0,
		     new_l, E.t);
	  ext_param F1(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift - alpha*mm, E.tau + sigma,
		     new_l, E.t);

	  int t_alpha=E.t.dot(alpha);
	  z_align(E,F0,flipped,m*t_alpha);
	  z_align(E,F1,flipped,mm*t_alpha);
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

	  ext_param F0(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift - alpha*m,
		     E.tau - alpha*((at+m)/2) - beta*((bt-m)/2),
		     E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t);
	  ext_param F1(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift - alpha*(1-m) - beta,
		     E.tau - alpha*((at-m)/2) - beta*((bt+m)/2),
		     F0.l,E.t);

	  int ta = E.t.dot(alpha), tb=E.t.dot(beta);
	  z_align(E,F0,flipped,ta*m);
	  z_align(E,F1,flipped,ta*(1-m)+tb);
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley
	} // end type 2i22 case
      }

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 2 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	flipped = not flipped; // because of wedge correction for 2i/2r cases

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return two_real_nonparity; // no link added here

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	WeightInvolution theta_1 = i_tab.matrix(theta)-1; // upstairs

	const RatWeight new_gamma_lambda = E.gamma_lambda + rho_r_shift
	  - alpha*(a_level/2) - beta*(b_level/2);

	int ta = E.t.dot(alpha); int tb = E.t.dot(beta);
	ext_param E1=E; // another modifiable copy, like |E0|

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

	  ext_param F0(E.ctxt, new_tw,
		     new_gamma_lambda,E.tau, E.l+alpha_v*m, E0.t);
	  ext_param F1(E.ctxt, new_tw,
		     new_gamma_lambda,E.tau, E.l+alpha_v*(1-m)+beta_v,E1.t);

	  z_align(E0,F0,flipped,m*((b_level-a_level)/2));
	  z_align(E1,F1,flipped,m*((a_level-b_level)/2));

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
	  ext_param F0(E.ctxt, new_tw,
		     new_gamma_lambda, E.tau, E.l+alpha_v*m, E0.t);
	  ext_param F1(E.ctxt, new_tw,
		     new_gamma_lambda, E.tau, E.l+alpha_v*mm, E1.t);

	  z_align(E0,F0,flipped,m *((b_level-a_level)/2));
	  z_align(E1,F1,flipped,mm*((b_level-a_level)/2));
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

	  E1.gamma_lambda -= alpha+beta;
	  E1.t = E0.t; // cross action, keeps adaptation of |t| to |F| below
	  assert(not same_standard_reps(E0,E1));

	  ext_param F(E.ctxt, new_tw, new_gamma_lambda, E.tau, E.l, E0.t);

	  z_align(E0,F,flipped); // no 4th arg, as |E.t.dot(alpha)==0| etc.
	  z_align(F,E1,flipped);
	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E1)); // cross link
	} // end of case 2r22
      }
      else // length 2 complex case
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha != (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // twisted non-commutation with |s0.s1|
	  result = ascent ? two_complex_ascent : two_complex_descent;
	  links.push_back(complex_cross(ctxt,p,E));
	}
	else if (ascent)
	{ // twisted commutation with |s0.s1|: 2Ci
	  result = two_semi_imaginary;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p = i_tab.nr(new_tw); // upstairs
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = ctxt.shift_flip(theta,theta_p,S);
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  // downstairs cross by |ww| only has imaginary and complex steps, so
	  // $\alpha_v.(\gamma-\lambda_\rho)$ is unchanged across |ww|
	  const int f = E.gamma_lambda.dot(alpha_v);
	  // number of times $\alpha$ is subtracted from $\gamma-\lambda$


	  const RatWeight new_gamma_lambda =
	    E.gamma_lambda - alpha*f - rho_r_shift;
	  // both $\gamma-\lambda$ and $\tau$ get $f*alpha$ subtracted by
	  // $\alpha$-reflection; adapt $\tau$ for vanishing $1-\delta$ image
	  const Weight new_tau = rd.reflection(n_alpha,E.tau) + alpha*f;

	  // but |dual_v| needs correction by |ell_shift|
	  const int dual_f = (E.ctxt.g_rho_check() - E.l).dot(alpha);

	  const Coweight new_l = E.l + alpha_v*dual_f;
	  const Coweight new_t =
	    rd.coreflection(E.t,n_alpha) - alpha_v*dual_f;
	  ext_param F (E.ctxt, new_tw, new_gamma_lambda, new_tau, new_l, new_t,
		     E.is_flipped()!=flipped);

	  // do extra conditional flip for 2Ci case
	  int ab_tau = (alpha_v+beta_v).dot(E.tau);
	  assert (ab_tau%2==0);
	  F.flip((ab_tau*dual_f)%4!=0);
	  links.push_back(std::move(F));  // "Cayley" link
	} // end of 2Ci case
	else // twisted commutation with |s0.s1|, and not |ascent|: 2Cr
	{ result = two_semi_real;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p=i_tab.nr(new_tw);
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = ctxt.shift_flip(theta,theta_p,S);
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  const int f = level_a(E,rho_r_shift,n_alpha);

	  const RatWeight new_gamma_lambda = // \emph{reflect} parallel to alpha
	    E.gamma_lambda + rho_r_shift - alpha*f;
	  const Weight new_tau = rd.reflection(n_alpha,E.tau) - alpha*f;

	  const int dual_f = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	  const Coweight new_l = E.l + alpha_v*dual_f;
	  const Coweight new_t =
	    rd.coreflection(E.t,n_alpha) + alpha_v*dual_f;

	  ext_param F (E.ctxt, new_tw, new_gamma_lambda, new_tau, new_l, new_t,
		   E.is_flipped()!=flipped);

	  // do extra conditional flip for 2Cr case
	  int t_ab = E.t.dot(beta-alpha);
	  assert(t_ab%2==0);
	  F.flip((t_ab * (f+alpha_v.dot(E.tau)))%4!=0);
	  links.push_back(std::move(F));  // "Cayley" link
	} // end of 2Cr case
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
      const Coweight& beta_v = integr_datum.simpleCoroot(p.s1);

      RootNbr n_kappa =integr_datum.simple_reflected_root
	 (p.s1, integr_datum.simpleRootNbr(p.s0));
      WeylWord s_kappa = subs.reflection(integr_datum.posRootIndex(n_kappa));

      const Weight& kappa = integr_datum.root(n_kappa);
      assert (kappa==alpha+beta);
      const Coweight& kappa_v = integr_datum.coroot(n_kappa);

      const Weight beta_alpha = beta - alpha;

      const TwistedInvolution new_tw = tW.prod(s_kappa,E.tw); // when applicable

      if (theta_alpha==n_alpha) // length 3 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	int tf_beta = (E.ctxt.g_rho_check() - E.l).dot(beta);
	assert((tf_alpha-tf_beta)%2==0); // same compactness
	if (tf_alpha%2!=0) // then $\alpha$ and $\beta$ are compact
	  return three_imaginary_compact;

	// noncompact case
	result = three_imaginary_semi; // this is the 3i case

	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 3

	E0.tau -= alpha*kappa_v.dot(E.tau); // make |kappa_v.dot(E.tau)==0|
	E0.l += alpha_v*(tf_alpha+tf_beta);
	E0.t += (beta_v-alpha_v)*((tf_alpha+tf_beta)/2);

	ext_param F(E.ctxt, new_tw,
		  E0.gamma_lambda-rho_r_shift, E0.tau, E0.l, E0.t);

	flipped = not flipped; // January unsurprise for 3i: delta acts by -1
	z_align(E0,F,flipped^not same_sign(E,E0));

	links.push_back(std::move(F)); // Cayley link
      } // end of 3i case
      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 3 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return three_real_nonparity; // no link added here

	// parity case
	result = three_real_semi; // this is the 3r case

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	const RatWeight new_gamma_lambda = // make level for |kappa| 0
	  E.gamma_lambda+rho_r_shift - alpha*(a_level+b_level); // even multiple

	E0.t -= alpha_v*kappa.dot(E.t); // makes |E.t.dot(kappa)==0|
	E0.gamma_lambda -= alpha*(a_level+b_level); // even multiple of |alpha|
	E0.tau += beta_alpha*((a_level+b_level)/2);
	assert(same_sign(E,E0)); // neither |t| change nor 2*real_root matter
	assert(E0.gamma_lambda+rho_r_shift==new_gamma_lambda);
	validate(E0);

	flipped = not flipped; // January unsurprise for 3r

	ext_param F(E.ctxt, new_tw, new_gamma_lambda,E0.tau,E0.l,E0.t);

	z_align(E0,F,flipped); // no 4th arg since |E.t.dot(kappa)==0|
	links.push_back(std::move(F)); // Cayley link
      } // end of 3r case
      else // length 3 complex case (one of 3Ci or 3Cr or 3C+/-)
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha == (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // reflection by |alpha+beta| twisted commutes with |E.tw|: 3Ci or 3Cr
	  result = ascent ? three_semi_imaginary : three_semi_real;

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p=i_tab.nr(new_tw);
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = ascent
	    ?  ctxt.to_simple_shift(theta,theta_p,S)
	    : -ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = ctxt.shift_flip(theta,theta_p,S);
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	  int dtf_alpha = E.gamma_lambda.dot(alpha_v);
	  RatWeight new_gamma_lambda = E.gamma_lambda - rho_r_shift; // for now

	  if (ascent) // 3Ci
	  {
	    if (dtf_alpha%2!=0)
	    { new_gamma_lambda -= beta_alpha;
	      E0.gamma_lambda -= beta_alpha;
	      E0.tau -= beta_alpha;
	    }
	    E0.l += kappa_v*tf_alpha;
	    E0.tau -= kappa*(kappa_v.dot(E.tau)/2);
	    validate(E0);
	    assert(E0.t.dot(kappa)==0);

	    ext_param F(E.ctxt,new_tw, new_gamma_lambda, E0.tau, E0.l, E0.t);

	    flipped = not flipped; // January unsurprise for 3Ci
	    z_align(E0,F, flipped^(not same_sign(E,E0)));
	    links.push_back(std::move(F)); // Cayley link
	  }
	  else // descent, so 3Cr
	  {
	    E0.gamma_lambda  -= kappa*dtf_alpha;
	    new_gamma_lambda -= kappa*dtf_alpha;

	    E0.t -= kappa_v*(kappa.dot(E.t)/2); // makes |E0.t.dot(kappa)==0|
	    if (tf_alpha%2!=0)
	    {
	      auto b_a = beta_v-alpha_v; // or |kappa_v-alpha_v*2|
	      E0.l += b_a;
	      E0.t -= b_a;
	    }
	    ext_param F(E.ctxt, new_tw, new_gamma_lambda, E0.tau, E0.l, E0.t);

	    flipped = not flipped; // January unsurprise for 3Cr
	    z_align(E0,F,flipped^not same_sign(E,E0));
	    // there was no 4th argument there since |E.t.dot(kappa)==0|
	    links.push_back(std::move(F)); // Cayley link
	  }

	} // end of 3ci and 3Cr cases
	else // twisted non-commutation: 3C+ or 3C-
	{
	  result = ascent ? three_complex_ascent : three_complex_descent;
	  links.push_back(complex_cross(ctxt,p,E));
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

bool ext_block::tune_signs(const blocks::common_block& block)
{
  repr::Ext_common_context ctxt
    (block.context().real_group(),delta,block.integral_subsystem());
  repr::Ext_rep_context param_ctxt(block.context(),delta);
  containers::sl_list<ext_param> links;
  for (BlockElt n=0; n<size(); ++n)
  { auto z=this->z(n);
    const ext_param E(param_ctxt,block.x(z),block.gamma_lambda(z));
    for (weyl::Generator s=0; s<rank(); ++s)
    { const ext_gen& p=orbit(s); links.clear(); // output arguments for |star|
      auto tp = star(ctxt,E,p,links);
      if (might_be_uncertain(descent_type(s,n)) and
	  data[s][n].links.first==UndefBlock) // then reset the uncertain type
      {
	data[s][n].type=tp; // replace type, leave possible links to |UndefBlock|
	continue; // so in this case there are no link signs to tune
      }
      else if (tp!=descent_type(s,n))
	return false; // something is wrong

      switch (tp)
      {
      case one_imaginary_pair_switched: case one_real_pair_switched:
      case one_real_nonparity: case one_imaginary_compact:
      case two_imaginary_single_double_switched:
      case two_real_single_double_switched:
      case two_real_nonparity: case two_imaginary_compact:
      case three_real_nonparity: case three_imaginary_compact:
	assert(links.empty());
	break;

      case one_complex_ascent: case one_complex_descent:
      case two_complex_ascent: case two_complex_descent:
      case three_complex_ascent: case three_complex_descent:
	{ assert(links.size()==1);
	  const ext_param q = *links.begin();
	  BlockElt m=cross(s,n); // cross neighbour as bare element of |*this|
	  if (m==UndefBlock)
	    break; // don't fall off the edge of a partial block
	  BlockElt cz = this->z(m); // corresponding element of (parent) |block|
	  ext_param F(param_ctxt,block.x(cz),block.gamma_lambda(cz)); // default
	  assert(same_standard_reps(q,F)); // must lie over same parameter
	  if (not same_sign(q,F))
	    flip_edge(s,n,m);
	}
	break;

      case one_imaginary_single: case one_real_single:
      case two_imaginary_single_single: case two_real_single_single:
	{ assert(links.size()==2);
	  const ext_param q0 = *links.begin();
	  const ext_param q1 = *std::next(links.begin());
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  if (m!=UndefBlock) // don't fall off the edge of a partial block
	  {
	    BlockElt Cz = this->z(m); // corresponding element of block
	    ext_param F(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	    assert(same_standard_reps(q0,F));
	    if (not same_sign(q0,F))
	      flip_edge(s,n,m);
	  }
	  if ((m=cross(s,n))!=UndefBlock) // cross link, don't fall off the edge
	  {
	    BlockElt cz = this->z(m);
	    ext_param Fc(param_ctxt,block.x(cz),block.gamma_lambda(cz));
	    assert(same_standard_reps(q1,Fc));
	    if (not same_sign(q1,Fc))
	      flip_edge(s,n,m);
	  }
	} break;
      case two_semi_imaginary: case two_semi_real:
      case three_semi_imaginary: case three_real_semi:
      case three_imaginary_semi: case three_semi_real:
	{ assert(links.size()==1);
	  const ext_param q = *links.begin();
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  if (m==UndefBlock)
	    break; // don't fall off the edge of a partial block
	  BlockElt Cz = this->z(m); // corresponding element of block
	  ext_param F(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  assert(same_standard_reps(q,F));
	  if (not same_sign(q,F))
	    flip_edge(s,n,m);
	}
	break;

      case one_imaginary_pair_fixed: case one_real_pair_fixed:
      case two_imaginary_double_double: case two_real_double_double:
	{ assert(links.size()==2);
	  const ext_param q0 = *links.begin();
	  const ext_param q1 = *std::next(links.begin());
	  BlockEltPair m=Cayleys(s,n);

	  if (m.first==UndefBlock)
	    break; // nothing to do if both are undefined

	  BlockElt Cz = this->z(m.first);
	  ext_param F0(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  bool straight = same_standard_reps(q0,F0);
	  const auto& node0 = straight ? q0 : q1;
	  assert(same_standard_reps(node0,F0));
	  if (not same_sign(node0,F0))
	    flip_edge(s,n,m.first);

	  if (m.second==UndefBlock)
	    break;

	  Cz = this->z(m.second);
	  ext_param F1(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  const auto& node1 = straight ? q1 : q0;
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node1,F1))
	    flip_edge(s,n,m.second);
	}
	break;

      case two_imaginary_single_double_fixed: case two_real_single_double_fixed:
	{ assert(links.size()==2);
	  const ext_param q0 = *links.begin();
	  const ext_param q1 = *std::next(links.begin());
	  BlockEltPair m=Cayleys(s,n);

	  if (m.first==UndefBlock)
	    break; // nothing to do if both are undefined

	  BlockElt Cz = this->z(m.first);
	  ext_param F0(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  bool straight=same_standard_reps(q0,F0);
	  const auto& node0 = straight ? q0 : q1;
	  assert(same_standard_reps(node0,F0));
	  if (not same_sign(node0,F0))
	    flip_edge(s,n,m.first);

	  if (m.second==UndefBlock)
	    break;

	  Cz= this->z(m.second);
	  ext_param F1(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  const auto& node1 = straight ? q1 : q0;
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node1,F1))
	    flip_edge(s,n,m.second);
	}
	break;
      } // |switch(tp)|
    } // |for(s)|
  } // |for(n)|
#ifndef NDEBUG // when debugging test braid relations for each extended block
  BitMap dummy(size());
  for (BlockElt x=0; x<size(); ++x)
    for (weyl::Generator s=0; s<rank(); ++s)
    {
      if (not check_quadratic(*this,s,x))
	throw std::runtime_error
	  ("Quadratic relation failure in extended block construction");
      for (weyl::Generator t=s+1; t<rank(); ++t)
	if (not check_braid(*this,s,t,x,dummy))
	  throw std::runtime_error
	    ("Braid relation failure in extended block construction");
    }
#endif
  return true; // report success if we get here
} // |tune_signs|


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


  // |ext_param| methods and functions

ext_param::ext_param
(const repr::Ext_rep_context& ec, const TwistedInvolution& tw,
   RatWeight gamma_lambda, Weight tau, Coweight l, Coweight t, bool flipped)
  : ctxt(ec), tw(tw)
  , l(std::move(l))
  , gamma_lambda(std::move(gamma_lambda))
  , tau(std::move(tau))
  , t(std::move(t))
  , flipped(flipped)
{
  validate(*this);
}


// contructor used for default extension once |x| and |gamma_lamba| are chosen
ext_param::ext_param
(const repr::Ext_rep_context& ec,
   KGBElt x, const RatWeight& gamma_lambda, bool flipped)
  : ctxt(ec)
  , tw(ec.real_group().kgb().involution(x)) // now computing |theta()| is valid
  , l(ell(ec.real_group().kgb(),x))
  , gamma_lambda(gamma_lambda)
  , tau(matreduc::find_solution
	(1-theta(), gamma_lambda.integer_diff<int>(delta()*gamma_lambda)))
  , t(matreduc::find_solution
	(theta().transposed()+1,(delta()-1).right_prod(l)))
  , flipped(flipped)
{
  validate(*this);
}

// build a default extended parameter for |sr| in the context |ec|
/*
  Importantly, this does not use |sr.gamma()| otherwise than for asserting its
  $\delta$-stability: the same default is used in |common_block| for an entire
  family of blocks, so dependence on |gamma| must be limited to dependence on
  its reduction modulo 1. Even though |gamma_lambda| is computed at the non
  $\delta$-fixed |srm.gamma_mod1()|, it is also (due to the way |mod_reduce|
  works) a proper value of |gamma_lambda| at |sr|, so |(1-delta)*gamma_lambda|,
  gives a valid value for the equation of which |tau| is a solution. However
  |gamma_lambda| may be a different representative than |rc.gamma_lambda(sr)|,
  so don't use that latter: it would give an undesired dependence on |gamma|.
*/
ext_param ext_param::default_extend
  (const repr::Ext_rep_context& ec, const repr::StandardRepr& sr)
{
  assert(((1-ec.delta())*sr.gamma().numerator()).isZero());

  auto srm =  repr::StandardReprMod::mod_reduce(ec,sr);
  // get default representative at |gamma%1|, normalised
  auto gamma_lambda=ec.gamma_lambda(srm);
  return ext_param(ec,sr.x(),gamma_lambda);
}

ext_param& ext_param::operator= (const ext_param& p)
{ assert(&ctxt==&p.ctxt); // assignment should remain in the same context
  tw=p.tw;
  l=p.l; gamma_lambda=p.gamma_lambda; tau=p.tau; t=p.t;
  flipped=p.flipped;
  return *this;
}
ext_param& ext_param::operator= (ext_param&& p)
{ assert(&ctxt==&p.ctxt); // assignment should remain in the same context
  tw=std::move(p.tw);
  l=std::move(p.l); gamma_lambda=std::move(p.gamma_lambda);
  tau=std::move(p.tau); t=std::move(p.t);
  flipped=p.flipped;
  return *this;
}


KGBElt ext_param::x() const
{ TitsElt a(ctxt.inner_class().titsGroup(),TorusPart(l),tw);
  return rc().kgb().lookup(a);
}

const WeightInvolution& ext_param::theta () const
  { return ctxt.inner_class().matrix(tw); }

repr::StandardRepr ext_param::restrict(const RatWeight& gamma) const
{
  const RatWeight gamma_rho = gamma-rho(rc().root_datum());
  const auto lambda_rho = gamma_rho.integer_diff<int>(gamma_lambda);
  return rc().sr_gamma(x(),lambda_rho,gamma);
} // |restrict|




/*
  This function serves to replace and circumvent |Rep_context::make_dominant|
  applied to a scaled parameter (as occurs in the ordinary deformation function
  by calling |finals_for| or |normalise|, both of which call |make_dominant|,
  after calling |scale|), where |make_dominant| maps any ordinary parameter to
  one with a dominant |gamma| component, and moreover descends through singular
  complex descents in the block to the lowest parameter equivalent to the
  initial parameter. The reason that this is necessary is that scaling only
  affects the |nu| component of the infinitesimal character, so it may make it
  traverse walls of Weyl chambers. Indeed the caller should make sure |sr|
  itself has dominant |gamma|, which moreover is assumed to be fixed by |delta|
  (if not, don't use this function).

  The difference with the functioning of |make_dominant| is that here we keep
  track of all extended parameter components inherited from |sr| (so before
  scaling its |nu| part by |factor|), transforming them from the default choices
  for |sr|, and at the end comparing the transformed values to the default
  choices at the final parameter reached, recording the sign in |flipped|.
 */
StandardRepr scaled_extended_dominant // result will have its |gamma()| dominant
(const Rep_context rc,
 const StandardRepr& sr, const WeightInvolution& delta,
 Rational factor, // |z.nu()| is scaled by |factor| first
 bool& flipped // records whether a net extended flip was computed
 )
{
  const RootDatum& rd=rc.root_datum();
  RealReductiveGroup& G=rc.real_group(); const KGB& kgb = rc.kgb();
  repr::Ext_rep_context ctxt(rc,delta);
  const ext_gens orbits = rootdata::fold_orbits(rd,delta);
  assert(is_dominant_ratweight(rd,sr.gamma())); // dominant
  assert(((delta-1)*sr.gamma().numerator()).isZero()); // $\delta$-fixed

  // first approximation to result is scaled input
  // importantly, $\lambda$ (or equivalently |lambda_rho|) is held fixed here
  auto scaled_sr = rc.sr(sr.x(),rc.lambda_rho(sr),sr.gamma()*factor);
  // it will be convenent to have a working (modifiable) copy of |gamma|
  RatWeight gamma = scaled_sr.gamma(); // a working copy
  KGBElt x = scaled_sr.x(); // another variable, for convenience

  ext_param E0 = ext_param::default_extend(ctxt,sr);

  E0.gamma_lambda += gamma-sr.gamma(); // shift |E0.gamma_lambda| by $\nu$ change

  int_Vector r_g_eval (rd.semisimpleRank()); // simple root evaluations at |-gr|
  { const RatCoweight& g_r=ctxt.g_rho_check();
    for (unsigned i=0; i<r_g_eval.size(); ++i)
      r_g_eval[i] = -g_r.dot(rd.simpleRoot(i));
  }

  { unsigned i; // index into |orbits|
    do // make |gamma_numer| dominant, uses only complex simple root reflections
      for (i=0; i<orbits.size(); ++i)
	if (kgb.status(x).isComplex(orbits[i].s0))
	{ const auto& s=orbits[i];
	  const auto& alpha_v = rd.simpleCoroot(s.s0);
	  if (alpha_v.dot(gamma.numerator())<0)
	  {
	    rd.act(s.w_kappa,gamma); // change infin.character representative
	    rd.act(s.w_kappa,E0.gamma_lambda);
	    rd.act(s.w_kappa,E0.tau);
	    rd.shifted_dual_act(E0.l,s.w_kappa,r_g_eval);
	    rd.dual_act(E0.t,s.w_kappa);
	    E0.flip(s.length()==2); // record flip for every 2C+/2C- done
	    x = kgb.cross(s.w_kappa,x);
	    break; // indicate we advanced; restart search for |s|
	  }
	} // |for(i)|, if |isComplex|
    while(i<orbits.size()); // continue until above |for| runs to completion
  } // end of transformation of extended parameter components

  // now ensure that |E| gets matching |gamma| and |theta| (for flipped test)
  ext_param E1(ctxt,kgb.involution(x),
	       E0.gamma_lambda,E0.tau,E0.l,E0.t,E0.flipped);

  { // descend through complex singular simple descents
    repr::Ext_common_context block_ctxt(G,delta, SubSystem::integral(rd,gamma));
    const auto& int_datum = block_ctxt.id();
    const ext_gens integral_orbits = rootdata::fold_orbits(int_datum,delta);
    const RankFlags singular_orbits = // flag singular among integral orbits
      reduce_to(integral_orbits,singular_generators(int_datum,gamma));
    /*
      |singular_orbits| are in fact orbits of simple roots, because we have
      ensured |gamma| is dominant, so the singular subsystem of |rd|
      is generated by simple (co)roots (evaluting to 0 on |ctxt.gamma()|)
    */
     // record the corresponding simple root indices in |rd|, in order
    containers::sl_list<unsigned> orbit_simple;
    for (auto it=singular_orbits.begin(); it(); ++it)
    { auto alpha=block_ctxt.subsys().parent_nr_simple(integral_orbits[*it].s0);
      orbit_simple.push_back(rd.simpleRootIndex(alpha));
    }

    while(true) // will |break| below if no singular complex descent exists
    {
      auto soit = singular_orbits.begin(); auto it=orbit_simple.begin();
      for (; not orbit_simple.at_end(it); ++it,++soit)
	if (kgb.isComplexDescent(*it,x))
	  break;

      assert((not soit())==orbit_simple.at_end(it));
      if (not soit()) // previous loop ran to completion
	break;

      // find orbit among |integral_orbits| corresponding to that ComplexDescent
      ext_gen p=integral_orbits[*soit];
      assert(block_ctxt.subsys().parent_nr_simple(p.s0)
	     ==rd.simpleRootNbr(*it)); // check that we located it

      containers::sl_list<ext_param> links;
      auto type = // compute neighbours in extended block
	star(block_ctxt,E1,p,links);
      assert(is_complex(type) or type==two_semi_real);
      E1 = *links.begin(); // replace |E| by descended parameter
      E1.flip(has_october_surprise(type)); // to undo extra flip |star|
      assert(x>E1.x()); // make sure we advance; we did simple complex descents
      x = E1.x(); // adapt |x| for complex descent test
    } // |while| a singular complex descent exists
  }

  // finally extract |StandardRepr| from |E|
  StandardRepr result = E1.restrict(gamma);

  // but the whole point of this function is to record the relative flip too!
  flipped = // compare |E1| to default
    not same_sign(E1,ext_param::default_extend(ctxt,result));
  return result;

} // |scaled_extended_dominant|

/*
  The following function determines whether an extended parameter has a descent
  for generator |kappa|, which is an orbit of singularly-simple roots, and since
  |gamma| is supposed dominant here these are actually simple roots. Rather than
  call |star| here to do the full analysis, we can do a simplified (there is no
  |to_simple_shift|) test of notably the parity condition in the real case.
 */
bool is_descent
  (const repr::Ext_common_context& ctxt, const ext_gen& kappa, const ext_param& E)
{ // easy solution would be to |return is_descent(star(E,kappa,dummy))|;
  const InvolutionTable& i_tab = E.rc().inner_class().involution_table();
  const InvolutionNbr theta = i_tab.nr(E.tw); // so use root action of |E.tw|
  const RootNbr n_alpha = ctxt.subsys().parent_nr_simple(kappa.s0);
  const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
  const RootDatum& rd = E.rc().root_datum();
  assert(rd.is_simple_root(n_alpha)); // as explained in the comment above

  // we don't need to inspect |kappa.type|, it does not affect descent status
  if (theta_alpha==n_alpha) // imaginary case, return whether compact
    return (E.ctxt.g_rho_check()-E.l).dot(rd.root(n_alpha)) %2!=0;
  if (theta_alpha==rd.rootMinus(n_alpha)) // real, return whether parity
    // whether $\<\alpha^\vee,\gamma-\lambda>$ even (since $\alpha$ will be
    // singular here, |gamma_lambda| gives the same result as |lambda| would)
    return E.gamma_lambda.dot(rd.coroot(n_alpha)) %2==0;
  else // complex
    return rd.is_negroot(theta_alpha);
} // |is_descent|

weyl::Generator first_descent_among
  (const repr::Ext_common_context& ctxt, RankFlags singular_orbits,
   const ext_gens& orbits, const ext_param& E)
{ for (auto it=singular_orbits.begin(); it(); ++it)
    if (is_descent(ctxt,orbits[*it],E))
      return *it;
  return orbits.size(); // no singular descents found
}

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
  assert(is_dominant_ratweight(rc.root_datum(),sr.gamma()));
  // we must assume |gamma| already dominant, DON'T call |make_dominant| here!

  repr::Ext_rep_context param_ctxt(rc,delta);
  repr::Ext_common_context ctxt
    (rc.real_group(),delta,SubSystem::integral(rc.root_datum(),sr.gamma()));

  const ext_gens orbits = rootdata::fold_orbits(ctxt.id(),delta);
  const RankFlags singular_orbits =
    reduce_to(orbits,singular_generators(ctxt.id(),sr.gamma()));

  containers::queue<ext_param> to_do { ext_param::default_extend(param_ctxt,sr) };
  containers::sl_list<std::pair<StandardRepr,bool> > result;

  do
  { const ext_param E= to_do.front();
    to_do.pop(); // we are done with |head|
    auto s = first_descent_among(ctxt,singular_orbits,orbits,E);
    if (s>=orbits.size()) // no singular descents, so append to result
      result.emplace_back
	(std::make_pair(E.restrict(sr.gamma()),not is_default(E)));
    else // |s| is a singular descent orbit
    { containers::sl_list<ext_param> links;
      auto type = star(ctxt,E,orbits[s],links);
      if (not is_like_compact(type)) // some descent, push to front of |to_do|
      { bool flip = has_october_surprise(type); // to undo extra flip |star|
	auto l_it=links.begin(); RootNbr witness;
	if (rc.is_nonzero(l_it->restrict(sr.gamma()),witness))
	{
	  l_it->flip(flip);
	  to_do.push(*l_it);
	}
	if (has_double_image(type)  // then maybe add second node after |head|
	    and rc.is_nonzero((++l_it)->restrict(sr.gamma()),witness))
	{ l_it->flip(flip);
	  to_do.push(*l_it);
	}
      }
    }
  }
  while(not to_do.empty());

  return result;
} // |extended_finalise|

template containers::simple_list<BlockElt> ext_block::condense
(matrix::Matrix<ext_kl::Pol>& M, RankFlags sing_orbs) const;

} // |namespace ext_block|

} // |namespace atlas|
