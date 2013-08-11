/*
  This is ext_block.cpp

  Copyright (C) 2013 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_block.h"

#include <cassert>
#include <vector>

#include "blocks.h"
#include "weyl.h"

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

bool is_proper_ascent(DescValue v)
{
  return not(is_descent(v) or is_like_nonparity(v));
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

    // single valued extended Cayleys use second link for cross action
  case one_imaginary_single:
  case one_real_single:
  case two_imaginary_single_single:
  case two_real_single_single:
    return data[s][n].links.second;

  }
  assert(false); return UndefBlock; // keep compiler happy
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
	? two_imaginary_single_single :  two_imaginary_single_double;
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

extended_block::extended_block
  (const Block_base& block,const TwistedWeylGroup& W)
  : parent(block)
  , tW(W)
  , info()
  , data(parent.folded_rank())
{
  unsigned int folded_rank = data.size();
  if (folded_rank==0 or block.Hermitian_dual(0)==UndefBlock)
    return; // block not stable under twist, so leave extended block empty

  std::vector<BlockElt> child_nr(block.size(),UndefBlock);
  std::vector<BlockElt> parent_nr;

  for (BlockElt z=0; z<block.size(); ++z)
    if (block.Hermitian_dual(z)==z)
    {
      child_nr[z]=parent_nr.size();
      parent_nr.push_back(z);
    }

  info.reserve(parent_nr.size());
  for (weyl::Generator s=0; s<folded_rank; ++s)
    data[s].reserve(parent_nr.size());

  for (BlockElt n=0; n<parent_nr.size(); ++n)
  {
    BlockElt z=parent_nr[n];
    info.push_back(elt_info(z));
    info.back().length = 0; // default value if no descents are found
    for (weyl::Generator s=0; s<folded_rank; ++s)
    {
      BlockElt link;
      DescValue type = extended_type(block,z,parent.orbit(s),link);
      data[s].push_back(block_fields(type));
      if (link==UndefBlock)
	continue; // |s| done for imaginary compact and real nonparity cases
      if (is_descent(type))
	info.back().length = length(child_nr[link])+1;
      data[s].back().links.first=child_nr[link];
      BlockElt second = UndefBlock;
      switch (type)
      {
      default: continue;

      case one_imaginary_single:
      case one_real_single: // in these cases find "cross neighbour" for |s|
	second = block.cross(parent.orbit(s).s0,z);
	break;

      case one_real_pair_fixed:
      case one_imaginary_pair_fixed: // in these cases find second Cayley image
	second = block.cross(parent.orbit(s).s0,link);
	break;

      case two_imaginary_single_single:
      case two_real_single_single: // now find "double cross neighbour" for |s|
	if (block.cross(parent.orbit(s).s0,z)==UndefBlock)
	  if (block.cross(parent.orbit(s).s1,z)==UndefBlock)
	    continue; // can't get there, let's hope it doesn't exist
	  else second = block.cross(parent.orbit(s).s0,
				    block.cross(parent.orbit(s).s1,z));
	else
	{
	  second = block.cross(parent.orbit(s).s1,
			       block.cross(parent.orbit(s).s0,z));
	  if (block.cross(parent.orbit(s).s1,z)!=UndefBlock)
	    assert(block.cross(parent.orbit(s).s0,
			       block.cross(parent.orbit(s).s1,z))==second);
	}
	break;

      case two_imaginary_single_double:
      case two_real_single_double: // find second Cayley image
	second = block.cross(parent.orbit(s).s0,link);
	assert(block.cross(parent.orbit(s).s1,link)==second);
	break;

      case two_imaginary_double_double:
      case two_real_double_double: // find second Cayley image
	if (block.cross(parent.orbit(s).s0,link)==UndefBlock)
	  if (block.cross(parent.orbit(s).s1,link)==UndefBlock)
	    continue; // can't get there, let's hope it doesn't exist
	  else second = block.cross(parent.orbit(s).s0,
				    block.cross(parent.orbit(s).s1,link));
	else
	{
	  second = block.cross(parent.orbit(s).s1,
			       block.cross(parent.orbit(s).s0,link));
	  if (block.cross(parent.orbit(s).s1,link)!=UndefBlock)
	    assert(block.cross(parent.orbit(s).s0,
			       block.cross(parent.orbit(s).s1,link))
		   ==second);
	}
	break;
      }
      if (second!=UndefBlock)
	data[s].back().links.second = child_nr[second];
    }
  } // |for(n)|
}

} // namespace ext_block

} // namespace atlas
