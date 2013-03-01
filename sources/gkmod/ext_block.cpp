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

std::vector<ext_gen> twist_orbits (const TwistedWeylGroup& W)
{
  const WeylGroup& Wg = W.weylGroup();
  unsigned int size=0;
  for (weyl::Generator s=0; s<W.rank(); ++s)
    if (W.twisted(s)>=s)
      ++size;

  std::vector<ext_gen> result; result.reserve(size);

  for (weyl::Generator s=0; s<W.rank(); ++s)
    if (W.twisted(s)==s)
      result.push_back(ext_gen(s));
    else if (W.twisted(s)>s)
      result.push_back(ext_gen(Wg.commutes(s,W.twisted(s)),s,W.twisted(s)));

  return result;
} // |twist_orbits|

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

bool has_double_image(DescValue v)
{
  switch (v)
  {
  case one_real_pair_fixed:
  case one_imaginary_pair_fixed:
  case two_real_double_double: case two_imaginary_double_double:
  case two_imaginary_single_double: case two_real_single_double:
    return true;
  case one_complex_ascent: case one_complex_descent:
  case two_complex_ascent: case two_complex_descent:
  case three_complex_ascent: case three_complex_descent:
  case one_real_nonparity: case one_imaginary_compact:
  case two_real_nonparity: case two_imaginary_compact:
  case three_real_nonparity: case three_imaginary_compact:
  case one_imaginary_single: case one_real_single:
  case one_real_pair_switched: case one_imaginary_pair_switched:
  case two_semi_imaginary: case two_semi_real:
  case two_imaginary_single_single: case two_real_single_single:
  case three_semi_imaginary: case three_real_semi:
  case three_imaginary_semi: case three_semi_real:
    return false;
  }
  assert(false); return false; // keep compiler happy
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
      if (block.Hermitian_dual(link)==link)
	return one_imaginary_pair_fixed;
      link=UndefBlock; return one_imaginary_pair_switched;
    case DescentStatus::RealTypeI:
      link=block.inverseCayley(p.s0,z).first;
      if (block.Hermitian_dual(link)==link)
	return one_imaginary_pair_fixed;
      link=UndefBlock; return one_imaginary_pair_switched;
    }
  case ext_gen::two:
    switch (block.descentValue(p.s0,z))
    {
       case DescentStatus::ComplexAscent:
	 link=block.cross(p.s0,z);
	 if (link==block.cross(p.s1,z))
	   return two_semi_imaginary;
	 link=block.cross(p.s1,link); return two_complex_ascent;
       case DescentStatus::ComplexDescent:
	 link=block.cross(p.s0,z);
	 if (link==block.cross(p.s1,z))
	   return two_semi_real;
	 link=block.cross(p.s1,link); return two_complex_descent;
       case DescentStatus::RealNonparity:
	 link=UndefBlock; return two_real_nonparity;
       case DescentStatus::ImaginaryCompact:
	 link=UndefBlock; return two_imaginary_compact;
       case DescentStatus::ImaginaryTypeI:
	 link=block.cayley(p.s1,block.cayley(p.s0,z).first).first;
	 assert(block.Hermitian_dual(link)==link);
	 return block.descentValue(p.s0,link)==DescentStatus::RealTypeI
	   ? two_imaginary_single_single :  two_imaginary_single_double;
       case DescentStatus::RealTypeII:
	 link=block.inverseCayley(p.s1,block.inverseCayley(p.s0,z).first).first;
	 assert(block.Hermitian_dual(link)==link);
	 return block.descentValue(p.s0,link)==DescentStatus::ImaginaryTypeII
	   ? two_real_single_single : two_real_single_double;
       case DescentStatus::ImaginaryTypeII:
	 link=block.cayley(p.s1,block.cayley(p.s0,z).first).first;
	 if (block.Hermitian_dual(link)!=link)
	 {
	   link=block.cross(p.s0,link);
	   assert(block.Hermitian_dual(link)==link);
	 }
	 return two_imaginary_double_double;
       case DescentStatus::RealTypeI:
	 link=block.inverseCayley(p.s1,block.inverseCayley(p.s0,z).first).first;
	 if (block.Hermitian_dual(link)!=link)
	 {
	   link=block.cross(p.s0,link);
	   assert(block.Hermitian_dual(link)==link);
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
	 if (link==block.cross(p.s1,link))
	 {
	   assert(block.descentValue(p.s1,link)==
		  DescentStatus::ImaginaryTypeII);
	   link=block.cayley(p.s1,link).first;
	   if (block.Hermitian_dual(link)!=link)
	   {
	     link=block.cross(p.s1,link);
	     assert(block.Hermitian_dual(link)==link);
	   }
	   return three_semi_imaginary;
	 }
	 link=block.cross(p.s0,block.cross(p.s1,link));
	 assert(block.Hermitian_dual(link)==link);
	 return three_complex_ascent;
       case DescentStatus::ComplexDescent:
	 link=block.cross(p.s0,z);
	 if (link==block.cross(p.s1,link))
	 {
	   assert(block.descentValue(p.s1,link)==DescentStatus::RealTypeI);
	   link=block.inverseCayley(p.s1,link).first;
	   if (block.Hermitian_dual(link)!=link)
	   {
	     link=block.cross(p.s1,link);
	     assert(block.Hermitian_dual(link)==link);
	   }
	   return three_semi_real;
	 }
	 link=block.cross(p.s0,block.cross(p.s1,link));
	 assert(block.Hermitian_dual(link)==link);
	 return three_complex_descent;
       case DescentStatus::ImaginaryTypeI:
	 link=block.cross(p.s1,block.cayley(p.s0,z).first);
	 assert(link==block.cross(p.s0,block.cayley(p.s1,z).first));
	 assert(block.Hermitian_dual(link)==link);
	 return three_imaginary_semi;
       case DescentStatus::RealTypeII:
	 link=block.cross(p.s1,block.inverseCayley(p.s0,z).first);
	 assert(link==block.cross(p.s0,block.inverseCayley(p.s1,z).first));
	 assert(block.Hermitian_dual(link)==link);
	 return three_real_semi;
       case DescentStatus::ImaginaryTypeII: case DescentStatus::RealTypeI:
	 assert(false);
    }
  } // |switch (p.type)|
  assert(false); return one_complex_ascent; // keep compiler happy
}

extended_block::extended_block
  (const Block_base& block,const TwistedWeylGroup& W)
  : parent(block)
  , tW(W)
  , orbit(twist_orbits(tW))
  , info()
  , data(orbit.size())
{
  unsigned int folded_rank = orbit.size();

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
      DescValue type = extended_type(block,z,orbit[s],link);
      data[s].push_back(block_fields(type));
      if (link==UndefBlock)
	continue; // |s| done for imaginary compact and real nonparity cases
      if (is_descent(type))
	info.back().length = length(child_nr[link])+1;
      data[s].back().links.first=child_nr[link];
      switch (type)
      {
      default: break;

      case one_imaginary_single:
      case one_real_single: // in these cases find "cross neighbour" for |s|
	data[s].back().links.second = child_nr[block.cross(orbit[s].s0,z)];
	break;

      case one_real_pair_fixed:
      case one_imaginary_pair_fixed: // in these cases find second Cayley image
	data[s].back().links.second =
	  child_nr[block.cross(orbit[s].s0,link)];
	break;

      case two_imaginary_single_single:
      case two_real_single_single: // now find "double cross neighbour" for |s|
	assert(block.cross(orbit[s].s0,block.cross(orbit[s].s1,z))==
	       block.cross(orbit[s].s1,block.cross(orbit[s].s0,z)));
	data[s].back().links.second =
	  child_nr[block.cross(orbit[s].s1,block.cross(orbit[s].s0,z))];
	break;

      case two_imaginary_single_double:
      case two_real_single_double: // find second Cayley image
	assert(block.cross(orbit[s].s0,link)==block.cross(orbit[s].s1,link));
	data[s].back().links.second =
	  child_nr[block.cross(orbit[s].s0,link)];
	break;

      case two_imaginary_double_double:
      case two_real_double_double: // find second Cayley image
	assert(block.cross(orbit[s].s0,block.cross(orbit[s].s1,link))==
	       block.cross(orbit[s].s1,block.cross(orbit[s].s0,link)));
	data[s].back().links.second =
	  child_nr[block.cross(orbit[s].s1,block.cross(orbit[s].s0,link))];
	break;
      }
    }
  } // |for(n)|
}

} // namespace ext_block

} // namespace atlas
