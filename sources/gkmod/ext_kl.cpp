/*
  This is ext_kl.cpp

  Copyright 2013, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_kl.h"

namespace atlas {
namespace ext_kl {

descent_table::descent_table(const ext_block::extended_block& eb)
  : descents(eb.size()), good_ascents(eb.size())
  , prim_index(1<<eb.rank(),std::vector<unsigned int>(eb.size(),0))
  , prim_count(1<<eb.rank(),0) // one for each descent set
  , block(eb)
{
  // following loop must decrease for primitivization calculation below
  for (BlockElt x = block.size(); x-->0; )
  {
    RankFlags desc, good_asc;
    for (weyl::Generator s=0; s<block.rank(); ++s)
    {
      ext_block::DescValue v = block.descent_type(s,x);
      if (is_descent(v))
	desc.set(s);
      else if (not has_double_image(v))
	good_asc.set(s); // good ascent
    }
    descents[x]=desc;
    good_ascents[x]=good_asc;

    // compute primitivizations of |x|, storing index among primitives for |D|
    // loop must be downwards so initially index is w.r.t. descending order
    for (unsigned long desc=prim_index.size(); desc-->0;)
    {
      RankFlags D(desc); // descent set for which to primitive
      D &= good_asc;
      if (D.none()) // then element |x| is primitive for the descent set
	prim_index[desc][x] = prim_count[desc]++; // self-ref; increment count
      else
      {
	weyl::Generator s = D.firstBit();
	ext_block::DescValue v = block.descent_type(s,x);
	BlockElt sx = is_complex(v) ? block.cross(s,x) : block.Cayley(s,x);
	assert(sx>x); // ascents go up in block; also |~0| exceeds everyobe
	prim_index[desc][x] = sx==UndefBlock ? ~0 : prim_index[desc][sx];
      }
    } // |for (desc)|
  } // |for(x)|

  // primitive lists will actually be stored increasing, so reverse indices
  for (unsigned long desc=prim_index.size(); desc-->0;)
  {
    BlockElt last=prim_count[desc]-1;
    for (BlockElt x = block.size(); x-->0; )
    {
      unsigned int& slot = prim_index[desc][x];
      if (slot != ~0u) // leave "cop out" indices as such
	slot = last-slot; // reverse all other indices
    }
  }

} // |descent_table|

// number of primimitive elements for descents(y) of length less than y
unsigned int descent_table::col_size(BlockElt y) const
{
  BlockElt x=length_floor(y);
  if (prim_back_up(x,y)) // find last primitive $x$ of length less than $y$
    return x_index(x,y)+1; //
  return 0; // no primitives below length of |y| at all
} // |descent_table::col_size|


bool descent_table::prim_back_up(BlockElt& x, BlockElt y) const
{
  RankFlags desc=descents[y];
  while (x-->0)
    if ((good_ascents[x]&desc).none())
      return true;
  return false;
} // |descent_table::prim_back_up|

bool descent_table::extr_back_up(BlockElt& x, BlockElt y) const
{
  RankFlags desc=descents[y];
  while (x-->0)
    if (descents[x].contains(desc))
      return true; // stop when no descents of |y| are (any) ascents of |x|
  return false;
} // |descent_table::prim_back_up|

kl::KLIndex KL_table::KL_pol_index(BlockElt x, BlockElt y) const
{
  const kl::KLRow& col_y = column[y];
  unsigned int xx=aux.x_index(x,y);
  return xx<col_y.size() ? col_y[xx] :
    x==aux.self_index(y) ? 1 : 0;
}

kl::KLCoeff KL_table::mu(int i,BlockElt x, BlockElt y) const
{
  unsigned d=aux.block.l(y,x)-i;
  if (d%2!=0)
    return 0;
  d/=2;
  kl::KLPolRef Pxy=P(x,y);
  return Pxy.degree()<d ? 0 : Pxy[d];
}

// Find descents |d| for |s| in interval $(x,y)$ with |mu(1,d,y)| nonzero
BlockEltList KL_table::mu1top(weyl::Generator s,BlockElt x, BlockElt y) const
{
  size_t ly=aux.block.length(y);
  size_t l0=aux.block.length(x)+1;
  assert(l0<=ly);
  BlockEltList result;
  result.reserve(aux.block.length_first(ly)-aux.block.length_first(l0));
  for (size_t l=l0+((ly+1-l0)&1); l<ly; l+=2)
    for (BlockElt d=aux.block.length_first(l);
	 d<aux.block.length_first(l+1); ++d)
      if (is_descent(type(s,d)))
	if (mu(1,d,y)!=0)
	  result.push_back(d);
  return result;
}

// Find descents |d| for |s| in interval $(x,y)$ with |mu(1,x,d)| nonzero
BlockEltList KL_table::mu1bot(weyl::Generator s,BlockElt x, BlockElt y) const
{
  size_t ly=aux.block.length(y);
  size_t l0=aux.block.length(x)+1;
  assert(l0<=ly);
  BlockEltList result;
  result.reserve(aux.block.length_first(ly)-aux.block.length_first(l0));
  for (size_t l=l0; l<ly; l+=2)
    for (BlockElt d=aux.block.length_first(l);
	 d<aux.block.length_first(l+1); ++d)
      if (is_descent(type(s,d)))
	if (mu(1,x,d)!=0)
	  result.push_back(d);
  return result;
}

// See theorem 9.3.10, whose case (1) $y\overset\kappa\to y$ does not apply
kl::KLPol KL_table::m(weyl::Generator s,BlockElt x, BlockElt y) const
{ // check that we are not being called in unexpected conditions
  assert(aux.block.length(y)>aux.block.length(x));
  assert(is_descent(type(s,x)));
  assert(not is_descent(type(s,y)));
  const unsigned d=aux.block.l(y,x)%2; // 0 or 1
  switch(aux.block.orbit(s).type)
  {
  case ext_gen::one: return KLPol(0,d==0 ? 0 : mu(1,x,y));
  case ext_gen::two:
    {
      KLPol result(d,mu(2-d,x,y)); // $\mu_{-2}(x,y)$ or $q\mu_{-1}(x,y)$
      if (d!=0)
	result[0]=result[1]; // change to $(q+1)\mu_{-1}(x,y)$
      else
      { // |d==0|
	BlockEltList interval=mu1top(s,x,y);
	KLCoeff sum=0;
	for (size_t i=interval.size(); i-->0; )
	{
	  BlockElt t=interval[i];
	  sum+=mu(1,x,t)*mu(1,t,y);
	}
	if (has_defect(type(s,y)))
	  sum+=mu(1,x,aux.block.Cayley(s,y));
	if (has_defect(type(s,x)))
	  result[0]+=mu(1,aux.block.Cayley(s,x),y); // positive contribution
	result[0]-=sum; // subtract negative contributions last
      }
      return result;
    }
  case ext_gen::three:
    {
      KLPol result(1+d,mu(2+d,x,y)); // $q\mu_{-2}(x,y)$ or $q^2\mu_{-3}(x,y)$
      if (d==0)
      {
	BlockEltList interval=mu1top(s,x,y);
	KLCoeff sum=0;
	for (size_t i=interval.size(); i-->0; )
	{
	  BlockElt t=interval[i];
	  sum+=mu(1,x,t)*mu(1,t,y);
	}
	result[1]-=sum;
	result[0]=result[1]; // change multiple of $q$ to multiple of $q+1$
      }
      else
      { // |d==1|
	result[0]=result[2]; // $\mu_{-3}(x,y)(q^2+1)$
	BlockEltList top=mu1top(s,x,y);
	KLCoeff sum=0;
	for (size_t i=top.size(); i-->0; )
	{
	  BlockElt t=top[i];
	  sum+=mu(2,x,t)*mu(1,t,y);
	}
	BlockEltList bot=mu1bot(s,x,y);
	for (size_t i=bot.size(); i-->0; )
	{
	  BlockElt t=bot[i];
	  KLCoeff m1xt=mu(1,x,t),acc=0;
	  sum+=mu(1,x,t)*mu(2,t,y);
	  for (size_t j=top.size(); j-->0 and top[j]>t; )
	  {
	    BlockElt u=top[j];
	    acc += mu(1,t,u)*mu(1,u,y);
	  }
	  result[1]=m1xt*acc; // start with positive contribution only
	} // fot |t|
	if (has_defect(type(s,y)))
	  sum+=mu(1,x,aux.block.Cayley(s,y));
	if (has_defect(type(s,x)))
	  result[1] += mu(1,aux.block.Cayley(s,x),y); // positive contribution
	result[1]-=sum; // subtract negative contributions last
      }
      return result;
    }
  } // |switch(kappa.type)|
  assert(false); return KLPol(); // dummy to keep compiler happy
} // |KL_table::m|

bool KL_table::direct_recursion(BlockElt y,
				weyl::Generator& s,
				std::vector<KLPol>& out) const
{
  ext_block::DescValue v; // make value survice loop
  for (s=0; s<rank(); ++s)
  {
    v=type(s,y);
    if (is_descent(v) and is_unique_image(v))
      break;
  }
  if (s==rank())
    return false; // none of the generators gives a direct recursion

  int k = aux.block.orbit(s).length();
  BlockElt sy = is_complex(v) ? aux.block.cross(s,y) : aux.block.Cayley(s,y);

  out.reserve(aux.col_size(y)); // but we'll do only extremal elements

  for (BlockElt x=aux.length_floor(y); aux.extr_back_up(x,y);)
  {
    assert(is_descent(type(s,x))); // or else |x| was certainly not extremal
    out.push_back(KLPol());
    KLPol& Q = out.back();
    switch(type(s,x))
    { default: assert(false); break; // list will only contain descent types
	// complex
    case ext_block::one_complex_descent:
    case ext_block::two_complex_descent:
    case ext_block::three_complex_descent:
      // contribute $P_{sx,sy}+q^kP_{x,sy}$
      Q = P(aux.block.cross(s,x),sy);
      Q.safeAdd(P(x,sy),k);
      break;
      // imaginary compact, real switched
    case ext_block::one_imaginary_compact:
    case ext_block::two_imaginary_compact:
    case ext_block::three_imaginary_compact:
    case ext_block::one_real_pair_switched:
      // contribute $(q^k+1)P_{x,sy}$
      Q = P(x,sy);
      Q.safeAdd(Q,k);
      break;
      // real type 1
    case ext_block::one_real_pair_fixed:
    case ext_block::two_real_double_double:
      { // contribute $P_{x',sy}+P_{x'',sy}+(q^k-1)P_{x,sy}$
	BlockEltPair sx = aux.block.inverse_Cayleys(s,x);
	Q = P(sx.first,sy);
	Q.safeAdd(P(sx.second,sy));
	kl::KLPolRef Pxsy = P(x,sy);
	Q.safeAdd(Pxsy,k);
	Q.safeSubtract(Pxsy);
      }
      break;
      // real type 2
    case ext_block::one_real_single:
    case ext_block::two_real_single_single:
      { // contribute $P_{x_s,sy}+q^kP_{x,sy}-P_{s*x,sy}$
	Q = P(aux.block.inverse_Cayley(s,x),sy);
	Q.safeAdd(P(x,sy),k);
	Q.safeSubtract(P(aux.block.cross(s,x),sy));
      }
      break;
      // defect type descents: 2Cr, 3Cr, 3r
      case ext_block::two_semi_real:
      case ext_block::three_semi_real:
      case ext_block::three_real_semi:
      { // contribute $(q+1)P_{x_s,sy}+(q^k-q)P_{x,sy}$
	Q = P(aux.block.inverse_Cayley(s,x),sy);
	Q.safeAdd(Q,1);
	kl::KLPolRef Pxsy = P(x,sy);
	Q.safeAdd(Pxsy,k);
	Q.safeSubtract(Pxsy,1);
      }
      break;
      // epsilon case: 2r21
    case ext_block::two_real_single_double:
      { // contribute $P_{x',sy}\pm P_{x'',sy}+(q^2-1)P_{x,sy}$
	BlockEltPair sx = aux.block.inverse_Cayleys(s,x);
	kl::KLPolRef Pxsy = P(x,sy);
	assert(aux.block.epsilon(s,x,0)>0); // for now this is assumed
	Q = P(sx.first,sy);
	Q.safeAdd(Pxsy,2); // since $k=2$ here
	Q.safeSubtract(Pxsy);
	if (aux.block.epsilon(s,x,1)>0)
	  Q.safeAdd(P(sx.second,sy));
	else
	  Q.safeSubtract(P(sx.second,sy));
      }
    } // |switch(type(s,x))|
  } // |for(x)|
  return true;
}

void KL_table::fill_next_column()
{
  BlockElt y = column.size();
  column.push_back(kl::KLRow()); column.back().resize(aux.col_size(y));
  weyl::Generator s;
  for (s=0; s<rank(); ++s)
  { ext_block::DescValue v=type(s,y);
    if (is_descent(v) and is_unique_image(v))
      break;
  }
  if (s<rank()) // now for all $x$ a recursion via $s$ is possible
    for (BlockElt x=aux.length_floor(y); aux.prim_back_up(x,y); )
    {
    }

  else // direct recursion was not possible
  {}
} // |KL_table::fill_next_column|

} // |namespace kl|
} // |namespace atlas|
