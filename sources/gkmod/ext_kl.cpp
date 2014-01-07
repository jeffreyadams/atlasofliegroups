/*
  This is ext_kl.cpp

  Copyright 2013, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_kl.h"

namespace atlas {
namespace ext_kl {

descent_table::descent_table(const ext_block::extended_block& block)
  : descents(block.size()), good_ascents(block.size())
{
  for (BlockElt n = 0; n < block.size(); ++n)
  {
    RankFlags desc, good_asc;
    for (weyl::Generator s=0; s<block.rank(); ++s)
    {
      ext_block::DescValue v = block.descent_type(s,n);
      if (ext_block::is_descent(v))
	desc.set(s);
      else if (not ext_block::has_double_image(v))
	good_asc.set(s); // good ascent
    }
    descents[n]=desc;
    good_ascents[n]=good_asc;
  }

} // |descent_table|

unsigned int descent_table::col_size(BlockElt y) const
{
  return 0; // to be implemented
} // |descent_table::col_size|

bool descent_table::back_up(BlockElt& x, BlockElt y) const
{
  return false; // to be implemented
} // |descent_table::back_up|



kl::KLIndex KL_table::KL_pol_index(BlockElt x, BlockElt y) const
{
  const kl::KLRow& col_y = column[y];
  unsigned int xx=aux.x_index(x,y);
  return xx<col_y.size() ? col_y[xx] :
    x==aux.self_index(y) ? 1 : 0;
}

kl::KLCoeff KL_table::mu(int i,BlockElt x, BlockElt y) const
{
  unsigned d=block.length(y)-block.length(x)-i;
  if (d%2!=0)
    return 0;
  d/=2;
  kl::KLPolRef Pxy=P(x,y);
  return Pxy.degree()<d ? 0 : Pxy[d];
}

// Find descents |d| for |s| in interval $(x,y)$ with |mu(1,d,y)| nonzero
BlockEltList KL_table::mu1top(weyl::Generator s,BlockElt x, BlockElt y) const
{
  size_t ly=block.length(y);
  size_t l0=block.length(x)+1;
  assert(l0<=ly);
  BlockEltList result;
  result.reserve(block.length_first(ly)-block.length_first(l0));
  for (size_t l=l0+((ly+1-l0)&1); l<ly; l+=2)
    for (BlockElt d=block.length_first(l); d<block.length_first(l+1); ++d)
      if (is_descent(block.descent_type(s,d)))
	if (mu(1,d,y)!=0)
	  result.push_back(d);
  return result;
}

// Find descents |d| for |s| in interval $(x,y)$ with |mu(1,x,d)| nonzero
BlockEltList KL_table::mu1bot(weyl::Generator s,BlockElt x, BlockElt y) const
{
  size_t ly=block.length(y);
  size_t l0=block.length(x)+1;
  assert(l0<=ly);
  BlockEltList result;
  result.reserve(block.length_first(ly)-block.length_first(l0));
  for (size_t l=l0; l<ly; l+=2)
    for (BlockElt d=block.length_first(l); d<block.length_first(l+1); ++d)
      if (is_descent(block.descent_type(s,d)))
	if (mu(1,x,d)!=0)
	  result.push_back(d);
  return result;
}

// See theorem 9.3.10, whose case (1) $y\overset\kappa\to y$ does not apply
kl::KLPol KL_table::m(weyl::Generator s,BlockElt x, BlockElt y) const
{ assert(block.length(y)>block.length(x)); // should not be called otherwise
  assert(ext_block::is_descent(block.descent_type(s,x)));
  assert(not ext_block::is_descent(block.descent_type(s,y)));
  const unsigned d=block.l(y,x)%2; // 0 or 1
  switch(block.orbit(s).type)
  {
  case ext_gen::one: return KLPol(0,d==0 ? 0 : mu(1,x,y));
  case ext_gen::two:
    {
      KLPol result(d,mu(2-d,x,y));
      if (d!=0)
	result[0]=result[1]; // $\mu_1(x,y)*(q+1)$
      else
      { // |d==0|
	BlockEltList interval=mu1top(s,x,y);
	KLCoeff sum=0;
	for (size_t i=interval.size(); i-->0; )
	{
	  BlockElt t=interval[i];
	  sum+=mu(1,x,t)*mu(1,t,y);
	}
	if (ext_block::has_defect(block.descent_type(s,y)))
	  sum+=mu(1,x,block.Cayleys(s,y).first);
	if (ext_block::has_defect(block.descent_type(s,x)))
	  sum+=mu(1,block.Cayleys(s,x).first,y);
	result[0]-=sum;
      }
      return result;
    }
  case ext_gen::three:
    {
      KLPol result(1+d,mu(2+d,x,y));
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
	result[0]=result[1]; // multiply by $q+1$
      }
      else
      { // |d==1|
	result[0]=result[2]; // $\mu(3,x,y)(q^2+1)$
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
	  KLCoeff sum0=mu(2,t,y);
	  for (size_t j=top.size(); j-->0 and top[j]>t; )
	  {
	    BlockElt u=top[j];
	    sum0+=mu(1,t,u)*mu(1,u,y);
	  }
	  sum+=mu(1,x,t)*sum0;
	} // fot |t|
	if (ext_block::has_defect(block.descent_type(s,y)))
	  sum+=mu(1,x,block.Cayleys(s,y).first);
	if (ext_block::has_defect(block.descent_type(s,x)))
	  sum+=mu(1,block.Cayleys(s,x).first,y);
	result[1]-=sum;
      }
      return result;
    }
  } // |switch(kappa.type)|
  assert(false); return KLPol(); // dummy to keep compiler happy
}

void KL_table::fill_next_column()
{
  BlockElt y = column.size();
  column.push_back(kl::KLRow()); column.back().resize(aux.col_size(y));
  weyl::Generator s;
  for (s=0; s<block.rank(); ++s)
  { ext_block::DescValue v=block.descent_type(s,v);
    if (is_descent(v) and is_unique_image(v))
      break;
  }
  if (s<block.rank()) // now for all $x$ a recursion via $s$ is possible
    for (BlockElt x=block.length_first(block.length(y)); aux.back_up(x,y); )
    {
    }

  else // direct recursion was not possible
  {}
} // |KL_table::fill_next_column|

} // |namespace kl|
} // |namespace atlas|
