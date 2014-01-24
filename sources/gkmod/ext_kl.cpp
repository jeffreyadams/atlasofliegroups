/*
  This is ext_kl.cpp

  Copyright 2013, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_kl.h"

namespace atlas {
namespace ext_kl {

class PolEntry : public Pol
{
public:
  // constructors
  PolEntry() : Pol() {} // default constructor builds zero polynomial
  PolEntry(const Pol& p) : Pol(p) {} // lift polynomial to this class

  // members required for an Entry parameter to the HashTable template
  typedef std::vector<Pol> Pooltype;  // associated storage type
  size_t hashCode(size_t modulus) const; // hash function

  // compare polynomial with one from storage
  bool operator!=(Pooltype::const_reference e) const;
}; // |class KLPolEntry|

inline size_t PolEntry::hashCode(size_t modulus) const
{ const Pol& P=*this;
  if (P.isZero()) return 0;
  polynomials::Degree i=P.degree();
  size_t h=P[i]; // start with leading coefficient
  while (i-->0) h= ((h<<21)+(h<<13)+(h<<8)+(h<<5)+h+P[i]) & (modulus-1);
  return h;
}

bool PolEntry::operator!=(PolEntry::Pooltype::const_reference e) const
{
  if (degree()!=e.degree()) return true;
  if (isZero()) return false; // since degrees match
  for (polynomials::Degree i=0; i<=degree(); ++i)
    if ((*this)[i]!=e[i]) return true;
  return false; // no difference found
}


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
	good_asc.set(s); // good ascent: at most one upward neighbour
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
  unsigned inx=aux.x_index(x,y);
  return inx<col_y.size() ? col_y[inx] : inx==aux.self_index(y) ? 1 : 0;
}

int KL_table::mu(int i,BlockElt x, BlockElt y) const
{
  unsigned d=aux.block.l(y,x)-i;
  if (d%2!=0)
    return 0;
  d/=2;
  PolRef Pxy=P(x,y);
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
Pol KL_table::m(weyl::Generator s,BlockElt x, BlockElt y, bool go_up) const
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
	int sum=0;
	for (size_t i=interval.size(); i-->0; )
	{
	  BlockElt t=interval[i];
	  sum -= mu(1,x,t)*mu(1,t,y);
	}
	if (go_up and has_defect(type(s,y)))
	  sum -= mu(1,x,aux.block.Cayley(s,y));
	if (has_defect(type(s,x)))
	  sum += mu(1,aux.block.Cayley(s,x),y); // positive contribution
	result[0] += sum;
      }
      return result;
    }
  case ext_gen::three:
    {
      KLPol result(1+d,mu(2+d,x,y)); // $q\mu_{-2}(x,y)$ or $q^2\mu_{-3}(x,y)$
      if (d==0)
      {
	BlockEltList interval=mu1top(s,x,y);
	int sum=0;
	for (size_t i=interval.size(); i-->0; )
	{
	  BlockElt t=interval[i];
	  sum -= mu(1,x,t)*mu(1,t,y);
	}
	result[1] += sum;
	result[0] = result[1]; // change multiple of $q$ to multiple of $q+1$
      }
      else
      { // |d==1|
	result[0]=result[2]; // $\mu_{-3}(x,y)(q^2+1)$

	BlockEltList top=mu1top(s,x,y);
	int sum=0;
	for (size_t i=top.size(); i-->0; )
	{
	  BlockElt t=top[i];
	  sum -= mu(2,x,t)*mu(1,t,y);
	}

	BlockEltList bot=mu1bot(s,x,y);
	for (size_t i=bot.size(); i-->0; )
	{
	  BlockElt t=bot[i];
	  int m1xt=mu(1,x,t),acc = -mu(2,t,y);
	  for (size_t j=top.size(); j-->0 and top[j]>t; )
	  {
	    BlockElt u=top[j];
	    acc += mu(1,t,u)*mu(1,u,y); // positive contribution
	  }
	  sum += m1xt*acc;
	} // fot |t|
	if (go_up and has_defect(type(s,y)))
	  sum -= mu(1,x,aux.block.Cayley(s,y));
	if (has_defect(type(s,x)))
	  sum += mu(1,aux.block.Cayley(s,x),y); // positive contribution
	result[1] += sum;
      }
      return result;
    }
  } // |switch(kappa.type)|
  assert(false); return KLPol(); // dummy to keep compiler happy
} // |KL_table::m|


/* Use recursive formula (in degree-shifted Laurent polynomials in $r$)
  $m(x)\cong r^k p_{x,sy} + r^def(s,x)p_{s_x,sy}-\sum_{x<u<y}m(u)p_{x,u}$
  where congruence is modulo $r^{-1+def(s,y)}\Z[r^{-1}]$, and use symmetry of
  $m(x)$ to complete. There is a complication when $def(s,y)=1$, since a
  congruence modulo $\Z[r^{-1}]$ cannot be used to determine the coefficient
  of $r^0$ in $m(x)$, and instead one must use that the difference between the
  members of above congruence should be a multiple of $(r+r^{-1})$.
 */
Pol KL_table::set_m(weyl::Generator s,BlockElt x, BlockElt sy,
		    std::vector<Pol>& m) const
{
  const BlockElt y = // unique ascent by |s| of |sy|
    is_complex(type(s,sy)) ? aux.block.cross(s,sy) : aux.block.Cayley(s,sy);
  const int k = aux.block.orbit(s).length();
  // in the following |l(sy,x)| might be $-1$, which is OK
  const unsigned d=(k+aux.block.l(sy,x))%2; // degree parity for $m(x)$
  Pol Q= Pol(k,P(x,sy));

  // INCOMPLETE
  return Q;
}


bool KL_table::direct_recursion(BlockElt y,
				weyl::Generator& s,
				BlockElt& sy,
				std::vector<Pol>& out) const
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

  const int k = aux.block.orbit(s).length();
  const Pol qk_plus_1 = Pol(k,1)+Pol(1);
  const Pol qk_minus_1 = Pol(k,1)-Pol(1);
  static const Pol q1 = Pol(1,1)+Pol(1);
  const Pol qk_minus_q = Pol(k,1)-Pol(1,1);
  sy = is_complex(v) ? aux.block.cross(s,y) : aux.block.inverse_Cayley(s,y);

  out.resize(aux.length_floor(y)); // vector indexed by block elements |x|

  for (BlockElt x=aux.length_floor(y); x-->0; )
    if (aux.descent_set(x)[s]) // compute contributions at all descent elts
    {
      Pol& Q = out[x];
      switch(type(s,x))
      { default: assert(false); break; // list will only contain descent types
	  // complex
      case ext_block::one_complex_descent:
      case ext_block::two_complex_descent:
      case ext_block::three_complex_descent:
	// contribute $P_{sx,sy}+q^kP_{x,sy}$
	Q = P(aux.block.cross(s,x),sy) + Pol(k,P(x,sy));
	break;
	// imaginary compact, real switched
      case ext_block::one_imaginary_compact:
      case ext_block::two_imaginary_compact:
      case ext_block::three_imaginary_compact:
      case ext_block::one_real_pair_switched:
	// contribute $(q^k+1)P_{x,sy}$
	Q = qk_plus_1 * P(x,sy);
	break;
	// real type 1
      case ext_block::one_real_pair_fixed:
      case ext_block::two_real_double_double:
	{ // contribute $P_{x',sy}+P_{x'',sy}+(q^k-1)P_{x,sy}$
	  BlockEltPair sx = aux.block.inverse_Cayleys(s,x);
	  Q = P(sx.first,sy) + P(sx.second,sy) + qk_minus_1 * P(x,sy);
	}
	break;
	// real type 2
      case ext_block::one_real_single:
      case ext_block::two_real_single_single:
	// contribute $P_{x_s,sy}+q^kP_{x,sy}-P_{s*x,sy}$
	Q = P(aux.block.inverse_Cayley(s,x),sy) - P(aux.block.cross(s,x),sy)
	  + Pol(k,P(x,sy));
	break;
	// defect type descents: 2Cr, 3Cr, 3r
      case ext_block::two_semi_real:
      case ext_block::three_semi_real:
      case ext_block::three_real_semi:
	// contribute $(q+1)P_{x_s,sy}+(q^k-q)P_{x,sy}$
	Q = q1 * P(aux.block.inverse_Cayley(s,x),sy) + qk_minus_q * P(x,sy);
	break;
	// epsilon case: 2r21
      case ext_block::two_real_single_double:
	{ // contribute $P_{x',sy}\pm P_{x'',sy}+(q^2-1)P_{x,sy}$
	  BlockEltPair sx = aux.block.inverse_Cayleys(s,x);
	  Q = P(sx.first,sy)*aux.block.epsilon(s,x,0)
	    + P(sx.second,sy)*aux.block.epsilon(s,x,1)
	    + qk_minus_1 * P(x,sy); // since $k=2$ here
	}
      } // |switch(type(s,x))|
    } // |for(x)|
  return true;
}

void KL_table::fill_columns(BlockElt y)
{
  PolHash hash(storage_pool); // (re)construct hash table for the polynomials
  if (y==0 or y>aux.block.size())
    y=aux.block.size(); // fill whole block if no explicit stop was indicated
  column.reserve(y);
  kl::KLIndex zero=hash.match(Pol(0)); // ensure 0 and 1 hold these constants
  kl::KLIndex one=hash.match(Pol(1)); // as |KL_pol_index| and |P| assume this
  assert (zero==0 and one==1); ndebug_use(zero); ndebug_use(one);
  while (column.size()<y)
    fill_next_column(hash);
}

void KL_table::fill_next_column(PolHash& hash)
{
  const BlockElt y = column.size();
  column.push_back(kl::KLRow());
  if (aux.col_size(y)==0)
    return; // there is just the non-recorded $P(y,y)=1$
  column.back().resize(aux.col_size(y));
  weyl::Generator s;
  BlockElt sy;
  std::vector<Pol> extr_contribs,Ms;
  Ms.resize(y);
  if (direct_recursion(y,s,sy,extr_contribs))
  {
    const int defect = has_defect(type(s,y)) ? 1 : 0;
    const int k = aux.block.orbit(s).length();

    for (BlockElt u=aux.length_floor(y); u-->0; )
      if (aux.descent_set(u)[s])
      {
	unsigned d=aux.block.l(y,u)+defect; // doubled implicit degree shift $Q$
	Pol& Q = extr_contribs[u];
	if (Q.isZero())
	  continue;
	assert(2*Q.degree()<d+k);
	unsigned M_deg = 2*Q.degree()-d; // might be negative; if so, unused
	Pol& M_u = Ms[u];
	if (defect==0) // then just pick up high degree terms from $Q$
	{
	  if (2*Q.degree()<d)
	    continue; // no correction needed

	  // compute $m_s(u,sy)$, the correction coefficient for $c_u$
	  M_u = Pol(M_deg,Q[Q.degree()]);
	  if (M_u.degree()==2)
	    M_u[1]=Q[Q.degree()-1];
	  if (M_deg>0)
	    M_u[0] = M_u[M_deg]; // symmetrise if non-constant
	  Q -= Pol(Q.degree()-M_deg,M_u);
	}
	else // |defect==1|
	{ // we must ensure that $q+1$ divides $Q-q^{(d-M_deg)/2}M_u$
	  if (2*Q.degree()<=d) // only possible "constant" term to consider
	  {
	    M_u=Pol(Q.factor_by(1)); // but divide by $q+1$ always
	    assert(M_deg==0 or M_u.isZero()); // correct only constant term
	    if (M_u.isZero())
	      continue; // if no correction was necessary, move to next |u|
	  }
	  else
	  {
	    M_u = Pol(M_deg,Q[Q.degree()]);
	    M_u[0]=M_u[M_deg]; // symmetrise
	    Q -= Pol(Q.degree()-M_deg,M_u); // subtract contribution
	    int r = Q.factor_by(1);
	    if (M_deg==1)
	    {
	      assert(r==0); ndebug_use(r);
	    }
	    else if (r!=0)
	    {
	      assert(d%2==0 and Q.degree()==d/2);
	      M_u[1]=r;
	    }
	  }
	} // |if(defect)|
	for (BlockElt x=aux.length_floor(u); x-->0; )
	  if (aux.descent_set(x)[s])
	  { // subtract $q^{(d-M_deg)/2}M_u*P_{x,u}$ from contribution for $x$
	    extr_contribs[x] -= Pol((d-M_deg)/2,M_u*P(x,u));
	  }
      } // |for(u)|

    // Now |extr_contribs|

    kl::KLRow::reverse_iterator it = column.back().rbegin();
    for (BlockElt x=aux.length_floor(y); aux.prim_back_up(x,y); it++)
      if (aux.descent_set(x)[s]) // then we computed $P(x,y)$ above
        *it = hash.match(extr_contribs[x]);
      else // |x| might not be descent for |s| if primitive but not extremal
      { // find and use a double-valued ascent for |x| that is decsent for |y|
	assert(has_double_image(type(s,x))); // since |s| non-good ascent
	BlockEltPair sx = aux.block.Cayleys(s,x);
	Pol Q = P(sx.first,y) + P(sx.second,y); // computed earlier in this loop
        *it = hash.match(Q);
      }
    assert(it==column.back().rend()); // check that we've traversed the column
  }
  else // direct recursion was not possible
  {
    std::cerr << "No direct recursion for element " << aux.block.z(y)
	      << std::endl;
  }
} // |KL_table::fill_next_column|

} // |namespace kl|
} // |namespace atlas|
