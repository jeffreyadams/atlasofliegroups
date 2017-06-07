/*
  This is ext_kl.cpp

  Copyright 2013-2016, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_kl.h"
#include "basic_io.h"
#include "wgraph.h"

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
  while (i-->0)
    h= (h<<21)+(h<<13)+(h<<8)+(h<<5)+h+P[i];
  return h & (modulus-1);
}

bool PolEntry::operator!=(PolEntry::Pooltype::const_reference e) const
{
  if (degree()!=e.degree()) return true;
  if (isZero()) return false; // since degrees match
  for (polynomials::Degree i=0; i<=degree(); ++i)
    if ((*this)[i]!=e[i]) return true;
  return false; // no difference found
}


descent_table::descent_table(const ext_block::ext_block& eb)
  : descents(eb.size()), good_ascents(eb.size())
  , prim_index(1<<eb.rank(),std::vector<unsigned int>(eb.size(),0))
  , prim_flip(eb.size(),BitMap(prim_index.size()))
  , block(eb)
{
  // counts of primitive block elements, one for each descent set
  std::vector<BlockElt> prim_count(1<<eb.rank(),0);

  // following loop must decrease for primitivisation calculation below
  for (BlockElt x = block.size(); x-->0; )
  {
    RankFlags desc, good_asc;
    for (weyl::Generator s=0; s<block.rank(); ++s)
    {
      ext_block::DescValue v = block.descent_type(s,x);
      if (ext_block::is_descent(v))
	desc.set(s);
      else if (not has_double_image(v))
	good_asc.set(s); // good ascent: at most one upward neighbour
    }
    descents[x]=desc;
    good_ascents[x]=good_asc;

    BitMap& flip = prim_flip[x]; // place to record primitivisation flips $x$

    // compute primitivisations of |x|, storing index among primitives for |D|
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
	if (is_like_nonparity(block.descent_type(s,x)))
	  prim_index[desc][x] = ~0; // stop primitivisation with zero result
	else
	{
	  BlockElt sx = block.some_scent(s,x);
	  assert(sx>x); // ascents go up in block
	  prim_index[desc][x] = prim_index[desc][sx];
	  if ((block.epsilon(s,x,sx)<0)!=prim_flip[sx].isMember(desc))
	    flip.insert(desc);
	}
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

} // |descent_table| constructor

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
} // |descent_table::extr_back_up|

KL_table::KL_table(const ext_block::ext_block& b, std::vector<Pol>& pool)
  : aux(b), storage_pool(pool), column()
{ // ensure first two pool entries are constant polynomials $0$, and $1$
  if (pool.empty())
    pool.push_back(Pol(0));
  else
    assert(pool[0]==Pol(0));
  if (pool.size()==1)
    pool.push_back(Pol(1));
  else
    assert(pool[1]==Pol(1));
}

std::pair<kl::KLIndex,bool>
KL_table::KL_pol_index(BlockElt x, BlockElt y) const
{ const kl::KLRow& col_y = column[y];
  unsigned inx=aux.x_index(x,y);
  if (inx<col_y.size())
    return std::make_pair(col_y[inx],aux.flips(x,y));
  else if (inx==aux.self_index(y)) // diagonal entries are unrecorded
    return std::make_pair(kl::KLIndex(1),aux.flips(x,y));
  else
    return std::make_pair(kl::KLIndex(0),false); // out of bounds implies zero
}

Pol KL_table::P(BlockElt x, BlockElt y) const
{
  auto index = KL_pol_index(x,y);
  return index.second ? -storage_pool[index.first] : storage_pool[index.first];
}

// coefficient of P_{x,y} of $q^{(l(y/x)-i)/2}$ (used with i=1,2,3 only)
int KL_table::mu(short unsigned int i,BlockElt x, BlockElt y) const
{
  unsigned int d=l(y,x);
  if (d<i or (d-=i)%2!=0) // or |l(y,x)<i|
    return 0; // coefficient would be at negative or non-integral degree
  PolRef Pxy=P(x,y);
  return Pxy.degree_less_than(d/=2) ? 0 : Pxy[d];
}

Pol qk_plus_1(int k)
{
  assert(k>=1 and k<=3);
  Pol res = Pol(k,1);
  res[0] = 1;
  return res;
}

Pol qk_minus_1(int k)
{
  assert(k>=1 and k<=3);
  Pol res = Pol(k,1);
  res[0] = -1;
  return res;
}

Pol qk_minus_q(int k)
{
  assert(k>1 and k<=3);
  Pol res = Pol(k,1);
  res[1] = -1;
  return res;
}


// component of element $a_x$ in product $(T_s+1)C_{sy}$
Pol KL_table::product_comp (BlockElt x, weyl::Generator s, BlockElt sy) const
{
  assert(is_descent(type(s,x)));  // otherwise don't call this function
  BlockEltList neighbours; neighbours.reserve(s);
  aux.block.add_neighbours(neighbours,s,x);
  Pol result=aux.block.T_coef(s,x,x)*P(x,sy); // start with term from diagonal
  for (auto it=neighbours.begin(); it!=neighbours.end(); ++it)
    result += aux.block.T_coef(s,x,*it)*P(*it,sy);
  return result;
} // |KL_table::product_comp|


/*
  In the descriptions below $r=\sqrt q$ and $p_{x,y}=r^{-l(y/x)}P_{x,y}$,
  which is a polynomial in $\Z[r^{-1}]$, without constant term unless $x=y$.

  our analogue of $\mu$ in ordinary KL computations takes the form of a
  symmetric Laurent polynomial $m$ in $r$. Since we have no data structure
  for Laurent polynomials, we shift exponents to get into $\Z[q]$
*/

// auxiliary to form $aq^{-1}+b+aq$, shifted to an ordinary polynomial
inline Pol m(int a,int b) { return a==0 ? Pol(b) : qk_plus_1(2)*a + Pol(1,b); }

/*
  Find $m_s(x,y)$ (symmetric Laurent polynomials in $r$) by recursive formula
  $m(x)\cong r^k p_{x,y} + def(s,x) r p_{s_x,y}-\sum_{x<u<y}p_{x,u}m(u)$ where
  $m(x)$ is $m_s(x,y)$, and congruence is modulo $r^{-1+def(s,y)}\Z[r^{-1}]$;
  then use symmetry of $m(x)$ to complete. The actual result returned is
  shifted minimally to an ordinary polynomial in $q=r^2$ (function |m| above).

  There is a complication when $def(s,y)=1$, since a congruence modulo
  $r^0\Z[r^{-1}]$ cannot be used to determine the coefficient of $r^0$ in
  $m(x)$. Instead one must use the fact that the difference between members of
  above congruence should be a multiple of $(r+r^{-1})$. To that end we use
  for our computations, instead of the coefficient of $r^0$ of polynomials,
  the appropriate |up_remainder(1,d)| values. As |up_remainder| is, with
  appropriate shifts, a ring morphism, it can be applied to separate factors.

  This function will be called when the values of $m(u)$ have already been
  computed for $u$ of length greater than $x$, and these values are passed in
  the final argument |M|, so as to avoid inefficient recursive calls.
*/
Pol KL_table::get_M(weyl::Generator s, BlockElt x, BlockElt y,
		    const std::vector<Pol>& M) const
{
  const ext_block::ext_block& bl=aux.block;
  const BlockElt z =  bl.some_scent(s,y); // unique ascent by |s| of |y|

  const unsigned defect = has_defect(type(s,z)) ? 1 : 0;
  const unsigned k = bl.orbit(s).length();

  if (k==1)
    return Pol(mu(1,x,y)); // will be zero if $l(y,x)$ is even

  if (k==2)
  {
    if (l(y,x)%2!=0)
      return q_plus_1() * Pol(mu(1,x,y));
    if (defect==0)
    {
      int acc = mu(2,x,y);
      if (has_defect(type(s,x)))
      {
	BlockElt sx=bl.Cayley(s,x);
	acc += mu(1,sx,y)*bl.epsilon(s,sx,x); // sign is |T_coef(s,x,sx)[1]|
      }
      for (unsigned l=bl.length(x)+1; l<bl.length(z); l+=2)
	for (BlockElt u=bl.length_first(l); u<bl.length_first(l+1); ++u)
	  if (aux.is_descent(s,u) and not M[u].isZero())
	    acc -= mu(1,x,u)*M[u][0];
      return Pol(acc);
    }
    // for |k==2| there remains the defect (for the |y| to |z| link) case
    int acc= product_comp(x,s,y).up_remainder(1,(l(z,x)+1)/2);
    for (unsigned lu=bl.length(x)+2; lu<bl.length(z); lu+=2)
      for (BlockElt u=bl.length_first(lu); u<bl.length_first(lu+1); ++u)
	if (aux.is_descent(s,u) and not M[u].isZero())
	  acc -= P(x,u).up_remainder(1,l(u,x)/2)*M[u][0];
    return Pol(acc);
  }
  if (k==3)
  {
    if (l(y,x)%2==0) // now we need a multiple of $1+q$
    {
      int acc = mu(2,x,y); // degree $1$ coefficient of |product_comp(x,s,y)|
      for (unsigned l=bl.length(x)+1; l<bl.length(z); l+=2)
	for (BlockElt u=bl.length_first(l); u<bl.length_first(l+1); ++u)
	  if (aux.is_descent(s,u) and M[u].degree()==2)
	    acc -= mu(1,x,u)*M[u][2];
      return q_plus_1() * acc;
    }

    // now we need a polynomial of the form $a+bq+aq^2$ for some $a,b$
    int a = mu(1,x,y); // coefficient $r^2$ in centered |product_comp(x,s,y)|
    if (defect==0)
    {
      int b = mu(3,x,y); // coefficient $r^0$ in same
      if (has_defect(type(s,x)))
      {
	BlockElt sx=bl.Cayley(s,x);
	b += mu(1,sx,y)*bl.epsilon(s,sx,x); // sign is |T_coef(s,x,sx)[1]|
      }
      for (BlockElt u=bl.length_first(bl.length(x)+1);
	   u<aux.length_floor(z); u++ )
	if (aux.is_descent(s,u) and M[u].degree()==2-l(u,x)%2)
	  b -= mu(M[u].degree(),x,u)*M[u][M[u].degree()];
      return m(a,b);
    }
    // there remains the |k==3|, even degree $m_s(x,y)$, defect $y$ case
    Pol Q = product_comp(x,s,y);
    if (a!=0)
      Q -= Pol((l(z,x)-1)/2,qk_plus_1(2)*a); // shaves top term
    assert(Q.degree_less_than((l(z,x)+3)/2)); // $\deg(Q)\leq(l(z,x)+1)/2$
    int b= Q.up_remainder(1,(l(z,x)+1)/2); // remainder by $q+1$

    /* since odd degree $m_s(u,y)$ do not contribute to $q+1$ remainder
       we restrict to those $u$ of the same length parity as $x$ */
    for (unsigned lu=bl.length(x)+2; lu<bl.length(z); lu+=2)
      for (BlockElt u=bl.length_first(lu); u<bl.length_first(lu+1); ++u)
	if (aux.is_descent(s,u) and not M[u].isZero())
	{ // extract $b-2a$ from |M[u]| representing $ar^{-2}+br^0+ar^2$
	  assert(M[u].degree()%2==0);
	  int mu_rem = M[u].degree()==0 ? M[u][0] : M[u][1]-2*M[u][0];
	  b -= P(x,u).up_remainder(1,l(u,x)/2)*mu_rem;
	}
    return m(a,b);
  }

  // this was the general case, now unused; it always works but not optimally
  Pol Q= product_comp(x,s,y);

  for (BlockElt u=bl.length_first(bl.length(x)+1); u<aux.length_floor(z); u++)
    if (aux.is_descent(s,u) and not M[u].isZero())
    { // subtract $q^{(d-deg(M))/2}M_u*P_{x,u}$ from contribution for $x$
      unsigned d=l(z,u)+defect; // doubled implicit degree shift
      assert(M[u].degree()<=d);
      Q -= Pol((d-M[u].degree())/2,P(x,u)*M[u]);
    }

  return extract_M(Q,l(z,x)+defect,defect);
} // |KL_table::get_M|

/*
  This is largely the same formula, but used in different context, which
  obliges to possibly leave out some term. Here one knows that $y$ is real
  nonparity for $s$, so in particular has no defect and there is no element
  $z$; also $x$ is known to be a descent for $s$ (unlike in the code above).
  On the other hand one is still busy computing the Hecke element $C_y$. It is
  a precondition that its coefficient $P_{x,y}$ has already been determined,
  and stored, but $P_{x',y}$ for $x'<x$ need not be; we must therefore refrain
  from (implicit) references to such polynomials. Again the vector $M$ can be
  used to safely access the (complete) values $M_s(u,y)$ for all $u>x$.

  Comparing with the formulas above, the terms to skip are those involving
  |Cayley(s,x)| (since |s| is a descent for |x|, in fact a downward Cayley).
 */
Pol KL_table::get_Mp(weyl::Generator s, BlockElt x, BlockElt y,
		     const std::vector<Pol>& M) const
{
  const ext_block::ext_block& bl=aux.block;
  const unsigned k = bl.orbit(s).length();
  if (k==1) // nothing changed for this case
    return  Pol(l(y,x)%2==0 ? 0 : mu(1,x,y));
  if (k==2)
  {
    if (l(y,x)%2!=0)
      return q_plus_1() * Pol(mu(1,x,y));
    int acc = mu(2,x,y);
    for (unsigned l=bl.length(x)+1; l<bl.length(y); l+=2)
      for (BlockElt u=bl.length_first(l); u<bl.length_first(l+1); ++u)
	if (aux.is_descent(s,u) and not M[u].isZero())
	{
	  assert(M[u].degree()==1 and M[u][0]==M[u][1]);
	  acc -= mu(1,x,u)*M[u][1];
	}
      return Pol(acc);
  }

  assert(k==3); // this case remains
  if (l(y,x)%2==0) // now we need a multiple of $1+q$
  {
    int acc = mu(2,x,y); // degree $1$ coefficient of |product_comp(x,s,y)|
    for (unsigned l=bl.length(x)+1; l<bl.length(y); l+=2)
      for (BlockElt u=bl.length_first(l); u<bl.length_first(l+1); ++u)
	if (aux.is_descent(s,u) and M[u].degree()==2)
	  acc -= mu(1,x,u)*M[u][2];
    return q_plus_1() * acc;
  }

  // now we need a polynomial of the form $a+bq+aq^2$ for some $a,b$
  int a = mu(1,x,y); // degree $2$ coefficient of |product_comp(x,s,y)|
  int b = mu(3,x,y); // degree $0$ coefficient of |product_comp(x,s,y)|
  for (BlockElt u=bl.length_first(bl.length(x)+1); u<aux.length_floor(y); u++)
    if (aux.is_descent(s,u) and M[u].degree()==2-l(u,x)%2)
      b -= mu(M[u].degree(),x,u)*M[u][M[u].degree()];
  return m(a,b);
} // |KL_table::get_Mp|



bool KL_table::has_direct_recursion(BlockElt y,
				    weyl::Generator& s, BlockElt& sy) const
{
  for (s=0; s<rank(); ++s)
  {
    const ext_block::DescValue v=type(s,y);
    if (is_descent(v) and is_unique_image(v))
    {
      sy = aux.block.some_scent(s,y); // some descent by $s$ of $y$
      return true;
    }
  }
  return false; // none of the generators gives a direct recursion
}


void KL_table::fill_columns(BlockElt y)
{
  PolHash hash(storage_pool); // (re)construct hash table for the polynomials
  if (y==0 or y>aux.block.size())
    y=aux.block.size(); // fill whole block if no explicit stop was indicated
  column.reserve(y);
  while (column.size()<y)
    fill_next_column(hash);
}

/*
  Clear terms of degree $\geq d/2$ in $Q$ by subtracting $r^d*m$ where $m$
  is a symmetric Laurent polynomial in $r=\sqrt q$, and if $defect>0$ dividing
  what remains by $q+1$, which division must be exact. Return $r^{deg(m)}m$.
  Implemented only under the hypothesis that $\deg(Q)<(d+3)/2$ initially,
  and $\deg(Q)\leq d$ (so if $d=0$ then $Q$ must be constant). $Q$ can be zero.
*/
Pol KL_table::extract_M(Pol& Q,unsigned d,unsigned defect) const
{
  assert(Q.degree_less_than(d/2+2) and Q.degree_less_than(d+1));
  unsigned M_deg = 2*Q.degree()-d; // might be negative; if so, unused
  Pol M(0); // result

  if (defect==0) // easy case; just pick up too high degree terms from $Q$
  {
    if (Q.degree_less_than((d+1)/2)) // that is, $deg(Q)<d/2$ mathematically
      return M; // no correction needed
    // now cases where |M_deg| was "negative" are gone, so we can safely use it

    // compute $m_s(u,sy)$, the correction coefficient for $c_u$
    assert(M_deg<3);
    M=Pol(M_deg,Q[Q.degree()]); // top term of |Q|, shifted to top term of |M|
    assert(M.degree()==M_deg); // in particular |M| is nonzero
    if (M_deg>0)
    {
      M[0] = M[M_deg]; // symmetrise if non-constant
      if (M_deg==2)
	M[1]=Q[Q.degree()-1]; // set sub-dominant coefficient here
    }

    assert(Q.degree()>=M_deg);
    // the need to have this "explains" the precondition $\deg(Q)\leq d$

    Q -= Pol(Q.degree()-M_deg,M); // subtract monomial multiple of |M|
    return M;
  } // |if(defect==0)|

  // now $defect=1$; we must ensure that $q+1$ divides $Q-q^{(d-M_deg)/2}M$
  if (not Q.degree_less_than(d/2+1))// that is, $deg(Q)>d/2$ mathematically
  {
    assert(M_deg!=0 and M_deg<3); // now |0<M_deg<3|
    M=Pol(M_deg,Q[Q.degree()]); // top term of |Q|, shifted to top term of |M|
    M[0]=M[M_deg]; // symmetrise (might leave middle or 3 terms as zero)
    assert(Q.degree()>=M_deg); // |Q| should have sufficient degree for:
    Q -= Pol(Q.degree()-M_deg,M); // subtract contribution of |M| from |Q|
    assert(Q.degree_less_than(d/2+1)); // terms conceptually of degree>0 are gone
  }
  // now divide by $1+q$, allowing remainder (degree $d/2$) from middle term |M|
  int c = Q.factor_by_1_plus_q(d/2);
  assert (c==0 or d%2==0); // if $d$ odd, there should be no such remainder
  if (c==0) // and in any case, if there was no remainder, leave |M| as is
    return M;

  // otherwise add constant $c$ to $m$, as $cr^d$ had to be subtracted from $Q$
  if (M.isZero())
    M=Pol(c); // if there were no terms, create one of degree $0$
  else
  {
    assert(M_deg==2); // now the constant term of |m| has index 1 in |M|
    M[1]=c;
  }
  return M;
} // |KL_table::extract_M|

void KL_table::fill_next_column(PolHash& hash)
{
  const BlockElt y = column.size();
  column.push_back(kl::KLRow());
  if (aux.col_size(y)==0)
    return; // there is just the non-recorded $P(y,y)=1$
  column.back().resize(aux.col_size(y));

  weyl::Generator s;
  BlockElt sy;
  if (has_direct_recursion(y,s,sy))
  {
    const unsigned defect = has_defect(type(s,y)) ? 1 : 0;
    const int sign = aux.block.epsilon(s,sy,y);
    const BlockElt floor_y =aux.length_floor(y);

    std::vector<PolEntry> cy(floor_y,(PolEntry()));
    std::vector<Pol>Ms(floor_y,(Pol()));

    // fill array |cy| with initial contributions from $(T_s+1)*c_{sy}$
    for (BlockElt x=floor_y; x-->0; )
      if (aux.is_descent(s,x)) // compute contributions at all descent elts
	cy[x] = product_comp(x,s,sy);

    // next loop downwards, refining coefficient polys to meet degree bound
    for (BlockElt u=floor_y; u-->0; )
      if (aux.is_descent(s,u))
      {
	unsigned d=aux.block.l(y,u)+defect; // doubled degree shift of |cy[u]|
	assert(u<cy.size());
	if (cy[u].isZero())
	  continue;

	assert(u<Ms.size());
	Pol gM = get_M(s,u,sy,Ms);
	Ms[u]=extract_M(cy[u],d,defect);
	assert(Ms[u]==gM);

	if (Ms[u].isZero())
	  continue;

	d -= Ms[u].degree();
	assert (d%2==0);
	d/=2;
	// now update contributions to all lower descents |x|
	for (BlockElt x=aux.length_floor(u); x-->0; )
	  if (aux.is_descent(s,x))
	  { // subtract $q^{(d-M_deg)/2}M_u*P_{x,u}$ from contribution for $x$
	    assert(x<cy.size());
	    cy[x] -= Pol(d,P(x,u)*Ms[u]);
	  }
      } // |for(u)|

    // finally copy relevant coefficients from |cy| array to |column[y]|
    kl::KLRow::reverse_iterator it = column.back().rbegin();
    for (BlockElt x=floor_y; aux.prim_back_up(x,y); it++)
      if (aux.is_descent(s,x)) // then we computed $P(x,y)$ above
        *it = hash.match(cy[x]*sign);
      else // |x| might not be descent for |s| if primitive but not extremal
      { // find and use a double-valued ascent for |x| that is decsent for |y|
	assert(has_double_image(type(s,x))); // since |s| non-good ascent
	BlockEltPair sx = aux.block.Cayleys(s,x);
	PolEntry Q = P(sx.first,y);  // computed earlier in this loop
	if (aux.block.epsilon(s,x,sx.first)<0)
	  Q *= -1;
	if (aux.block.epsilon(s,x,sx.second)>0)
	  Q += P(sx.second,y);
	else
	  Q -= P(sx.second,y);
        *it = hash.match(Q);
      }
    assert(it==column.back().rend()); // check that we've traversed the column
  }
  else // direct recursion was not possible
    do_new_recursion(y,hash);

  assert(check_polys(y));
 } // |KL_table::fill_next_column|

/*
  Basic idea for new recursion: if some $s$ is real nonparity for $y$ and a
  proper ascent for $x$ (with some restriction in case of type 1) then one has

  $$
  0 = [T_x](T_s+1).C_y - [T_x]\sum_u [s\in\tau(u)]r^{l(y/u)+k} P_{x,u}m_s(u,y)
  $$

  where in the first term occurs $P_{x,y}$ (which we are after) with a nonzero
  coefficient, and $P_{x^s,y}$ (known by descending induction on $x$), while
  the terms in the summation (where $u=x$ does not contribute since
  $s\notin\tau(x)$) are also known by descending induction on $x$
 */
void KL_table::do_new_recursion(BlockElt y,PolHash& hash)
{
  const BlockElt floor_y =aux.length_floor(y);
  std::vector<PolEntry> cy(floor_y,(PolEntry()));
  kl::KLRow::iterator out_it = column.back().end();
  std::vector<weyl::Generator> rn_s; rn_s.reserve(rank());
  std::vector<std::vector<Pol> > M_s; M_s.reserve(rank());
  for (weyl::Generator s=0; s<rank(); ++s)
    if (is_like_nonparity(type(s,y)))
    {
      rn_s.push_back(s);
      M_s.push_back(std::vector<Pol>(floor_y,Pol()));
    }

  { // check the absence of elements for which recursion would not work
    BlockEltList downs = aux.block.down_set(y);
    for (BlockEltList::const_iterator it=downs.begin(); it!=downs.end(); ++it)
    {
      BlockElt u = *it; weyl::Generator s;
      for (unsigned i=0; i<rn_s.size(); ++i)
      {
	const ext_block::DescValue tsu=type(s=rn_s[i],u);
	if (not is_descent(tsu) and has_defect(tsu)) // defect ascent: not-good
	  std::cerr << "Bad element " << aux.block.z(u)
		    << "in down-set for " << aux.block.z(y)
		    << "; would need M_" << s+1 << '['
		    << aux.block.z(aux.block.Cayley(s,u)) << ']' << std::endl;
      } // |for(i)|, loop over |s|
    } // |for(it)|, check downset
  }

  for (BlockElt x=floor_y; x-->0; )
  {
    if (aux.easy_set(x,y).any())
    {
      weyl::Generator s=aux.very_easy_set(x,y).firstBit();
      if (s<rank()) // non primitive case; equate to a previous polynomial
	cy[x] = P(x,y); // primitivisation returns copy of such a polynomial
      else // do primitive but not extremal case
      { // must compute by hand here, since |P(x,y)| would be yet undefined
	s = aux.easy_set(x,y).firstBit();
	assert(has_double_image(type(s,x))); // since |s| non-good ascent
	BlockEltPair sx = aux.block.Cayleys(s,x);
	cy[x] = P(sx.first,y);  // computed earlier this loop
	if (aux.block.epsilon(s,x,sx.first)<0)
	  cy[x] *= -1;
	if (aux.block.epsilon(s,x,sx.second)>0)
	  cy[x] += P(sx.second,y);
	else
	  cy[x] -= P(sx.second,y);
      }
    }
    else // |x| is extremal for |y|, so we must do real computation
    { // first seek proper |s|
      unsigned i; ext_block::DescValue tsx; weyl::Generator s;
      for (i=0; i<rn_s.size(); ++i)
      { tsx=type(s=rn_s[i],x); // save values for in case we |break| below
	if (is_like_type_1(tsx) // take cases like i1 only if (endgame case):
	   ? aux.easy_set(aux.block.cross(s,x),y).any() // another |t| helps;
	   : is_proper_ascent(tsx) // all other ascents except rn for x are OK
	   )
	  break;
	if (is_like_compact(tsx))
	  break; // having types (ic,rn) for $x,y$ allows new recursion too
      }

      if (i<rn_s.size()) // that is, we did |break| above
      {	// so we still have |s==rn_s[i]| and |tsx==type(s,x)|
	std::vector<Pol>& M = M_s[i];
	const unsigned k = aux.block.orbit(s).length();
	PolEntry& Q=cy[x];
	const BlockElt last_u=aux.block.length_first(aux.block.length(x)+1);

	// initialise $Q=\sum_{x<u<y}[s\in\tau(u)]r^{l(y/u)+k}P_{x,u}m_s(u,y)$
	for (BlockElt u=floor_y; u-->last_u; )
	  if (is_descent(type(s,u)) and not M[u].isZero())
	    Q += Pol((aux.block.l(y,u)+k-M[u].degree())/2, P(x,u)*M[u]);

	// subtract terms for ascent(s) of |x|; divide by its own coefficient
	switch(tsx)
	{
 	case ext_block::one_complex_ascent:
	case ext_block::two_complex_ascent:
	case ext_block::three_complex_ascent:
	  { // |(is_complex(tsx))|
	    BlockElt sx=aux.block.cross(s,x);
	    if (sx<floor_y) // subtract contrib. from $[T_x](T_s+1).T_{sx}=q^k$
	      Q -= aux.block.T_coef(s,x,sx)*cy[sx]; // coef is $\pm q^k$
	  } // implicit division of $Q$ here is by |T_coef(s,x,x)==1|
	  break;
	case ext_block::two_semi_imaginary:
	case ext_block::three_semi_imaginary:
	case ext_block::three_imaginary_semi:
	  { // |has_defect(tsx)|
	    BlockElt sx=aux.block.Cayley(s,x);
	    if (sx<floor_y)
	      Q -= aux.block.T_coef(s,x,sx)*cy[sx]; // coef is $\pm(q^k-q)$

	    /* divide by |T_coef(s,x,x)==1+q|, knowing that $Q$ may be missing
	       a term in effective degree $r^1$, from |mu(1,x,y)| that is not
	       included in $M[sx]$ when used to fill |Q| in loop above */
	    int c = // remainder in upward division by $[T_x](T_s+1).T_x=1+q$
	      // the remainder being taken in degree $\ceil l(y/x)/2$
	      Q.factor_by_1_plus_q((aux.block.l(y,x)+1)/2);
	    if (aux.block.l(y,x)%2!=0)
	      assert(-c==Q.coef(aux.block.l(y,x)/2));
	    else
	      assert(c==0); // division should be exact, leaving no $r^0$ term
 	    ndebug_use(c);
	  }
	  break;
	case ext_block::one_imaginary_pair_fixed:
	case ext_block::two_imaginary_double_double:
	  { // |is_like_type_2(tsx)|
	    assert(has_double_image(tsx)); // since it is a type 2 ascent
	    BlockEltPair sx=aux.block.Cayleys(s,x);
	    if (sx.first<floor_y)
	      Q -= aux.block.T_coef(s,x,sx.first)*cy[sx.first]; // $\pm(q^k-1)$
	    if (sx.second<floor_y)
	      Q -= aux.block.T_coef(s,x,sx.second)*cy[sx.second]; // idem
	    Q/=2;                     // divide by |T_coef(s,x,x)==2|
	  }
	  break;
	case ext_block::one_imaginary_single:
	case ext_block::two_imaginary_single_single:
	  { // |is_like_type_1(tsx)|, this used to be called the endgame case
	    BlockElt x_prime=aux.block.Cayley(s,x);
	    if (x_prime<floor_y)
	      Q -= aux.block.T_coef(s,x,x_prime)*cy[x_prime]; // $\pm(q^k-1)$
	    // implict division of $Q$ here is by |T_coef(s,x,x)==1|

	    // now subtract off $P_{s\times x,y}$, easily computed on the fly
	    BlockElt s_cross_x = aux.block.cross(s,x);
	    assert(aux.easy_set(s_cross_x,y).any()); // this was tested above
	    const weyl::Generator t=aux.easy_set(s_cross_x,y).firstBit();
	    auto ttscx=type(t,s_cross_x);
	    if (has_double_image(ttscx))
	    {
	      BlockEltPair sx_up_t = aux.block.Cayleys(t,s_cross_x);
	      if (sx_up_t.first<floor_y)
		Q -= cy[sx_up_t.first];
	      if (sx_up_t.second<floor_y)
		Q -= cy[sx_up_t.second];
	    }
	    else
	    {
	      BlockElt sx_up_t=aux.block.Cayley(t,s_cross_x);
	      if (sx_up_t<floor_y)
		Q -= cy[sx_up_t];
	    }
	  }
	  break;
	case ext_block::two_imaginary_single_double_fixed:
	  {
	    BlockEltPair sx=aux.block.Cayleys(s,x);
	    if (sx.first<floor_y)
	      Q -= aux.block.T_coef(s,x,sx.first)*cy[sx.first]; // $\pm(q^2-1)$
	    if (sx.second<floor_y)
	      Q -= aux.block.T_coef(s,x,sx.second)*cy[sx.second]; // idem
 	    Q/=2; // divide by |T_coef(s,x,x)==2|
	  }
	  break;
	case ext_block::one_imaginary_compact:
	case ext_block::one_real_pair_switched:
	case ext_block::two_imaginary_compact:
	case ext_block::two_real_single_double_switched:
	case ext_block::three_imaginary_compact:
	  // these cases require no additional terms to be substracted
	  Q.factor_by_1_plus_q_to_the(k,(aux.block.l(y,x)-1)/2+k); // degree
	  assert(Q.degree_less_than((aux.block.l(y,x)+1)/2));
	    // that was the condition $\deg(Q) \leq l(y,x)-1/2|, computed safely
	  break;
	default: assert(false); // other cases should not have selected |s|
	} // |switch(tsx)|
	// now |Q| is stored in |cy[x]|
      } // end of |if(i<rn_s.size())|; no |else|, just leave |cy[x]==Pol(0)|
    } // end of easy/hard condition

    if (aux.very_easy_set(x,y).none())
      *--out_it = hash.match(cy[x]); // store result whenever primitive

    // now if there is a defect ascent from x, update |M_s| for |mu(1,x,y)|
    if (aux.block.l(y,x)%2!=0 and cy[x].degree()==aux.block.l(y,x)/2)
      for (unsigned j=0; j<rn_s.size(); ++j)
      {
	const weyl::Generator s=rn_s[j];
	const ext_block::DescValue tsx=type(s,x);
	if (not is_descent(tsx) and has_defect(tsx))
	{
	  const BlockElt sx = aux.block.Cayley(s,x);
	  assert(sx<floor_y); // could only fail if |x| in downset tested above
	  Pol& dst = M_s[j][sx];
	  int mu = cy[x].coef(aux.block.l(y,x)/2) * aux.block.epsilon(s,x,sx);
	  dst += Pol (dst.degree()==2 ? 1 : 0, mu);
	}
      }
    // and update the entries |M_s[j][x]|
    for (unsigned j=0; j<rn_s.size(); ++j)
      if (is_descent(type(rn_s[j],x)))
	M_s[j][x] = get_Mp(rn_s[j],x,y,M_s[j]);
  } // |for(x)|

  assert(out_it==column[y].begin()); // check that we've traversed the column
} // |KL_table::do_new_recursion|


bool check(const Pol& P_sigma, const KLPol& P)
{
  if (P_sigma.isZero())
  { if (P.isZero())
      return true;
    for (unsigned i=0; i<=P.degree(); ++i)
      if (P[i]%2!=0)
	return false;
    return true;
  }
  if (P.degree_less_than(P_sigma.degree()))
    return false; // there are terms of |P_sigma| outside the degree bound
  for (polynomials::Degree i=0; i<=P.degree(); ++i)
  {
    KLCoeff d = P[i]+KLCoeff(P_sigma.coef(i)); // unsigned addition
    if (d%2!=0 or d>2*P[i]) // unequal parity, or |abs(P_sigma[i])>abs(P[i])|
      return false;
  }
  return true;
}

bool KL_table::check_polys(BlockElt y) const
{
  bool result = true;
#ifndef NDEBUG
  kl::KLContext untwisted(aux.block.untwisted());
  untwisted.fill();
  for (BlockElt x=y; x-->0; )
    if (not check(P(x,y),untwisted.klPol(aux.block.z(x),aux.block.z(y))))
    {
      std::cerr << "Mismatch at (" << aux.block.z(x) << ',' << aux.block.z(y)
		<< "): ";
      std::cerr << P(x,y) << " and "
		<< untwisted.klPol(aux.block.z(x),aux.block.z(y)) << std::endl;
      result=false;
    }
#endif
  return result;
}

void ext_KL_matrix (const StandardRepr p, const int_Matrix& delta,
		    const Rep_context& rc, // the rest is output
		    std::vector<StandardRepr>& block_list,
		    int_Matrix& P_mat,
		    int_Vector& lengths)
{ BlockElt entry_element;
  if (not ((delta-1)*p.gamma().numerator()).isZero())
  {
    std::cout << "Delta does not fix gamma=" << p.gamma() << "." << std::endl;
    throw std::runtime_error("No valid extended bock");
  }
  param_block B(rc,p,entry_element);
  ext_block::ext_block eblock (rc.innerClass(), B, delta);

  BlockElt size= // size of extended block we shall use; before compression
    eblock.element(entry_element+1);

  std::vector<ext_kl::Pol> pool;
  KL_table twisted_KLV(eblock,pool);
  twisted_KLV.fill_columns(size); // fill up to and including |p|

  int_Vector pol_value (pool.size()); // polynomials evaluated as $q=-1$
  for (auto it=pool.begin(); it!=pool.end(); ++it)
    if (it->isZero())
      pol_value[it-pool.begin()]=0;
    else
    { auto d=it->degree();
      int sum=(*it)[d];
      while (d-->0)
	sum=(*it)[d]-sum;
      pol_value[it-pool.begin()]=sum;
    }

  P_mat = int_Matrix(size);
  for (BlockElt x=0; x<size; ++x) // could stop at |size-1|
    for (BlockElt y=x+1; y<size; ++y)
    { auto pair= twisted_KLV.KL_pol_index(x,y);
      P_mat(x,y) = // |pol_value| at index of |it|, negated if |pair.second|
	pair.second ? -pol_value[pair.first] : pol_value[pair.first];
    }

/*
  singular blocks must be condensed by pushing down rows with singular
  descents to those descents (with negative sign as is implicit in |P_mat|),
  recursively, until reaching a "survivor" row without singular descents.
  Then afterwards non surviving rows and columns (should be 0) are removed.
*/

  containers::simple_list<BlockElt> survivors = eblock.condense(P_mat,B);

  BlockEltList compressed (survivors.wcbegin(), survivors.wcend());
  if (compressed.size()<size) // there were non survivors, so compress |P_mat|
  { size=compressed.size(); // henceforth this is our size
    int_Matrix M (size,size);
    for (BlockElt i=0; i<M.numRows(); ++i)
    { auto comp_i = compressed[i];
      for (BlockElt j=0; j<M.numColumns(); ++j)
	M(i,j)=P_mat(comp_i,compressed[j]);
    }
    P_mat = std::move(M); // replace |P_mat| by its expunged version
  }

  // flip signs for odd length distance, since that is what deformation wants
  for (BlockElt i=0; i<P_mat.numRows(); ++i)
  { auto parity = eblock.length(compressed[i])%2;
    for (BlockElt j=0; j<P_mat.numColumns(); ++j)
      if (eblock.length(compressed[j])%2!=parity)
	P_mat(i,j) *= -1;
  }

  block_list.clear(); block_list.reserve(size);
  lengths = int_Vector(0); lengths.reserve(size);

  const auto gamma = B.gamma();
  assert(is_dominant_ratweight(rc.rootDatum(),gamma)); // from |param_block|
#ifndef incompletecpp11
  for (auto ez : compressed)
  {
   auto z = eblock.z(ez);
 #else
  for (auto it=compressed.begin(); it!=compressed.end(); ++it)
  {
    auto z= eblock.z(*it);
#endif
    block_list.push_back(rc.sr_gamma(B.x(z),B.lambda_rho(z),gamma));
    lengths.push_back(B.length(z));
  }


} // |ext_KL_matrix|

} // |namespace kl|
} // |namespace atlas|
