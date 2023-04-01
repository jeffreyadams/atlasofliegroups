/*
  This is ext_kl.cpp

  Copyright 2013-2016, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ext_kl.h"
#include "kl.h" // for presence of |kl::KL_table| in |check_polys|

#include "basic_io.h"
#include "wgraph.h"

namespace atlas {
namespace ext_kl {


descent_table::descent_table(const ext_block::ext_block& eb)
  : info()
  , prim_index(1<<eb.rank(),std::vector<unsigned int>(eb.size(),0))
  , prim_flip(eb.size(),BitMap(prim_index.size()))
  , block(eb)
{
  info.reserve(block.size());
  for (BlockElt x=0; x<block.size(); ++x)
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
    info.emplace_back(desc,good_asc);
  }

  for (auto& prindex_vec : prim_index)
  {
    auto descs = &prindex_vec - &prim_index[0]; // position within |prim_index|
    constexpr unsigned int dead_end = -1;

    BlockElt count = 0;
    // following loop must decrease for primitivisation calculation below
    for (BlockElt x = block.size(); x-->0; )
    {
    // store index of primitivized |x| among primitives for |RankFlags(descs)|
    // since |x| is decreasing, initially count _larger_ primitive elements
      BitMap& flip_x = prim_flip[x]; // place to record primitivisation flips $x$
      unsigned int& dest = prindex_vec[x]; // the slot to fill

      const RankFlags D = good_ascent_set(x) & RankFlags(descs);

      if (D.none()) // then element |x| is primitive for the descent set
      {
	dest = count++; // store a self-reference, then increment count
	continue;
      }

      const weyl::Generator s = D.firstBit();
      if (is_like_nonparity(block.descent_type(s,x)))
      {
	dest = dead_end; // stop primitivisation with zero result
	continue;
      }

      BlockElt sx = block.some_scent(s,x);
      if (sx==UndefBlock) // primitivization would cross partial block edge
      {
	dest = dead_end; // stop primitivisation with zero result
	continue;
      }

      assert(sx>x); // ascents go up in block
      dest = prindex_vec[sx]; // |x| has the same primitivization as |sx|
      flip_x.set_to(descs, // and a flip that is relative to that of |sx|
		    (block.epsilon(s,x,sx)<0)!=prim_flip[sx].isMember(descs));
    } // |for(x)|, decreasing
    // primitive lists are actually to be stored increasing, so reverse indices

    const BlockElt last=count-1;
    for (unsigned int& slot : prindex_vec)
      if (slot != dead_end) // leave dead end indices as such
	slot = last-slot; // reverse all other indices

  } // |for (desc)|

} // |descent_table| constructor

// number of primitive elements for |descent_set(y)| of length less than |y|
unsigned int descent_table::col_size(BlockElt y) const
{
  BlockElt x=length_floor(y);
  if (prim_back_up(x,y)) // find last primitive |x| of length less than |y|
    return x_index(x,y)+1; // size in one more than index of that |x| for |y|
  return 0; // no primitives below length of |y| at all
} // |descent_table::col_size|

bool descent_table::prim_back_up(BlockElt& x, BlockElt y) const
{
  RankFlags desc=descent_set(y);
  while (x-->0)
    if ((good_ascent_set(x) & desc).none()) // disjoint sets
      return true;
  return false; // in which case |x| has crashed through 0 and should be ignored
} // |descent_table::prim_back_up|

bool descent_table::extr_back_up(BlockElt& x, BlockElt y) const
{
  RankFlags desc=descent_set(y);
  while (x-->0)
    if (descent_set(x).contains(desc)) // ascent set of |x| disjoint from |desc|
      return true; // stop when no descents of |y| are (any) ascents of |x|
  return false; // in which case |x| has crashed through 0 and should be ignored
} // |descent_table::extr_back_up|

KL_table::KL_table(const ext_block::ext_block& b, ext_KL_hash_Table* pol_hash)
  : aux(b)
  , pol_hash(pol_hash)
  , own(pol_hash!=nullptr ? nullptr : new IntPolEntry::Pooltype {Pol(0),Pol(1)})
  , storage_pool(pol_hash!=nullptr ? pol_hash->pool() : *own )
  , column(b.size(),KLColumn()) // start with empty columns
{ // ensure first two pool entries are constant polynomials $0$, and $1$
  if (pol_hash!=nullptr and pol_hash->size()<2)
  {
    assert(pol_hash->size()==0);
    pol_hash->match(Pol(0));
    pol_hash->match(Pol(1));
  }
  assert(storage_pool[zero]==Pol(0));
  assert(storage_pool[one] ==Pol(1));
}

std::pair<ext_kl::KLIndex,bool>
  KL_table::KL_pol_index(BlockElt x, BlockElt y) const
{ const KLColumn& col_y = column[y];
  unsigned inx=aux.x_index(x,y);
  if (inx<col_y.size())
    return std::make_pair(col_y[inx],aux.flips(x,y));
  else if (inx==aux.self_index(y)) // diagonal entries are unrecorded
    return {one,aux.flips(x,y)};
  else
    return {zero,false}; // out of bounds implies zero
}

Pol KL_table::P(BlockElt x, BlockElt y) const
{
  auto index = KL_pol_index(x,y);
  assert(index.first<storage_pool.size());
  return index.second ? -storage_pool[index.first] : storage_pool[index.first];
}

containers::sl_list<BlockElt> KL_table::nonzero_column(BlockElt y) const
{
  const KLColumn& col_y = column[y];
  containers::sl_list<BlockElt> result({y});
  for (BlockElt x=aux.length_floor(y); x-->0;)
  {
    unsigned inx=aux.x_index(x,y);
    if (inx<col_y.size() ? col_y[inx]!=zero : inx==aux.self_index(y))
      result.push_back(x);
  }
  return result;
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
  containers::sl_list<BlockElt> neighbours;
  aux.block.add_neighbours(neighbours,s,x);
  Pol result=aux.block.T_coef(s,x,x)*P(x,sy); // start with term from diagonal
  for (BlockElt sx : neighbours)
    result += aux.block.T_coef(s,x,sx)*P(sx,sy);
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

#ifndef NDEBUG // use |get_M| only for double-checking the result of |extract_M|
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
  const BlockElt z =  bl.some_scent(s,y); // ascent by |s| of |y|
  assert(z!=UndefBlock); // since in caller |y| was obtained as an |s|-descent

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
#endif

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

Poly_hash_export KL_table::polynomial_hash_table ()
{
  return pol_hash!=nullptr ? Poly_hash_export(pol_hash) : Poly_hash_export(*own);
}

// ensure all columns |y<limit| are computed
void KL_table::fill_columns(BlockElt limit)
{
  auto hash_object = polynomial_hash_table();
  if (limit==0)
    limit=aux.block.size(); // fill whole block if no explicit stop was indicated
  for (BlockElt y=aux.block.length_first(1); y<limit; ++y)
    if (column[y].size()!=aux.col_size(y))
    { assert(column[y].empty()); // there should not be partially filled columns
      try
      {
	fill_column(y,hash_object.ref);
      }
      catch(...)
      {
	column[y].clear(); // ensure partially filled columns are removed
      }
    }

  assert(check_polys(limit));
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

  void KL_table::fill_column(BlockElt y,PolHash& hash)
{
  // initialise column with dummy zero values; necessary for backwards filling
  column[y].assign(aux.col_size(y),ext_kl::KLIndex(0));

  weyl::Generator s;
  BlockElt sy; // gets set to unique descent for |s| of |y|, if one can be found
  if (has_direct_recursion(y,s,sy))
  {
    const unsigned defect = has_defect(type(s,y)) ? 1 : 0;
    const int sign = aux.block.epsilon(s,sy,y);
    const BlockElt floor_y =aux.length_floor(y);

    std::vector<IntPolEntry> cy(floor_y,(IntPolEntry()));
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
	Ms[u]=extract_M(cy[u],d,defect);
	assert(Ms[u]==get_M(s,u,sy,Ms));

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
    KLColumn::reverse_iterator it = column[y].rbegin();
    for (BlockElt x=floor_y; aux.prim_back_up(x,y); it++)
      if (aux.is_descent(s,x)) // then we computed $P(x,y)$ above
        *it = hash.match(cy[x]*sign);
      else // |s| might not be descent for |x| if it's primitive but not extremal
      { // use the double-valued ascent |s| for |x| that is descent for |y|
	assert(has_double_image(type(s,x))); // since |s| non-good ascent
	BlockEltPair sx = aux.block.Cayleys(s,x);
	if (sx.first==UndefBlock) // if we cross the edge of a partial block
	  *it=0; // then there can be no contribution form the Cayleys
	else
	{
	  IntPolEntry Q = P(sx.first,y);  // computed earlier in this loop
	  if (aux.block.epsilon(s,x,sx.first)<0)
	    Q *= -1;
	  if (sx.second!=UndefBlock)
	  {
	    if (aux.block.epsilon(s,x,sx.second)>0)
	      Q += P(sx.second,y);
	    else
	      Q -= P(sx.second,y);
	  }
	  *it = hash.match(Q);
	}
      }
    assert(it==column[y].rend()); // check that we've traversed the column
  } // end of |if (has_direct_recursion(y,s,sy))|
  else // direct recursion was not possible
    do_new_recursion(y,hash);
 } // |KL_table::fill_column|

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

  struct non_parity_info { weyl::Generator s; std::vector<Pol> M; };
  containers::sl_list<non_parity_info> rn_for_y;
  for (weyl::Generator s=0; s<rank(); ++s)
    if (is_like_nonparity(type(s,y)))
      rn_for_y.push_back(non_parity_info{s,std::vector<Pol>(floor_y,Pol())});

#ifndef NDEBUG
  { // check the absence of elements for which new recursion would not work
    BlockEltList downs = aux.block.down_set(y);
    for (const BlockElt u : downs)
      for (const auto& info : rn_for_y)
      {
	const auto tsu=type(info.s,u);
	if (not is_descent(tsu) and has_defect(tsu)) // defect ascent: not-good
	{
	  auto Csu = aux.block.Cayley(info.s,u);
	  std::cerr << "Bad element " << aux.block.z(u)
		    << "in down-set for " << aux.block.z(y)
		    << "; would need M_" << info.s+1 << '['
		    << (Csu<aux.block.size() ? aux.block.z(Csu) : UndefBlock)
		    << "] at defect ascent\n";
	  assert(false);
	}
      } // |for(i)|, loop over |s| and |for(u:downs)|, check downset
  }
#endif

  // for the primitive |x| we transfer to |column.back()| so |P(xx,y)| works
  auto out_it = column[y].rbegin();
  for (BlockElt x=floor_y; x-->0; )
  {
    if (is_primitive(x,y))
    { // compute $P_{x,y}$ for all |x| primitive for |y|, and store at |*out_it|
      if (not is_extremal(x,y)) // primitive though not extremal
      { // combine contributions from (at most) two Cayley ascents by hand
	weyl::Generator s = aux.easy_set(x,y).firstBit();
	assert(s<rank() and not aux.is_descent(s,x)); // since |x| not extremal
	assert(has_double_image(type(s,x))); // as ascent |s| for |x| not good
	BlockEltPair sx = aux.block.Cayleys(s,x);
	auto Pxy = Pol(0);
	if (sx.first!=UndefBlock)
	{ // although |P(x,y)| cannot be called yet, |P(x',y)| for |x'>x| is OK
	  Pxy = P(sx.first,y); // computed earlier this loop
	  if (aux.block.epsilon(s,x,sx.first)<0)
	    Pxy *= -1;
	  if (sx.second!=UndefBlock)
	  {
	    if (aux.block.epsilon(s,x,sx.second)>0)
	      Pxy += P(sx.second,y);
	    else
	      Pxy -= P(sx.second,y);
	  }
	}
	*out_it = hash.match(Pxy); // store result in primitive (only) case
      } // end of "then" branch for |if (not is_extremal(x,y))|
      else // |x| is extremal for |y|, so we must do real computation
      { // first seek proper |s| that is real non-parity for |y|
	const non_parity_info* info_ptr=nullptr; ext_block::DescValue tsx;
	for (const auto& info : rn_for_y)
	{ tsx=type(info.s,x); // this will also be reused in case of |break|
	  if (is_proper_ascent(tsx))
	  { // consider only |s| that are proper (not rn) ascents for |x|
	    if (not is_like_type_1(tsx)) // then |s| is certainly good
	    { info_ptr = &info; break; }
	    BlockElt csx =
	      aux.block.cross(info.s,x); // cross neighbour might help
	    if (csx!=UndefBlock and not is_extremal(csx,y)) // if so, endgame
	    { info_ptr = &info; break; } // which we consider good as well
	  }
	  else if (is_like_compact(tsx)) // also accept imaginary compact |x|
	  { info_ptr = &info; break; }
	}

	auto Q = Pol(0); // default value for when no good |s| was found
	if (info_ptr!=nullptr) // |break| above, found |s| with good type |tsx|
	{ // we still have |info_ptr| and |tsx==type(info_ptr->s,x)| at hand
	  const auto s = info_ptr->s;
	  const auto& M = info_ptr->M;
	  const unsigned k = aux.block.orbit(s).length();
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
	      // if |sx==UndefBlock| the next condition always fails
	      if (sx<floor_y) // subtract contr. from $[T_x](T_s+1).T_{sx}=q^k$
		Q -= aux.block.T_coef(s,x,sx)*P(sx,y); // coef is $\pm q^k$
	    } // implicit division of $Q$ here is by |T_coef(s,x,x)==1|
	    break;
	  case ext_block::two_semi_imaginary:
	  case ext_block::three_semi_imaginary:
	  case ext_block::three_imaginary_semi:
	    { // |has_defect(tsx)|
	      BlockElt sx=aux.block.Cayley(s,x);
	      if (sx<floor_y) // then in particular |sx!=UndefBlock|
		Q -= aux.block.T_coef(s,x,sx)*P(sx,y); // coef is $\pm(q^k-q)$

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
	      BlockEltPair sx=aux.block.Cayleys(s,x); // |UndefBlock|s are OK
	      if (sx.first<floor_y)
		Q -= aux.block.T_coef(s,x,sx.first) // $\pm(q^k-1)$
		  *P(sx.first,y);
	      if (sx.second<floor_y)
		Q -= aux.block.T_coef(s,x,sx.second)*P(sx.second,y); // idem
	      Q/=2; // divide by |T_coef(s,x,x)==2|
	    }
	    break;
	  case ext_block::one_imaginary_single:
	  case ext_block::two_imaginary_single_single:
	    { // |is_like_type_1(tsx)|, this used to be called the endgame case
	      BlockElt x_prime=aux.block.Cayley(s,x);
	      if (x_prime<floor_y)
		Q -= aux.block.T_coef(s,x,x_prime)*P(x_prime,y); // $\pm(q^k-1)$
	      // implicit division of $Q$ here is by |T_coef(s,x,x)==1|

	      // now subtract off $P_{s\times x,y}$, easily computed on the fly
	      BlockElt s_cross_x = aux.block.cross(s,x);
	      assert(not is_extremal(s_cross_x,y)); // this was tested above
	      const weyl::Generator t=aux.easy_set(s_cross_x,y).firstBit();
	      auto ttscx=type(t,s_cross_x);
	      if (has_double_image(ttscx))
	      {
		BlockEltPair sx_up_t = aux.block.Cayleys(t,s_cross_x);
		if (sx_up_t.first<floor_y)
		  Q -= P(sx_up_t.first,y)
		      *(aux.block.epsilon(s,x,s_cross_x)
		        *aux.block.epsilon(t,s_cross_x,sx_up_t.first));
		if (sx_up_t.second<floor_y)
		  Q -= P(sx_up_t.second,y)
		      *(aux.block.epsilon(s,x,s_cross_x)
		        *aux.block.epsilon(t,s_cross_x,sx_up_t.second));
	      }
	      else
	      {
		BlockElt sx_up_t=aux.block.Cayley(t,s_cross_x);
		if (sx_up_t<floor_y)
		  Q -= P(sx_up_t,y)
		      *(aux.block.epsilon(s,x,s_cross_x)
		        *aux.block.epsilon(t,s_cross_x,sx_up_t));
	      }
	    }
	    break;
	  case ext_block::two_imaginary_single_double_fixed:
	    {
	      BlockEltPair sx=aux.block.Cayleys(s,x);
	      if (sx.first<floor_y)
		Q -=
		  aux.block.T_coef(s,x,sx.first)*P(sx.first,y); // $\pm(q^2-1)$
	      if (sx.second<floor_y)
		Q -= aux.block.T_coef(s,x,sx.second)*P(sx.second,y); // idem
	      Q/=2; // divide by |T_coef(s,x,x)==2|
	    }
	    break;
	  case ext_block::one_imaginary_compact:
	  case ext_block::one_real_pair_switched:
	  case ext_block::two_imaginary_compact:
	  case ext_block::two_real_single_double_switched:
	  case ext_block::three_imaginary_compact:
	    // these cases require no additional terms to be subtracted
	    Q.factor_by_1_plus_q_to_the(k,(aux.block.l(y,x)-1)/2+k); // degree
	    assert(Q.degree_less_than((aux.block.l(y,x)+1)/2));
	    // that was the condition $\deg(Q) \leq l(y,x)-1/2|, computed safely
	    break;
	  default: assert(false); // other cases should not have selected |s|
	  } // |switch(tsx)|
	  // now |Q| is completely computed
	} // |if (info_ptr!=nullptr)|
	*out_it = hash.match(Q); // extremal implies primitive: store result
      } // end of |else| of |if (not is_extremal(x,y))|
      ++out_it; // output iterator advance only for primitive elements
    } // end of |if (is_primitive(x,y)|; the remainder is done is for all |x|

    // now if there is a defect ascent from |x|, update |M| for |mu(1,x,y)|
    if (aux.block.l(y,x)==1+2*P(x,y).degree())
      for (auto& info : rn_for_y)
      {
	const weyl::Generator s=info.s;	auto& M = info.M;
	const ext_block::DescValue tsx=type(s,x);
	if (not is_descent(tsx) and has_defect(tsx))
	{
	  const BlockElt sx = aux.block.Cayley(s,x);
	  assert(sx<floor_y); // could only fail if |x| in downset, tested above
	  int mu = P(x,y).coef(aux.block.l(y,x)/2) * aux.block.epsilon(s,x,sx);
	  M[sx] += Pol (M[sx].degree()==2 ? 1 : 0, mu);
	}
      }
    // and update the entries |M[x]| for every |M| in |rn_for_y|
    for (auto& info : rn_for_y)
      if (is_descent(type(info.s,x)))
	info.M[x] = get_Mp(info.s,x,y,info.M);
  } // |for(x)|

  assert(out_it==column[y].rend()); // check that we've traversed the column
} // |KL_table::do_new_recursion|


void KL_table::swallow(KL_table&& sub, const BlockEltList& embed)
{
  if (pol_hash!=nullptr and sub.pol_hash!=nullptr and
      &pol_hash->pool()==&sub.pol_hash->pool()) // case of shared hash tables
  {
    for (BlockElt y=0; y<sub.aux.block.size(); ++y)
      if (sub.column[y].size()==sub.aux.col_size(y) and column[embed[y]].empty())
      { // then transfer |sub.column[y]| to new block
	auto& cur_col=column[embed[y]];
	cur_col.assign(aux.col_size(embed[y]),zero); // default to 0
	BlockElt x=sub.aux.length_floor(y);
	const auto desc = sub.descent_set(y);
	assert(desc==descent_set(embed[y]));
	for (auto it=sub.column[y].crbegin(); sub.aux.prim_back_up(x,desc); ++it)
	  cur_col.at(aux.x_index(embed[x],embed[y])) = *it;
      }
    return;
  }
  // distict polynomial hash tables requires setting up polynomial translation
  std::vector<ext_kl::KLIndex> poly_trans(sub.storage_pool.size(),KLIndex(-1));
  {
    auto hash_object = polynomial_hash_table ();
    auto& hash = hash_object.ref;
    for (const auto& c : sub.column)
      for (auto ind : c)
	poly_trans[ind] = hash.match(sub.storage_pool[ind]);
    // besides filling |poly_trans| this also extends |storage_pool| as needed
  }

  // it remains to do the same as above but passing values through |poly_trans|
  for (BlockElt y=0; y<sub.aux.block.size(); ++y)
    if (sub.column[y].size()==sub.aux.col_size(y) and column[embed[y]].empty())
    { // then transfer |sub.column[y]| to new block
      auto& cur_col=column[embed[y]];
      cur_col.assign(aux.col_size(embed[y]),zero); // default to 0
      BlockElt x=sub.aux.length_floor(y);
      const auto desc = sub.descent_set(y);
      assert(desc==descent_set(embed[y]));
      for (auto it=sub.column[y].crbegin(); sub.aux.prim_back_up(x,desc); ++it)
	cur_col.at(aux.x_index(embed[x],embed[y])) = poly_trans[*it];
    }
}

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

bool KL_table::check_polys(BlockElt limit) const
{
  bool result = true;
#ifndef NDEBUG
  kl::KL_table untwisted(aux.block.untwisted());
  untwisted.fill(); // fill KL table silently
  for (BlockElt y=aux.block.length_first(1); y<limit; ++y)
    for (BlockElt x=y; x-->0; )
      if (not check(P(x,y),untwisted.KL_pol(aux.block.z(x),aux.block.z(y))))
      {
	std::cerr << "Mismatch at (" << aux.block.z(x) << ',' << aux.block.z(y)
		  << "): ";
	std::cerr << P(x,y) << " and "
		  << untwisted.KL_pol(aux.block.z(x),aux.block.z(y))
		  << std::endl;
	result=false;
      }
#endif
  return result;
}

// this function serves uniquely to implement the built-in function of that name
void ext_KL_matrix (const StandardRepr p, const int_Matrix& delta,
		    const Rep_context& rc, // the rest is output
		    std::vector<StandardRepr>& block_list,
		    int_Matrix& P_index_mat,
		    IntPolEntry::Pooltype& polys)
{ BlockElt entry_element;
  if (not ((delta-1)*p.gamma().numerator()).isZero())
  {
    std::cout << "Delta does not fix gamma=" << p.gamma() << "." << std::endl;
    throw std::runtime_error("No valid extended bock");
  }
  auto srm = repr::StandardReprMod::mod_reduce(rc,p); // modular |z|
  common_context ctxt(rc,srm.gamma_lambda());
  blocks::common_block B(ctxt,srm,entry_element); // build full block
  const auto& gamma = p.gamma();
  RatWeight diff(gamma.size()); // a custom made |StandardReprMod| gives zero offset
  assert(is_dominant_ratweight(rc.root_datum(),gamma)); // from |common_block|
  const RankFlags singular = B.singular(gamma);
  ext_block::ext_block eblock(B,delta,nullptr);

  BlockElt size= // size of extended block we shall use; before compression
    eblock.element(entry_element+1);

  IntPolEntry::Pooltype pool; // keep a separate pool, |polys| is for later
  ext_KL_hash_Table hash(pool,4);
  KL_table twisted_KLV(eblock,&hash);
  twisted_KLV.fill_columns(size); // fill up to and including |p|

  matrix::Matrix<Pol> P_mat(size); // will hold expanded KL matrix
  for (BlockElt x=0; x<size; ++x) // could stop at |size-1|
    for (BlockElt y=x+1; y<size; ++y)
      P_mat(x,y) = twisted_KLV.P(x,y);

/*
  singular blocks must be condensed by pushing down rows with singular
  descents to those descents (with negative sign as is implicit in |P_mat|),
  recursively, until reaching a "survivor" row without singular descents.
  Then afterwards non surviving rows and columns (should be 0) are removed.
*/

  containers::sl_list<BlockElt> survivors // convert from |simple_list|
    (eblock.condense(P_mat,eblock.singular_orbits(singular)));

  if (survivors.size()<size) // if any non-survivors, we need to compress |P_mat|
  { size=survivors.size(); // henceforth this is our size
    matrix::Matrix<Pol> M (size,size,Pol(0));
    BlockElt i=0;
    for (auto it=survivors.wbegin(); not survivors.at_end(it); ++it,++i)
    { BlockElt j=i; // only upper triangular part can be nonzero
      for (auto jt=it; not survivors.at_end(jt); ++jt,++j)
	M(i,j) = std::move(P_mat(*it,*jt));
    }
    P_mat = std::move(M); // replace |P_mat| by its expunged version
  }

  { // flip signs for odd length distance, since that is what deformation wants
    auto jt = survivors.wcbegin();
    for (BlockElt j=0; j<P_mat.n_rows(); ++j,++jt)
    { auto parity = eblock.length(*jt)%2;
      auto it = survivors.wcbegin();
      for (BlockElt i=0; i<j; ++i,++it)
	if (eblock.length(*it)%2!=parity)
	  P_mat(i,j) *= -1;
    }
  }

  block_list.clear(); block_list.reserve(size);
  for (auto ez : survivors)
    block_list.push_back(rc.sr(B.representative(eblock.z(ez)),diff,gamma));

  polys.assign({Pol(),Pol(1)}); // set up initial values of table
  P_index_mat = int_Matrix(size);
  {
    ext_KL_hash_Table pol_hash(polys,4);
    for (BlockElt j=1; j<size; ++j)
      for (BlockElt i=0; i<j; ++i)
	P_index_mat(i,j) = pol_hash.match(P_mat(i,j));
  }

} // |ext_KL_matrix|

} // |namespace kl|
} // |namespace atlas|
