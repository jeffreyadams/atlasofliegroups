/*
  This is tori.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


// implementation of real torus related functions, and |class RealTorus|


#include "tori.h"

#include <algorithm>
#include <cassert>

#include "matreduc.h"
#include "bitmap.h" // for destructing the result of |matreduc::column_echelon|
#include "lattice.h"

namespace atlas {
namespace tori {

namespace {

  WeightList fullSignBasis(const WeightInvolution& tau,int sign,int_Matrix& M);

} // |namespace|


/*****************************************************************************

        Chapter I -- The RealTorus class

  The real torus class represents the datum of a torus $H=H(\R)$ defined over
  $\R$. This is equivalent to the datum of a (character) lattice equipped with
  an involution. To be consistent with the rest of the program, we use the
  Cartan involution $\tau$ which is the negative of the Galois involution on the
  character lattice. A fundamental fact is that the component group of $H(R)$
  (the fixed points of the Galois involution) is canonically isomorphic to the
  component group of the complex group $H(C)^\tau$ fixed under the Cartan
  involution (even though the initial groups can be of quite different
  dimensions), the latter of which component groups is easy to determine in
  terms of the action of $\tau$ on the character lattice.

  The fundamental data are the rank stored in |d_rank| (allowing us to represent
  the character lattice $X=X^*$ by $Z^n$ for $n=d_rank)$, and the action of
  $\tau$ on $X$ by the integral $n*n$ involution matrix |d_involution|. It
  determines two sublattices $X_+$ and $X_-$, cut out by its eigenspaces for
  $+1$ and $-1$ respectively. Both are sublattices are saturated, and therefore
  supplementable in $X$, but in general $X$ is not equal to their direct sum
  $X_+ + X_-$, since that sum need not be saturated (for instance in $Z^2$ when
  $\tau$ interchanges the vectors of a basis).

  The most delicate invariant we shall have to deal with is the component group
  of the real torus $H(R)$. This is an elementary abelian 2-group; we shall
  rather consider its dual group $dpi0(H)$, which like $X$ lives on the weight
  side. To describe its rank is fairly easy. Indeed, $H$ may be decomposed as a
  product of compact, split and complex factors; this corresponds to a
  decomposition of $X$ into a direct sum of sublattices on which $\tau$ acts
  respectively as $1$ (identity), $-1$, and by a permutation of some basis with
  only 2-cycles. Such a decomposition is not canonical, but the number of
  factors of each type is uniquely determined; call them $r_u$, $r_s$ and $r_c$
  (they give the ranks of the first two sublattices, and half the rank of the
  third, $n = r_u + r_s + 2r_c$). Then the rank of the component group is $r_s$.

  We refrain from seeking any non-canonical decomposition of $X$, and instead
  work with the subspaces $X_+$ and $X_-$. One has $rk(X_+) = r_u+r_c$ and
  $rk(X_-) = r_s+r_c$. Put $V = X/2X$, a vector space over the two-element field
  $Z/2Z$, and denote by $V_+,V_-$ the respective images of $X_+, X_-$ in $V$.
  Then factors $r_u$ and $r_s$ contribute one-dimensional subspaces to $V_+$ and
  $V_-$, respectively, while a factor $r_c$ contributes the \emph{same}
  one-dimensional subspace both to $V_+$ and to $V_-$ (since after reduction
  modulo $2$ the is no distinction between the diagonal and anti-diagonal
  subpaces). Therefore one has $r_c = \dim(V_+ \cap V_-)$, and it follows that
  $r_s = rk(X_-) - r_c=\dim(V_-/(V_+ \cap V_-))$. It can be seen that $V_+ \cap
  V_-$ is the image of the endomorphism of $V$ induced by $\tau-1$ (or by
  $\tau+1$), while $V_+ + V_-$ is the kernel of that endomorphism.

  Having used only $X_+$ and $X_-$ we almost have a functorial description of
  the component group of $H(R)$, which is an elementary 2-group isomorphic to
  $V_-/(V_ \cap V_-)$, or equivalently to $(V_+ + V_-)/V_+$. However, having
  started with the weight lattice $X=X^*$, the latter groups are actually dual
  objects: $X/2X$ can be interpreted as $\Hom(H(2),\C^*)$ where $H(2)$ denotes
  the set of elements of order 2 (or 1) in $H$. Therefore each of the mentioned
  subquotients naturally models the dual group $dpi0(H)$ of the component group
  of $H(\R)$ (or of $H(\C)^\tau$). In practice we shall work only with these
  dual groups, and use the subquotient $(V_+ + V_-)/V_+$ to represent it (the
  latter point is relevant, because even when two subquotients are canonically
  isomorphic, they cannot be used interchangeably in all situations).

******************************************************************************/

/*
  Construct the torus with involution |i| (the rank is recovered
  from the size of the matrix.)
*/
RealTorus::RealTorus(const WeightInvolution& theta)
  : d_rank(theta.numColumns())
  , d_complexRank(d_rank) // this value is too large; it is modified below
  , d_involution(theta)
  , d_toPlus()
  , d_toMinus()
  , d_plus (fullSignBasis(theta,+1,d_toPlus ))
  , d_minus(fullSignBasis(theta,-1,d_toMinus))
  , d_topology(dualPi0(theta))
{
  /* one has $n=r_u+r_s+2r_c$ while $\dim(V_+ + V_-)=r_u+r_s+r_c$, so $r_c$
     can be recovered by subtracting the latter from the former: */
  d_complexRank -= d_topology.space().dimension();
}

RealTorus RealTorus::dual() const
{ return RealTorus(d_involution.negative_transposed()); }


/******** accessors *********************************************************/


/*
  Construct the map |cm| induced, at the level of dual component groups of
  real tori, by a map |m| from the weight lattice of our own |RealTorus| to
  the one of |T_dest|.

  The idea is simple enough: look at the columns of m modulo 2 and project
  them onto dpi0, and then restrict to our dual component representatives.
  But in fact this is hidden inside |subquotient::subquotientMap|.
*/
BinaryMap RealTorus::componentMap
  (const LatticeMatrix& m,
   const RealTorus& T_dest)
  const
{
  BinaryMap m2(m); // reduce mod2

  return subquotient::subquotientMap
    (dual_component_group(),T_dest.dual_component_group(),m2);
}


/*****************************************************************************

        Chapter II -- Functions declared in tori.h

******************************************************************************/

/*
  Compute $(V_+ + V_-) / V_+$ where $V_+$ is the image mod 2 of the $q$-fixed
  sublattice ($+1$-eigenspace) and $V_-$ of the $-q$-fixed sublattice.

  When $q$ describes a real torus, this is the dual of the component group of
  the torus, whence the curious name.

  Actually the numerator subgroup is computed as the kernel of $q \mod 2$.

  In practice we also want this at the dual side with $V_-$ in the denominator
  (to compute fiber groups): to that end call with |q.negative_transposed()|
*/
SmallSubquotient dualPi0(const WeightInvolution& q)
{
  assert(q.numRows()==q.numColumns());
  auto rank = q.numRows();

  BinaryMap i2(q);  // reduce modulo 2
  i2 += BinaryMap::identity(rank);
  // now |i2| is mod-2 image of |theta-1| (and also of |theta+1|)
  // whose kernel is the sum $V_+ + V_-$, which will be our numerator subgroup

  SmallSubquotient cs(i2.kernel(),plusBasis(q),rank);
  assert(cs.rank()==rank);

  return cs;
}

// Return a basis for the +1 eigenlattice of |theta|
WeightList plusBasis(const WeightInvolution& theta)
{
  return lattice::eigen_lattice(theta,1).columns();
}

WeightList minusBasis(const WeightInvolution& theta)
{
  return lattice::eigen_lattice(theta,-1).columns();
}

std::tuple<unsigned,unsigned,unsigned> classify(const WeightInvolution& tau)
{
  int_Matrix tau1=tau+1,dummy; bool flip;
  matreduc::column_echelon(tau1,dummy,flip);
  const auto r=tau1.numRows(),plus_rank = tau1.numColumns();
  const auto complex_rank = BinaryMap(tau1).image().size();
  const auto compact_rank = plus_rank-complex_rank;
  return std::make_tuple(compact_rank,complex_rank,r-plus_rank-complex_rank);
}

/*****************************************************************************

        Chapter III -- Auxiliary functions

******************************************************************************/

namespace {


/*
  Return a basis for the |sign| ($\in\{-1,+1\}$) eigenlattice of the involution,
  while storing in |M| the matrix of a coordinate transformation from the
  standard basis to the resturned basis, valid for weights in the eigenlattice.

  Algorithm: a maximal rank sublattice of the eigenlattice in question is
  spanned by the columns of the matrix $\tau\mp1$. An adapted basis for this
  sublattice will start with a basis for the eigenlattice (whose cardinality is
  given by the number of invariant factors), followed by a basis for a
  complement. The inverse of the matrix with as columns the adapted basis gives
  a coordinate transform from the standard basis to the adapted basis; retaining
  only the rows that give the initial coordinates we obtain the matrix |M|.
*/
WeightList fullSignBasis(const WeightInvolution& tau,int sign,int_Matrix& M)
{
  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(tau+sign,factor);

  size_t r=factor.size();

  // copy significant part of basis in mb
  WeightList result;  result.reserve(r);
  for (size_t j=0; j<r; ++j)
    result.push_back(b.column(j));

  // make projection matrix; rows are the first rows of the inverse of |b|
  M = b.inverse().block(0,0,r,b.numRows());
  return result;
}

} // |namespace|

} // |namespace tori|

} // |namespace atlas|
