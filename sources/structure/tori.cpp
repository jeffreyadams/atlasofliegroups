/*!
\file
\brief
Implementation of the class RealTorus

*/
/*
  This is tori.cpp

  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "tori.h"

#include <algorithm>
#include <cassert>

#include "matreduc.h"

namespace atlas {

namespace tori {

namespace {

void makeTopology(SmallSubquotient&, const RealTorus&);

void fullPlusBasis(WeightList&,
		   int_Matrix&,
		   const WeightInvolution&);

void fullMinusBasis(WeightList&,
		    int_Matrix&,
		    const WeightInvolution&);

} // |namespace|

} // |namespace tori|

/*****************************************************************************

        Chapter I -- The RealTorus class

  The real torus class represents the datum of a torus $H$ defined over $R$.
  This is equivalent to the datum of a (character) lattice equipped with an
  involution. To be consistent with the rest of the program, we use the Cartan
  involution $\tau$ which is the negative of the Galois involution on the
  character lattice. A fundamental fact is that the component group of $H(R)$
  (the fixed points of the Galois involution) is canonically isomorphic to the
  component group of the complex group $H(C)^\tau$ fixed under the Cartan
  involution (even thought the initial groups can be of quite different
  dimensions), the latter of which component groups is easy to determine in
  terms of the action of $\tau$ on the character lattice.

  The fundamental data are the rank stored in |d_rank| (allowing us to
  represent the character lattice $X$ as $Z^n$ for $n=d_rank)$, and the
  integral $n*n$ involution matrix |d_involution| giving the action of
  $\tau$ on $X$. It determines two sublattices $X_+$ and $X_-$, cut out by
  its eigenspaces for $+1$ and $-1$ respectively. Both are sublattices are
  saturated, and therefore supplementable in $X$, but in general $X$ is not
  equal to the direct sum $X_+ + X_-$, since that sum need not be saturated
  (for instance in $Z^2$ when $\tau$ interchanges the vectors of a basis).

  The most delicate invariant we shall have to deal with is the component
  group of the real torus $H(R)$. This is an elementary abelian 2-group; we
  shall rather consider its dual group $dpi0(H)$. To describe its rank is
  fairly easy. Indeed, $H$ may be decomposed as a product of compact, split
  and complex factors; this corresponds to a decomposition of $Z^n$ into a
  direct sum of sublattices on which $\tau$ acts respectively as $1$
  (identity), $-1$, and by a permutation of some basis with only 2-cycles.
  While such a decomposition is not canonical, the number of factors of each
  type is uniquely determined; denote them by $r_u$, $r_s$ and $r_c$ (they
  give the ranks of the first two sublattices, and half the rank of the
  third); the rank $n$ of $H$ is equal to $r_u + r_s + 2 r_c$. Then the rank
  of the component group is $r_s$.

  We shall refrain from seeking such a non-canonical decomposition of $X$, and
  work with the subspaces $X_+$ and $X_-$ instead. One has $rk(X_+) = r_u+r_c$
  and $rk(X_-) = r_s+r_c$. Denote V = $X/2X$, a vector space over the
  two-element field $Z/2Z$, and denote by $V_+,V_-$ the respective images of
  $X_+, X_-$ in $V$. Then factors $r_u$ and $r_s$ contribute one-dimensional
  subspaces to $V_+$ and $V_-$, respectively, while a factor $r_c$ contributes
  the _same_ one-dimensional subspace both to $V_+$ and to $V_-$ (since after
  reduction modulo $2$ the is no distinction between the diagonal and
  anti-diagonal subpaces). Therefore on has $r_c = \dim(V_ +\cap V_-)$,
  and it follows that $r_s = rk(X_-) - r_c=\dim(V_-/(V_ \cap V_-))$. It
  can be seen that $V_+ \cap V_-$ is the image of the endomorphism of $V$
  induced by $\tau-1$ (or by $\tau+1$), while $V_+ + V_-$ is the
  kernel of that endomorphism.

  Having used only $X_+$ and $X_-$ we almost have a functorial description of
  the component group of $H(R)$, which is an elementary 2-group isomorphic to
  $V_-/(V_ \cap V_-)$, or equivalently to $(V_+ + V_-)/V_+$. However,
  having started with the weight lattice $X$, the latter groups are actually
  dual objects: $X/2X$ can be interpreted as $Hom(H(2),C^*)$ where $H(2)$
  denotes the set of elements of order 2 (or 1) in $H$. Therefore each of the
  mentioned subquotients naturally models the dual group $dpi0(H)$ of the
  component group of $H(R)$ (or of $H(C)^\tau$). In practice we shall work
  only with these dual groups, and use the subquotient $(V_+ + V_-)/V_+$ to
  represent it (the latter precision is relevant, because even if the two
  subquotients are canonically isomorphic, they cannot be used interchangeably
  in all situations).

******************************************************************************/

namespace tori {

/*!
  \brief Constructs the torus with involution |i| (the rank is recovered
  from the size of the matrix.)
*/
RealTorus::RealTorus(const WeightInvolution& i)
  : d_rank(i.numColumns())
  , d_complexRank(d_rank) // this value is too large; it is modified below
  , d_involution(i)
  , d_plus()
  , d_minus()
  , d_toPlus()
  , d_topology(d_rank) // set dimension of ambient space for subquotient
{
  // make bases of +1 and -1 eigenlattices

  fullPlusBasis(d_plus,d_toPlus,d_involution);
  fullMinusBasis(d_minus,d_toMinus,d_involution);

  // find component group data

  makeTopology(d_topology,*this);

  /* one has $n=r_u+r_s+2r_c$ while $\dim(V_+ + V_-)=r_u+r_s+r_c$, so $r_c$
     can be recovered by subtracting the latter from the former: */
  d_complexRank -= d_topology.space().dimension();
}


/*!
  Synopsis: constructs the torus "dual" to T.

  This is the torus whose weight lattice is dual to that of T, equipped with
  the involution -tau^t, where tau is the involution for T. The minus sign is
  there because it is the way we wish to dualize our real forms.
*/
RealTorus::RealTorus(const RealTorus& T, tags::DualTag)
  : d_rank(T.d_rank)
  , d_complexRank(d_rank) // this value is too large; it is modified below
  , d_involution(d_rank,d_rank)
  , d_plus()
  , d_minus()
  , d_toPlus()
  , d_topology(d_rank) // set dimension of ambient space for subquotient
{
  // set the involution matrix

  for (size_t j = 0; j < d_rank; ++j)
    for (size_t i = 0; i < d_rank; ++i)
      d_involution(i,j) = -T.d_involution(j,i);

  // make bases of +1 and -1 eigenlattices

  fullPlusBasis(d_plus,d_toPlus,d_involution);
  fullMinusBasis(d_minus,d_toMinus,d_involution);

  // find component group data

  makeTopology(d_topology,*this);
  d_complexRank -= d_topology.space().dimension();
}

/******** accessors *********************************************************/


/*!
  \brief Constructs the map |cm| induced, at the level of dual component
  groups of real tori, by a map |m| from the weight lattice of our own
  |RealTorus| to the one of |T_dest|.

  The idea is simple enough: look at the columns of m modulo 2 and project
  them onto dpi0, and then restrict to our dual component representatives.
  But in fact this is hidden inside |subquotient::subquotientMap|.
*/
BinaryMap RealTorus::componentMap
  (const LatticeMatrix& m,
   const RealTorus& T_dest) const
{
  BinaryMap m2(m); // reduce mod2

  return subquotient::subquotientMap(d_topology,T_dest.topology(),m2);
}

/******** manipulators ******************************************************/

}

/*****************************************************************************

        Chapter II -- Functions declared in tori.h

******************************************************************************/

namespace tori {


/*!
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

  SmallBitVectorList plus2(plusBasis(q)); // mod 2: denominator subgroup

  BinaryMap i2(q);  // reduce modulo 2
  auto id = BinaryMap::identity(q.numRows());

  i2 += id; // now |i2| is mod-2 image of |theta-1| (and also of |theta+1|)

  SmallBitVectorList b=i2.kernel();
  // the kernel of |i2| is the sum $V_+ + V_-$: numerator subgroup

  SmallSubquotient cs(b,plus2,q.numRows());
  assert(cs.rank()==q.numRows());

  return cs;
}

/*!
  \brief Puts in mb a basis for the +1 eigenspace of the involution;

  Algorithm: the vectors e+i(e), when e runs through b, generate a lattice
  commensurate with the eigenspace.
*/
void plusBasis(WeightList& pb,
	       const WeightInvolution& i)

/*!
  Synopsis: puts in mb a basis for the +1 eigenspace of the involution;

  Algorithm: the vectors e+i(e), when e runs through b, generate a lattice
  commensurate with the eigenspace.
*/

{
  size_t n = i.numRows();

  // put in q the matrix of vectors i(e)+e in the basis b

  WeightInvolution q(i);
  for (size_t i=0; i<n; ++i)
    q(i,i) += 1;

  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(q,factor);

  size_t r=factor.size();

  // copy significant part of basis in pb
  pb.reserve(r);
  for (size_t j=0; j<r; ++j)
    pb.push_back(b.column(j));
}

WeightList plusBasis(const WeightInvolution& i)
{
  WeightList result; plusBasis(result,i); return result;
}

/*!
  Synopsis: puts in mb a basis for the -1 eigenspace of the involution.

  Algorithm: the vectors e-i(e), when e runs through b, generate a lattice
  commensurate with the eigenspace.
*/
void minusBasis(WeightList& mb,
		const WeightInvolution& i)
{
  size_t n = i.numRows();

  // put in q the matrix of vectors i(e)-e in the basis b

  WeightInvolution q(i);
  for (size_t i=0; i<n; ++i)
    q(i,i) -= 1;

  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(q,factor);

  size_t r=factor.size();

  // copy significant part of basis in pb
  mb.reserve(r);
  for (size_t j=0; j<r; ++j)
    mb.push_back(b.column(j));
}

WeightList minusBasis(const WeightInvolution& i)
{
  WeightList result; minusBasis(result,i); return result;
}

/*!
  Synopsis: writes in qm the matrix of the restriction of q to X_-.

  Precondition: q commutes with the involution;
*/
void minusMatrix(int_Matrix& qm,
		 const int_Matrix& q, const RealTorus& t)
{
  const WeightList& bm = t.minusLattice();
  int_Matrix(bm.size(),bm.size()).swap(qm);

  for (size_t j = 0; j < bm.size(); ++j)
  {
    Weight v= q*bm[j];
    Weight vm(bm.size());
    t.toMinus(vm,v);
    for (size_t i = 0; i < bm.size(); ++i)
      qm(i,j) = vm[i];
  }
}


/*!
  Synopsis: writes in qp the matrix of the restriction of q to X_+.

  Precondition: q commutes with the involution;
*/
void plusMatrix(WeightInvolution& qp,
		const WeightInvolution& q, const RealTorus& t)
{
  const WeightList& bp = t.plusLattice();
  WeightInvolution(bp.size(),bp.size()).swap(qp); // create

  for (size_t j = 0; j < bp.size(); ++j)
  {
    Weight v= q*bp[j];
    Weight vp(bp.size());
    t.toPlus(vp,v);
    for (size_t i = 0; i < bp.size(); ++i)
      qp(i,j) = vp[i];
  }
}

}  // |namespace tori|

/*****************************************************************************

        Chapter III -- Auxiliary functions

******************************************************************************/

namespace tori {

namespace {


/*!
  Synopsis: puts the subquotient V^tau/V_+ in cs.

  Explanation: V is X/2X; V_+ is the image of X_+ in V, similarly for V_-,
  and V^tau=V_+ + V_- is the kernel of the map induced by tau-1 on V.
  So the subquotient V^tau/V_+ is (V_+ + V_-)/V_+ which is canonically
  isomorphic to V_-/(V_+  cap V_-), which we know is the _dual_ of the
  component group of the torus.
*/
void makeTopology(SmallSubquotient& cs, const RealTorus& T)
{
  SmallBitVectorList plus2(T.plusLattice()); // reduce mod 2

  BinaryMap i2(T.involution()); // reduce mod 2
  auto id = BinaryMap::identity(T.rank());
  i2 += id;

  SmallBitVectorList b=i2.kernel(); // kernel of |theta+1| mod 2, contains V_+

  cs = SmallSubquotient(b,plus2,T.rank());
}


/*!
  Synopsis: puts in mb a basis for the -1 eigenspace of the involution; puts
  in tm the matrix of a projection onto X_-.

  The matrix tm is in terms of the standard basis of X, and the chosen basis
  of X_-. There is no significance to the chosen complement; the intention is
  that it should be used only for vectors already in X_-.

  Algorithm: the vectors e-i(e), when e runs through b, generate a lattice
  commensurate with the eigenspace. The basis adapted to this sublattice also
  yields the projection matrix, after inversion.
*/
void fullMinusBasis(WeightList& mb,
		    int_Matrix& tm,
		    const WeightInvolution& i)
{
  size_t n = i.numRows();

  // put in q the matrix of vectors i(e)-e in the basis b

  WeightInvolution q(i);
  for (size_t i=0; i<n; ++i)
    q(i,i) -= 1;

  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(q,factor);

  size_t r=factor.size();

  // copy significant part of basis in mb
  mb.reserve(r);
  for (size_t j=0; j<r; ++j)
    mb.push_back(b.column(j));

  // make projection matrix
  // rows of projection matrix are the first rows of the inverse of |b|
  tm = b.inverse().block(0,0,r,n);
}


/*!
  \brief puts in |pb| a basis for the +1 eigenlattice of the involution;
  puts in |tp| the matrix of a coordinate transformation from the standard
  basis to |pb|.

  Although |tp| can be applied to any element of $X$, the result is well
  defined only when applied to an element already in $X_+$; otherwise the
  image will have been projected parallel to an unspecified complement.

  Algorithm: a maximal rank sublattice of $X_+$ is spanned by the columns of
  the matrix \f$\tau+1\f$. A Smith basis for this sublattice will start with a
  basis for $X$ (whose cardinality is given by the number of invariant
  factors), followed by a basis for a complement. The inveres of the matrix
  with as columns the full Smith basis provides a coordinate transform from
  the standard basis to the Smith basis; retaining only the rows that give the
  initial coordinates we obtain the matrix |tp|.
*/
void fullPlusBasis(WeightList& pb,
		   int_Matrix& tp,
		   const WeightInvolution& tau)
{
  size_t n = tau.numRows();

  // put in q the matrix of vectors tau(e)+e in the basis b

  WeightInvolution q(tau);
  for (size_t i=0; i<n; ++i)
    q(i,i) += 1;

  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(q,factor);

  size_t r=factor.size();

  // copy significant part of basis in pb
  pb.reserve(r);
  for (size_t j=0; j<r; ++j)
    pb.push_back(b.column(j));

  // make coordinate transformation matrix
  tp=b.inverse().block(0,0,r,n);
}

} // |namespace|

} // |namespace tori|

} // |namespace atlas|
