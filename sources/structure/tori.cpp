/*!
\file
\brief
Implementation of the class RealTorus

*/
/*
  This is tori.cpp

  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "tori.h"

#include <algorithm>

#include "lattice.h"
#include "smithnormal.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace tori {

  using namespace latticetypes;

}

namespace {

using namespace tori;

void makeTopology(ComponentSubquotient&, const RealTorus&);

void fullPlusBasis(WeightList&, LatticeMatrix&, const LatticeMatrix&);

void fullMinusBasis(WeightList&, LatticeMatrix&, const LatticeMatrix&);

}

/*****************************************************************************

        Chapter I -- The RealTorus class

  ... the pi0 part is outdated ...

  The real torus class represents the datum of a torus defined over R. This
  is equivalent to the datum of a lattice with an involution; to be consistent
  with the rest of the program, we use the involution theta which is the
  negative of the Galois involution on the character lattice.

  We always assume that the character lattice X is the standard lattice Z^n,
  where n >= 0 is the rank. Then we have two sublattices X_+ and X_-, the
  +1 and -1 eigenspaces for the involution. Both are supplementable in X, but
  in general X is not equal to the direct sum X_+ + X_-.

  The most delicate invariant we shall have to deal with is the component
  group of the group of real points of T. This is an elementary abelian
  2-group; we shall rather consider its dual dpi0(T). To describe its rank
  is fairly easy : indeed, T may be decomposed as a product of compact, split
  and complex factors; denote r_u, r_s and r_c the number of factors of each
  type, so that the rank n of T is r_u + r_s + 2 r_c. Then we have : rk(X_+)
  = r_u + r_c; rk(X_-) = r_s + r_c. Denote V = X/2X, a vector space over the
  two-element field F_2, and denote V_+,V_- the images of X_+, X_- in V.
  Then it is not hard to show that r_c = dim(V_+\cap V_-); one may prove
  that this is also the image of theta - 1 in V.

  It is a little bit harder to describe dpi0(T) functorially as a vector space.
  Denote V_+- the intersection V_+ \cap V_-. Then one may show that the group
  T(2)(R) of real points of the group of elements of order 2 in T is in natural
  duality with V/V_+-. There is a natural surjection from T(2)(R) to pi0(T),
  so dpi0(T) is a sub-vector space of V/V_+-, which may in fact be described
  as the image of V_-. This is how we will consider it.

******************************************************************************/

namespace tori {

RealTorus::RealTorus(const LatticeMatrix& i)
  :d_rank(i.numColumns()),
   d_involution(i),
   d_topology(d_rank)

/*!
  This function constructs the torus with involution i (the rank is recovered
  from the size of the matrix.)
*/

{
  using namespace lattice;

  // make bases of +1 and -1 eigenlattices

  fullPlusBasis(d_plus,d_toPlus,d_involution);
  fullMinusBasis(d_minus,d_toMinus,d_involution);

  // find component group data

  makeTopology(d_topology,*this);
  d_complexRank = d_rank - d_topology.space().dimension();
}

RealTorus::RealTorus(const RealTorus& T, tags::DualTag)
  :d_rank(T.d_rank),
   d_involution(d_rank)

/*!
  Synopsis: constructs the torus dual to T.

  This is the torus in the dual lattice (but the lattice and the dual lattice
  are identical for us!) with involution -tau^t, where tau is the involution
  for T. The minus sign is there because it is the way we wish to dualize
  our real forms.
*/

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
  d_complexRank = d_rank - d_topology.space().dimension();
}

/******** accessors *********************************************************/

void RealTorus::componentMap(ComponentMap& cm, const LatticeMatrix& m, 
			     const RealTorus& T_dest) const

/*!
  Constructs the map induced at the level of dual component groups by a
  map from X_source to X_dest, given by its matrix m. The idea is simple
  enough : look at the columns of m modulo 2 and project them onto dpi0,
  and then restrict to component representatives.
*/

{
  using namespace lattice;
  using namespace subquotient;

  ComponentMap m2;
  mod2(m2,m);

  subquotientMap(cm,d_topology,T_dest.topology(),m2);

  return;
}

/******** manipulators ******************************************************/

}

/*****************************************************************************

        Chapter II -- Functions declared in tori.h

  ... explain here when it is stable ...

******************************************************************************/

namespace tori {

void dualPi0(LT::ComponentSubquotient& dpi0, const LT::LatticeMatrix& q)

/*!
  Synopsis: puts in dpi0 what would have been the topology field for the
  corresponding torus.

  Precondition: q is an involution matrix; its size doesn't exceed RankMax;

  Explanation: this is canonically the dual of the component group; see
  makeTopology.
*/

{
  using namespace lattice;

  WeightList plus;
  plusBasis(plus,q);

  ComponentList plus2;
  mod2(plus2,plus);

  ComponentMap i2;
  mod2(i2,q);

  ComponentMap id;
  identityMatrix(id,q.numRows());
  i2 += id;

  ComponentList b;
  i2.kernel(b);

  ComponentSubquotient cs(b,plus2,q.numRows());
  dpi0.swap(cs);

  return;
}

void minusBasis(latticetypes::WeightList& mb, 
		const latticetypes::LatticeMatrix& i)

/*!
  Synopsis: puts in mb a basis for the -1 eigenspace of the involution.

  Algorithm: the vectors e-i(e), when e runs through b, generate a lattice 
  commensurate with the eigenspace. Thus a smith basis for this lattice will 
  do the trick.
*/

{
  using namespace latticetypes;
  using namespace matrix;
  using namespace smithnormal;

  size_t n = i.numColumns();

  // put in q the matrix of vectors i(e)-e in the basis b

  LatticeMatrix q(i);
  for (size_t j = 0; j < n; ++j)
    q(j,j) -= 1;

  // find smith normal form

  CoeffList invf;
  WeightList bs;
  initBasis(bs,n);

  smithNormal(invf,bs.begin(),q);

  // copy significant part of basis in mb

  mb.reserve(invf.size());

  for (size_t j = 0; j < invf.size(); ++j)
    mb.push_back(bs[j]);

  return;
}

void minusMatrix(latticetypes::LatticeMatrix& qm, 
		 const latticetypes::LatticeMatrix& q, const RealTorus& t)

/*!
  Synopsis: writes in qm the matrix of the restriction of q to X_-.

  Precondition: q commutes with the involution;
*/

{
  const WeightList& bm = t.minusLattice();
  qm.resize(bm.size(),bm.size());

  for (size_t j = 0; j < bm.size(); ++j) {
    Weight v(t.rank());
    q.apply(v,bm[j]);
    Weight vm(bm.size());
    t.toMinus(vm,v);
    for (size_t i = 0; i < bm.size(); ++i)
      qm(i,j) = vm[i];
  }

  return;
}

void plusBasis(latticetypes::WeightList& pb,
	       const latticetypes::LatticeMatrix& i)

/*!
  Synopsis: puts in mb a basis for the +1 eigenspace of the involution;

  Algorithm: the vectors e+i(e), when e runs through b, generate a lattice 
  commensurate with the eigenspace.
*/

{
  using namespace latticetypes;
  using namespace matrix;
  using namespace smithnormal;

  size_t n = i.numColumns();

  // put in q the matrix of vectors i(e)+e in the basis b

  LatticeMatrix q(i);
  for (size_t j = 0; j < n; ++j)
    q(j,j) += 1;

  // find smith normal form

  CoeffList invf;
  WeightList bs;
  initBasis(bs,n);

  smithNormal(invf,bs.begin(),q);

  // copy significant part of basis in pb

  pb.reserve(invf.size());

  for (size_t j = 0; j < invf.size(); ++j)
    pb.push_back(bs[j]);

  return;
}

void plusMatrix(latticetypes::LatticeMatrix& qp, 
		const latticetypes::LatticeMatrix& q, const RealTorus& t)

/*!
  Synopsis: writes in qp the matrix of the restriction of q to X_+.

  Precondition: q commutes with the involution;
*/

{
  const WeightList& bp = t.plusLattice();
  qp.resize(bp.size(),bp.size());

  for (size_t j = 0; j < bp.size(); ++j) {
    Weight v(t.rank());
    q.apply(v,bp[j]);
    Weight vp(bp.size());
    t.toPlus(vp,v);
    for (size_t i = 0; i < bp.size(); ++i)
      qp(i,j) = vp[i];
  }

  return;
}

}

/*****************************************************************************

        Chapter III -- Auxiliary functions

  ... explain here when it is stable ...

******************************************************************************/

namespace {

void makeTopology(ComponentSubquotient& cs, const RealTorus& T)

/*!
  Synopsis: puts the subquotient V^tau/V_+ in cs.

  Explanation: V is X/2X; V_+ is the image of X_+ in V. So the subquotient
  is the _dual_ of the component group of the torus.
*/
{
  using namespace lattice;

  ComponentList plus2;
  mod2(plus2,T.plusLattice());

  ComponentMap i2;
  mod2(i2,T.involution());

  ComponentMap id;
  identityMatrix(id,T.rank());
  i2 += id;

  ComponentList b;
  i2.kernel(b);

  ComponentSubquotient cs1(b,plus2,T.rank());
  cs.swap(cs1);

  return;
}

void fullMinusBasis(latticetypes::WeightList& mb, 
		    latticetypes::LatticeMatrix& tm,
		    const latticetypes::LatticeMatrix& i)

/*!
  Synopsis: puts in mb a basis for the -1 eigenspace of the involution; puts
  in tm the matrix of a projection onto X_-.

  The matrix tm is in terms of the standard basis of X, and the chosen basis
  of X_-. There is no significance to the chosen complement; the intention is
  that it should be used only for vectors already in X_-.

  Algorithm: the vectors e-i(e), when e runs through b, generate a lattice 
  commensurate with the eigenspace. Thus a smith basis for this lattice will 
  do the trick. Such a basis also yields the projection matrix.
*/

{
  using namespace latticetypes;
  using namespace matrix;
  using namespace smithnormal;

  size_t n = i.numColumns();

  // put in q the matrix of vectors i(e)-e in the basis b

  LatticeMatrix q(i);
  for (size_t j = 0; j < n; ++j)
    q(j,j) -= 1;

  // find smith normal form

  CoeffList invf;
  WeightList bs;
  initBasis(bs,n);

  smithNormal(invf,bs.begin(),q);

  // copy significant part of basis in mb

  mb.reserve(invf.size());

  for (size_t j = 0; j < invf.size(); ++j)
    mb.push_back(bs[j]);

  // make projection matrix
  // rows of projection matrix are the first rows of the inverse of the
  // matrix of bs

  LatticeMatrix p(bs);
  LatticeCoeff d;
  p.invert(d);

  tm.resize(invf.size(),n);

  for (size_t i = 0; i < invf.size(); ++i)
    for (size_t j = 0; j < n; ++j)
      tm(i,j) = p(i,j);

  return;
}

void fullPlusBasis(latticetypes::WeightList& pb, 
		   latticetypes::LatticeMatrix& tp,
		   const latticetypes::LatticeMatrix& i)

/*!
  Synopsis: puts in mb a basis for the +1 eigenspace of the involution; puts
  in tm the matrix of a projection onto X_+.

  The matrix tm is in terms of the standard basis of X, and the chosen basis
  of X_+. There is no significance to the chosen complement; the intention is
  that it should be used only for vectors already in X_+.

  Algorithm: the vectors e+i(e), when e runs through b, generate a lattice 
  commensurate with the eigenspace. Thus a smith basis for this lattice will 
  do the trick. Such a basis also yields the projection matrix.
*/

{
  using namespace latticetypes;
  using namespace matrix;
  using namespace smithnormal;

  size_t n = i.numColumns();

  // put in q the matrix of vectors i(e)+e in the basis b

  LatticeMatrix q(i);
  for (size_t j = 0; j < n; ++j)
    q(j,j) += 1;

  // find smith normal form

  CoeffList invf;
  WeightList bs;
  initBasis(bs,n);

  smithNormal(invf,bs.begin(),q);

  // copy significant part of basis in pb

  pb.reserve(invf.size());

  for (size_t j = 0; j < invf.size(); ++j)
    pb.push_back(bs[j]);

  // make projection matrix
  // rows of projection matrix are the first rows of the inverse of the
  // matrix of bs

  LatticeMatrix p(bs);
  LatticeCoeff d;
  p.invert(d);

  tp.resize(invf.size(),n);

  for (size_t i = 0; i < invf.size(); ++i)
    for (size_t j = 0; j < n; ++j)
      tp(i,j) = p(i,j);

  return;
}

}

}
