/*!
\file
\brief Definitions for the type of a real or complex reductive Lie algebra.

  A complex reductive Lie algebra is simply the product of a number of simple
  complex Lie algebras, and a torus; the simple factors can be of types A-G,
  with ranks >=1 in type A, >=2 in type B, >=3 in type C, >=4 in type D, 6,7 or
  8 in type E, 4 in type F, and 2 in type G. A torus component will be
  designated by the letter T. So we express the type as a vector, each
  component of which is a pair (letter,rank).

  As is well-known, a real reductive Lie algebra is a product of factors which
  are (a) real forms of simple complex Lie algebras (b) simple complex Lie
  algebras where we forget the complex structure (these are real forms of the
  product of two isomorphic complex simple Lie algebras) (c) real forms of
  one-dimensional tori and (d) one-dimensional complex tori seen as
  two-dimensional real tori. In fact for the torus factors we will also allow
  higher-dimensional factors; and of course at the Lie algebra level case
  (d) above might be reduced to the sum of two instances of type (c), but
  the distinction is necessary at the group level, so we keep it here as well.

  The classification of real forms of simple complex Lie algebras is
  well-known;  we will follow the notation from Helgason, Differential
  Geometry, Lie Groups, and Symmetric Spaces, Academic Press, New York, 1978,
  Table VI in Chapter X, page 532, where the Satake diagrams of all non-compact
  and non-complex simple real reductive Lie algebras appear. Accordingly, we
  represent a real form of our complex Lie algebra by a string of symbols
  which may be of the form Split, Compact, Complex or I-IX, depending on the
  type; a symbol of type Complex really consumes two consecutive isomorphic
  complex factors. This representation is only for input/output purposes;
  internally, we work only with the Cartan matrix and the Cartan involution.

*/
/*
  This is lietype.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <cassert>

#include "lietype.h"

#include "constants.h"
#include "latticetypes.h"
#include "layout.h"
#include "smithnormal.h"

/*****************************************************************************

  This file contains the definitions for the type of a real or complex
  reductive Lie algebra.

  [Original comment by Fokko, more about the design than the actual code.]
  A complex reductive Lie algebra is simply the product of a number of simple
  complex Lie algebras, and a torus; the simple factors can be of types A-G,
  with ranks >=1 in type A, >=2 in type B, >=3 in type C, >=4 in type D, 6,7 or
  8 in type E, 4 in type F, and 2 in type G. A torus component will be
  designated by the letter T. So we express the type as a vector, each
  component of which is a pair (letter,rank).

  As is well-known, a real reductive Lie algebra is a product of factors which
  are (a) real forms of simple complex Lie algebras (b) simple complex Lie
  algebras where we forget the complex structure (these are real forms of the
  product of two isomorphic complex simple Lie algebras) (c) real forms of
  one-dimensional tori and (d) one-dimensional complex tori seen as
  two-dimensional real tori. In fact for the torus factors we will also allow
  higher-dimensional factors; and of course at the Lie algebra level case
  (d) above might be reduced to the sum of two instances of type (c), but
  the distinction is necessary at the group level, so we keep it here as well.

  The classification of real forms of simple complex Lie algebras is
  well-known; we will follow the notation from Helgason, Differential
  Geometry, Lie Groups, and Symmetric Spaces, Academic Press, New York, 1978,
  Table VI in Chapter X, page 532, where the Satake diagrams of all
  non-compact and non-complex simple real reductive Lie algebras appear.
  Accordingly, we represent a real form of our complex Lie algebra by a string
  of symbols which may be of the form Split, Compact, Complex or I-IX,
  depending on the type. [The types I-IX are not used; they do not seem to
  have survived more that the initial design of this module. MvL] A symbol of
  type Complex really consumes two consecutive isomorphic complex factors.
  This representation is only for input/output purposes; internally, we work
  only with the Cartan matrix and the Cartan involution.

******************************************************************************/

/* In fact almost all that happens below is computing involution matrices */

namespace atlas {

namespace lietype {

  void addCompactInvolution(latticetypes::LatticeMatrix&, size_t, size_t,
			    const setutils::Permutation& pi);

  void addDInvolution(latticetypes::LatticeMatrix&, size_t, size_t,
		      const setutils::Permutation& pi);

  void addMinusIdentity(latticetypes::LatticeMatrix&, size_t, size_t,
			const setutils::Permutation& pi);

  void addSimpleInvolution(latticetypes::LatticeMatrix&, size_t,
			   const SimpleLieType&, TypeLetter,
			   const setutils::Permutation& pi);
}

/*****************************************************************************

        Chapter I -- Functions declared in lietype.h

******************************************************************************/

namespace lietype {

// Cartan matrix info, simple type |tp|, distance $d\leq2$ off diagonal
inline int dispatch(TypeLetter tp, size_t r,size_t min,size_t d, bool lower)
{
  if (d==0) return 2;
  if (tp=='D') return (min<r-3 ? d==1 : min==r-3) ? -1 : 0;
  if (tp=='E') return d==(min<2 ? 2 : 1) ? -1 : 0;
  // now diagram is linear
  if (d==2) return 0;
  if (tp<'D') return tp=='A' or min<r-2 or lower==(tp=='C') ? -1 : -2 ;
  if (tp=='F' or tp=='f') return min!=1 or lower==(tp=='f') ? -1 : -2 ;
  return lower==(tp=='G') ? -1 : -3;
}

int SimpleLieType::Cartan_entry(size_t i,size_t j) const
{
  if (type()=='T') return 0;

  size_t min,d;
  if (i<j) min=i,d=j-i;
  else     min=j,d=i-j;

  return d>2 ? 0 : dispatch(type(),rank(),min,d,i<j);
}

// implicitly define the Cartan matrix corresponding to the type
int LieType::Cartan_entry(size_t i,size_t j) const
{ size_t min,d;
  if (i<j)
    min=i,d=j-i;
  else
    min=j,d=i-j;

  if (d>2) return 0;

  for (base::const_iterator it=begin(); it!=end(); ++it)
  {
    size_t r=it->rank();
    if (r<=min)
      min-=r; // the only case that continues the loop
    else
    { TypeLetter tp=it->type();
      if (r<=min+d or tp=='T') // distinct simple factors or torus
	return 0;
      else return dispatch(tp,r,min,d,i<j);
    }
  }
  assert(false); // indices out of bounds
  return 0;
}

latticetypes::LatticeMatrix SimpleLieType::Cartan_matrix() const
{ size_t r=rank();
  latticetypes::LatticeMatrix result(r,r);
  for (size_t i=0; i<r; ++i)
    for (size_t j=0; j<r; ++j)
      result(i,j)=Cartan_entry(i,j);

  return result;
}

latticetypes::LatticeMatrix SimpleLieType::transpose_Cartan_matrix() const
{ size_t r=rank();
  latticetypes::LatticeMatrix result(r,r);
  for (size_t i=0; i<r; ++i)
    for (size_t j=0; j<r; ++j)
      result(i,j)=Cartan_entry(j,i);

  return result;
}

latticetypes::LatticeMatrix LieType::Cartan_matrix() const
{ size_t r=rank();
  latticetypes::LatticeMatrix result(r,r);
  for (size_t i=0; i<r; ++i)
    for (size_t j=0; j<r; ++j)
      result(i,j)=Cartan_entry(i,j);

  return result;
}

latticetypes::LatticeMatrix LieType::transpose_Cartan_matrix() const
{ size_t r=rank();
  latticetypes::LatticeMatrix result(r,r);
  for (size_t i=0; i<r; ++i)
    for (size_t j=0; j<r; ++j)
      result(i,j)=Cartan_entry(j,i);

  return result;
}

size_t LieType::rank() const
{
  size_t r = 0;
  for (base::const_iterator it=begin(); it!=end(); ++it)
    r += it->rank();
  return r;
}


size_t LieType::semisimple_rank() const
{
  size_t r = 0;
  for (base::const_iterator it=begin(); it!=end(); ++it)
    r += it->semisimple_rank();
  return r;
}

/*
  This function constructs a "blockwise Smith normal" basis for the root
  lattice inside the weight lattice (i.e., it does just that for each
  simple factor, and returns the canonical basis for the torus factors.)

  The purpose of doing this blockwise instead of globally is to permit a
  better reading of the quotient group: this will be presented as a sequence
  of factors, corresponding to each simple block.
*/
latticetypes::WeightList
LieType::Smith_basis(latticetypes::CoeffList& invf) const
{

  // Smith-normalize for each simple factor
  latticetypes::WeightList result;
  matrix::initBasis(result,rank());
  latticetypes::WeightList::iterator rp = result.begin();

  for (const_iterator it=begin(); it!=end(); ++it)
  {
    size_t r =it->rank();

    if (it->type() == 'T') // torus type T_r
      invf.insert(invf.end(),r,0); // add |r| factors 0
    else
    {
      latticetypes::LatticeMatrix tC=it->transpose_Cartan_matrix();
      smithnormal::smithNormal(invf,rp,tC); // adapt next |r| vectors to lattice

      //make a small adjustment for types $D_{2n}$
      if (it->type() == 'D' and it->rank()%2 == 0)
	rp[r-2] += rp[r-1];

    }
    rp += r;

  }
  return result;
}

/*! brief Returns the dual Lie type of lt. In fact this applies to ordered
  Dynkin diagrams, whence the distinction between (B2,C2), (F4,f4) and (G2,g2)
*/
LieType dual_type(LieType lt)
{
  for (size_t i=0; i<lt.size(); ++i)
    switch (lt[i].first) {
    case 'B':
      lt[i].first = 'C';
      break;
    case 'C':
      lt[i].first = 'B';
      break;
    case 'F':
      lt[i].first = 'f';
      break;
    case 'f':
      lt[i].first = 'F';
      break;
    case 'G':
      lt[i].first = 'g';
      break;
    case 'g':
      lt[i].first = 'G';
      break;
    }

  return lt;
}


/*!
  \brief Returns dual inner class type of ict.
*/
InnerClassType dual_type(InnerClassType ict, const LieType& lt)
{

  size_t ltj = 0;

  for (size_t i=0; i<ict.size(); ++i)
  {
    if (ict[i] == 'C') // dual type remains complex
    {
      ltj += 2;
      continue;
    }
    SimpleLieType slt = lt[ltj];
    switch(slt.first)
    {
    case 'B':
    case 'C':
    case 'F':
    case 'f':
    case 'G':
    case 'g':
      break;
    case 'A':
    case 'E':
    case 'T':
      // Interchange split and compact inner classes
      if (ict[i] == 's')
	ict[i] = 'c';
      else
	ict[i] = 's';
      break;
    case 'D':
      if (slt.second%2 !=0) {
	// interchange split and compact inner classes for D_{2n+1}
	if (ict[i] == 's' or ict[i] == 'u')
	  ict[i] = 'c';
	else
	  ict[i] = 's';
      }
      break;
    }
    ++ltj;
  }

  return ict;
}


/*!
  Synopsis: checks if the rank l is in the valid range for x.
*/
bool checkRank(const TypeLetter& x, size_t l)
{
  switch (x)
  {
  case 'A': // rank must be >= 1
    if ((l < 1) or (l > constants::RANK_MAX))
      return false;
    break;
  case 'B': // rank must be >= 2
    if ((l < 2) or (l > constants::RANK_MAX))
      return false;
    break;
  case 'C': // rank must be >= 2
    if ((l < 2) or (l > constants::RANK_MAX))
      return false;
    break;
  case 'D': // rank must be >= 4
    if ((l < 4) or (l > constants::RANK_MAX))
      return false;
    break;
  case 'E': // rank must be 6, 7 or 8
    if ((l < 6) or (l > 8))
      return false;
    break;
  case 'F':
  case 'f': // rank must be 4
    if (l != 4)
      return false;
    break;
  case 'G':
  case 'g': // rank must be 2
    if (l != 2)
      return false;
    break;
  case 'T': // rank must be >= 1
    if ((l < 1) or (l > constants::RANK_MAX))
      return false;
    break;
  default: // this cannot happen!
    assert(false && "unexpected type in checkRank");
    break;
  }

  return true;
}


/*!
  Synopsis: constructs the fundamental involution for the Lie type lt and the
  inner class ic, in the weight basis for the simply connected group.

  Precondition: it has already been checked that ic holds a valid inner class
  type for lt.
*/
latticetypes::LatticeMatrix involution(const lietype::LieType& lt,
				       const lietype::InnerClassType& ic)
{ // make |Layout| structure with empty basis and trivial permutation
  layout::Layout lo(lt,ic);
  return involution(lo); // and call the more general procedure
}

/* The layout structure provides arguments; |d_basis| ignored if empty */
latticetypes::LatticeMatrix involution(const layout::Layout& lo)
  throw (std::runtime_error,std::bad_alloc)
{
  const lietype::LieType& lt = lo.d_type;
  const lietype::InnerClassType& ic = lo.d_inner;
  const setutils::Permutation& pi = lo.d_perm;

  latticetypes::LatticeMatrix result(rank(lt),rank(lt),0);

  size_t r = 0;   // position in flattened Dynkin diagram; index to |pi|
  size_t pos = 0; // position in |lt|

  for (size_t j=0; j<ic.size(); ++j) // |r|,|pos| are also advanced, near end
  {
    SimpleLieType slt = lt[pos];
    size_t rs = rank(slt);

    switch (ic[j])
    {
    case 'c': // add the identity
      addCompactInvolution(result,r,rs,pi);
      break;
    case 's': // add split involution
      switch (type(slt))
      {
      case 'A': // antidiagonal matrix
	for (size_t i=0; i<rs; ++i)
	  result(pi[r+i],pi[r+rs-1-i]) = 1;
	break;
      case 'D':
	if (rank(slt)%2 != 0)
	  addDInvolution(result,r,rs,pi);
	else
	  addCompactInvolution(result,r,rs,pi);
	break;
      case 'E':
	if (rs == 6)
	{
	  result(pi[r+1],pi[r+1]) = 1;
	  result(pi[r+3],pi[r+3]) = 1;
	  result(pi[r],pi[r+5]) = 1;
	  result(pi[r+5],pi[r]) = 1;
	  result(pi[r+2],pi[r+4]) = 1;
	  result(pi[r+4],pi[r+2]) = 1;
	}
	else
	  addCompactInvolution(result,r,rs,pi);
	break;
      case 'T':
	addMinusIdentity(result,r,rs,pi);
	break;
      default: // identity involution for types B,C,E7,E8,F,f,G,g
	addCompactInvolution(result,r,rs,pi);
	break;
      } // |switch (type(slt))|
      break;
    case 'C': // Compact: parallel interchange of |rs| vertices with next |rs|
      for (size_t i=0; i<rs; ++i)
      {
	result(pi[r+i],pi[r+rs+i]) = 1;
	result(pi[r+rs+i],pi[r+i]) = 1;
      }
      ++pos; r += rs; // account for consumption of extra simple factor
      break;
    case 'u': // flip the last two vectors
      addDInvolution(result,r,rs,pi);
      break;
    default: // this should not happen!
      assert(false and "wrong inner class letter");
      break;
    }
    ++pos; r += rs; // consume simple factor
  } // |for (j)|

  if (lo.d_basis.size()!=0)
    invConjugate(result,latticetypes::LatticeMatrix(lo.d_basis));

  return result;
}

} // |namespace lietype|

/*****************************************************************************

        Chapter II -- Private functions

******************************************************************************/

namespace lietype {


/*!
  Synopsis: sets the block of size (rs,rs) starting from (r,r) to the
  identity

  Precondition: the block is set to zero.
*/
void addCompactInvolution(latticetypes::LatticeMatrix& m, size_t r,
			  size_t rs,
			  const setutils::Permutation& pi)
{
  for (size_t i=0; i<rs; ++i)
    m(pi[r+i],pi[r+i]) = 1;
}


/*!
  Synopsis: flips the last two vectors in the block of size rs starting
  from (r,r).

  Precondition: the block is set to zero.
*/
void addDInvolution(latticetypes::LatticeMatrix& m, size_t r, size_t rs,
		    const setutils::Permutation& pi)
{
  for (size_t i=0; i<rs-2; ++i)
    m(pi[r+i],pi[r+i]) = 1;

  m(pi[r+rs-2],pi[r+rs-1]) = 1;
  m(pi[r+rs-1],pi[r+rs-2]) = 1;
}


/*!
  Synopsis: sets the block of size (rs,rs) starting from (r,r) to minus the
  identity

  Precondition: the block is set to zero.
*/
void addMinusIdentity(latticetypes::LatticeMatrix& m, size_t r, size_t rs,
		      const setutils::Permutation& pi)
{
  for (size_t i=0; i<rs; ++i)
    m(pi[r+i],pi[r+i]) = -1;
}


/*!
  Synopsis: appends to m, from position (r,r), the fundamental involution
  corresponding to x in size rs.
*/
void addSimpleInvolution(latticetypes::LatticeMatrix& m, size_t r,
			 const SimpleLieType& slt, TypeLetter x,
			 const setutils::Permutation& pi)
{
  size_t rs = rank(slt);

  switch (x) {
  case 'c': // add the identity
    addCompactInvolution(m,r,rs,pi);
    break;
  case 's': // add split involution
    switch (type(slt)) {
    case 'A': // antidiagonal matrix
      for (size_t i=0; i<rs; ++i)
	m(pi[r+i],pi[r+rs-1-i]) = 1;
      break;
    case 'D':
      if (rank(slt)%2 != 0)
	addDInvolution(m,r,rs,pi);
      else
	addCompactInvolution(m,r,rs,pi);
      break;
    case 'E':
      if (rank(slt) == 6) {
	m(pi[r+1],pi[r+1]) = 1;
	m(pi[r+3],pi[r+3]) = 1;
	m(pi[r],pi[r+5]) = 1;
	m(pi[r+5],pi[r]) = 1;
	m(pi[r+2],pi[r+4]) = 1;
	m(pi[r+4],pi[r+2]) = 1;
      }
      else
	addCompactInvolution(m,r,rs,pi);
      break;
    case 'T':
      addMinusIdentity(m,r,rs,pi);
      break;
    default: // identity involution for types B,C,E7,E8,F,f,G,g
      addCompactInvolution(m,r,rs,pi);
      break;
    }
    break;
  case 'C': // rs-dimensional flip
    for (size_t i=0; i<rs; ++i)
    {
      m(pi[r+i],pi[r+rs+i]) = 1;
      m(pi[r+rs+i],pi[r+i]) = 1;
    }
    break;
  case 'u': // flip the last two vectors
    addDInvolution(m,r,rs,pi);
    break;
  default: // this should not happen!
    assert(false && "wrong inner class letter in addSimpleInvolution");
    break;
  }
}

} // |namespace lietype|

} // |namespace atlas|
