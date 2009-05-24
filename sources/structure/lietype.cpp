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
  \brief Returns dual inner class type of |ict| with respect to |lt|

  The result is independent of whether |lt| is the original or dual type
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
    switch(slt.type())
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
      if (slt.rank()%2 !=0)// for D_{2n+1}
      { //  interchange split and compact inner classes
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

/* compute dual of Lie type (including diagram numbering) and inner class.
   For the dual of types $F_4$ end $G_2$ we reverse the numbering of the
   nodes, thus avoiding the introduction of the variant types $f_4$ and $g_2$.
   We just interchange types $B_n$ and $C_n$, even for $n=2$ where we could
   have instead interchanged the numbering of the two roots as for $G_2$. This
   is because users like the (fictive) distinction between $B_2$ and $C_2$.

   The dual inner class contains (as non-based root datum involution) the
   negated transpose of the distinguished involution of the inner class. This
   amounts to interchanging inner class types 'c' and 's'; for complex ('C')
   inner classes we currently do nothing, although the Lie type is XX where -1
   is not in the Weyl group of X (i.e., such that 's' is not 'c' for X), the
   negated transpose distinguished involution is not in the "same" inner
   class, but rather in "another" complex class of Lie type XX; the limited
   use currently made of the dual Layout (in realform_io) justifies this.
 */
Layout dual(const Layout& lo)
{
  Layout result=lo;
  for (size_t i=0,k=0; i<lo.d_type.size(); k+=lo.d_type[i].rank(),++i)
    switch(lo.d_type[i].type())
    {
    case 'B': result.d_type[i].type()='C'; break;
    case 'C': result.d_type[i].type()='B'; break;
    case 'F': // reverse order of simple factor in permutation
      std::swap(result.d_perm[k],result.d_perm[k+3]);
      std::swap(result.d_perm[k+1],result.d_perm[k+2]); break;
    case 'G': // reverse order of simple factor in permutation
      std::swap(result.d_perm[k],result.d_perm[k+1]); break;
    default: break;
    }

  size_t i=0; // index into |lo.d_type|
  for (size_t j=0; j<lo.d_inner.size(); ++i,++j)
  {

    if (lo.d_inner[j]=='C') // dual type remains complex
      ++i; // skip additional simple factor, 'C' remains unchanged
    else
    { bool swap_sc;
      switch(lo.d_type[i].type())
      {
      case 'A': swap_sc = lo.d_type[i].rank()>1; break;
      case 'D': swap_sc = lo.d_type[i].rank()%2!=0; break;
      case 'E': swap_sc = lo.d_type[i].rank()==6; break;
      case 'T': swap_sc = true;
      default: swap_sc=false; break;
      }
      if (swap_sc) // then interchange 'c' and 's'; note that 'u'->'s' here
      {
	if (lo.d_inner[j]=='c')
	  result.d_inner[j]='s';
	else if (lo.d_inner[j]=='s')
	  result.d_inner[j]='c';
      }
    } // |if(lo.d_inner[j]...)|
  } // |for(j)|

  return result;
}


/*!
  Synopsis: checks if the rank l is in the valid range for x.
*/
bool checkRank(const TypeLetter& x, size_t l)
{
  if (l>constants::RANK_MAX) return false;
  switch (x)
  {
  case 'A': return l>=1;
  case 'B': return l>=2;
  case 'C': return l>=2;
  case 'D': return l>=4;
  case 'E': return l>=6 and l<=8;
  case 'F':
  case 'f': return l==4;
  case 'G':
  case 'g': return l==2;
  case 'T': return l>=1;
  default: // this cannot happen!
    assert(false && "unexpected type in checkRank");
    return false;
  }
}


/*!
  Synopsis: constructs the fundamental involution for the Lie type lt and the
  inner class ic, in the weight basis for the simply connected group.

  Precondition: validity if |lo.d_inner| has been checked
*/
latticetypes::LatticeMatrix involution(const Layout& lo)
  throw (std::runtime_error,std::bad_alloc)
{
  const lietype::LieType& lt = lo.d_type;
  const lietype::InnerClassType& ic = lo.d_inner;
  const setutils::Permutation& pi = lo.d_perm;

  latticetypes::LatticeMatrix result(lt.rank(),lt.rank(),0);

  size_t r = 0;   // position in flattened Dynkin diagram; index to |pi|
  size_t pos = 0; // position in |lt|

  for (size_t j=0; j<ic.size(); ++j) // |r|,|pos| are also advanced, near end
  {
    SimpleLieType slt = lt[pos];
    size_t rs = slt.rank();

    switch (ic[j])
    {
    case 'c': // add the identity
      addCompactInvolution(result,r,rs,pi);
      break;
    case 's': // add split involution
      switch (slt.type())
      {
      case 'A': // antidiagonal matrix
	for (size_t i=0; i<rs; ++i)
	  result(pi[r+i],pi[r+rs-1-i]) = 1;
	break;
      case 'D':
	if (slt.rank()%2 != 0)
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
  size_t rs = slt.rank();

  switch (x) {
  case 'c': // add the identity
    addCompactInvolution(m,r,rs,pi);
    break;
  case 's': // add split involution
    switch (slt.type()) {
    case 'A': // antidiagonal matrix
      for (size_t i=0; i<rs; ++i)
	m(pi[r+i],pi[r+rs-1-i]) = 1;
      break;
    case 'D':
      if (slt.rank()%2 != 0)
	addDInvolution(m,r,rs,pi);
      else
	addCompactInvolution(m,r,rs,pi);
      break;
    case 'E':
      if (slt.rank() == 6) {
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
