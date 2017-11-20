/*
  This is lietype.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  Definitions for the type of a real or complex reductive Lie algebra.

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

#include <cassert>

#include "lietype.h"

#include "constants.h"
#include "matreduc.h"
#include "matrix.h"

#include "../Atlas.h"


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

  void addCompactInvolution(WeightInvolution&, unsigned int, unsigned int,
			    const Permutation& pi);

  void addDInvolution(WeightInvolution&, unsigned int, unsigned int,
		      const Permutation& pi);

  void addMinusIdentity(WeightInvolution&, unsigned int, unsigned int,
			const Permutation& pi);

  void addSimpleInvolution(WeightInvolution&, unsigned int,
			   const SimpleLieType&, TypeLetter,
			   const Permutation& pi);
}

/*****************************************************************************

        Chapter I -- Functions declared in lietype.h

******************************************************************************/

namespace lietype {

// Cartan matrix info, simple type |tp|, distance $d\leq2$ off diagonal
inline int dispatch
  (TypeLetter tp, unsigned int r,unsigned int min,unsigned int d, bool lower)
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

int SimpleLieType::Cartan_entry(unsigned int i,unsigned int j) const
{
  if (type()=='T') return 0;

  unsigned int min,d;
  if (i<j) min=i,d=j-i;
  else     min=j,d=i-j;

  return d>2 ? 0 : dispatch(type(),rank(),min,d,i<j);
}

// implicitly define the Cartan matrix corresponding to the type
int LieType::Cartan_entry(unsigned int i,unsigned int j) const
{ unsigned int min,d;
  if (i<j)
    min=i,d=j-i;
  else
    min=j,d=i-j;

  if (d>2) return 0;

  for (base::const_iterator it=begin(); it!=end(); ++it)
  {
    auto r=it->rank();
    if (min>=r)
      min-=r; // the only case that continues the loop
    else // |min<r|, so least index belongs to current simple factor
    { TypeLetter tp=it->type();
      if (min+d>=r or tp=='T') // distinct simple factors or torus
	return 0;
      else return dispatch(tp,r,min,d,i<j);
    }
  }
  assert(false); // indices out of bounds
  return 0;
}

int_Matrix SimpleLieType::Cartan_matrix() const
{ const auto r=rank();
  int_Matrix result(r,r);
  for (unsigned int i=0; i<r; ++i)
    for (unsigned int j=0; j<r; ++j)
      result(i,j)=Cartan_entry(i,j);

  return result;
}

int_Matrix SimpleLieType::transpose_Cartan_matrix() const
{ const auto r=rank();
  int_Matrix result(r,r);
  for (unsigned int i=0; i<r; ++i)
    for (unsigned int j=0; j<r; ++j)
      result(i,j)=Cartan_entry(j,i);

  return result;
}

int_Matrix LieType::Cartan_matrix() const
{ const auto r=rank();
  int_Matrix result(r,r);
  for (unsigned int i=0; i<r; ++i)
    for (unsigned int j=0; j<r; ++j)
      result(i,j)=Cartan_entry(i,j);

  return result;
}

int_Matrix LieType::transpose_Cartan_matrix() const
{ const auto r=rank();
  int_Matrix result(r,r);
  for (unsigned int i=0; i<r; ++i)
    for (unsigned int j=0; j<r; ++j)
      result(i,j)=Cartan_entry(j,i);

  return result;
}

unsigned int LieType::rank() const
{
  unsigned int r = 0;
  for (base::const_iterator it=begin(); it!=end(); ++it)
    r += it->rank();
  return r;
}


unsigned int LieType::semisimple_rank() const
{
  unsigned int r = 0;
  for (base::const_iterator it=begin(); it!=end(); ++it)
    r += it->semisimple_rank();
  return r;
}

/*
  This function constructs a "blockwise adapted" basis for the root
  lattice inside the weight lattice (i.e., it does just that for each
  simple factor, and returns the canonical basis for the torus factors.)

  The purpose of doing this blockwise instead of finding a global Smuth normal
  basis is to permit a better reading of the quotient group: this will be
  presented as a sequence of factors, corresponding to each simple block.

  In fact, |matreduc::adapted_basis| applied to a block matrix like the full
  Cartan matrix would return a similar block result, since it does not try to
  get to the true Smith normal form (which might involve combining invariant
  factors). It is clearer and more efficient though to do this blockwise.
*/
int_Matrix LieType::Smith_basis(CoeffList& invf) const
{

  const auto R=rank();
  int_Matrix result(R,R,0);
  invf.reserve(R); // every element of |result| has its invariant factor

  // get adapted basis for each simple factor
  unsigned int s=0; //offset
  for (const_iterator it=begin(); it!=end(); ++it)
  {
    const auto r =it->rank();

    if (it->type() == 'T') // torus type T_r
    {
      for (unsigned int i=s; i<s+r; ++i)
	result(i,i)=1;
      invf.insert(invf.end(),r,0); // add |r| factors 0
    }
    else
    {
      int_Matrix tC=it->transpose_Cartan_matrix();
      CoeffList new_invf;
      int_Matrix Sb = matreduc::adapted_basis(tC,new_invf);

      //make a small adjustment for types $D_{2n}$
      if (it->type() == 'D' and it->rank()%2 == 0)
      {
	assert(new_invf[r-2]==2 and new_invf[r-1]==2);
	Sb.columnOperation(r-2,r-1,1); // this makes Sb[r-2] a unit vector
      }

      // copy matrix |Sb| into block of result
      for (unsigned int j=0; j<r; ++j)
	for (unsigned int i=0; i<r; ++i)
	  result(s+i,s+j)=Sb(i,j);

      invf.insert(invf.end(),new_invf.begin(),new_invf.end()); // append |invf|
    }
    s += r;
  } // |for it|

  return result;
}

/*
  Return the dual Lie type of lt. In fact this applies to ordered Dynkin
  diagrams, whence the distinction between (B2,C2), (F4,f4) and (G2,g2)
*/
LieType dual_type(LieType lt)
{
  for (unsigned int i=0; i<lt.size(); ++i)
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


/*
  Return dual inner class type of |ict| with respect to |lt|
  The result is independent of whether |lt| is the original or dual type
*/
InnerClassType dual_type(InnerClassType ict, const LieType& lt)
{

  unsigned int ltj = 0;

  for (unsigned int i=0; i<ict.size(); ++i)
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
  amounts to interchanging inner class types 'c' and 's'. For complex ('C')
  inner classes we currently do nothing, although if the Lie type is XX where
  -1 is not in the Weyl group of X (i.e., such that 's' is not 'c' for X),
  then the negated transpose distinguished involution is not in the "same"
  inner class, but rather in "another" complex class of Lie type XX. The very
  limited use made (in output) of the dual Layout justifies this imprecision.
*/
Layout dual(const Layout& lo)
{
  Layout result=lo;
  for (unsigned int i=0,k=0; i<lo.d_type.size(); k+=lo.d_type[i].rank(),++i)
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

  unsigned int i=0; // index into |lo.d_type|
  for (unsigned int j=0; j<lo.d_inner.size(); ++i,++j)
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
      case 'T': swap_sc = true; break;
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
bool checkRank(const TypeLetter& x, unsigned int l)
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

WeightInvolution simple_involution(const SimpleLieType& slt, simple_ict tp)
{
  unsigned r=slt.rank();
  WeightInvolution result
    (tp==simple_ict::complex ? 2*r : r); // start with identity
  if (tp==simple_ict::complex)
    for (unsigned i=0; i<r; ++i)
    {
      std::swap(result(i,i),result(i,i+r));
      std::swap(result(i+r,i),result(i+r,i+r));
    }
  else if (tp==simple_ict::unequal_rank)
    switch (slt.type())
    {
    case 'A': // antidiagonal matrix
      assert(r>1);
      for (unsigned i=0; 2*i<r-1; ++i)
      {
	std::swap(result(i,i),result(i,r-1-i));
	std::swap(result(r-1-i,i),result(r-1-i,r-1-i));
      }
      break;
    case 'B': case 'C': case 'F': case 'G': default:
      assert(false); break;
    case 'D':
      std::swap(result(r-2,r-2),result(r-2,r-1));
      std::swap(result(r-1,r-2),result(r-1,r-1));
      break;
    case 'E':
      assert(r==6);
      std::swap(result(0,0),result(0,5));
      std::swap(result(5,0),result(5,5));
      std::swap(result(2,2),result(2,4));
      std::swap(result(4,2),result(4,4));
      break;
    case 'T':
      assert(r==1);
      result(0,0)=-1;
    }
  return result;
}

/*
  Construct the fundamental involution for the Lie type lt and the
  inner class ic, in the weight basis for the simply connected group.

  Precondition: validity if |lo.d_inner| has been checked
*/
WeightInvolution involution(const Layout& lo)
{
  const LieType& lt = lo.d_type;
  const InnerClassType& ic = lo.d_inner;
  const Permutation& pi = lo.d_perm;

  WeightInvolution result(lt.rank(),lt.rank(),0);

  unsigned int r = 0;   // position in flattened Dynkin diagram; an index into |pi|
  unsigned int pos = 0; // position in |lt|

  for (unsigned int j=0; j<ic.size(); ++j) // |r|,|pos| are also advanced, near end
  {
    SimpleLieType slt = lt[pos];
    const auto rs = slt.rank();

    switch (ic[j])
    {
    case 'c': // add the identity
      addCompactInvolution(result,r,rs,pi);
      break;
    case 's': // add split involution
      switch (slt.type())
      {
      case 'A': // antidiagonal matrix
	for (unsigned int i=0; i<rs; ++i)
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
      for (unsigned int i=0; i<rs; ++i)
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

WeightInvolution involution(const LieType& lt,
			    const InnerClassType& ict)
{ return involution(Layout(lt,ict)); }

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
void addCompactInvolution(WeightInvolution& m, unsigned int r,
			  unsigned int rs,
			  const Permutation& pi)
{
  for (unsigned int i=0; i<rs; ++i)
    m(pi[r+i],pi[r+i]) = 1;
}


/*!
  Synopsis: flips the last two vectors in the block of size rs starting
  from (r,r).

  Precondition: the block is set to zero.
*/
void addDInvolution(WeightInvolution& m, unsigned int r, unsigned int rs,
		    const Permutation& pi)
{
  for (unsigned int i=0; i<rs-2; ++i)
    m(pi[r+i],pi[r+i]) = 1;

  m(pi[r+rs-2],pi[r+rs-1]) = 1;
  m(pi[r+rs-1],pi[r+rs-2]) = 1;
}


/*!
  Synopsis: sets the block of size (rs,rs) starting from (r,r) to minus the
  identity

  Precondition: the block is set to zero.
*/
void addMinusIdentity(WeightInvolution& m, unsigned int r, unsigned int rs,
		      const Permutation& pi)
{
  for (unsigned int i=0; i<rs; ++i)
    m(pi[r+i],pi[r+i]) = -1;
}


/*!
  Synopsis: appends to m, from position (r,r), the fundamental involution
  corresponding to x in size rs.
*/
void addSimpleInvolution(WeightInvolution& m, unsigned int r,
			 const SimpleLieType& slt, TypeLetter x,
			 const Permutation& pi)
{
  const auto rs = slt.rank();

  switch (x) {
  case 'c': // add the identity
    addCompactInvolution(m,r,rs,pi);
    break;
  case 's': // add split involution
    switch (slt.type()) {
    case 'A': // antidiagonal matrix
      for (unsigned int i=0; i<rs; ++i)
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
    for (unsigned int i=0; i<rs; ++i)
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
