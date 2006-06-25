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
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "lietype.h"

#include <cassert>

#include "constants.h"
#include "latticetypes.h"

/*****************************************************************************

  This file contains the definitions for the type of a real or complex
  reductive Lie algebra.

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

******************************************************************************/

namespace atlas {

namespace lietype {

  void addCompactInvolution(latticetypes::LatticeMatrix&, size_t, size_t);

  void addDInvolution(latticetypes::LatticeMatrix&, size_t, size_t);

  void addMinusIdentity(latticetypes::LatticeMatrix&, size_t, size_t);

  void addSimpleInvolution(latticetypes::LatticeMatrix&, size_t, 
			   const SimpleLieType&, TypeLetter);
}

/*****************************************************************************

        Chapter I -- Functions declared in lietype.h

  ... explain here when stable ...

******************************************************************************/

namespace lietype {

void dualLieType(LieType& dlt, const LieType& lt)

/*!
  Synopsis: puts in dlt the dual Lie type of lt.
*/

{
  dlt = lt;

  for (size_t j = 0; j < dlt.size(); ++j)
    switch (dlt[j].first) {
    case 'B':
      dlt[j].first = 'C';
      break;
    case 'C':
      dlt[j].first = 'B';
      break;
    case 'F':
      dlt[j].first = 'f';
      break;
    case 'f':
      dlt[j].first = 'F';
      break;
    case 'G':
      dlt[j].first = 'g';
      break;
    case 'g':
      dlt[j].first = 'G';
      break;
    }

  return;
}

void dualInnerClassType(InnerClassType& dict, const InnerClassType& ict,
			const LieType& lt)

/*!
  Synopsis: puts in dict the dual inner class type of ict.
*/

{
  dict = ict;

  size_t ltj = 0;

  for (size_t j = 0; j < dict.size(); ++j) {
    if (dict[j] == 'C') { // dual type is complex
      ltj += 2;
      continue;
    }
    SimpleLieType slt = lt[ltj];
    switch(slt.first) {
    case 'B':
    case 'C':
    case 'F':
    case 'f':
    case 'G':
    case 'g':
      break;
    case 'A':
    case 'E': 
      // interchange split and compact inner classes
      if (dict[j] == 's')
	dict[j] = 'c';
      else      
	dict[j] = 's';
      break;
    case 'D':
      if (slt.second & 1ul) {
	// interchange split and compact inner classes
	if (dict[j] == 's' or dict[j] == 'u')
	  dict[j] = 'c';
	else
	  dict[j] = 's';
      }
      break;
    }
    ++ltj;
  }

  return;
}

bool checkRank(const TypeLetter& x, size_t l)

/*!
  Synopsis: checks if the rank l is in the valid range for x.
*/

{
  using namespace constants;

  switch (x) {
  case 'A': // rank must be >= 1
    if ((l < 1) or (l > RANK_MAX))
      return false;
    break;
  case 'B': // rank must be >= 2
    if ((l < 2) or (l > RANK_MAX))
      return false;
    break;
  case 'C': // rank must be >= 2
    if ((l < 2) or (l > RANK_MAX))
      return false;
    break;
  case 'D': // rank must be >= 4
    if ((l < 4) or (l > RANK_MAX))
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
    if ((l < 1) or (l > RANK_MAX))
      return false;
    break;
  default: // this cannot happen!
    assert(false && "unexpected type in checkRank");
    break;
  }

  return true;
}

void involution(latticetypes::LatticeMatrix& i, const lietype::LieType& lt, 
		const lietype::InnerClassType& ic)

/*!
  Synopsis: constructs the fundamental involution for the Lie type lt and the 
  inner class ic, in the weight basis for the simply connected group.

  Precondition: it has already been checked that ic holds a valid inner class
  type for lt.
*/

{
  using namespace lietype;

  size_t n = rank(lt);
  i.resize(n,n,0);

  size_t r = 0;
  size_t pos = 0;

  for (size_t j = 0; j < ic.size(); ++j) {
    SimpleLieType slt = lt[pos];
    addSimpleInvolution(i,r,slt,ic[j]);
    if (ic[j] == 'C') { // consume an extra Lie type
      pos += 2;
      r += 2*rank(slt);
    } else {
      ++pos;
      r += rank(slt);
    }
  }

}

size_t rank(const LieType& lt)

/*!
  Synopsis: returns the rank of the group.
*/

{
  size_t r = 0;

  for (size_t j = 0; j < lt.size(); ++j)
    r += rank(lt[j]);

  return r;
}

size_t semisimpleRank(const LieType& lt)

/*!
  Synopsis: returns the semisimple rank of the group.
*/

{
  size_t r = 0;

  for (size_t j = 0; j < lt.size(); ++j)
    r += semisimpleRank(lt[j]);

  return r;
}

}

/*****************************************************************************

        Chapter II -- Private functions

  ... explain here when stable ...

******************************************************************************/

namespace lietype {

void addCompactInvolution(latticetypes::LatticeMatrix& m, size_t r, 
			  size_t rs)

/*!
  Synopsis: sets the block of size (rs,rs) starting from (r,r) to the
  identity

  Precondition: the block is set to zero.
*/

{  
  for (size_t j = 0; j < rs; ++j)
    m(r+j,r+j) = 1;

  return;
}

void addDInvolution(latticetypes::LatticeMatrix& m, size_t r, size_t rs)

/*!
  Synopsis: flips the last two vectors in the block of size rs starting
  from (r,r).

  Precondition: the block is set to zero.
*/

{  
  for (size_t j = 0; j < rs-2; ++j)
    m(r+j,r+j) = 1;

  m(r+rs-2,r+rs-1) = 1;
  m(r+rs-1,r+rs-2) = 1;

  return;
}

void addMinusIdentity(latticetypes::LatticeMatrix& m, size_t r, size_t rs)

/*!
  Synopsis: sets the block of size (rs,rs) starting from (r,r) to minus the
  identity

  Precondition: the block is set to zero.
*/

{  
  for (size_t j = 0; j < rs; ++j)
    m(r+j,r+j) = -1;

  return;
}

void addSimpleInvolution(latticetypes::LatticeMatrix& m, size_t r, 
			 const SimpleLieType& slt, TypeLetter x)

/*!
  Synopsis: appends to m, from position (r,r), the fundamental involution 
  corresponding to x in size rs.
*/

{
  size_t rs = rank(slt);

  switch (x) {
  case 'c': // add the identity
    addCompactInvolution(m,r,rs);
    break;
  case 's': // add split involution
    switch (type(slt)) {
    case 'A': // antidiagonal matrix
      for (size_t j = 0; j < rs; ++j)
	m(r+j,r+rs-1-j) = 1;
      break;
    case 'D':
      if (rank(slt)&1UL)
	addDInvolution(m,r,rs);
      else
	addCompactInvolution(m,r,rs);
      break;
    case 'E':
      if (rank(slt) == 6) {
	m(r+1,r+1) = 1;
	m(r+3,r+3) = 1;
	m(r,r+5) = 1;
	m(r+5,r) = 1;
	m(r+2,r+4) = 1;
	m(r+4,r+2) = 1;
      }
      else
	addCompactInvolution(m,r,rs);
      break;
    case 'T':
      addMinusIdentity(m,r,rs);
      break;
    default: // identity involution for types B,C,E7,E8,F,f,G,g
      addCompactInvolution(m,r,rs);
      break;
    };
    break;
  case 'C': // rs-dimensional flip
    for (size_t j = 0; j < rs; ++j) {
      m(r+j,r+rs+j) = 1;
      m(r+rs+j,r+j) = 1;
    }
    break;
  case 'u': // flip the first two vectors
    addDInvolution(m,r,rs);
    break;
  default: // this should not happen!
    assert(false && "wrong inner class letter in addSimpleInvolution");
    break;
  };

  return;
}

}

}

