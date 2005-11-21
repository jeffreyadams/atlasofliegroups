/*
  This is abelian.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "abelian.h"

#include <cassert>

#include "arithmetic.h"
#include "bitmap.h"
#include "latticetypes.h"
#include "matrix.h"
#include "smithnormal.h"

/*****************************************************************************

  This is a partial and tentative implementation of the concept of a finite
  abelian group. Typically we have in mind the center of a reductive semisimple
  group.

******************************************************************************/

/*****************************************************************************

        Chapter I -- The FiniteAbelianGroup class

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace abelian {

/******** constructors and destructors ***************************************/

FiniteAbelianGroup::FiniteAbelianGroup(const std::vector<unsigned long>& t)
  :d_type(t),d_cotype(t.size())

{
  d_size = 1;

  for (size_t j = 0; j < t.size(); ++j)
    d_size *= t[j];

  for (size_t j = 0; j < t.size(); ++j)
    d_cotype[j] = t.back()/t[j];
}

/******** type conversions ***************************************************/

void FiniteAbelianGroup::toWeight(latticetypes::Weight& v, const GrpArr& a) 
  const

/*
  Synopsis: put in v a representative of a.

  NOTE: sloppy implementation; we don't check for overflow, which may happen
  in all cases, as the coefficients of v will be signed quantities.
*/

{
  v.assign(a.begin(),a.end());

  return;
}

void FiniteAbelianGroup::toWeight(latticetypes::Weight& v, GrpNbr x) const

/*
  Synopsis: put in v a representative of a.

  NOTE: sloppy implementation; we don't check for overflow, which may happen
  in all cases, as the coefficients of v will be signed quantities.
*/

{    
  const GroupType& t = type();

  for (size_t j = 0; j < rank(); ++j) {
    v[j] = x % t[j];
    x /= t[j];
  }

  return;
}

/******** accessors **********************************************************/

GrpArr& FiniteAbelianGroup::add(GrpArr& a, const GrpArr& b) const

/*
  Synopsis: a += b.

  Precondition: a and b hold valid arrays for the group;

  The only difficulty is doing it without triggering overflow.
*/

{
  for (size_t j = 0; j < a.size(); ++j) {
    if (a[j] < d_type[j] - b[j])
      a[j] += b[j];
    else // overflow
      a[j] -= d_type[j] - b[j];
  }

  return a;
}

GrpArr& FiniteAbelianGroup::add(GrpArr& a, GrpNbr x) const

/*
  Synopsis: a += x.
*/

{
  GrpArr b(rank());
  toArray(b,x);

  return add(a,b);
}

GrpNbr FiniteAbelianGroup::add(GrpNbr x, const GrpArr& b) const

/*
  Synopsis: x += b.

  The main problem is to deal with the modular addition, so that the overflow 
  is handled correctly.
*/

{
  GrpArr a(d_type.size());
  toArray(a,x);
  add(a,b);

  return toGrpNbr(a);
}

GrpNbr FiniteAbelianGroup::add(GrpNbr x, GrpNbr y) const

/*
  Synopsis: x += y.

  The main problem is to deal with the modular addition, so that the overflow 
  is handled correctly.
*/

{
  GrpArr a(rank());
  toArray(a,x);

  for (size_t j = 0; j < a.size(); ++j) {
    unsigned long yj = y % d_type[j];
    if (a[j] < d_type[j] - yj)
      a[j] += yj;
    else // overflow
      a[j] -= d_type[j] - yj;
    y /= d_type[j];
  }

  return toGrpNbr(a);
}

unsigned long FiniteAbelianGroup::annihilator() const

/*
  Synopsis: returns the s.c.m. of the orders of the elements of the group.
*/

{
  if (d_type.size())
    return d_type.back();
  else
    return 1;
}

GrpNbr FiniteAbelianGroup::leftApply(GrpNbr x, const Endomorphism& q) const

/*
  Synopsis: applies the matrix q to the element x, on the left.

  The idea is that the matrix defines an endomorphism of the group in terms
  of the array representation.
*/

{
  GrpArr a(rank());
  toArray(a,x);
  leftApply(a,q);

  return toGrpNbr(a);
}

GrpArr& FiniteAbelianGroup::leftApply(GrpArr& a, const Endomorphism& q) const

/*
  Synopsis: applies the matrix q to the array a, on the left.

  The idea is that the matrix defines an endomorphism of the group in terms
  of the array representation.
*/

{
  using namespace arithmetic;

  GrpArr tmp = a;

  for (size_t i = 0; i < q.numRows(); ++i) {
    unsigned long r = 0;
    for (size_t j = 0; j < q.numColumns(); ++j) {
      unsigned long s = tmp[j];
      modProd(s,q(i,j),d_type[i]);
      modAdd(r,s,d_type[i]);
    }
    a[i] = r;
  }

  return a;
}

unsigned long FiniteAbelianGroup::order(GrpNbr x) const

/*
  Synopsis: computes the order of x in the group.
*/

{
  using namespace arithmetic;

  unsigned long n = 1;
  GrpArr v(rank());
  toArray(v,x);

  for (size_t j = 0; j < d_type.size(); ++j) {
    unsigned long m = gcd(static_cast<long>(v[j]),d_type[j]);
    n = lcm(n,d_type[j]/m);
  }

  return n;
}

unsigned long FiniteAbelianGroup::order(const bitmap::BitMap& B, 
					GrpNbr x) const

/*
  Synopsis: computes the order of x modulo B.

  NOTE : we have not tried to be smart at all here. For large groups this
  would be very unsatisfactory.
*/

{
  if (B.isMember(x))
    return 1;

  unsigned long n = 2;
  GrpNbr xn = x;

  for (;; ++n) {
    xn = add(xn,x);
    if (B.isMember(xn))
      break;
  }

  return n;
}

unsigned long FiniteAbelianGroup::pairing(const GrpArr& a, const GrpArr& b) 
  const

/*
  Synopsis: computes the m in [0,n[ s.t. a(b) = e^{2i pi m/n}, where a is 
  interpreted as an element of the dual group, b as an element of the group, 
  and n is the annihilator of the group (the last entry in d_type).

  NOTE: this is a sloppy implementation, that doesn't deal carefully with
  overflow. It is expected to be used only for very small groups.
*/

{
  unsigned long m = 0;

  for (size_t j = 0; j < rank(); ++j)
    m += a[j]*b[j]*d_cotype[j];

  m %= annihilator();

  return m;
}

unsigned long FiniteAbelianGroup::pairing(const GrpArr& a, const GrpArr& b,
					  unsigned long t) 
  const

/*
  Synopsis: computes the m in [0,t[ s.t. a(b) = e^{2i pi m/t}, where a is 
  interpreted as an element of the dual group, b as an element of the group.

  Precondition: pairing(a,b) is an element of t-torsion in Z/n, where n is
  the annihilator of the group;

  The required m is just pairing(a,b)/(n/t).

  NOTE: this is a sloppy implementation, that doesn't deal carefully with
  overflow. It is expected to be used only for very small groups.

  NOTE: we put in an assertion for safety.
*/

{
  unsigned long n = annihilator();
  unsigned long r = pairing(a,b)*t;

  assert(r%n == 0);

  return r/n;
}

GrpNbr FiniteAbelianGroup::prod(GrpNbr x, unsigned long n) const

/*
  Replaces x with n.x. We use the classic logarithmic algorithm, where the
  result is obtained through a sequence of additions and multiplications by
  two.
*/

{
  using namespace constants;

  if (n == 0)
    return 0;

  unsigned long p = n;

  for (; ~p & hiBit; p <<= 1)  /* shift n up to high powers */
    ;

  GrpNbr y = x; // save the original value of x
    
  for (unsigned long j = n >> 1; j; j >>= 1) 
    {
      p <<= 1;
      x = add(x,x);
      if (p & hiBit)
	x = add(x,y);
    }

  return x;
}

GrpArr& FiniteAbelianGroup::subtract(GrpArr& a, const GrpArr& b) const

/*
  Synopsis: a -= b.

  Precondition: a and b hold valid arrays for the group;

  The only difficulty is doing it without triggering overflow.
*/

{
  for (size_t j = 0; j < a.size(); ++j) {
    if (b[j] <= a[j])
      a[j] -= b[j];
    else // underflow
      a[j] += d_type[j] - b[j];
  }

  return a;
}

}

/*****************************************************************************

        Chapter II -- The Homomorphism class

  A Homomorhism object represents a homomorphism from a group of type d_source
  to a group of type d_dest, where the elements are represented as arrays.

  As we know, for such a homomorphism to be well defined, the matrix entry
  (i,j) should be a multiple of d_dest[i]/g, with g = 
  gcd(d_dest[i],d_source[j]).

  Here we take a different approach. We set M a common multiple of the
  annihilators of d_source and d_dest. Then
  we allow arbitrary coefficients in Z/M. The corresponding matrix is
  interpreted as follows. It is defined only for those tuples (x_1,...,x_m)
  s.t. for all i, sum_j a_ij q_j x_j is of index t_i in Z/M, where (q_j)
  is the cotype of d_source w.r.t. M; and then the corresponding value y_i is
  (sum_j a_ij q_j x_j)/(M/t_i). It is easy to see that this definition is
  in fact independent of the choice of M (of course we will take the lcm
  of the annihilators.)

******************************************************************************/

namespace abelian {

/******** constructors and destructors ***************************************/

Homomorphism::Homomorphism(const std::vector<GrpArr>& al, 
			   const GroupType& source, const GroupType& dest)
  :d_source(source),
   d_dest(dest),
   d_cosource(source.size()),
   d_codest(dest.size()),
   d_matrix(dest.size(),source.size())

/*
  Synopsis: constructs the homomorphism with matrix al, from a group of
  type source to a group of type dest.

  Precondition: the elements of al are of size source.size(), and their
  number is dest.size();
*/

{
  using namespace arithmetic;

  for (size_t i = 0; i < dest.size(); ++i)
    for (size_t j = 0; j < source.size(); ++j)
      d_matrix(i,j) = al[i][j];

  unsigned long a = 1;
  if (source.size())
    a = source.back(); // annihilator of source

  unsigned long b = 1;
  if (dest.size())
    b = dest.back();   // annihilator of dest

  d_annihilator = lcm(a,b);

  for (size_t j = 0; j < d_cosource.size(); ++j)
    d_cosource[j] = d_annihilator/d_source[j];

  for (size_t j = 0; j < d_codest.size(); ++j)
    d_codest[j] = d_annihilator/d_dest[j];
}

/******** accessors **********************************************************/

void Homomorphism::apply(GrpArr& dest, const GrpArr& source) const

/*
  Synopsis: applies the homomorphism to source according to the rules explained
  in the introduction to this section, and puts the result in dest.

  Precondition: source is in the subgroup for which this makes sense; dest
  is already allocated to the correct size.

  NOTE: we put in an assertion for safety.

  NOTE: sloppy implementation; we don't worry about overflow.
*/

{
  for (size_t i = 0; i < dest.size(); ++i) {
    dest[i] = 0;
    for (size_t j = 0; j < source.size(); ++j) {
      dest[i] += d_matrix(i,j)*d_cosource[j]*source[j];
    }
    assert((dest[i]*d_dest[i])%d_annihilator == 0);
    dest[i] /= (d_annihilator/d_dest[i]);
    dest[i] %= d_dest[i];
  }

  return;
}

GrpNbr Homomorphism::apply(GrpNbr source) const

/*
  Synopsis: return h(source).

  Precondition: source is in the subgroup for which this makes sense;

  Forwarded to the GrpArr form.
*/

{
  GrpArr a(d_source.size());
  toArray(a,source,d_source);
  GrpArr dest(d_dest.size());
  apply(dest,a);

  return toGrpNbr(dest,d_dest);
}

bool Homomorphism::defined(const GrpArr& a) const

/*
  Synopsis: tells whether a is in the domain.

  This means that a satisfies the congruences stated in the apply function.
*/

{  
  for (size_t i = 0; i < d_dest.size(); ++i) {
    unsigned long b = 0;
    for (size_t j = 0; j < a.size(); ++j) {
      b += d_matrix(i,j)*d_cosource[j]*a[j];
    }
    if ((b*d_dest[i])%d_annihilator != 0)
      return false;
  }

  return true;
}

bool Homomorphism::defined(GrpNbr x) const

/*
  Synopsis: tells whether x is in the domain.

  Forwarded to the array-version.
*/

{
  GrpArr a(d_source.size());
  toArray(a,x,d_source);
  
  return defined(a);
}

}

/*****************************************************************************

        Chapter III -- Functions declared in abelian.h

  ... explain here when it is stable ...

******************************************************************************/

namespace abelian {

void basis(latticetypes::WeightList& b, const bitmap::BitMap& B,
	   const FiniteAbelianGroup& A) 

/*
  Synopsis: writes A/B in canonical form.

  Explanation: we see the current group A as a quotient of Z^d, where d is 
  the rank of A, and the generators of the kernel are given by d_type. Then we 
  wish to write the quotient A/B in canonical form (i.e., as a product of
  cyclic groups with cardinalities dividing each other.) For this, we put in b 
  a scaled Smith normal basis for the inverse image of B in Z^d.

  Note that A is not necessarily in canonical form, so even when B is the 
  trivial subgroup this might yield a basis rather different from the kernel 
  basis.
*/

{  
  using namespace bitmap;
  using namespace latticetypes;
  using namespace matrix;
  using namespace smithnormal;

  WeightList gl;

  // put in the kernel basis

  for (size_t j = 0; j < A.rank(); ++j) {
    gl.push_back(Weight(A.rank(),0));
    gl.back()[j] = A.type()[j];
  }

  // add generators of the subgroup

  GrpNbrList gen;
  generators(gen,B,A);

  for (size_t j = 0; j < gen.size(); ++j) {
    Weight v(A.rank());
    A.toWeight(v,gen[j]);
    gl.push_back(v);
  }

  // write matrix corresponding to gl

  LatticeMatrix m(gl);
  CoeffList invf;

  // get smith normal basis

  initBasis(b,A.rank());
  smithNormal(invf,b.begin(),m);

  // scale

  for (size_t j = 0; j < invf.size(); ++j)
    b[j] *= invf[j];

  return;
}

void coset(bitmap::BitMap& C, const bitmap::BitMap& B, GrpNbr x,
	   const FiniteAbelianGroup& A)

/*
  Synopsis: puts in C the coset x+B in A.

  Precondition: C.size() == A.size();

  NOTE : this is a straightforward implementation, shifting elements of B
  individually.
*/

{
  using namespace bitmap;

  C.reset();
  BitMap::iterator B_end = B.end();

  for (BitMap::iterator i = B.begin(); i != B_end; ++i) {
    GrpNbr y = *i;
    C.insert(A.add(y,x));
  }

  return;
}

const bitmap::BitMap& cycGenerators(const FiniteAbelianGroup& A)

/*
  Synopsis: returns a reference to a bitmap containing exactly one generator
  for each cyclic subgroup of A.

  NOTE: it is expected that this function will typically be called repeatedly
  for the same group. The bitmap is constructed on the first call for the
  given group, following the principle of lazy evaluation.
*/

{
  using namespace arithmetic;
  using namespace bitmap;

  static FiniteAbelianGroup Astat;
  static bitmap::BitMap cyc;

  if (Astat.type() != A.type()) { // update cyc

    cyc.resize(A.size());
    cyc.fill();
    BitMap::iterator cyc_end = cyc.end();

    for (BitMap::iterator i = ++cyc.begin(); i != cyc_end; ++i) {
      GrpNbr x = *i;
      unsigned long n = A.order(x);
      for (unsigned long j = 2; j < n; ++j) {
	if (gcd(n,j) != 1) // j is not prime to n
	  continue;
	GrpNbr xj = A.prod(x,j);
	cyc.remove(xj);
      }
    }

    Astat = A;

  }

  return cyc;
}

void generateSubgroup(bitmap::BitMap& B, GrpNbr x, const FiniteAbelianGroup& A)

/*
  Synopsis: transforms B into the subgroup generated by B and x in A.

  NOTE : this is a simple-minded implementation; we do not aim for speed.
*/

{
  using namespace bitmap;

  unsigned long n = A.order(x);
  BitMap b(B);
  BitMap c(A.size());

  for (unsigned long j = 1; j < n; ++j) {
    GrpNbr xj = A.prod(x,j);
    coset(c,B,xj,A);
    b |= c;
  }

  B.swap(b);

  return;
}

void generators(GrpNbrList& gen, const bitmap::BitMap& B,
		const FiniteAbelianGroup& A)

/*
  Synopsis: puts in gen a list of generators of the subgroup B.
*/

{
  gen.assign(B.begin(),B.end());

  return;
}

bool isElementaryAbelian(const std::vector<unsigned long>& c)

/*
  Tells if c is all twos.
*/

{
  for (size_t j = 0; j < c.size(); ++j)
    if (c[j] != 2)
      return false;

  return true;
}

void quotReps(bitmap::BitMap& qr, const bitmap::BitMap& B,
	      const FiniteAbelianGroup& A) 

/*
  Synopsis: puts in qr a list of representatives of the cosets modulo B. 

  Precondition: qr.size() = A.size();

  Algorithm: fill up qr; traverse qr, and for each x that is reached, remove
  from qr the elements of the coset of x other than x.

  The efficiency depends on the efficiency of coset computations.
*/

{
  using namespace bitmap;

  qr.fill();
  BitMap::iterator qr_end = qr.end();

  for (BitMap::iterator i = qr.begin(); i != qr_end; ++i) {
    BitMap c(A.size());
    coset(c,B,*i,A);
    c.remove(*i);
    qr.andnot(c);
  }

  return;
}

void toArray(GrpArr& a, GrpNbr x, const GroupType& t)

/*
  Synopsis: puts the array-form of x in a.

  Precondition: x is representable for t (x < prod t[j]); otherwise the
  result is the expression of x modulo that product.

  Explanation: the array-form (relative to type) is the unique expression
  of x in the variable base defined by type.
*/

{    
  for (size_t j = 0; j < a.size(); ++j) {
    a[j] = x % t[j];
    x /= t[j];
  }

  return;
}

void toArray(GrpArr& a, const latticetypes::Weight& v, const GroupType& t)

/*
  Synopsis: reduces v mod d_type.

  The only difficulty is to make sure that negative values are reduced in
  the way that we want.

  Precondition: a is set to v.size();
*/

{
  using namespace arithmetic;

  for (size_t j = 0; j < v.size(); ++j)
    a[j] = remainder(v[j],t[j]);

  return;
}

void toEndomorphism(Endomorphism& e, const latticetypes::LatticeMatrix& q,
		    const FiniteAbelianGroup& A)

/*
  Synopsis: transforms q into the corresponding Endomorphism.

  This just involves rewriting the coeficients as unsigned longs modulo
  the type factors. Just in case the conventions for modular arithmetic
  with negative numbers are not what we want, we make sure that we
  deal with nonnegative numbers all along.

  NOTE: there should be some better way to do this!
*/

{
  using namespace arithmetic;

  e.resize(q.numRows(),q.numColumns());

  for (size_t j = 0; j < q.numColumns(); ++j)
    for (size_t i = 0; i < q.numRows(); ++i)
      e(i,j) = remainder(q(i,j),A.type()[i]);
     
  return;
}

GrpNbr toGrpNbr(const GrpArr& a, const GroupType& t)

/*
  Synopsis: returns the number-form of a.

  Precondition: a is representable as a GrpNbr;
*/

{    
  GrpNbr x = 0;

  for (size_t j = a.size(); j;) {
    --j;
    x *= t[j];
    x += a[j];
  }

  return x;
}

void transpose(Endomorphism& e, const FiniteAbelianGroup& A)

/*
  Synopsis: transposes the endomorphism e.

  This is more subtle than one might think; the main point is that when
  m|n, the transpose of the canonical map Z/nZ -> Z/mZ is the injection
  Z/mZ -> Z/nZ which takes 1 to n/m. In general, any map Z/nZ -> Z/mZ
  factors thru Z/dZ, where d = gcd(m,n), so it is a multiple of the map
  defined by 1 -> m/d; the transpose is the same multiple of 1 -> n/d.

  NOTE: this implementation works or homomorphisms between diferent groups
  just as well (except that one needs two grouptypes then.)
*/

{
  using namespace arithmetic;

  const GroupType& type = A.type();
  Endomorphism tmp = e;
  e.resize(tmp.numColumns(),tmp.numRows());

  for (size_t j = 0; j < e.numColumns(); ++j)
    for (size_t i = 0; i < e.numRows(); ++i) {
      unsigned long d = gcd(type[i],type[j]);
      // tmp(i,j) has to be divisible by type[i]/d
      unsigned long ni = type[i]/d;
      unsigned long nj = type[j]/d;
      e(j,i) = nj*(tmp(i,j)/ni);
    }

  return;
}

}

}
