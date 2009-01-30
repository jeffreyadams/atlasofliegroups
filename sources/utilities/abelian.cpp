/*
  This is abelian.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "abelian.h"

#include <cassert>

#include "arithmetic.h"
#include "bitmap.h"
#include "latticetypes.h"
#include "matrix.h"
#include "smithnormal.h"

/*!****************************************************************************
\file
  This is a partial and tentative implementation of the concept of a finite
  abelian group. Typically we have in mind the center of a reductive
  semisimple group.

******************************************************************************/

/*****************************************************************************

        Chapter I -- The FiniteAbelianGroup class

******************************************************************************/

namespace atlas {

namespace abelian {

/******** constructors and destructors ***************************************/

FiniteAbelianGroup::FiniteAbelianGroup(const std::vector<unsigned long>& t)
  :d_type(t),d_cotype(t.size())
{
  d_size = 1;

  for (size_t j = 0; j < t.size(); ++j)
  {
    d_size *= t[j];
    assert(t.back()%t[j]==0);
    d_cotype[j] = t.back()/t[j];
  }
}

/******** type conversions ***************************************************/

/* What is a |GrpArr|, and what is a representative weight for it?
   Why is this a method of |FiniteAbelianGroup|, without using |*this|? MvL
*/

/*!
  Synopsis: put in |v| a representative of |a|.

  NOTE: sloppy implementation; we don't check for overflow, which may happen
  in all cases, as the coefficients of v will be signed quantities.
*/
void FiniteAbelianGroup::toWeight(latticetypes::Weight& v, const GrpArr& a)
  const
{
  v.assign(a.begin(),a.end()); // copy all elements; |v=a| would mismatch type
}


/*!
  Synopsis: put in |v| a representative of |x|.
*/
void FiniteAbelianGroup::toWeight(latticetypes::Weight& v, GrpNbr x) const
{
  const GroupType& t = type();

  for (size_t j = 0; j < rank(); ++j) {
    v[j] = x % t[j]; // extract coefficient in cyclic factor |j|
    x /= t[j];
  }
}

/******** accessors **********************************************************/


/*!
  Synopsis: a += b.

  Precondition: a and b hold valid arrays for the group;

  The only difficulty is doing it without triggering overflow.
*/
GrpArr& FiniteAbelianGroup::add(GrpArr& a, const GrpArr& b) const
{
  for (size_t j = 0; j < a.size(); ++j) {
    if (a[j] < d_type[j] - b[j])
      a[j] += b[j];
    else // overflow
      a[j] -= d_type[j] - b[j];
  }

  return a;
}


/*!
  Synopsis: a += x.
*/
GrpArr& FiniteAbelianGroup::add(GrpArr& a, GrpNbr x) const
{
  return add(a,toArray(x));
}


/*!
  Synopsis: x += b.

  The main problem is to deal with the modular addition, so that the overflow
  is handled correctly.
*/
GrpNbr FiniteAbelianGroup::add(GrpNbr x, const GrpArr& b) const
{
  GrpArr a=toArray(x);
  return toGrpNbr(add(a,b));
}


/*!
  Synopsis: x += y.

  The main problem is to deal with the modular addition, so that the overflow
  is handled correctly; |arithmetic::modAdd| does just that
*/
GrpNbr FiniteAbelianGroup::add(GrpNbr x, GrpNbr y) const
{
  GrpArr a=toArray(x);

  for (size_t j = 0; j < a.size(); ++j)
  {
    arithmetic::modAdd(a[j], y%d_type[j], d_type[j]);
    y /= d_type[j];
  }

  return toGrpNbr(a);
}


/*!
  Synopsis: returns the l.c.m. of the orders of the elements of the group.
*/
unsigned long FiniteAbelianGroup::annihilator() const
{
  return d_type.size()==0 ? 1 : d_type.back();
}


/*!
  Synopsis: applies the matrix q to the element x, on the left.

  The idea is that the matrix defines an endomorphism of the group in terms
  of the array representation.
*/
GrpNbr FiniteAbelianGroup::leftApply(GrpNbr x, const Endomorphism& q) const
{
  GrpArr a=toArray(x);
  leftApply(a,q);

  return toGrpNbr(a);
}


/*!
  Synopsis: applies the matrix q to the array a, on the left.

  The idea is that the matrix defines an endomorphism of the group in terms
  of the array representation.
*/
GrpArr& FiniteAbelianGroup::leftApply(GrpArr& a, const Endomorphism& q) const
{
  GrpArr tmp = a;

  for (size_t i = 0; i < q.numRows(); ++i) {
    unsigned long r = 0;
    for (size_t j = 0; j < q.numColumns(); ++j) {
      unsigned long s = tmp[j];
      arithmetic::modProd(s,q(i,j),d_type[i]);
      arithmetic::modAdd(r,s,d_type[i]);
    }
    a[i] = r;
  }

  return a;
}


/*!
  Synopsis: computes the order of x in the group.
*/
unsigned long FiniteAbelianGroup::order(GrpNbr x) const
{
  unsigned long n = 1;
  GrpArr v=toArray(x);

  for (size_t i=0; i<d_type.size(); ++i)
  {
    unsigned long m = arithmetic::unsigned_gcd(v[i],d_type[i]);
    n = arithmetic::lcm(n,d_type[i]/m);
  }

  return n;
}

unsigned long FiniteAbelianGroup::order(const bitmap::BitMap& B,
					GrpNbr x) const

/*!
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


/*!
  Synopsis: computes the m in [0,n[ s.t. a(b) = e^{2i pi m/n}, where a is
  interpreted as an element of the dual group, b as an element of the group,
  and n is the annihilator of the group (the last entry in d_type).

  NOTE: this is a sloppy implementation, that doesn't deal carefully with
  overflow. It is expected to be used only for very small groups.
*/
unsigned long FiniteAbelianGroup::pairing(const GrpArr& a, const GrpArr& b)
  const
{
  unsigned long m = 0;

  for (size_t i = 0; i < rank(); ++i)
    m += a[i]*b[i]*d_cotype[i];

  m %= annihilator();

  return m;
}

unsigned long FiniteAbelianGroup::pairing(const GrpArr& a, const GrpArr& b,
					  unsigned long t)
  const

/*!
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


/*!
  Replaces x with n.x. We use the classic logarithmic algorithm, where the
  result is obtained through a sequence of additions and multiplications by
  two.
*/
GrpNbr FiniteAbelianGroup::prod(GrpNbr x, unsigned long n) const
{
  if (n == 0)
    return 0;

  unsigned long p = n;

  while((p & constants::hiBit)==0)
    p <<= 1;  /* shift n up to high powers */

  GrpNbr y = x; // save the original value of x

  for (unsigned long j = n >> 1; j!=0; j >>= 1)
    {
      p <<= 1;
      x = add(x,x);
      if ((p & constants::hiBit)!=0)
	x = add(x,y);
    }

  return x;
}


/*!
  Synopsis: a -= b.

  Precondition: a and b hold valid arrays for the group;

  The only difficulty is doing it without triggering overflow.
*/
GrpArr& FiniteAbelianGroup::subtract(GrpArr& a, const GrpArr& b) const
{
  for (size_t j = 0; j < a.size(); ++j)
    if (b[j] <= a[j])
      a[j] -= b[j];
    else // underflow
      a[j] += d_type[j] - b[j];

  return a;
}

}

/*****************************************************************************

        Chapter II -- The Homomorphism class

  A |Homomorphism| object represents a homomorphism from a group of type
  |d_source| to a group of type |d_dest|. The essential data is a matrix
  |d_matrix| with integral coefficients. As usual the matrix entry (i,j)
  expresses component $i$ of the image of generator $j$.

  For such a homomorphism to be well defined, the matrix entry (i,j) should be
  a multiple of |d_dest[i]/g|, with |g == gcd(d_dest[i],d_source[j])|.

  Computationally we shall proceed as follows. We set $M$ to a common multiple
  of the annihilators of |d_source| and |d_dest|, and think of all cyclic
  groups as embedded in $Z/M$; this means that before acting upon a source
  element $(x_1,...,x_m)$ we should multiply each component $x_j$ by the
  cotype $q_j=M/s_j$ of that component with respect to $M$, where $s_j$ is
  |d_source[j]|. After forming $y_j=\sum_j a_{i,j} q_j x_j$ for each $i$, the
  component $y_i$ must be in the subgroup of index $t_i$ in $Z/M$, where $t_i$
  is |d_dest[i]|, in other words it must be divisible by $M/t_i$, and the
  compoentn $i$ of the final result will be $y_i/(M/t_i)$. This definition is
  independent of the choice of $M$; we will take the lcm of the annihilators.

  In fact this definition might work for some group elements of the source but
  not for all. In that case we don't really have a homomorphism, but we allow
  applying to those group elements for which it works nontheless. The
  predicates |defined| serve to find out whether a source element is OK.

******************************************************************************/

namespace abelian {

/******** constructors and destructors ***************************************/


/*!
  Constructs the homomorphism with matrix al, from groups |source| to |dest|.

  Precondition: the elements of |al| are of size source.rank(), and their
  number is dest.rank(). Note that the matrix is therefore given by rows
  (linear forms on |source|), not by columns (images in |dest|); unusal!
*/
Homomorphism::Homomorphism(const std::vector<GrpArr>& al,
			   const FiniteAbelianGroup& source,
			   const FiniteAbelianGroup& dest)
  :d_source(source.type()),
   d_dest(dest.type()),
   d_cosource(d_source.size()),
   d_codest(d_dest.size()),
   d_matrix(d_dest.size(),d_source.size())
{
  for (size_t i = 0; i < dest.size(); ++i)
    for (size_t j = 0; j < source.size(); ++j)
      d_matrix(i,j) = al[i][j];

  d_annihilator = arithmetic::lcm(source.annihilator(),dest.annihilator());

  for (size_t j = 0; j < d_cosource.size(); ++j)
    d_cosource[j] = d_annihilator/d_source[j];

  for (size_t j = 0; j < d_codest.size(); ++j)
    d_codest[j] = d_annihilator/d_dest[j];
}

/******** accessors **********************************************************/

/*!
  Synopsis: applies the homomorphism to source according to the rules explained
  in the introduction to this section, and puts the result in dest.

  Precondition: source is in the subgroup for which this makes sense.

  NOTE: sloppy implementation; we don't worry about overflow.
*/
GrpArr Homomorphism::apply(const GrpArr& source) const
{
  GrpArr dest(d_dest.size(),0);
  for (size_t i = 0; i < dest.size(); ++i)
  {
    for (size_t j = 0; j < source.size(); ++j)
      dest[i] += d_matrix(i,j)*d_cosource[j]*source[j];

    assert(dest[i]%(d_annihilator/d_dest[i]) == 0);
    dest[i] /= (d_annihilator/d_dest[i]);
    dest[i] %= d_dest[i];
  }

  return dest;
}


/*!
  Synopsis: return h(source).

  Precondition: source is in the subgroup for which this makes sense;

  Forwarded to the GrpArr form.
*/
GrpNbr Homomorphism::apply(GrpNbr source) const
{
  GrpArr a(d_source.size());
  to_array(a,source,d_source);

  return to_GrpNbr(apply(a),d_dest);
}


/*!
  Synopsis: tells whether a is in the domain.

  This means that a satisfies the congruences stated in the apply function.
*/
bool Homomorphism::defined(const GrpArr& a) const
{
  for (size_t i = 0; i < d_dest.size(); ++i)
  {
    unsigned long b = 0;
    for (size_t j = 0; j < a.size(); ++j)
      b += d_matrix(i,j)*d_cosource[j]*a[j];

    if (b%(d_annihilator/d_dest[i]) != 0)
      return false;
  }

  return true;
}


/*!
  Synopsis: tells whether x is in the domain.

  Forwarded to the array-version.
*/
bool Homomorphism::defined(GrpNbr x) const
{
  GrpArr a(d_source.size());
  to_array(a,x,d_source);

  return defined(a);
}

} // |namespace abelian|

/*****************************************************************************

        Chapter III -- Functions declared in abelian.h

******************************************************************************/

namespace abelian {

/*!
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
void basis(latticetypes::WeightList& b, const bitmap::BitMap& B,
	   const FiniteAbelianGroup& A)
{
  latticetypes::WeightList gl;

  // put in the kernel basis

  for (size_t j = 0; j < A.rank(); ++j)
  {
    gl.push_back(latticetypes::Weight(A.rank(),0));
    gl.back()[j] = A.type()[j];
  }

  // add generators of the subgroup

  GrpNbrList gen;
  generators(gen,B,A);

  for (size_t j = 0; j < gen.size(); ++j)
  {
    latticetypes::Weight v(A.rank());
    A.toWeight(v,gen[j]);
    gl.push_back(v);
  }

  // write matrix corresponding to gl

  latticetypes::LatticeMatrix m(gl);
  latticetypes::CoeffList invf;

  // get smith normal basis

  matrix::initBasis(b,A.rank());
  smithnormal::smithNormal(invf,b.begin(),m);

  // scale

  for (size_t j = 0; j < invf.size(); ++j)
    b[j] *= invf[j];

  return;
}


/*!
  Synopsis: puts in C the coset x+B in A.

  Precondition: C.size() == A.size();

  NOTE : this is a straightforward implementation, shifting elements of B
  individually.
*/
void coset(bitmap::BitMap& C, const bitmap::BitMap& B, GrpNbr x,
	   const FiniteAbelianGroup& A)
{
  C.reset();
  bitmap::BitMap::iterator B_end = B.end();

  for (bitmap::BitMap::iterator i = B.begin(); i != B_end; ++i)
  {
    GrpNbr y = *i;
    C.insert(A.add(y,x));
  }
}

const bitmap::BitMap& cycGenerators(const FiniteAbelianGroup& A)

/*!
  Synopsis: returns a reference to a bitmap containing exactly one generator
  for each cyclic subgroup of A.

  NOTE: it is expected that this function will typically be called repeatedly
  for the same group. The bitmap is constructed on the first call for the
  given group, following the principle of lazy evaluation.
*/

{
  static FiniteAbelianGroup Astat;
  static bitmap::BitMap cyc;

  if (Astat.type() != A.type()) { // update cyc

    cyc.set_capacity(A.size());
    cyc.fill();
    bitmap::BitMap::iterator cyc_end = cyc.end();

    for (bitmap::BitMap::iterator i = ++cyc.begin(); i != cyc_end; ++i)
    {
      GrpNbr x = *i;
      unsigned long n = A.order(x);
      for (unsigned long j = 2; j < n; ++j)
      {
	if (arithmetic::unsigned_gcd(n,j) != 1) // j is not prime to n
	  continue;
	GrpNbr xj = A.prod(x,j);
	cyc.remove(xj);
      }
    }

    Astat = A;

  }

  return cyc;
}


/*!
  Synopsis: transforms B into the subgroup generated by B and x in A.

  NOTE : this is a simple-minded implementation; we do not aim for speed.
*/
void generateSubgroup(bitmap::BitMap& B, GrpNbr x, const FiniteAbelianGroup& A)
{
  unsigned long n = A.order(x);
  bitmap::BitMap b(B);
  bitmap::BitMap c(A.size());

  for (unsigned long j = 1; j < n; ++j)
  {
    GrpNbr xj = A.prod(x,j);
    coset(c,B,xj,A);
    b |= c;
  }

  B.swap(b);
}


/*!
  Synopsis: puts in gen a list of generators of the subgroup B.
*/
void generators(GrpNbrList& gen, const bitmap::BitMap& B,
		const FiniteAbelianGroup& A)
{
  gen.assign(B.begin(),B.end());
}


/*!
  Tells if c is all twos.
*/
bool isElementaryAbelian(const std::vector<unsigned long>& c)
{
  for (size_t j = 0; j < c.size(); ++j)
    if (c[j] != 2)
      return false;

  return true;
}


/*!
  Synopsis: puts in qr a list of representatives of the cosets modulo B.

  Precondition: qr.size() = A.size();

  Algorithm: fill up qr; traverse qr, and for each x that is reached, remove
  from qr the elements of the coset of x other than x.

  The efficiency depends on the efficiency of coset computations.
*/
void quotReps(bitmap::BitMap& qr, const bitmap::BitMap& B,
	      const FiniteAbelianGroup& A)
{
  qr.fill();
  bitmap::BitMap::iterator qr_end = qr.end();

  for (bitmap::BitMap::iterator i = qr.begin(); i != qr_end; ++i)
  {
    bitmap::BitMap c(A.size());
    coset(c,B,*i,A);
    c.remove(*i);
    qr.andnot(c);
  }
}


/*!
  Synopsis: puts the array-form of x in a.

  Precondition: x is representable for t (x < prod t[j]); otherwise the
  result is the expression of x modulo that product.

  Explanation: the array-form (relative to type) is the unique expression
  of x in the variable-radix base defined by type.
*/
void to_array(GrpArr& a, GrpNbr x, const GroupType& t)
{
  for (size_t j = 0; j < a.size(); ++j) {
    a[j] = x % t[j];
    x /= t[j];
  }
}


/*!
  Synopsis: reduces v mod d_type.

  The only difficulty is to make sure that negative values are reduced in
  the way that we want.

  Precondition: a is set to v.size();
*/
void to_array(GrpArr& a, const latticetypes::Weight& v, const GroupType& t)
{
  for (size_t j = 0; j < v.size(); ++j)
    a[j] = arithmetic::remainder(v[j],t[j]);
}


/*!
  Synopsis: transforms q into the corresponding Endomorphism.

  This just involves rewriting the coeficients as unsigned longs modulo
  the type factors.
*/
void toEndomorphism(Endomorphism& e, const latticetypes::LatticeMatrix& q,
		    const FiniteAbelianGroup& A)
{
  e.resize(q.numRows(),q.numColumns());

  for (size_t j = 0; j < q.numColumns(); ++j)
    for (size_t i = 0; i < q.numRows(); ++i)
      e(i,j) = arithmetic::remainder(q(i,j),A.type()[i]);
}


/*!
  Synopsis: returns the number-form of a.

  Precondition: a is representable as a GrpNbr;
*/
GrpNbr to_GrpNbr(const GrpArr& a, const GroupType& t)
{
  GrpNbr x = 0;

  for (size_t j = a.size(); j-->0;)
    (x *= t[j]) += a[j];

  return x;
}


/*!
  Synopsis: transposes the endomorphism e.

  This is more subtle than one might think; the main point is that when
  m|n, the transpose of the canonical map Z/nZ -> Z/mZ is the injection
  Z/mZ -> Z/nZ which takes 1 to n/m. In general, any map Z/nZ -> Z/mZ
  factors thru Z/dZ, where d = gcd(m,n), so it is a multiple of the map
  defined by 1 -> m/d; the transpose is the same multiple of 1 -> n/d.

  NOTE: this implementation works for homomorphisms between different groups
  just as well (except that one needs two grouptypes then.)
*/
void transpose(Endomorphism& e, const FiniteAbelianGroup& A)
{
  const GroupType& type = A.type(); // vector of |unsigned long|
  Endomorphism tmp = e;
  e.resize(tmp.numColumns(),tmp.numRows());

  for (size_t j = 0; j < e.numColumns(); ++j)
    for (size_t i = 0; i < e.numRows(); ++i)
    {
      unsigned long d = arithmetic::unsigned_gcd(type[i],type[j]);
      // tmp(i,j) has to be divisible by type[i]/d
      unsigned long ni = type[i]/d;
      unsigned long nj = type[j]/d;
      e(j,i) = nj*(tmp(i,j)/ni);
    }
}

} // |namespace abelian|

} // |namespace atlas|
