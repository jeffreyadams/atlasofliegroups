/*
  This is abelian.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "abelian.h"

#include <cassert>

#include "arithmetic.h"
#include "bitmap.h"
#include "arithmetic.h"
#include "matrix.h"
#include "matreduc.h"

/*****************************************************************************

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

FiniteAbelianGroup::FiniteAbelianGroup(const GroupType& t)
  : d_size(1)
  , d_type(t)
  , d_cotype(t.size())
{
  for (size_t i=0; i<t.size(); ++i)
  {
    d_size *= t[i];
    assert(i+1 == t.size() or t[i+1]%t[i]==0);
    d_cotype[i] = t.back()/t[i];
  }
}

/******** type conversions ***************************************************/

/* What is a |GrpArr|, and what is a representative weight for it?
   Why is this a method of |FiniteAbelianGroup|, without using |*this|? MvL
*/

/*
  Put in |v| a representative of |a|.

  NOTE: sloppy implementation; we don't check for overflow, which may happen
  in all cases, as the coefficients of v will be signed quantities.
*/
void FiniteAbelianGroup::toWeight(matrix::Vector<int>& v, const GrpArr& a)
  const
{
  v.assign(a.begin(),a.end()); // copy all elements; |v=a| would mismatch type
}


// Put in |v| a representative of |x|.
void FiniteAbelianGroup::toWeight(matrix::Vector<int>& v, GrpNbr x) const
{
  const GroupType& t = type();

  for (size_t i=0; i<rank(); ++i)
  {
    v[i] = x % t[i]; // extract coefficient in cyclic factor |i|
    x /= t[i];
  }
}

/******** accessors **********************************************************/


/*
  Perform |a += b|.

  Precondition: |a| and |b| hold valid arrays for the group;

  The only difficulty is doing it without triggering overflow.
*/
GrpArr& FiniteAbelianGroup::add(GrpArr& a, const GrpArr& b) const
{
  for (size_t i=0; i<a.size(); ++i)
    if (a[i] < d_type[i] - b[i])
      a[i] += b[i];
    else // overflow
      a[i] -= d_type[i] - b[i];

  return a;
}

// Perform |a -= b|.
GrpArr& FiniteAbelianGroup::subtract(GrpArr& a, const GrpArr& b) const
{
  for (size_t i=0; i<a.size(); ++i)
    if (b[i] <= a[i])
      a[i] -= b[i];
    else // underflow
      a[i] += d_type[i] - b[i];

  return a;
}

// Perform |a += x|.
GrpArr& FiniteAbelianGroup::add(GrpArr& a, GrpNbr x) const
{
  return add(a,toArray(x)); // expand in components, then add each
}


// Compute |x + b|. Note that |x| is by value, so here |add| cannot mean |+=|
GrpNbr FiniteAbelianGroup::add(GrpNbr x, const GrpArr& b) const
{
  GrpArr a=toArray(x); // we need an lvalue, so give it a name
  return toGrpNbr(add(a,b));
}


/*
  Compute |x + y|.

  Of course adding as integers does not do the right thing; decompose
*/
GrpNbr FiniteAbelianGroup::add(GrpNbr x, GrpNbr y) const
{
  return add(x,toArray(y));
}


// Return the least common multiple of the orders of the elements of the group.
unsigned long FiniteAbelianGroup::annihilator() const
{
  return d_type.size()==0 ? 1 : d_type.back();
}


/*
  Apply the matrix |q| to the element |x|, on the left.

  The idea is that the matrix defines an endomorphism of the group in terms
  of the array representation.
*/
GrpNbr FiniteAbelianGroup::leftApply(GrpNbr x, const Endomorphism& q) const
{
  GrpArr a=toArray(x);
  leftApply(a,q);

  return toGrpNbr(a);
}


/*
  Apply the matrix |q| to the array |a|, on the left.

  The idea is that the matrix defines an endomorphism of the group in terms
  of the array representation.
*/
GrpArr& FiniteAbelianGroup::leftApply(GrpArr& a, const Endomorphism& q) const
{
  GrpArr tmp = a; // need a copy, because reading/writing are not in step

  for (size_t i=0; i<q.n_rows(); ++i)
  {
    unsigned long r = 0;
    for (size_t j=0; j<q.n_columns(); ++j)
      r = arithmetic::modAdd(r,arithmetic::modProd(tmp[j],q(i,j),d_type[i]),
			     d_type[i]);
    a[i] = r; // $\sum_j q_{i,j}a_j$
  }

  return a;
}


//  Compute the order of |x| in the group.
unsigned long FiniteAbelianGroup::order(GrpNbr x) const
{
  unsigned long n = 1;
  GrpArr v=toArray(x);

  // compute least common multiple of order of each |v[i]|
  for (size_t i=0; i<d_type.size(); ++i)
    n = arithmetic::lcm(n,arithmetic::div_gcd(d_type[i],v[i]));

  return n;
}


/*
  Compute the order of |x| modulo |B|, where |B| is a subgroup
  represented by flagging the |GrpNbr| values of its elements in a bitmap

  We just keep adding |x| until we end up at an element flagged in |B|
*/
unsigned long FiniteAbelianGroup::order(const bitmap::BitMap& B,
					      GrpNbr x) const
{
  if (B.isMember(x))
    return 1;

  GrpArr a = toArray(x), b=a;
  for (unsigned long n=2; ; ++n)
    if (B.isMember(toGrpNbr(add(a,b)))) // get next multiple and test
      return n;
}


/*
  Compute the $m$ in $[0,n[$ such that $a(b) = e^{m.2i\pi/n}$, where $a$ is
  interpreted as an element of the dual group ($\Hom(G,\C^\times)$), which in
  our representation can be identified with the group itself, $b$ as an element
  of the group, and $n$ is the annihilator of the group (the last entry in
  |d_type|). This is a sort of scalar product, weighted by cotype, modulo $n$.
*/
unsigned long FiniteAbelianGroup::pairing(const GrpArr& a, const GrpArr& b)
  const
{
  unsigned long m = 0;
  unsigned long n = annihilator();

  for (size_t i=0; i<rank(); ++i)
    m = arithmetic::modAdd(m,arithmetic::modProd(a[i],b[i]*d_cotype[i],n),n);

  return m;
}


/*
  Compute the $m$ in $[0,t[$ such that $a(b) = e^{2i pi m/t}$, where $a$ is
  interpreted as an element of the dual group, $b$ as an element of the group.

  Precondition: |pairing(a,b)| is an element of t-torsion in $\Z/n$, where
  $n$ is the annihilator of the group;

  The required $m$ is just |pairing(a,b)/(n/t)|.

  NOTE: this is a sloppy implementation, that doesn't deal carefully with
  overflow. It is expected to be used only for very small groups.

  NOTE: we put in an assertion for safety.
*/
unsigned long FiniteAbelianGroup::pairing(const GrpArr& a, const GrpArr& b,
					  unsigned long t) const
{
  unsigned long n = annihilator();
  unsigned long r = pairing(a,b)*t;

  assert(r%n == 0);

  return r/n;
}


// Compute $n.x$ in the group
GrpNbr FiniteAbelianGroup::prod(GrpNbr x, unsigned long n) const
{
  GrpArr a = toArray(x);
  for (size_t i=0; i<rank(); ++i)
    a[i]=arithmetic::modProd(a[i],n,d_type[i]);

  return toGrpNbr(a);
}



}

/*****************************************************************************

        Chapter II -- The Homomorphism class

  A |Homomorphism| object represents a homomorphism from a group of type
  |d_source| to a group of type |d_dest|. The essential data is a matrix
  |d_matrix| with integral coefficients. As usual the matrix entry (i,j)
  expresses component $i$ of the image of generator $j$.

  For such a homomorphism to be well defined, the matrix entry (i,j) should be
  a multiple of |arithmetic::div_gcd(d_dest[i],d_source[j])|, since the
  maximal order the image of a generator of the cyclic subgroup $j$ of
  |d_source| can have in the cyclic subgroup $i$ of |d_dest| is
  |gcd(d_dest[i],d_source[j])|.

  Computationally we shall proceed as follows. We set $M$ to a common multiple
  of the annihilators of |d_source| and |d_dest|, and think of all cyclic
  groups as embedded in $Z/M$; this means that before acting upon a source
  element $(x_1,...,x_m)$ we should multiply each component $x_j$ by the
  cotype $q_j=M/s_j$ of that component with respect to $M$, where $s_j$ is
  |d_source[j]|. After forming $y_j=\sum_j a_{i,j} q_j x_j$ for each $i$, the
  component $y_i$ must be in the subgroup of index $t_i$ in $Z/M$, where $t_i$
  is |d_dest[i]|, in other words it must be divisible by $M/t_i$, and the
  component $i$ of the final result will be $y_i/(M/t_i)$. This definition is
  independent of the choice of $M$; we will take the lcm of the annihilators.


******************************************************************************/

namespace abelian {

/******** constructors and destructors ***************************************/


/*
  Construct the homomorphism with matrix |al|, from groups |source| to |dest|.

  Precondition: the elements of |al| are of size source.rank(), and their
  number is dest.rank(). Note that the matrix is therefore given by rows
  (linear forms on |source|), not by columns (images in |dest|); unusual!
*/
Homomorphism::Homomorphism(const std::vector<GrpArr>& al,
			   const FiniteAbelianGroup& source,
			   const FiniteAbelianGroup& dest)
  : d_source(source.type())
  , d_dest(dest.type())
  , d_cosource(d_source.size())
  , d_codest(d_dest.size())
  , d_annihilator(arithmetic::lcm(source.annihilator(),dest.annihilator()))
  , d_matrix(d_dest.size(),d_source.size())
{
  for (size_t i=0; i<d_dest.size(); ++i)
    for (size_t j=0; j<d_source.size(); ++j)
      d_matrix(i,j) = al[i][j];

  for (size_t i=0; i<d_cosource.size(); ++i)
    d_cosource[i] = d_annihilator/d_source[i]; // stride of source $i$

  for (size_t i=0; i<d_codest.size(); ++i)
    d_codest[i] = d_annihilator/d_dest[i]; // stride of destination $i$
}

/******** accessors **********************************************************/

/*
  applies the homomorphism to x according to the rules explained
  in the introduction to this section, and puts the result in dest.

  Precondition: x is in the subgroup for which this makes sense.
*/
GrpArr Homomorphism::operator*(const GrpArr& a) const
{
  GrpArr dest(d_dest.size(),0);
  for (size_t i=0; i<dest.size(); ++i)
  {
    for (size_t j=0; j<a.size(); ++j)
      dest[i] =
	arithmetic::modAdd(dest[i],
			   arithmetic::modProd(d_matrix(i,j),d_cosource[j]*a[j],
					       d_annihilator),
			   d_annihilator);

    assert(dest[i]%d_codest[i] == 0);
    dest[i] /= d_codest[i];
  }

  return dest;
}


/*
  Return $h(x)$.

  Precondition: |x| is in the subgroup for which this makes sense;

  Forwarded to the GrpArr form.
*/
GrpNbr Homomorphism::operator*(GrpNbr x) const
{
  GrpArr a; to_array(a,x,d_source);
  return to_GrpNbr(operator*(a),d_dest);
}


} // |namespace abelian|

/*****************************************************************************

        Chapter III -- Functions declared in abelian.h

******************************************************************************/

namespace abelian {

/*
  Write $A/B$ in canonical form.

  We see the current group $A$ as a quotient of $\Z^d$, where $d$ is the
  rank of $A$, and the generators of the kernel are multiples of the standard
  basis vectors given by |d_type|. The bitmap |B| specifies a subgroup by
  generators through the |GrpNbr| encoding. Then we wish to write the quotient
  $A/B$ in canonical form (i.e., as a product of cyclic groups with
  cardinalities dividing each other.) For this, we put in $b$ a scaled Smith
  normal basis for the inverse image of $B$ in $\Z^d$.

  Note that $A$ is not necessarily in canonical form, so even when $B$ is the
  trivial subgroup this might yield a basis rather different from the kernel
  basis.
*/
  void basis(std::vector<matrix::Vector<int> >& b, const bitmap::BitMap& B,
	   const FiniteAbelianGroup& A)
{
  GrpNbrList gen;
  generators(gen,B,A); // transform bitset to list

  matrix::PID_Matrix<int> M(A.rank(),A.rank()+gen.size(),0);

  // put in the kernel basis

  for (size_t i=0; i<A.rank(); ++i)
    M(i,i) = A.type()[i];

  // add generators of the subgroup
  for (size_t i=0; i<gen.size(); ++i)
  {
    matrix::Vector<int> v(A.rank());
    A.toWeight(v,gen[i]);
    M.set_column(A.rank()+i,v);
  }

// get Smith normal basis for columns span of |M|, and invariant factors

  std::vector<int> inv_factors;
  matrix::PID_Matrix<int> basis = matreduc::Smith_basis(M,inv_factors);

  // export scaled columns
  b.resize(inv_factors.size());
  for (size_t j=0; j<inv_factors.size(); ++j)
    b[j] = basis.column(j)*inv_factors[j];
}


/*
  Put in $C$ the coset $x+B$ in $A$.

  Precondition: |.capacity() == A.order()|;

  NOTE : this is a straightforward implementation, shifting elements of $B$
  individually.
*/
void coset(bitmap::BitMap& C, const bitmap::BitMap& B, GrpNbr x,
	   const FiniteAbelianGroup& A)
{
  C.reset();

  for (bitmap::BitMap::iterator it = B.begin(); it(); ++it)
  {
    GrpNbr y = *it;
    C.insert(A.add(y,x));
  }
}


/*
  Return a reference to a bitmap containing exactly one generator
  for each cyclic subgroup of |A|.

  We simply set all bits, then traverse all elements and whenever a bit is set
  we clear the bits of all elements that generate the same cyclic subgroup as
  it, namely its multiples by a number relatively prime to its order in |A|.
  This uses that iterators over a |bitmap::BitMap| adapt to changes to the
  underlying bitmap during traversal, unlike |bitset::BitSet::iterator|s.

  NOTE: it is expected that this function will typically be called repeatedly
  for the same group. The bitmap is constructed on the first call for the
  given group, following the principle of lazy evaluation.
*/
const bitmap::BitMap& cycGenerators(const FiniteAbelianGroup& A)
{
  static FiniteAbelianGroup Astat;
  static bitmap::BitMap cyc;

  if (Astat.type() != A.type())  // then we need to update cyc
  {
    cyc.set_capacity(A.order());
    cyc.fill();

    for (bitmap::BitMap::iterator it = ++cyc.begin(); it(); ++it)
    {
      GrpNbr x = *it;
      unsigned long n = A.order(x);
      for (unsigned long i=2; i<n; ++i)
	if (arithmetic::unsigned_gcd(n,i) == 1) // i is prime relative to n
	{
	  GrpNbr xj = A.prod(x,i);
	  cyc.remove(xj);
	}
    }

    Astat = A;

  }

  return cyc;
}


/*
  Transform |B| into the subgroup generated by |B| and |x| in |A|.

  NOTE : this is a simple-minded implementation; we do not aim for speed.
*/
void generateSubgroup(bitmap::BitMap& B, GrpNbr x, const FiniteAbelianGroup& A)
{
  unsigned long n = A.order(x);
  bitmap::BitMap b(B);
  bitmap::BitMap c(A.order());

  for (unsigned long i=1; i<n; ++i)
  {
    GrpNbr xj = A.prod(x,i);
    coset(c,B,xj,A);
    b |= c;
  }

  B.swap(b);
}


/*
  puts in gen a list of generators of the subgroup B.
*/
void generators(GrpNbrList& gen, const bitmap::BitMap& B,
		const FiniteAbelianGroup& A)
{
  gen.assign(B.begin(),B.end());
}


// Whether |c| is all twos.
bool isElementaryAbelian(const std::vector<unsigned long>& c)
{
  for (size_t i=0; i<c.size(); ++i)
    if (c[i] != 2)
      return false;

  return true;
}


/*
  Return a list of representatives of the cosets modulo |B|.

  Precondition: |qr.capacity() = A.order()|;

  Algorithm: fill up |qr|; traverse |qr|, and for each |x| that is reached,
  remove from |qr| the elements of the coset of |x| other than |x|.

  The efficiency depends on the efficiency of coset computations.
*/
bitmap::BitMap quotReps(const bitmap::BitMap& B, const FiniteAbelianGroup& A)
{
  bitmap::BitMap qr(A.order());
  qr.fill();

  for (bitmap::BitMap::iterator it = qr.begin(); it(); ++it)
  {
    bitmap::BitMap c(A.order());
    coset(c,B,*it,A);
    c.remove(*it);
    qr.andnot(c);
  }
  return qr;
}


/*
  Put the array-form of |x| into |a|.

  Precondition: |x| is representable for |t| (|x < prod t[i]|); otherwise the
  result is the expression of |x| modulo that product.

  Explanation: the array-form (relative to type) is the unique expression
  of |x| in the variable-radix base defined by type.
*/
void to_array(GrpArr& a, GrpNbr x, const GroupType& t)
{
  for (size_t i=0; i<a.size(); ++i)
  {
    a[i] = x % t[i];
    x /= t[i];
  }
}


/*
  Reduce |v| modulo |d_type|.

  The only difficulty is to make sure that negative values are reduced in
  the way that we want, which |arithmetic::remainder| does for us>

  Precondition: |a| is set to |v.size()|;
*/
void to_array(GrpArr& a, const matrix::Vector<int>& v, const GroupType& t)
{
  for (size_t i=0; i<v.size(); ++i)
    a[i] = arithmetic::remainder<long int>(v[i],t[i]);
}


/*
  Transform |q| into the corresponding |Endomorphism|.

  This just involves rewriting the coefficients as unsigned longs modulo
  the type factors.
*/
void toEndomorphism(Endomorphism& e, const matrix::PID_Matrix<int>& q,
		    const FiniteAbelianGroup& A)
{
  Endomorphism(q.n_rows(),q.n_columns()).swap(e);

  for (size_t j=0; j<q.n_columns(); ++j)
    for (size_t i=0; i<q.n_rows(); ++i)
      e(i,j) = arithmetic::remainder<long int>(q(i,j),A.type()[i]);
}


/*
  Return the number-form of |a|.

  Precondition: |a| is representable as a |GrpNbr|;
*/
GrpNbr to_GrpNbr(const GrpArr& a, const GroupType& t)
{
  GrpNbr x = 0;

  for (size_t i = a.size(); i-->0;)
    (x *= t[i]) += a[i];

  return x;
}


/*
  Transpose the endomorphism |e|.

  This is more subtle than one might think; the main point is that when
  $m\div n$, the transpose of the canonical map $\Z/n\Z to \Z/m\Z$ is the
  injection $\Z/m\Z to \Z/n\Z$ which takes 1 to $n/m$. In general, any map
  $\Z/n\Z \to \Z/m\Z$ factors through $\Z/d\Z$, where $d =gcd(m,n)$, so it is
  a multiple of the map defined by $1 \to m/d$; the transpose is the same
  multiple of $1 \to n/d$.

  NOTE: this implementation works for homomorphisms between different groups
  just as well (except that one needs two group types then.)
*/
void transpose(Endomorphism& e, const FiniteAbelianGroup& A)
{
  const GroupType& type = A.type(); // vector of |unsigned long|
  Endomorphism tr(e.n_columns(),e.n_rows());

  for (size_t i=0; i<e.n_rows(); ++i)
    for (size_t j=0; j<e.n_columns(); ++j)
    {
      unsigned long d = arithmetic::unsigned_gcd(type[i],type[j]);
      // e(i,j) has to be divisible by type[i]/d
      unsigned long ni = type[i]/d;
      unsigned long nj = type[j]/d;
      tr(j,i) = nj*(e(i,j)/ni);
    }

  e.swap(tr);
}

} // |namespace abelian|

} // |namespace atlas|
