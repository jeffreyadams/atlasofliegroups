/*!
\file
\brief Declarations and definitions for templates for BitVector.

This is a vector in the first d_size coordinates of the vector space
(Z/2Z)^dim over the two-element field.  The software envisions dim
between 0 and four times the machine word length (precisely, four
times the constant longBits, which is the number of bits in an
unsigned long integer).

*/
/*
  This is bitvector.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2008,2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <algorithm>
#include <cassert>

#include "bitvector.h"

#include "comparison.h"
#include "constants.h"
#include "bitset.h"
#include "matrix.h"

/*****************************************************************************

        Chapter I -- The BitVector class

******************************************************************************/

namespace atlas {

namespace bitvector {

template<size_t dim>
BitVector<dim>::BitVector(const matrix::Vector<int>& v) // reduce mod 2
  : d_data()
  , d_size(v.size())
{
  assert(d_size<=dim);
  for (size_t j = 0; j < d_size; ++j)
    d_data.set(j, v[j]%2!=0 ); // WARNING: not |v[j]%2==1|, think |v[j]<0| !
}


/*!
\brief Adds b to the bitvector (as the last coordinate), increasing the
  size by one.

  It is the user's responsibility to make sure that the |size| does not exceed
  |dim|. This is typically applied with |d_size<=RANK_MAX| initially, and
  |dim>=RANK_MAX+1|.
*/
template<size_t dim> BitVector<dim>& BitVector<dim>::pushBack(bool b)

{
  assert(d_size<dim); // inequality must be strict here
  d_data.set(d_size,b);

  ++d_size;

  return *this;
}

/*! \brief Extracts the bits flagged by |t|, and packs their value into
  consecutive positions; resets the size to the number of bits so packed.

  This is value is useful to transform an element known to lie in a subspace
  to the canonical basis of the subspace; that notion is so defined as to make
  this base change just a selection of coordinate values at increasing
  positions. This is so because the canonical basis $b$ comes with a set of
  coordinate positions $j$ such that only one $b_i$ has $b_i[j]==1; for a
  linear combination of basis vectors to be equal to $v$, the coefficient of
  $b_i$ must then be $v[j]$.

  Note that |slice| changes (weakly decreases) the size of its |BitVector|.
*/
template<size_t dim>
void BitVector<dim>::slice(const bitset::BitSet<dim>& t)

{
  size_t c = 0;

  // in the following loop |c<=j|, so we can inspect bit |j| and then set |c|
  for (size_t j = 0; j < size(); ++j)
    if (t[j])
    { // insert bit j at position c
      d_data.set(c,d_data[j]);
      ++c;
    }
  /* actually we could have used a bitset iterator over |t| */

  d_data&=bitset::BitSet<dim>(constants::lMask[c]); // clear remainder

  d_size = c; // new size equals (counted) number of sliced out bits
}

/*! \brief Undoes the effect of |slice|, inserting zero bits where needed.

  This cannot be used to undo the use of |slice| to express a vector on a
  subspace basis (you need to form a linear combination of the basis for
  that), but when |slice| is used during a projection modulo a subspace (as
  happens in subquotients), then |unslice| can reconstruct a preimage.

  Note that |unslice| changes (weakly increases) the size of its |BitVector|.
*/
template<size_t dim>
void BitVector<dim>::unslice(bitset::BitSet<dim> t,size_t new_size)
{
  bitset::BitSet<dim> result;
  for (typename bitset::BitSet<dim>::iterator it=t.begin(); it();
       d_data>>=1,++it)
    result.set(*it,d_data[0]);
  d_data=result;
  d_size=new_size;
}

} // namespace bitvector

/*****************************************************************************

        Chapter II -- The BitMatrix class

******************************************************************************/

namespace bitvector {

/*!
  Constructs the matrix whose columns are given by the vectors in |b|.

  NOTE : it is assumed that all the vectors in |b| have the same size.
*/
template<size_t dim>
BitMatrix<dim>::BitMatrix(const std::vector<BitVector<dim> >& b,
			  unsigned short int num_rows)  // size of columns
  : d_data() // start out empty
  , d_rows(num_rows)
  , d_columns(b.size())
{
  assert(d_rows<=dim);

  d_data.reserve(b.size());

  // this loop is necessary to reduce each column to a |bitset::BitSet<dim>|
  for (size_t j = 0; j < b.size(); ++j)
  {
    assert(b[j].size()==num_rows);
    d_data.push_back(b[j].data());
  }

}

template<size_t dim>
BitMatrix<dim>::BitMatrix(const matrix::Matrix<int>& m) // set modulo 2
  : d_data(m.numColumns(),bitset::BitSet<dim>())
  , d_rows(m.numRows())
  , d_columns(m.numColumns())
{
  assert(m.numRows()<=dim);
  for (size_t i=0; i<d_rows; ++i)
    for (size_t j=0; j<d_columns; ++j)
      d_data[j].set(i, (m(i,j)&1)!=0 );
}

/******** accessors *********************************************************/

/*!
  \brief Applies our |BitMatrix| to |source| and returns the result

  It is assumed that |d_columns| is equal to |source.size()|. The result size
  will be set from |d_rows|.
*/
template<size_t dim>
BitVector<dim> BitMatrix<dim>::operator*(const BitVector<dim>& source) const
{
  assert(d_columns==source.size());

  BitVector<dim> result(combination(d_data,source.data()),d_rows);
  return result;
}

/*!
  \brief A pipe-dream version of apply.

  We assume that |I| is an InputIterator with value-type |BitVector<dim>|, and
  |O| an OutputIterator with the same value-type. Then we apply our matrix to
  each vector in [first,last[ and output it to out.
*/
template<size_t dim> template<typename I, typename O>
void BitMatrix<dim>::apply(const I& first, const I& last, O out) const
{
  for (I i = first; i != last; ++i)
  {
    *out = apply(*i);
    ++out;
  }
}

/*!
  Puts in r the i-th row of the matrix.
*/

template<size_t dim> void BitMatrix<dim>::get_row(BitVector<dim>& r, size_t i)
  const
{
  assert(d_columns<=dim);
  r.resize(d_columns);

  for (size_t j = 0; j < d_columns; ++j)
    r.set(j,test(i,j));
}


/*!
  \brief Puts in |b| a basis of the image of the matrix.
*/
template<size_t dim>
BitVectorList<dim> BitMatrix<dim>::image() const
{
  BitVectorList<dim> b;
  std::vector<size_t> f; // auxialiry to record leading bit positions

  for (size_t j = 0; j < d_columns; ++j)
    spanAdd(b,f,column(j));  // add column |j| to span of |b|

  return b;
}


/*! \brief Puts in |b| an echelon basis for the kernel of the matrix (but not
  the canonical basis of that subspace; it is echelon from high to low bit
  positions).

  This is based on making the transpose matrix echelon, and then solving the
  resulting trivialised system of equations.
*/
template<size_t dim>
void BitMatrix<dim>::kernel(std::vector<BitVector<dim> >& b) const
{
  b.clear();

  if (isEmpty())
    return; // no unknowns, no nontrivial solutions

  assert(d_columns<=dim);
  std::vector<BitVector<dim> > eqn(d_rows);

  // get rows of the matrix into |eqn|

  for (size_t i = 0; i < d_rows; ++i)
    get_row(eqn[i],i);

  // normalize |eqn|

  bitset::BitSet<dim> t; // will flag a subset of [0,d_columns[

  normalize(t,eqn);

/* Now the coordinates in the positions not flagged by |t| are free parameters
   of the kernel, and each element of |eqn| allows one of the remaining
   coordinates to be expressed in terms of them. As for generators of the
   kernel, there is one such generator $v$ for each position $j$ not flagged
   by |t|, and apart from $v[j]=1$ it can only have nonzero bits on positions
   flagged by |t|; the bits at those positions are given by coordinates $j$
   taken from the succesive elements of the normalised basis |eqn|.

   Note that this is already a normal basis for the complement of |t| as
   explained below under |normalSpanAdd|, but that complement not lex-minimal.
*/


  for (size_t j = 0; j < d_columns; ++j)
    if (not t[j]) {
      BitVector<dim> v(d_columns,j); // start with unit vector $e_j$
      size_t c = 0;
      for (typename bitset::BitSet<dim>::iterator it=t.begin(); it(); ++it,++c)
	v.set(*it,eqn[c][j]); // use |eqn[c]| to find |v[*it]|
      b.push_back(v);
    }
}

/******** manipulators *******************************************************/


/*!
  \brief Increment *this by adding m.

  Precondition: m has the same size as the current matrix.
*/
template<size_t dim>
BitMatrix<dim>& BitMatrix<dim>::operator+= (const BitMatrix<dim>& m)
{
  assert(d_rows==m.d_rows);
  assert(d_columns==m.d_columns);
  for (size_t j = 0; j < d_columns; ++j)
    d_data[j] ^= m.d_data[j]; // not |+|, since these are |BitSet|s

  return *this;
}


/*!
  \brief Right multiply our BitMatrix by |m|.

  As in apply, this can be done rather efficently by computing a whole column
  at a time.

  NOTE : of course |m.d_rows| must be equal to |d_columns|.
*/
template<size_t dim>
BitMatrix<dim>& BitMatrix<dim>::operator*= (const BitMatrix<dim>& m)
{
  assert(d_columns==m.d_rows);
  BitMatrix<dim> res(d_rows,m.d_columns);

  for (size_t j = 0; j < m.d_columns; ++j) {
    /* set column |j| of result as linear combination of columns of left
       factor |*this| determined by column |j| of right factor |m| */
    res.d_data[j]=combination(d_data,m.d_data[j]);
  }

  swap(res); // put result in |*this|

  return *this;
}


/*!
  \brief Replaces the current matrix by its inverse.

  [This code is untested by me, and may be removed if no use is found. MvL]

  It is the caller's responsibility to make sure that m is in fact
  invertible (in particular, that it is square).  If necessary, this
  may be done by a call to isInvertible().

  For the algorithm, we use the normalSpanAdd function. We start out with
  a matrix of size (2r,c) (if r,c is the size of our original matrix)
  containing the identity matrix below the given one. Then we apply
  normalSpanAdd; at the end, the upper part of our matrix is (some permutattion
  of) the identity. Setting it right will finish the job.

  Actually, since this would require working with BitMatrix<2*dim>, and since
  conversion functions don't seem to be forthcoming, we work with a pair
  of matrices, and just copy the code from normalSpanAdd.
*/
template<size_t dim> BitMatrix<dim>& BitMatrix<dim>::invert()
{
  BitMatrix<dim> i(d_rows);

  for (size_t j = 0; j < d_rows; ++j)
    i.set(j,j);

  std::vector<size_t> f;

  for (size_t k = 0; k < d_rows; ++k) {// add column k

    for (size_t j = 0; j < k; ++j)
      if (test(k,f[j])) { // set bit f[j] of data[k] to zero
	d_data[k] ^= d_data[j];
	i.d_data[k] ^= i.d_data[j];
      }

    // adjust the basis a

    size_t n = d_data[k].firstBit();

    for (size_t j = 0; j < k; ++j)
      if (d_data[j][n]) {
	d_data[j] ^= d_data[k];
	i.d_data[j] ^= i.d_data[k];
      }

    f.push_back(n);
  }

  // write the appropriate column-permutation of i in the current matrix

  for (size_t j = 0; j < d_rows; ++j)
    d_data[f[j]] = i.d_data[j];

  return *this;
}


//! \brief Resets the matrix to zero.
template<size_t dim> void BitMatrix<dim>::reset()
{
  for (unsigned long j = 0; j < d_data.size(); ++j)
    d_data[j].reset();
}

//! \brief Resizes the matrix to |m| rows, |n| columns, leaving data around
template<size_t dim> void BitMatrix<dim>::resize(size_t m, size_t n)
{
  assert(m<=dim);
  d_data.resize(n);

  d_rows = m;
  d_columns = n;
}

template<size_t dim> void BitMatrix<dim>::swap(BitMatrix<dim>& m)
{
  d_data.swap(m.d_data);
  std::swap(d_rows,m.d_rows);
  std::swap(d_columns,m.d_columns);
}

//! \brief Transposes the matrix.
template<size_t dim> BitMatrix<dim>& BitMatrix<dim>::transpose()
{
  BitMatrix<dim> result(d_columns,d_rows);

  for (size_t i = 0; i < d_rows; ++i)
    for (size_t j = 0; j < d_columns; ++j)
      result.set(j,i,test(i,j));

  swap(result);
  return *this;
}

}

/*****************************************************************************

        Chapter III -- Functions defined in bitvector.h

******************************************************************************/

namespace bitvector {

/*!
  \brief Puts in |v| the linear combination of the elements of |b| given by
  |e|.

  NOTE : it is the caller's responsibility to check that |v| already has the
  correct size. (The right size cannot be determined here if |b.size()=0|.)
*/
template<size_t dim>
  BitVector<dim> combination
   (const std::vector<BitVector<dim> >& b,
    size_t n,
    const bitset::BitSet<dim>& e)
{
  BitVector<dim> result(n);

  for (size_t i = 0; i<b.size(); ++i)
    if (e[i])
      result += b[i];

  return result;
}

/*!
  \brief Returns the linear combination of the elements of |b| given by |coef|.

  Contrary to the previous case, there is no notion of size for |v|, and
  there is no need for any particular preparation.
*/
template<size_t dim>
  bitset::BitSet<dim> combination(const std::vector<bitset::BitSet<dim> >& b,
				  const bitset::BitSet<dim>& coef)
  {
    bitset::BitSet<dim> result(0);
    for (size_t i = 0; i < b.size(); ++i)
      if (coef[i])
	result ^= b[i]; // not |+| here, these are |BitSet|s.
    return result;
  }

/*!
  \brief Put into |c| a solution of the system with as left hand sides the
  rows of a matrix whose columns are given by |b|, and as right hand sides the
  bits of |rhs|, and return |true|; if no solution exists just return |false|.
*/
template<size_t dim>
  bool firstSolution(bitset::BitSet<dim>& c,
		     const std::vector<BitVector<dim> >& b,
		     const BitVector<dim>& rhs)


{
  if (b.size() == 0) // then there are no unknowns to solve
  {
    if (rhs.isZero()) {
      c.reset(); // clear any previously set bits
      return true; // list of 0 "solved" unknowns
    }
    else
      return false; // no unknowns, no solutions
  }

  size_t n = b[0].size();
  assert(n==rhs.size());

  std::vector<bitset::BitSet<dim> > a; // list of normalised eqns, lhs part
  bitset::BitSet<dim> rh;      // corresponding right hand sides
  std::vector<size_t> f;       // list indicating "pivot" positions in |a|

  for (size_t i = 0; i < n; ++i) {
    bitset::BitSet<dim> r;
    // set r to i-th row of matrix whose columns are the |b[j]|
    for (size_t j = 0; j < b.size(); ++j)
      r.set(j,b[j][i]);
    bool x = rhs[i]; // now $(r,x)$ is one of the equations to solve

    // normalize |r| with respect to |a|: clear coefficients at previous pivots
    for (size_t j = 0; j < f.size(); ++j)
      if (r[f[j]]) {
	r ^= a[j];
	x ^= rh[j];
      }

    // if |r| is independent, normalize the elements of |a| and add |r| to it
    if (r.any()) // there are nonzero coefficients in the remaining equation
    {
      size_t m = r.firstBit(); // pivot position; this equation solves bit |m|
      f.push_back(m);          // record pivot

      // update previous equations, clearing their coefficient at position |m|
      for (size_t j = 0; j < a.size(); ++j)
	if (a[j][m]) {
	  a[j] ^= r;
	  rh.set(j,rh[j]^x);
	}

      // add new equation (could precede loop above, but it changes |a.size()|)
      rh.set(a.size(),x);
      a.push_back(r);
    }
    else // trivial equation
      if (x)
	return false;  // if right hand side non-null system is contradictary
      else {}          // otherwise just drop the null equation
  }

  // now we know the system is solvable

  // the solution has values of |rh| at positions indicated by the pivots

  c.reset(); // clear all bits, then define those at positions |f[j]|
  for (size_t j = 0; j < f.size(); ++j)
    c.set(f[j],rh[j]);

  return true;
}

/*!
  \brief Either find a solution of the system of equations |eqn|, putting it
  into |sol| and returning |true|, or return |false| if no solition exists.

  Here |eqn| holds a system of equations, the last bit of each being
  interpreted as the right hand side.

  However, to solve it we may introduce an extra indeterminate for this final
  position, solve the homogenous system (which is like finding the kernel of
  the transpose matrix), and look for a kernel element with final coordinate
  equal to $-1$ (which working over $Z/2Z$ is of course the same as $1$).
*/
template<size_t dimsol, size_t dimeq>
bool firstSolution(BitVector<dimsol>& sol,
		   const std::vector<BitVector<dimeq> >& eqns)


{
  std::vector<BitVector<dimeq> > eqn(eqns); // local copy
  bitset::BitSet<dimeq> t;

  normalize(t,eqn); // normalize matrix, treating the last bit like all others

  /* now the system is contradictory if and only if the last bit in |t| is
     set, since that means some equation sets the corresponding indeterminate
     to 0 (because it has its \emph{leading} bit at the final position), while
     absence of any equation for that indeterminate means we can make it $-1$.
   */
  if (eqn.size()>0 and t[eqn[0].size()-1])
    return false;

  sol.reset(); // only now can we clear the solution

  if (eqn.size() == 0) // do nothing, then zero solution solves the null system
    return true; // note that we leave |sol| at its original size (only) here

  sol.resize(eqn[0].size()-1); // otherwise |sol| is one shorter than equations

  // we set the bits flagged by |t| to the corresponding right hand side

  size_t c = 0;
  const size_t rhs = sol.size();

  for (size_t j = 0; j < rhs; ++j)
    if (t[j]) {
      sol.set(j,eqn[c][rhs]);
      ++c;
    }

  return true;
}


/*!
  \brief Puts in m the identity matrix in rank n.

  Precondition: n <= dim;
*/
template<size_t dim> void identityMatrix(BitMatrix<dim>& m, size_t n)
{
  m.resize(n,n);
  m.reset();

  for (size_t j = 0; j < n; ++j)
    m.set(j,j);
}


/*!
  \brief Initializes b to the canonical basis in dimension n.
*/
template<size_t dim> void initBasis(std::vector<BitVector<dim> >& b, size_t n)
{
  assert(n<=dim);
  b.assign(n,BitVector<dim>(n)); // set to |n| null vectors of size |n|

  for (size_t j = 0; j < n; ++j)
    b[j].set(j);
}

/*
   This auxiliary unary function class is used as comparison object, to
   partially sort bit vectors by the position of their leading bit
*/
template<size_t dim> struct FirstBit
{
  typedef const BitVector<dim>& argument_type;
  typedef size_t result_type;

  result_type operator() (argument_type v) const { return v.firstBit(); }
};

/*!
  \brief Replaces |b| by the ordered canonical basis of the vector space $V$
  it spans. Flags in |t| the set of coordinate positions associated to |b|.

  What is flagged in |t| is the set $J$ in described in the comment for
  |normalSpanAdd| below. For any $j$, only |b[j]| has a nonzero bit at the
  position $j'$ of set bit number |j| of |t|; consequently, for any \f$v\in
  V\f$, the coordinate of |b[j]| in $v$ is $v[j']$.

  This function works essentially be repeatedly calling |normalSpanAdd| for
  the vectors of |b|, replacing |b| by the resulting canonical basis at the
  end. However, the selected coordiante positions |f| do not come out
  increasingly this way, so we have to sort the canonical basis by leading bit
  position. This amounts to setting $a'[k]=a[p(k)]$ where
  \f$p(0)\ldots,p(l-1)\f$ is the result of sorting \f$f[0]\ldots,f[l-1]\f$
  with $l=f.size()=a.size()$.
*/
template<size_t dim>
  void normalize(bitset::BitSet<dim>& t, std::vector<BitVector<dim> >& b)

{
  std::vector<BitVector<dim> > a;
  std::vector<size_t> f;

  for (size_t j = 0; j < b.size(); ++j)
    normalSpanAdd(a,f,b[j]);

  // convert |f| to a |BitSet|
  t.reset();
  for (size_t j = 0; j < f.size(); ++j)
    t.set(f[j]);

  // reorder the basis elements according to their first bit
  std::sort(a.begin(),a.end(),comparison::compare(FirstBit<dim>()));

  // commit
  b.swap(a);
}


/*!
  \brief Transforms the normal basis defined by the unordered list $a$ into
  one for the span of $a$ and $v$

  Also updates the list |f| of the same length $l$ as |a| such that
  |a[i][f[j]]==(i==j?1:0)| for all $i,j<l$.

  For each subvectorspace $V$ of $k^d$, let $I$ be a subset of
  \f$\{0,\ldots,d-1\}\f$ such that the standard basis vectors $e_i$ for
  \f$i\in I\f$ generate a complementary subspace $e_I$ to $V$ (one can find
  such an $I$ by repeatedly throwing in $e_i$s linearly independent to $V$ and
  previously chosen ones). The normal basis of $V$ corresponding to $I$ is
  obtained by projecting the $e_j$ for $j$ in the complement $J$ of $I$ onto
  $V$ along $e_I$ (i.e., according to the direct sum decompostion
  \f$k^d=V\oplus e_I\f$). This can be visualised by viewing $V$ as the
  function-graph of a linear map from $k^J$ to $k^I$; then the normal basis is
  the lift to $V$ of the standard basis of $k^J$. We define the canonical
  basis of $V$ be the normal basis for the complement $I$ of the
  lexicographically minimal possible set $J$ (lexicographic for the increasing
  sequences representing the subsets; in fact $I$ is lexicographically maximal
  since complementation reverses this ordering one fixed-size subsets). One
  can find this $J$ by repeatedly choosing the smallest index such that the
  projection from $V$ defined by extracting the coordinates at the selected
  indices remains surjective.

  This function assumes that $a$ already contains the canonical basis of some
  subspace, and that the elements of |f| describe the corresponding set $J$.
  We add |v| to the subspace and extend |f|. Then |a| nor |f| are modified if
  |v| lies in the subspace generated by |a|. Otherwise a new element is added
  to |a|, a new coordinate index |n| is added to |f|, and the existing vectors
  in |a| are modified to clear their coordinate |n|.
*/

template<size_t dim>
  void normalSpanAdd(std::vector<BitVector<dim> >& a, std::vector<size_t>& f,
		     const BitVector<dim>& v)

{
  assert(a.size()==0 or a[0].size()==v.size());
  assert(a.size()==f.size());
  if (v.isZero()) // |v| is the zero vector do nothing (test needed?)
    return;

  // reduce |v| modulo |a|

  BitVector<dim> w = v;

  // substract canonical projection of |v| onto span of |a|
  for (size_t j = 0; j < a.size(); ++j)
    if (w[f[j]]) // if coordinate $v[f[j]]$ is $1$, subtract |a[j]|
      w -= a[j];

  if (w.isZero()) // that is, if |v| is in span of |a|
    return;

  // now |w| will be new basis vector; it already has its bits at $J$ cleared

  // determine coordinate index associated to |w|
  size_t n = w.firstBit();

  // clear bit |n| in previous basis vectors
  for (size_t j = 0; j < a.size(); ++j)
    if (a[j][n])
      a[j] -= w;

  f.push_back(n);
  a.push_back(w);
}


/*!
  \brief Enlarges the basis a to span v.

  This is a simplified version of |normalSpanAdd|

  It is assumed that a contains a list of independent bitvectors all of size
  |v.size()|; then a reduction of |v| modulo the vectors of |a| is added to
  the list if |v| was independent, and if it was dependent nothing happens.

  Here we still assume that the first set bits of the elements in |a| are all
  distinct, their positions are indicated in |f|, and bit number |f[i]| is
  cleared in |a[j]| whenever |i<j|; under these conditions reduction
  modulo~|a| can be performed by subtracting, for all |i| in increasing order,
  the vector |a[i]| if bit |f[i]| is currently set. We do not however assume
  that bit number |f[i]| is cleared in |a[j]| for all |j<i|, and as a
  consequence we do not need to modify previous vectors |a[i]| to maintain the
  condition for calling |spanAdd| again; this is the (only) difference with
  |normalSpanAdd|.
*/
template<size_t dim> void spanAdd(std::vector<BitVector<dim> >& a,
				  std::vector<size_t>& f,
				  const BitVector<dim>& v)
{
  assert(a.size()==0 or a[0].size()==v.size());
  assert(a.size()==f.size());
  if (v.isZero()) // v is the zero vector
    return;

  // reduce |v| modulo |a|

  BitVector<dim> w = v;

  for (unsigned long j = 0; j < a.size(); ++j)
    if (w[f[j]]) // if coordinate $v[f[j]]$ is $1$, subtract |a[j]|
      w -= a[j];

  if (w.isZero()) // that is, if |v| is in span of |a|
    return;

  f.push_back(w.firstBit());
  a.push_back(w);
}

/* functions never called have been grouped here, MvL */
#if 0
template<size_t dim>
  void complement(bitset::BitSet<dim>& c,
		  const std::vector<BitVector<dim> >& b,
		  size_t d)

/*!
  \brief Flags into |c| a subset of the standard basis that spans a
  complementary subspace to the subspace spanned by |b|.

  It is assumed that the vectors in b are all of the same size, but not
  necessarily independent.

  NOTE : we need to pass the dimension in case |b| is empty; we don't want
  to set bits beyond that (so that |c.count()|, for instance, yields the
  correct dimension of the complement.)
*/

{
  std::vector<size_t> f;
  std::vector<BitVector<dim> > a;

  for (unsigned long j = 0; j < b.size(); ++j)
    spanAdd(a,f,b[j]);

  // now f contains the indices we _don't_ want

  for (size_t j = 0; j < d; ++j)
    c.set(j);

  for (size_t j = 0; j < f.size(); ++j)
    c.reset(f[j]);
}

template<size_t dim> bool isIndependent(const std::vector<BitVector<dim> >& b)

/*!
  \brief Tells whether the system of bitvectors is independent.
*/

{
  std::vector<BitVector<dim> > a;
  std::vector<size_t> f;

  for (size_t j = 0; j < b.size(); ++j) {
    spanAdd(a,f,b[j]);
    if (a.size() == j) // first time a dependent vector is found
      return false;
  }

  return true;
}


/*!
  \brief Puts in p the matrix of the projection on the canonical
  complement to the span of b.

  We assume that b holds the normal basis for the
  subspace V that it spans (if not, this can be obtained by a call to
  normalize.) This function then puts in p the matrix of the projection
  to the canonical complement of V, parallel to V.

  NOTE : we need to pass the dimension in case b is empty.
*/
template<size_t dim>
  void projection(BitMatrix<dim>& p, const std::vector<BitVector<dim> >& b,
		  size_t d)
{
  // flag the indices for the canonical complement

  bitset::BitSet<dim> c;
  std::vector<size_t> f;

  bitset::set(c,d);

  for (size_t j = 0; j < b.size(); ++j) {
    size_t fb = b[j].firstBit();
    f.push_back(fb);
    c.reset(f[j]);
  }

  // the coordinates of the basis vectors corresponding to indices in s
  // give the projection matrix

  p.resize(d-b.size(),d);
  p.reset();

  size_t k = 0;

  for (size_t j = 0; j < d; ++j)
    if (c[j]) {
      p.set(k,j);
      ++k;
    }

  for (size_t j = 0; j < b.size(); ++j) {
    BitVector<dim> v = b[j];
    v.slice(c);
    const bitset::BitSet<dim>& vd = v.data();
    p.setColumn(f[j],vd);
  }
}

// this function does not appear to be used anywhere [MvL]
template<size_t dim>
  void reflectionMatrix(BitMatrix<dim>& m, const BitVector<dim>& a,
			const BitVector<dim>& a_check)

/*!
  \brief Puts in m the matrix of the reflection defined by a and a_check.

  Precondition: a and a_check have same size; <a,a_check> = 0;

  This is the operator x -> x + <x,a_check>a (we can write + because we are
  in characteristic two.)
*/

{
  identityMatrix(m,a.size());

  for (size_t j = 0; j < a.size(); ++j)
    if (a_check[j])
      m.addToColumn(j,a);
}

// this function does not appear to be used anywhere [MvL]
template<size_t dim>
  void relations(std::vector<BitVector<dim> >& rel,
		 const std::vector<BitVector<dim> >& b)

/*!
  \brief Writes in r the relations among the elements in b.

  In other words, it solves the system of equations defined by the _rows_ in
  the matrix whose columns are given by b.
*/

{
  rel.resize(0);

  if (b.size() == 0) // do nothing
    return;

  size_t r = b.size();       // dimension of the source space
  size_t d = b[0].size();    // dimension of the target space
  std::vector<BitVector<dim> > eqn(d);

  // "transpose" b in e

  for (size_t i = 0; i < d; ++i) {
    eqn[i].resize(b.size());
    for (size_t j = 0; j < b.size(); ++j)
      if (b[j][i])
	eqn[i].set(j);
  }

  // normalize e

  bitset::BitSet<dim> t; // will flag a subset of [0,r[
  normalize(t,eqn);

  // now we get a relation for each j in [0,d[ not flagged by t
  // the relations are of the form e_j = sum a_{i,j}e_i, j not in t,
  // i in t, where the a_{i,j} are given by the columns in eqn

  for (size_t j = 0; j < r; ++j)
    if (not t[j]) {
      BitVector<dim> v(r,j);
      size_t c;
      for (size_t i = 0; i < r; ++i)
	if (t[i]) {
	  if (eqn[c][j])
	    v.set(i);
	  ++c;
	}
      rel.push_back(v);
    }
}
#endif // end of unused functions

template
  BitVector<constants::RANK_MAX> combination
   (const std::vector<BitVector<constants::RANK_MAX> >& b,
    size_t n,
    const bitset::BitSet<constants::RANK_MAX>& e);

template
  bitset::BitSet<constants::RANK_MAX> combination
  (const std::vector<bitset::BitSet<constants::RANK_MAX> >&,
   const bitset::BitSet<constants::RANK_MAX>&);

template
  bool firstSolution(bitset::BitSet<constants::RANK_MAX>& c,
		     const std::vector<BitVector<constants::RANK_MAX> >& b,
		     const BitVector<constants::RANK_MAX>& rhs);
template
  bool firstSolution
   (BitVector<constants::RANK_MAX>& sol,
    const std::vector<BitVector<constants::RANK_MAX+1> >& eqns);

template void identityMatrix(BitMatrix<constants::RANK_MAX>&, size_t);
template void initBasis(std::vector<BitVector<constants::RANK_MAX> >&, size_t);

template class BitVector<constants::RANK_MAX>;
template class BitVector<constants::RANK_MAX+1>;//|BinaryEquation|
template class BitMatrix<constants::RANK_MAX>;

} // |namespace bitvector|

} // |namespace atlas|
