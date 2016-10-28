/*
  This is bitvector.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2008-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/* Declarations and definitions for templates for BitVector.

  A vector in the first d_size coordinates of the vector space $(Z/2Z)^dim$
*/

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "bitvector.h"

#include "../Atlas.h"

#include "comparison.h"
#include "constants.h"
#include "bitset.h"
#include "matrix.h"

/*****************************************************************************

        Chapter I -- The BitVector class

******************************************************************************/

namespace atlas {

namespace bitvector {

template<size_t dim> template<typename C>
BitVector<dim>::BitVector(const matrix::Vector<C>& v) // reduce mod 2
  : d_data()
  , d_size(v.size())
{
  assert(d_size<=dim);
  for (size_t j = 0; j < d_size; ++j)
    d_data.set(j, v[j]%2!=0 ); // WARNING: not |v[j]%2==1|, think |v[j]<0| !
}


/*
  Add |b| to the bitvector (as the last coordinate), increasing the
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

/* Extracts the bits flagged by |t|, and packs their value into
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
void BitVector<dim>::slice(const BitSet<dim>& t)

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

  d_data&=BitSet<dim>(constants::lMask[c]); // clear remainder

  d_size = c; // new size equals (counted) number of sliced out bits
}

/* Undoes the effect of |slice|, inserting zero bits where needed.

  This cannot be used to undo the use of |slice| to express a vector on a
  subspace basis (you need to form a linear combination of the basis for
  that), but when |slice| is used during a projection modulo a subspace (as
  happens in subquotients), then |unslice| can reconstruct a preimage.

  Note that |unslice| changes (weakly increases) the size of its |BitVector|.
*/
template<size_t dim>
void BitVector<dim>::unslice(BitSet<dim> t,size_t new_size)
{
  BitSet<dim> result;
  for (typename BitSet<dim>::iterator it=t.begin(); it();
       d_data>>=1,++it)
    result.set(*it,d_data[0]);
  d_data=result;
  d_size=new_size;
}

} // |namespace bitvector|

/*****************************************************************************

        Chapter II -- The BitMatrix class

******************************************************************************/

namespace bitvector {

/*
  Construct the matrix whose columns are given by the vectors in |b|.

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

  // this loop is necessary to reduce each column to a |BitSet<dim>|
  for (size_t j = 0; j < b.size(); ++j)
  {
    assert(b[j].size()==num_rows);
    d_data.push_back(b[j].data());
  }

}

template<size_t dim>
BitMatrix<dim>::BitMatrix(const matrix::Matrix<int>& m) // set modulo 2
  : d_data(m.numColumns(),BitSet<dim>())
  , d_rows(m.numRows())
  , d_columns(m.numColumns())
{
  assert(m.numRows()<=dim);
  for (size_t i=0; i<d_rows; ++i)
    for (size_t j=0; j<d_columns; ++j)
      d_data[j].set(i, (m(i,j)&1)!=0 );
}

/******** accessors *********************************************************/

/*
  Applies our |BitMatrix| to |source| and returns the result

  It is assumed that |d_columns| is equal to |source.size()|. The result size
  will be set from |d_rows|.
*/
template<size_t dim>
BitVector<dim> BitMatrix<dim>::operator*(const BitVector<dim>& source) const
{
  assert(d_columns==source.size());
  // use |BitSet| vertion of |combination|, then size the result to |d_rows|
  return BitVector<dim> (combination(d_data,source.data()),d_rows);
}

// the same, but acting on the right
template<size_t dim>
BitVector<dim> BitMatrix<dim>::right_act(const BitVector<dim>& source) const
{
  assert(d_rows==source.size());
  BitVector<dim> result(d_columns);
  for (unsigned j=0; j<numColumns(); ++j)
    result.set(j,source.dot(column(j)));
  return result;
}

/*
  A pipe-dream version of binary matrix * vector multiplication.

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

// The i-th row of the matrix.
template<size_t dim>
BitVector<dim> BitMatrix<dim>::row(size_t i) const
{
  assert(d_columns<=dim);
  BitVector<dim> r(d_columns);

  for (size_t j = 0; j < d_columns; ++j)
    r.set(j,test(i,j));
  return r;
}


// Put in |b| a basis of the image of the matrix.
template<size_t dim>
BitVectorList<dim> BitMatrix<dim>::image() const
{
  BitVectorList<dim> b;
  std::vector<size_t> f; // auxialiry to record leading bit positions

  for (size_t j = 0; j < d_columns; ++j)
    spanAdd(b,f,column(j));  // add column |j| to span of |b|

  return b;
}


/* Puts in |b| an echelon basis for the kernel of the matrix (but not
  the canonical basis of that subspace; it is echelon from high to low bit
  positions).

  This is based on making the transpose matrix echelon, and then solving the
  resulting trivialised system of equations.
*/
template<size_t dim> BitVectorList<dim> BitMatrix<dim>::kernel() const
{
  BitVectorList<dim> result;

  if (numColumns()==0)
    return result; // no unknowns, so there are no nontrivial solutions

  assert(numColumns()<=dim);
  std::vector<BitVector<dim> > eqn;
  eqn.reserve(d_rows);

  // get rows of the matrix into |eqn|

  for (size_t i = 0; i < d_rows; ++i)
    eqn.push_back(row(i));

  // normalize |eqn|

  BitSet<dim> t; // will flag a subset of [0,numColumns()[

  Gauss_Jordan(t,eqn);

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

  result.reserve(numColumns()-t.count());
  for (size_t j = 0; j < numColumns(); ++j)
    if (not t[j])
    {
      BitVector<dim> v(numColumns(),j); // start with unit vector $e_j$
      size_t c = 0;
      for (typename BitSet<dim>::iterator it=t.begin(); it(); ++it, ++c)
	v.set(*it,eqn[c][j]); // use |eqn[c]| to find |v[*it]|
      result.push_back(v);
    }

  return result;
}

/******** manipulators *******************************************************/


/*
  Increment |*this| by adding |m|.

  Precondition: |m| has the same size as the current matrix.
*/
template<size_t dim>
BitMatrix<dim>& BitMatrix<dim>::operator+= (const BitMatrix<dim>& m)
{
  assert(d_rows==m.d_rows);
  assert(d_columns==m.d_columns);
  for (size_t j = 0; j < numColumns(); ++j)
    d_data[j] ^= m.d_data[j]; // not |+|, since these are |BitSet|s

  return *this;
}


/*
  Right multiply our BitMatrix by |m|.

  As in apply, this can be done rather efficently by computing a whole column
  at a time.

  NOTE : of course |m.numRows()| must be equal to |numColumns()|.
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


/*
  For a given matrix $A$, return a solution to $ABA=A$. This always exists,
  and will be the inverse of $A$ in case $A$ is invertible. In case $A$ is
  surjective $B$ will be a left inverse (or section) of $A$ (whence the name)
  while if $A$ is injective it will be a right inverse.

  Let $A$ be a $r\times c$ matrix. We initialise a new matrix $B$ to the
  identity matrix of size $c$, and let another matrix $M$ to $A$; the
  invariant will be that $M=A.B$. We perform in parallel column operations to
  $M$ and to $B$, transforming the latter into one in which each nonzero
  column has a "pivot" entry that is the unique nonzero entry of its row, as
  follows. For each column of $B$ in order a nonzero pivot is chosen if
  possible and, then cleared out of the remainder of its row.

  At the end the nonzero columns of $M$ have become standard basis vectors of
  $(Z/2Z)^r$; the corresponding column of $B$ (its preimage by $A$) is moved
  to the appropriate column of $B$ (its remaining columns are zero).
*/
template<size_t dim> BitMatrix<dim> BitMatrix<dim>::section() const
{

  std::vector<BitSet<dim> > basis (d_columns,BitSet<dim>()); // square matrix
    for (unsigned int i = 0; i<d_columns; ++i)
      basis[i].set(i); // init: $n\times n$ identity matrix, for |n==d_columns|

  std::vector<BitSet<dim> > col(d_data); // copy columns of our matrix

  BitSet<dim> pivots; // $r$ will be set if some column has its pivot in row $r$
  unsigned int pivot_col[dim]; // column number having pivot in row $r$

  for (unsigned int k=0; k<d_columns; ++k)
  {
    const BitSet<dim> col_k = col[k];
    if (col_k.none()) // then |k| will not be stored in |pivot_col| at all
      continue; // we'll forget about |col_k|, and |basis[k]| remains as it is

    const unsigned cur_pivot = col_k.firstBit(); // row in which pivot is found
    pivots.set(cur_pivot); // mark this row as a pivot row
    pivot_col[cur_pivot]=k; // associate current column |k| to row |pivot|

    const auto b_k = basis[k]; // basis vector to be added to some others

    // ensure that existing pivot columns get zero entry at inr |cur_pivot|
    for (auto it=pivots.begin(); // traverse previous pivot rows
	 it() and *it<cur_pivot; // consider only rows with pivot above ours
	 ++it)
    {
      const unsigned int j = pivot_col[*it]; // column where that row has pivot
      if (col[j].test(cur_pivot)) // see if that column nonzero at |cur_pivot|
      {
	col[j] ^= col_k; // if so clear out the entry using our column |col_k|
	basis[j] ^= b_k; // and let the basis matrix follow suit
      }
    }

    // also clear row |cur_pivot| in (yet) non-pivot columns: those beyond |k|
    for (unsigned int j=k+1; j<d_columns; ++j)
      if (col[j].test(cur_pivot))
      {
	col[j] ^= col_k;
	basis[j] ^= b_k;
      }
    // now $j=k$ is the unique index for which |col[j].test(cur_pivot)| holds
  }

/* at this point, the number of pivots equals the rank $r$ of $A$, and if $P$
   is the row-selection matrix for |pivots|, ten $A'=P.A$ is surjective. We
   shall take a right-inverse $B'$ of $A'$, and put $B=B'.P$. From the fact
   that $A'.B'=I_r$ it follws that $B.A.B=B$, and from $\ker(A)=\ker(A')$ it
   can be deduced that $A.B.A=A$. By construction, column $j$ of $B$ is
   nonzero only if |pivots.test(j)| holds, and |basis[pivot_col[j]]| will do.
 */
  BitMatrix<dim> B(d_columns,d_rows); // transpose shaped zero bitmatrix

  for (auto it=pivots.begin(); it(); ++it) // for |r| member of |pivots|
    B.setColumn(*it,basis[pivot_col[*it]]); // store pre-image of $e_r$

  // remainder of |B| remains zero

  return B;
}


// Reset the matrix to zero.
template<size_t dim> void BitMatrix<dim>::reset()
{
  for (unsigned long j = 0; j < d_data.size(); ++j)
    d_data[j].reset();
}

// Resize the matrix to |m| rows, |n| columns, leaving data around
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

// Transpose the matrix.
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
    const BitSet<dim>& e)
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
  BitSet<dim> combination(const std::vector<BitSet<dim> >& b,
				  const BitSet<dim>& coef)
  {
    BitSet<dim> result(0);
    for (size_t i = 0; i < b.size(); ++i)
      if (coef[i])
	result ^= b[i]; // not |+| here, these are |BitSet|s.
    return result;
  }

/*!
  Find out whether any combination of the vectors in |b| adds to |rhs|, and if
  so flag such a combination in the bits of |c|. Nothing changes when |false|.
*/
template<size_t dim>
  bool combination_exists(const std::vector<BitVector<dim> >& b,
			  const BitVector<dim>& rhs,
			  BitSet<dim>& c)
{
  if (b.size() == 0) // then there are no unknowns to solve
  {
    if (rhs.isZero())
    {
      c.reset(); // clear any previously set bits
      return true; // list of 0 "solved" unknowns
    }
    else
      return false; // no unknowns, no solutions
  }

  size_t n = b[0].size();
  assert(n==rhs.size());

  std::vector<BitSet<dim> > a; // list of normalised eqns, lhs part
  BitSet<dim> rh;      // corresponding right hand sides
  std::vector<size_t> f;       // list indicating "pivot" positions in |a|

  for (size_t i = 0; i < n; ++i)
  {
    BitSet<dim> r;
    // set r to i-th row of matrix whose columns are the |b[j]|
    for (size_t j = 0; j < b.size(); ++j)
      r.set(j,b[j][i]);
    bool x = rhs[i]; // now $(r,x)$ is one of the equations to solve

    // normalize |r| with respect to |a|: clear coefficients at previous pivots
    for (size_t j = 0; j < f.size(); ++j)
      if (r[f[j]])
      {
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
	if (a[j][m])
	{
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
  into |sol| and returning |true|, or return |false| if no solution exists.

  Here |eqn| holds a system of equations, the last bit of each being
  interpreted as the right hand side.
*/
template<size_t dimsol, size_t dimeq>
bool solvable(const std::vector<BitVector<dimeq> >& eqns,
	      BitVector<dimsol>& sol)
{
  if (eqns.empty()) // nothing to do, and we cannot even dimension |sol|
    return true; // just say "anything goes", in particular the original |sol|

  std::vector<BitVector<dimeq> > eqn(eqns); // make a local copy
  BitSet<dimeq> pivots;
  /*
    We solve by imagining an extra unknown for the final position,
    solve the homogenous system (find the kernel of transpose matrix), and
    look for a kernel element with final coordinate equal to $-1$
    (which working over $Z/2Z$ is of course the same as $1$).
  */
  Gauss_Jordan(pivots,eqn); // get equation matrix into reduced echelon form

  /* now the system is contradictory if and only if the last bit in |t| is
     set, since that means some equation sets the corresponding indeterminate
     to 0 (because it has its \emph{leading} bit at the final position), while
     absence of any equation for that indeterminate means we can make it $-1$.
   */
  unsigned int last = eqns[0].size()-1;
  if (pivots[last]) // some equation sets the extra unknown to $0$
    return false; // note that we leave |sol| at its original size here

  sol.resize(last); // otherwise |sol| is one shorter than the equations
  sol.reset(); // zero out all unknowns: all free unknowns are taken to be $0$

  // we set any non-free unknown, flagged by |pivots|, to the right hand side
  // of the corresponding equation (which involves no other non-free unknowns)
  size_t c = 0;
  for (typename BitSet<dimeq>::iterator it=pivots.begin(); it(); ++it, ++c)
    sol.set(*it,eqn[c][last]); // set equated unknown |c| (at |*it|)

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
  position $j'$ where number |j| among the raised bits of |t| is found;
  consequently, for any $v\in V$, the coordinate of |b[j]| in $v$ is $v[j']$.

  This function works essentially by repeatedly calling |normalSpanAdd| for
  the vectors of |b|, replacing |b| by the resulting canonical basis at the
  end. However, the selected coordiante positions |f| do not come out
  increasingly this way, so we have to sort the canonical basis by leading bit
  position. This amounts to setting $a'[k]=a[p(k)]$ where
  $p(0)\ldots,p(l-1)$ is the result of sorting $f[0]\ldots,f[l-1]$
  with $l=f.size()=a.size()$.
*/
template<size_t dim>
  void Gauss_Jordan(BitSet<dim>& t, std::vector<BitVector<dim> >& b)

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
  $\{0,\ldots,d-1\}$ such that the standard basis vectors $e_i$ for $i\in I$
  generate a complementary subspace $e_I$ to $V$ (one can find such an $I$ by
  repeatedly throwing in $e_i$s linearly independent to $V$ and previously
  chosen ones). The normal basis of $V$ corresponding to $I$ is obtained by
  projecting the $e_j$ for $j$ in the complement $J$ of $I$ onto $V$ along
  $e_I$ (i.e., according to the direct sum decompostion $k^d=V\oplus e_I$).
  This can be visualised by viewing $V$ as the function-graph of a linear map
  from $k^J$ to $k^I$ (with the coordinates of domain and codomain interwoven
  at the positions $J$ and $I$, respectively); then the normal basis is the
  lift to the graph $V$ of the standard basis of $k^J$. We define the
  canonical basis of $V$ to be the normal basis for the complement $I$ of the
  lexicographically minimal possible set $J$ (lexicographic for the increasing
  sequences representing the subsets; in fact $I$ is lexicographically maximal
  since complementation reverses this ordering on fixed-size subsets). One can
  find this $J$ by repeatedly choosing the smallest index such that the
  projection from $V$ defined by extracting the coordinates at the selected
  indices remains surjective to the set of all possible coordinate tuples.
  (Intersect $V$ with the subspace defined by setting previously selected
  coordinates to zero; choose the smallest coordinate that can be nonzero.)

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
  \brief Enlarges the basis |a| so as to span |v|.

  This is a simplified version of |normalSpanAdd|

  It is assumed that |a| contains a list of independent bitvectors all of size
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

template<size_t dim> int_Vector lift(const BitVector<dim>& v)
{ int_Vector result(v.size(),0);
  for (auto it=v.data().begin(); it(); ++it)
    result[*it]=1;
  return result;
}


/* functions never called have been grouped here, MvL */
#if 0

/*
  Flag into |c| a subset of the standard basis that spans a complementary
  subspace to the subspace spanned by |b|.

  It is assumed that the vectors in b are all of the same size, but not
  necessarily independent.

  NOTE : we need to pass the dimension in case |b| is empty; we don't want
  to set bits beyond that (so that |c.count()|, for instance, yields the
  correct dimension of the complement.)
*/
template<size_t dim>
  void complement(BitSet<dim>& c,
		  const std::vector<BitVector<dim> >& b,
		  size_t d)
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

// Tell whether the system of bitvectors is independent.
template<size_t dim> bool isIndependent(const std::vector<BitVector<dim> >& b)
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


/*
  Put in p the matrix of the projection on the canonical complement to the
  span of b.

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
  BitSet<dim> c;
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
    const BitSet<dim>& vd = v.data();
    p.setColumn(f[j],vd);
  }
}

/*
  Put into m the matrix of the reflection defined by a and a_check.

  Precondition: a and a_check have same size; <a,a_check> = 0;

  This is the operator x -> x + <x,a_check>a (we can write + because we are
  in characteristic two.)
*/
template<size_t dim>
  void reflectionMatrix(BitMatrix<dim>& m, const BitVector<dim>& a,
			const BitVector<dim>& a_check)
{
  identityMatrix(m,a.size());

  for (size_t j = 0; j < a.size(); ++j)
    if (a_check[j])
      m.addToColumn(j,a);
}

/*
  Write into r the relations among the elements in b.

  In other words, it solves the system of equations defined by the _rows_ in
  the matrix whose columns are given by b.
*/
template<size_t dim>
  void relations(std::vector<BitVector<dim> >& rel,
		 const std::vector<BitVector<dim> >& b)
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

  BitSet<dim> t; // will flag a subset of [0,r[
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
  SmallBitVector combination
   (const std::vector<SmallBitVector>& b,
    size_t n,
    const BitSet<constants::RANK_MAX>& e);

template
  BitSet<constants::RANK_MAX> combination
  (const std::vector<BitSet<constants::RANK_MAX> >&,
   const BitSet<constants::RANK_MAX>&);

template
  void Gauss_Jordan(BitSet<constants::RANK_MAX>& t,
		    std::vector<SmallBitVector>& b);

template
  bool combination_exists(const std::vector<SmallBitVector>& b,
			  const SmallBitVector& rhs,
			  BitSet<constants::RANK_MAX>& c);
template
  bool solvable(const std::vector<BitVector<constants::RANK_MAX+1> >& eqns,
		SmallBitVector& sol);

template void identityMatrix(BitMatrix<constants::RANK_MAX>&, size_t);
template void initBasis(std::vector<SmallBitVector>&, size_t);

template int_Vector lift(const SmallBitVector& v);

template class BitVector<constants::RANK_MAX>;   // |SmallBitVector|
template class BitVector<constants::RANK_MAX+1>; // |BinaryEquation|
template class BitMatrix<constants::RANK_MAX>;   // |BinaryMap|

template class BitVector<64ul>; // used in atlas function |subspace_normal|
template
   void initBasis<64ul>(std::vector<BitVector<64ul> >& b, size_t r); // idem
template class BitMatrix<64ul>; // used in realex function |binary_invert|

template
  BitVector<constants::RANK_MAX>::BitVector
    (const matrix::Vector<int>& weight);
template
  BitVector<constants::RANK_MAX>::BitVector
    (const matrix::Vector<long long int>& weight);
template
  BitVector<constants::RANK_MAX+1>::BitVector
    (const matrix::Vector<int>& weight);
template
  BitVector<64ul>::BitVector (const matrix::Vector<int>& weight);

} // |namespace bitvector|

} // |namespace atlas|
