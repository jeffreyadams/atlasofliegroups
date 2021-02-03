/*
  matreduc.cpp

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
#include "matreduc.h"

#include <cstdlib>

#include "bitmap.h"
#include "permutations.h"
#include "matrix.h"
#include "arithmetic.h"
#include "bigint.h"
#include "sl_list.h"

namespace atlas {

namespace matreduc {

// make |M(i,j)==0| and |M(i,k)>0| by operations with columns |j| and |k|
// precondition |M(i,k)>0|. Returns whether determinant -1 operation applied.
template<typename C>
bool column_clear(matrix::PID_Matrix<C>& M, size_t i, size_t j, size_t k)
{
  bool neg=M(i,j)<C(0);
  if (neg) // happens often in Cartan matrices; result will be more pleasant
    M.columnMultiply(j,C(-1)); // if we ensure non-negativity the easy way
  do // no initial test needed because of precondition |M(i,k)>0|
  {
    M.columnOperation(j,k,-(M(i,j)/M(i,k))); // makes |M(i,j)>=0| smaller
    if (M(i,j)==C(0))
      return neg;
    M.columnOperation(k,j,-(M(i,k)/M(i,j))); // makes |M(i,k)>=0| smaller
  }
  while (M(i,k)!=C(0)); // perform final swap if leaving at this "normal" exit
  M.swapColumns(j,k);
  return not neg; // take into account additional sign change here
}

// make |M(i,j)==0| and |M(k,j)>0| by operations with rows |i| and |k|
// precondition |M(k,j)>0|. Returns whether determinant -1 operation applied.
template<typename C>
bool row_clear(matrix::PID_Matrix<C>& M, size_t i, size_t j, size_t k)
{
  bool neg=M(i,j)<C(0);
  if (neg)
    M.rowMultiply(i,C(-1));
  do // no initial test needed because of precondition |M(k,j)>0|
  {
    M.rowOperation(i,k,-(M(i,j)/M(k,j))); // makes |M(i,j)>=0| smaller
    if (M(i,j)==C(0))
      return neg;
    M.rowOperation(k,i,-(M(k,j)/M(i,j))); // makes |M(k,j)>=0| smaller
  }
  while (M(k,j)!=C(0)); // perform final swap if leaving at this "normal" exit
  M.swapRows(i,k);
  return not neg; // take into account additional sign change here
}

/*
   The same operations, but with "recording" matrix |rec| on which same ops
   are done. The number of operations needed is quite small on average, so
   the cheapest solution is probably to do recording operations immediately,
   rather than to compute an intermediate 2x2 matrix to be applied once.
*/
// Make |M(i,j)==0| and keep |M(i,k)>0| by operations with columns |j| and |k|
// precondition |M(i,k)>0|
template<typename C>
bool column_clear(matrix::PID_Matrix<C>& M, size_t i, size_t j, size_t k,
		  matrix::PID_Matrix<C>& rec)
{
  /* Rather than forcing the sign of |M(i,j)|, as above, we shall use
     |arithmetic::divide| the first time here. The reason is heuristic; this
     may sometimes avoid unnecessary negative entries in |lattice::kernel|
  */
  C q = -arithmetic::divide(M(i,j),M(i,k));
  while (true) // not entering loop at an exit point, though it has 2 of them
  {
    M.columnOperation(j,k,q); // now |M(i,j)>=0|
    rec.columnOperation(j,k,q);
    if (M(i,j)==C(0))
      return false; // determinant $1$ here
    q=-(M(i,k)/M(i,j));
    M.columnOperation(k,j,q); // we still have |M(i,k)>=0|
    rec.columnOperation(k,j,q);
    if (M(i,k)==C(0))
      break; // perform final swap at this exit
    q = -(M(i,j)/M(i,k)); // next time around prepare |q| without |divide|
  }
  M.swapColumns(j,k);
  rec.swapColumns(j,k);
  return true;  // determinant $-1$ here
}

/*
   Make |M(i,j)==0| and keep |M(k,j)>0| by operations with rows |i| and |k|,
   precondition |M(k,j)>0|. Here we also define a variant where row operations
   are recorded inversely, in other words by column operations with negated
   coefficient; template parameter |direct| tells whether recording is direct
*/
template<typename C, bool direct>
bool row_clear(matrix::PID_Matrix<C>& M, size_t i, size_t j, size_t k,
	       matrix::PID_Matrix<C>& rec)
{
  // Rather than forcing the sign of |M(i,j)|, we shall use |arithmetic::divide|
  // For |direct==false|, this often avoids touching column |i| of |rec| at all
  C q = arithmetic::divide(M(i,j),M(k,j));
  while (true) // not entering loop at an exit point, though it has 2 of them
  {
    M.rowOperation(i,k,-q); // now |M(i,j)>=0|, even if |M(i,j)| was negative
    if (direct)
      rec.rowOperation(i,k,-q);
    else // do inverse column operation
      rec.columnOperation(k,i,q);
    if (M(i,j)==C(0))
      return false;
    q=M(k,j)/M(i,j);
    M.rowOperation(k,i,-q);
    if (direct)
      rec.rowOperation(k,i,-q);
    else // inverse column operation
      rec.columnOperation(i,k,q);
    if (M(k,j)==C(0))
      break; // perform final swap at this exit
    q = M(i,j)/M(k,j); // next time around prepare |q| without |divide|
  }
  M.swapRows(i,k);
  if (direct)
    rec.swapRows(i,k);
  else
    rec.swapColumns(i,k);
  return true;
}


/*
  Find |row|, |col| of determinant $1$ such that $row*M*col$ is diagonal, and
  return diagonal entries (positive except maybe first). Result is not unique.
  Also coefficient growth in |row| and |col| can be extreme for large matrices.
*/
template<typename C>
std::vector<C> diagonalise(matrix::PID_Matrix<C> M, // by value
			   matrix::PID_Matrix<C>& row,
			   matrix::PID_Matrix<C>& col)
{
  const size_t m=M.numRows();
  const size_t n=M.numColumns();

  row=matrix::PID_Matrix<C>(m); // initialise |row| to identity matrix
  col=matrix::PID_Matrix<C>(n);
  std::vector<C> diagonal;
  if (n==0 or m==0) return diagonal; // take out these trivial cases
  diagonal.reserve(std::min(m,n));

  bool row_minus=false, col_minus=false; // whether $\det(row//col)=-1$
  bitmap::BitMap pivot_columns(n); // columns of |M| where a pivot was found
  matrix::PID_Matrix<C> ops;

  for (size_t k=0,l=0; l<n; ++l) // |k| may or may not be increased inside loop
  {
    bool flip=false;
    C d = gcd(M.partial_column(l,k,m),&ops,flip,0);
    if (d==0) // if partial column was already zero
      continue; // advance in loop on |l|, do not increment |k|

    pivot_columns.insert(l); // there will be a pivot at |M(k,l)|

    row_minus=flip;
    ops.transpose(); // because recorded column operations applied as row ops
    row_apply(M,ops,k);
    row_apply(row,ops,k);
    assert(M(k,l)==d);

    C old_d; // used in final condition of next loop

    do // exit when row and column from |M(k,l)| onwards are zero
    {
      old_d=d; flip=false;
      d=gcd(M.partial_row(k,l,n),&ops,flip,0);
      col_minus^=flip;
      column_apply(M,ops,l);
      column_apply(col,ops,l);
      assert(M(k,l)==d);
      if (d==old_d)
	break; // no improvement to |d| means row was cleared beyond |l|

      old_d=d; flip=false;
      d=gcd(M.partial_column(l,k,m),&ops,flip,0);
      row_minus^=flip;
      ops.transpose(); // because recorded column operations applied as row ops
      row_apply(M,ops,k);
      row_apply(row,ops,k);
      assert(M(k,l)==d);
    }
    while(d<old_d);

    row_minus^=flip;

    diagonal.push_back(d); // record positive gcd that finally remains
    ++k; // record that this row contains a pivot, so has been dealt with
  } // |for(l)|

  // adapt |col| by stable-sorting columns, moving |zero_columns| to the end
  { const auto n_piv=pivot_columns.size();
    if (pivot_columns.position(n_piv)<n_piv) // pivot columns not left-adjusted
    { permutations::Permutation pi(pivot_columns.begin(),pivot_columns.end());
      pivot_columns.take_complement(); // now get the non-pivot columns
      pi.insert(pi.end(),pivot_columns.begin(),pivot_columns.end()); // add them
      permute_columns(col,pi);
      col_minus ^= (sign(pi)<0);
    }
  }

  if (diagonal.size()>0 and row_minus!=col_minus)
    diagonal[0] = -diagonal[0];
  if (row_minus)
    row.rowMultiply(0,C(-1)); // ensure determinant of |row| is |1|
  if (col_minus)
    col.columnMultiply(0,C(-1)); // ensure determinant of |col| is |1|

  return diagonal;
}

// auxiliary function for |adapted_basis|
// find |k| with |k%row[0]| small but nonzero, return whether found
template<typename C>
  bool find_small_remainder (matrix::Vector<C> row, size_t& k)
{
  C a=row[0];
  assert(a>C(0));
  C min=a;
  for (auto it=std::next(row.begin()); it!=row.end(); ++it)
  { C r = arithmetic::remainder(*it,a);
    if (r>0 and r<min)
      min=r, k=it-row.begin();
  }
  return min<a; // at least one positive remainder found
}

/*
  The following is a variation of |diagonalise|, used in cases in which we are
  mostly interested in the matrix |B=row.inverse()| and possibly also the
  vector |diagonal|, because they give a transparent expression for the image
  (column span) of |M|: that image is the same as that of $B*D$ where $D$ is
  the diagonal matrix corresponding to |diagonal| (extended with null rows to
  match the height of |M| = size of |B|), in other words it is spanned by the
  multiples of the columns of $B$ by their |diagonal| factors.

  The procedure, which returns |B|, follows mostly the same steps as
  |diagonalise|, but the differences in handling the matrix |row| (where we
  apply instead of row operations the inverse column operations) are
  sufficiently important to justify the code duplication. We take advantage of
  this duplication to arrange the algorithm for minimal use of row operations,
  which besides being more efficient tends to give a basis more closely
  related to the original matrix. All diagonal entries are positive.
*/
template<typename C>
matrix::PID_Matrix<C> adapted_basis(matrix::PID_Matrix<C> M, // by value
				    std::vector<C>& diagonal)
{
  size_t m=M.numRows();
  size_t n=M.numColumns(); // in fact start of known null columns

  matrix::PID_Matrix<C> result (m); // initialise |result| to identity matrix
  matrix::PID_Matrix<C> ops; bool flip=false;
  diagonal.clear(); diagonal.reserve(n); // maximum, maybe not needed

  bitmap::BitMap kept_rows(m);
  size_t j=0;
  for (size_t i=0; i<m; ++i)
  {
    int d = gcd(M.partial_row(i,j,n),&ops,flip,0);
    if (d==0)
      continue; // without advancing |j| or recoring |i|
    kept_rows.insert(i);
    column_apply(M,ops,j);
    assert(M(i,j)==d);
    size_t k; // row number
    while (find_small_remainder(M.partial_column(j,i,m),k)) // sets |k| when true
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    { k += i; // convert relative |k| to actual row number
#pragma GCC diagnostic pop
      auto q=arithmetic::divide(M(k,j),M(i,j));
      for (size_t l=j; l<n; ++l) // columns |M| can be nonzero starting from |j|
      {
	C tmp=M(i,l);
	M(i,l)=M(k,l)-q*tmp;
	M(k,l) = tmp;
      }
      assert(M(i,j)>0 and M(i,j)<M(k,j));

      for (size_t l=0; l<m; ++l) // apply inverse to columns |i,k| of |result|
      {
	C tmp=result(l,k);
	result(l,k)=result(l,i)+q*tmp;
	result(l,i) = tmp;
      }

      d = gcd(M.partial_row(i,j,n),&ops,flip,0);
      column_apply(M,ops,j); // ensure remainder of row |i| is zero again
      assert(M(i,j)==d);
    }
 /* once divisibility by the pivot is attained, we need not actually clear out
    the remainder of column |j| of |M|: it will subsequently be ignored. But
    we do need to adapt |result| corresponding to what would would be done here
 */
    for (k=i+1; k<m; ++k)
    { assert(M(k,j) % M(i,j) == 0); // since |find_small_remainder| said so
      auto q=M(k,j)/M(i,j); // exact division, so safe even when |M(k,j)<0|
      // |M.rowOperation(k,i,-q);| is what we conceptually do
      result.columnOperation(i,k,q); // inverse multiplication on the right
    }
    assert(M(i,j)==d);
    diagonal.push_back(d);
    ++j;
  } // |for(i)|

  if (not kept_rows.full())
  {
    containers::sl_list<matrix::Vector<C> > kept_columns,dropped_columns;
    for (size_t i=0; i<m; ++i)
      (kept_rows.isMember(i) ? kept_columns : dropped_columns).
	push_back(result.column(i));
    kept_columns.append(std::move(dropped_columns));
    auto it=kept_columns.begin();
    for (size_t i=0; i<m; ++i,++it)
      result.set_column(i,*it);
  }

  return result;
} // |adapted_basis|

/*
  For a true Smith basis, we must assure divisibility of successive elements
  in |diagonal|. Since small factors tend to be extracted first, there
  probably remains little work to do, and there is not much against assuring
  divisibilty for adjacent pairs, and repeating this in a bubble-sort like
  manner. Rather than iterating matrix operations for a given pair, we use the
  following formula valid whenever gcd(a,b)=d=pa+qb:

    [   1     1   ]   [ a  0 ]   [ p  -b/d ]   [ d   0   ]
    [             ] * [      ] * [         ] = [         ],
    [ pa/d-1 pa/d ]   [ 0  b ]   [ q   a/d ]   [ 0  ab/d ]

  where the outer matrices have deteminant $1$. The inverse of the first one:

    [  pa/d  -1 ]    [ 1  -1 ]   [  1     0 ]
    [           ] =  [       ] * [          ]
    [ 1-pa/d  1 ]    [ 0   1 ]   [ 1-pa/d 1 ]

  is applied to the right to the initial |adapted_basis| for correction.
  The numbers |lcm=ab/d|, |d| and |pa| are obtained from |arithmetic::lcm|
*/
template<typename C>
matrix::PID_Matrix<C> Smith_basis(const matrix::PID_Matrix<C>& M,
				  std::vector<C>& diagonal)
{
  matrix::PID_Matrix<C> result = adapted_basis(M,diagonal);
  size_t start=0, stop=diagonal.size()>0 ? diagonal.size()-1 : 0, new_stop=stop;
  while (start<stop)
  {
    size_t new_start=new_stop; // if unchanged, this will be last of |while|
    for (size_t i=start; i<stop; ++i)
      if (diagonal[i+1]%diagonal[i]!=0) // failure of divisibility condition
      {
	arithmetic::Denom_t d, pa;
        diagonal[i+1]=arithmetic::lcm(diagonal[i],diagonal[i+1],d,pa);
	diagonal[i]=d;

	// correct the basis according to the row operations implicilty applied
	result.columnOperation(i+1,i,-1); result.columnOperation(i,i+1,1-pa/d);

	if (new_start==stop) // only change |new_start| on first swap
	  new_start = i>0 ? i-1 : 0; // need to back up one place next time
        new_stop=i; // if this is last switch, diagonal[i+1],... are final
      }
    start=new_start; stop=new_stop;
  }
  return result;
}


template<typename C> // find a solution |x| for |A*x==b|
bool has_solution(const matrix::PID_Matrix<C>& A, matrix::Vector<C> b)
{
  matrix::PID_Matrix<C> row,col;
  std::vector<C> diagonal = diagonalise(A,row,col); // $R*A*C=D$ diagonal
  row.apply_to(b); // left multiply equation by $R$, giving $D*C^{-1}*x=R*b$

  // now solve for the value of $C^{-1}*x$
  for (unsigned int i=b.size(); i-->0; )
    if ((i<diagonal.size() ? b[i] % diagonal[i] : b[i]) != 0)
      return false;
  return true;
}

template<typename C> // find a solution |x| for |A*x==b|
matrix::Vector<C> find_solution(const matrix::PID_Matrix<C>& A,
				matrix::Vector<C> b)
{
  matrix::PID_Matrix<C> row,col;
  std::vector<C> diagonal = diagonalise(A,row,col); // $R*A*C=D$ diagonal
  row.apply_to(b); // left multiply equation by $R$, giving $D*C^{-1}*x=R*b$

  // now solve for the value of $C^{-1}*x$
  for (unsigned int i=0; i<diagonal.size(); ++i)
  {
    C d = diagonal[i];
    if (b[i] % d != C(0))
      throw std::runtime_error("unsolvable integral system");
    b[i] /= d;
  }
  for (unsigned int i=diagonal.size(); i<b.size(); ++i)
    if (b[i] != C(0))
      throw std::runtime_error("unsolvable system");

  b.resize(col.numRows(),0); // adapt size in opposite sense to $A$
  col.apply_to(b); // finally reconstruct value of |x|
  return b;
}

 // instantiations

using Num = arithmetic::Numer_t;

template
int gcd (matrix::Vector<int> row, matrix::PID_Matrix<int>* col, bool& flip,
	 size_t dest);
template
bool column_clear(matrix::PID_Matrix<int>& M, size_t i, size_t j, size_t k);
template
bool row_clear(matrix::PID_Matrix<int>& M, size_t i, size_t j, size_t k);

template
std::vector<int> diagonalise(matrix::PID_Matrix<int> M,
			     matrix::PID_Matrix<int>& row,
			     matrix::PID_Matrix<int>& col);

template
matrix::PID_Matrix<int> Smith_basis(const matrix::PID_Matrix<int>& M,
				    std::vector<int>& diagonal);

template
matrix::PID_Matrix<int> adapted_basis(const matrix::PID_Matrix<int> M,
				      std::vector<int>& diagonal);

template
bool has_solution(const matrix::PID_Matrix<int>& A, matrix::Vector<int> b);

template
matrix::Vector<int> find_solution(const matrix::PID_Matrix<int>& A,
				  matrix::Vector<int> b);

} // |namespace matreduc|
} // |namespace atlas|
