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

/* transform |M| to column echelon form using only PID operations
   postcondition: |M| has unchanged column span, |result.size==M.numColumns|,
   and for |j| and |i=result.n_th(j)|: |M(i,j)>0| and |M(ii,j)==0| for |ii>i|
*/
template<typename C>
  bitmap::BitMap column_echelon(matrix::PID_Matrix<C>& M)
{
  bitmap::BitMap result(M.numRows()); // set of pivot rows found so far
  for (size_t j=M.numColumns(); j-->0;)
  { // all columns beyond |j| have pivot, column |j+1+p| in row |result.n_th(p)|
    size_t i=M.numRows();
    while (i-->0)
      if (M(i,j)!=C(0))
      {
	if (result.isMember(i)) // there is a previous pivot in row |i|
	{ // so use it to clear |M(i,j)| (and then continue loop decreasing |i|)
	  size_t k=j+1+result.position(i); // index of column with pivot at |i|
	  column_clear(M,i,j,k); // now column |j| is empty in row |i| and below
	}
	else // now column |j| will have pivot in row |i|
	{
          if (M(i,j)<0)
	    M.columnMultiply(j,-1); // ensure positive pivot entry
	  result.insert(i); // mark row |i| as a pivot row
	  size_t p=result.position(i);
	  if (p>0) // then move column |j| to the right |p| places
	  {
	    matrix::Vector<C> c=M.column(j); // save this column while shifting
	    for (size_t d=j; d<j+p; ++d)
	      M.set_column(d,M.column(d+1));
            M.set_column(j+p,c); // re-insert column at its destination place
	  }
	  break; // from |while| loop; done with decreasing |i|
	}
      } // |if (M(i,j)!=0)| and |while (i-->0)|
    if (i==size_t(-1)) // if no pivot found for column |j|; forget about it
      M.eraseColumn(j); // and no bit is set in |result| for |j| now
  } // |for(j-->0)|

  return result;
}

/* find |row|, |col| of determinant $1$ such that $row*M*col$ is diagonal, and
   return diagonal entries (positive except maybe first). Result is not unique.
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

  int row_sign = 1, col_sign=1; // determinants of row/column operations
  bitmap::BitMap pivot_columns(n); // columns of |M| where a pivot was found

  for (size_t k=0,l=0; l<n; ++l) // |k| may or may not be increased inside loop
  {
    size_t i=k;
    for (; i<m; ++i)
      if (M(i,l)!=C(0))
	break;
    if (i==m)
      continue; // advance in loop on |l|, do not increment |k|

    pivot_columns.insert(l); // there will be a pivot at |M(k,l)|

    if (i>k) // first non-zero entry below diagonal; add it to 0 on diagonal
    {
      C u = M(i,l)>0 ? 1 : -1;  // ensure diagonal entry becomes positive
      M.rowOperation(k,i,u);
      row.rowOperation(k,i,u);
    }
    else // we have |i==k|
    { ++i; // search for (further) nonzero column entries will start here
      if (M(k,l)<0) // negative pivot entry, make it positive
      {
	M.rowMultiply(k,-1);
	row.rowMultiply(k,-1);
	row_sign = -row_sign;
      }
    }
    // Now we have ensured |M(k,l)>0| and |M(ii,l)==0| for $k<ii<i$

    // sweep row and column alternatively with |M(k,l)| until both cleared
    while (true) // exit when row and column from |M(d,d)| onwards are zero
    {

      // |M(i,l)| is first (potentially) nonzero entry below |M(k,l)|
      for ( ; i<m; ++i) // ensure column |d| is null below diagonal
	if (row_clear<C,true>(M,i,l,k,row)) // makes |M(k,l)==gcd>0|,|M(i,l)==0|
	  row_sign = -row_sign;

      size_t j=l+1; // needs to survive next loop
      for ( ; j<n; ++j) // check if row |k| is zero beyond diagonal
	if (M(k,j)!=C(0))
	  break; // row not zero, work needs to be done
      if (j==n) // loop just terminated, so arm and leg are zero
	break; // so terminate outer loop

      // now |M(k,j)| is first nonzero entry right of |M(k,l)|
      for ( ; j<n; ++j) // ensure row |d| is null beyond diagonal
	if (column_clear(M,k,j,l,col)) // makes |M(k,l)==gcd>0| and |M(k,j)==0|
	  col_sign = -col_sign;

      for (i=k+1; i<m; ++i) // check if column |d| is still zero below diagonal
	if (M(i,l)!=C(0))
	  break; // column not zero, do not terminate outer loop
      if (i==m) // then whole column below |M(k,l)| is zero, and done for |k|
	break; // so terminate outer loop
    }

    diagonal.push_back(M(k,l)); // record positive gcd that finally remains

    ++k; // record that this row contains a pivot, so has been dealt with
  } // |for(l)|

  // adapt |col| by stable-sorting columns, moving |zero_columns| to the end
  { const auto n_piv=pivot_columns.size();
    if (pivot_columns.position(n_piv)<n_piv) // pivot columns not left-adjusted
    { permutations::Permutation pi(pivot_columns.begin(),pivot_columns.end());
      ~pivot_columns; // now get the non-pivot columns, and add them at end
      pi.insert(pi.end(),pivot_columns.begin(),pivot_columns.end());
      permute_columns(col,pi);
      col_sign *= sign(pi);
    }
  }

  if (diagonal.size()>0 and row_sign!=col_sign)
    diagonal[0] = -diagonal[0];
  if (row_sign<0)
    row.rowMultiply(0,-1); // ensure determinant of |row| is |1|
  if (col_sign<0)
    col.columnMultiply(0,-1); // ensure determinant of |col| is |1|

  return diagonal;
}

/* The following is a variation of |diagonalise|, used in cases in which we
   are mostly interested in the matrix |B=row.inverse()| and possibly also the
   vector |diagonal|, because they give a transparent expression for the image
   (column span) of |M|: that image is the same as that of $B*D$ where $D$ is
   the diagonal matrix corresponding to |diagonal| (extended with null rows to
   match the height of |M| = size of |B|), in other words it is spanned by the
   multiples of the columns of $B$ by their |diagonal| factors.

   The procedure, which returns |B|, follows mostly the same steps as
   |diagonalise|, but the differences in handling the matrix |row| (where we
   apply instead of row operations the inverse column operations) are
   sufficiently important to justify the code duplication. We take advantage
   of this duplication to arrange the algorithm for minimal use of row
   operations, which besides being more efficient tends to give a basis more
   closely related to the original matrix. All diagonal entries are positive.
 */
template<typename C>
matrix::PID_Matrix<C> adapted_basis(matrix::PID_Matrix<C> M, // by value
				    std::vector<C>& diagonal)
{
  size_t m=M.numRows();
  size_t n=M.numColumns(); // in fact start of known null columns

  matrix::PID_Matrix<C> result (m); // initialise |result| to identity matrix
  diagonal.clear(); diagonal.reserve(n); // maximum, maybe not needed

  for (size_t d=0; d<m and d<n; ++d)
  {
    while (M.column(d).isZero()) // ensure column |d| is nonzero
    {
      --n;
      if (d==n) // then this was the last column, quit
	return result;
      M.eraseColumn(d);
    }

    { // get nonzero entry from column |d| at (d,d), and make it positive
      size_t i=d;
      while (M(i,d)==C(0)) // guaranteed to terminate
	++i;
      if (i>d)
      {
	C u = M(i,d)>0 ? 1 : -1;
	M.rowOperation(d,i,u); // "copy" entry, making |M(d,d)==abs(M(i,d))|
	result.columnOperation(i,d,-u); // inverse operation on basis
      }
      else if (M(d,d)<0)
	M.columnMultiply(d,-1); // prefer a column operation here
    }

    // we prefer to start with column operations here, which need no recording
    for (size_t j=d+1; j<n; ++j)
      column_clear(M,d,j,d); // makes |M(d,d)==gcd>0| and |M(d,j)==0|
    // initial sweep is unlikely to be sufficient, so no termination test here

    size_t i=d+1,j;
    do // sweep column and row alternatively with |M(d,d)| until both cleared
    {
      for ( ; i<m; ++i)
	row_clear<C,false>(M,i,d,d,result);

      for (j=d+1; j<n; ++j)
	if (M(d,j)!=C(0))
	  break;

      if (j==n) // most likely because |row_clear| left column |d| unchanged
	break;

      for ( ; j<n; ++j)
	column_clear(M,d,j,d);

      for (size_t k=d+1; k<m; ++k)
	if (M(k,d)!=C(0))
	  break;
    }
    while(i<m); // then apparently some |column_clear| changed row |d|

    diagonal.push_back(M(d,d));
  } // |for d|

  return result;
}

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
  [ 1-pa/d  1 ]	   [ 0   1 ]   [ 1-pa/d 1 ]

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
template
bool column_clear(matrix::PID_Matrix<int>& M, size_t i, size_t j, size_t k);
template
bool row_clear(matrix::PID_Matrix<int>& M, size_t i, size_t j, size_t k);

template
bitmap::BitMap column_echelon<int>(matrix::PID_Matrix<int>& M);

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
matrix::Vector<int> find_solution(const matrix::PID_Matrix<int>& A,
				  matrix::Vector<int> b);

} // |namespace matreduc|
} // |namespace atlas|
