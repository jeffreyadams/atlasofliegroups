/*
  matreduc.cpp

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/
#include "matreduc.h"

#include "latticetypes_fwd.h"
#include "bitmap.h"
#include "matrix.h"
#include "intutils.h"

namespace atlas {

namespace matreduc {

// make |M(i,j)==0| and |M(i,k)>0| by operations with columns |j| and |k|
// precondition |M(i,k)>0|
template<typename C>
void column_clear(matrix::Matrix<C>& M, size_t i, size_t j, size_t k)
{
  do
  {
    long q = -intutils::divide(M(i,j),(unsigned long)(M(i,k)));
    M.columnOperation(j,k,q); // now |M(i,j)>=0|
    if (M(i,j)==C(0))
      return;
    M.columnOperation(k,j,-(M(i,k)/M(i,j))); // |operator/| is safe now
  }
  while (M(i,k)!=C(0));
  M.swapColumns(j,k); // only upon "normal" exit of loop swapping is needed
}

// make |M(i,j)==0| and |M(k,j)>0| by operations with rows |i| and |k|
// precondition |M(k,j)>0|
template<typename C>
void row_clear(matrix::Matrix<C>& M, size_t i, size_t j, size_t k)
{
  do
  {
    long q = -intutils::divide(M(i,j),(unsigned long)(M(k,j)));
    M.rowOperation(i,k,q); // now |M(i,j)>=0|
    if (M(i,j)==C(0))
      return;
    M.rowOperation(k,i,-(M(k,j)/M(i,j))); // |operator/| is safe now
  }
  while (M(k,j)!=C(0));
  M.swapRows(i,k); // only upon "normal" exit of loop swapping is needed
}

/* Same operations with recording matrix |rec| on which identical ops are done
   With the number of operations needed probably quite small on average, tying
   to delay and combine successive ones would in fact be more expensive
*/
// make |M(i,j)==0| and |M(i,k)>0| by operations with columns |j| and |k|
// precondition |M(i,k)>0|
template<typename C>
void column_clear(matrix::Matrix<C>& M, size_t i, size_t j, size_t k,
		  matrix::Matrix<C>& rec)
{
  do
  {
    long q = -intutils::divide(M(i,j),(unsigned long)(M(i,k)));
    M.columnOperation(j,k,q); // now |M(i,j)>=0|
    rec.columnOperation(j,k,q);
    if (M(i,j)==C(0))
      return;
    q=M(i,k)/M(i,j);  // |operator/| is safe now
    M.columnOperation(k,j,-q);
    rec.columnOperation(k,j,-q);
  }
  while (M(i,k)!=C(0));
  M.swapColumns(j,k);
  rec.swapColumns(j,k);
}

// make |M(i,j)==0| and |M(k,j)>0| by operations with rows |i| and |k|
// precondition |M(k,j)>0|
template<typename C>
void row_clear(matrix::Matrix<C>& M, size_t i, size_t j, size_t k,
	       matrix::Matrix<C>& rec)
{
  do
  {
    long q = -intutils::divide(M(i,j),(unsigned long)(M(k,j)));
    M.rowOperation(i,k,q); // now |M(i,j)>=0|
    rec.rowOperation(i,k,q);
    if (M(i,j)==C(0))
      return;
    q=-(M(k,j)/M(i,j));
    M.rowOperation(k,i,q); // |operator/| is safe now
    rec.rowOperation(k,i,q);
  }
  while (M(k,j)!=C(0));
  M.swapRows(i,k);
  rec.swapRows(i,k);
}

/* transform |M| to column echelon form
   postcondition: |M| unchanged column span, |result.size==M.numColumns|, and
   for all |j| and |i=result.n_th(j)|, |M(i,j)>0| and |M(ii,j)==0| for |ii>i|
*/
template<typename C>
  bitmap::BitMap column_echelon(matrix::Matrix<C>& M)
{
  bitmap::BitMap result(M.numRows());
  for (size_t j=M.numColumns(); j-->0;)
  { size_t i=M.numRows();
    while (i-->0)
      if (M(i,j)!=C(0))
      {
	if (result.isMember(i))
	{
	  size_t k=j+1+result.position(i); // index of column with pivot at |i|
	  column_clear(M,i,j,k); // now column |j| is empty in row |i| and below
	}
	else // now column |j| will have pivot in row |i|
	{
          if (M(i,j)<0)
	    M.changeColumnSign(j);
	  result.insert(i);
	  size_t p=result.position(i);
	  if (p>0) // then move column |j| to the right |p| places
	  {
	    matrix::Vector<C> c=M.column(j);
	    for (size_t d=j; d<j+p; ++d)
	      M.set_column(d,M.column(d+1));
            M.set_column(j+p,c);
	  }
	  break; // from |while| loop; done with decreasing |i|
	}
      } // |if (M(i,j)!=0)| and |while (i-->0)|
    if (i==size_t(~0)) // no pivot found for column |j|; forget about it
      M.eraseColumn(j); // and no bit is set in |result| for |j| now
  }

  return result;
}

/* find invertible |row|, |col| such that $row*M*col$ is diagonal, and
   return diagonal entries.
*/
template<typename C>
matrix::Vector<C> diagonalise(matrix::Matrix<C> M, // by value
			      matrix::Matrix<C>& row,
			      matrix::Matrix<C>& col)
{
  size_t m=M.numRows();
  size_t n=M.numColumns(); // in fact start of known null columns

  matrix::Matrix<C>(m).swap(row); // intialise |row| to identity matrix
  matrix::Matrix<C>(n).swap(col);
  matrix::Vector<C> diagonal(intutils::min(m,n),C(0));

  for (size_t d=0; d<m and d<n; ++d)
  {
    while (M.column(d).isZero()) // ensure column |d| is nonzero
    {
      --n;
      if (d==n) // then remainder of matrix is zero
	return diagonal;
      M.swapColumns(d,n); // now matrix is zero from column |n| on
      col.swapColumns(d,n);
    }

    { // get nonzero entry from column |d| at (d,d), and make it positive
      size_t i=d;
      while (M(i,d)==C(0)) // guaranteed to terminate
	++i;
      if (i>d)
      {
	C u = M(i,d)>0 ? 1 : -1;
	M.rowOperation(d,i,u);
	row.rowOperation(d,i,u);
      }
      else if (M(d,d)<0)
      {
	M.changeRowSign(d);
	row.changeRowSign(d);
      }
    }

    for (size_t i=d+1; i<m; ++i)
      row_clear(M,i,d,d,row); // makes |M(d,d)==gcd>0| and |M(i,d)==0|

    // initial sweep is unlikely to be sufficient, so no termination test here

    bool clear;
    do // sweep row and column alternatively with |M(d,d)| until both cleared
    {
     for (size_t j=d+1; j<n; ++j)
       column_clear(M,d,j,d,col); // makes |M(d,d)==gcd>0| and |M(d,j)==0|

     clear=true;
     for (size_t i=d+1; i<m; ++i)
       if (M(i,d)!=C(0))
       { clear=false; break; }
     if (clear) // most likely because |column_clear| left column |d| unchanged
       break;

     for (size_t i=d+1; i<m; ++i)
       row_clear(M,i,d,d,row); // makes |M(d,d)==gcd>0| and |M(i,d)==0|

     clear=true;
     for (size_t j=d+1; j<n; ++j)
       if (M(d,j)!=C(0))
       { clear=false; break; }
    }
    while(not clear); // then apparently some |row_clear| changed row |d|

    diagonal[d] = M(d,d);
  } // |for d|

  return diagonal;
}

 // instantiations
typedef latticetypes::LatticeCoeff T;

template
void column_clear(matrix::Matrix<T>& M, size_t i, size_t j, size_t k);
template
void row_clear(matrix::Matrix<T>& M, size_t i, size_t j, size_t k);

template
void column_clear(matrix::Matrix<T>& M, size_t i, size_t j, size_t k,
		  matrix::Matrix<T>& rec);
template
void row_clear(matrix::Matrix<T>& M, size_t i, size_t j, size_t k,
	       matrix::Matrix<T>& rec);
template
bitmap::BitMap column_echelon<T>(matrix::Matrix<T>& M);

template
matrix::Vector<T> diagonalise(matrix::Matrix<T> M, // by value
			      matrix::Matrix<T>& row,
			      matrix::Matrix<T>& col);

template // an abomination due to |abelian::Endomorphism|
matrix::Vector<unsigned long>
   diagonalise(matrix::Matrix<unsigned long> M, // by value
	       matrix::Matrix<unsigned long>& row,
	       matrix::Matrix<unsigned long>& col);

} // |namespace matreduc|
} // |namespace atlas|
