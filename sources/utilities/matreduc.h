/*
  matreduc.h

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef   	MATREDUC_H
# define   	MATREDUC_H

#include <vector>
#include "matrix_fwd.h"
#include "matrix.h"
#include "bitmap.h"
#include "bigint.h"
#include "sl_list.h"

#include <stdexcept>

namespace atlas {

namespace matreduc {

// compute |d=gcd| of all entries of |row|
// also if |col!=nullptr| assign matrix of col operations for |row->[d,0,0,...]|
template<typename C>
C gcd (matrix::Vector<C> row, matrix::PID_Matrix<C>* col, bool& flip,
       size_t dest=0); // place where |d| must end up, ignored if |col==nullptr|

// find reduced basis for span of |vectors| and column operations to get there
template<typename C>
  bitmap::BitMap column_echelon(matrix::PID_Matrix<C>& vectors,
				matrix::PID_Matrix<C>& col,
				bool& flip);

// solve $E*x=b*f$, with $f>0$ minimal where |E| is echelon with |pivots|
// non-pivot rows may cause solving to fail, then throw |std::runtime_error|
template<typename C>
  matrix::Vector<C> echelon_solve(const matrix::PID_Matrix<C>& E,
				  const bitmap::BitMap& pivots,
				  matrix::Vector<C> b,
				  arithmetic::big_int& f);

// find |row,col|, with $\det=1$ making |row*M*col| (the returned) diagonal
// WARNING first entry (only) of result might be negative (rest is positive)
template<typename C>
  std::vector<C> diagonalise(matrix::PID_Matrix<C> M, // by value
			     matrix::PID_Matrix<C>& row,
			     matrix::PID_Matrix<C>& col);
template<typename C>
  matrix::PID_Matrix<C> adapted_basis(matrix::PID_Matrix<C> M, // by value
				      std::vector<C>& diagonal);
template<typename C>
  matrix::PID_Matrix<C> Smith_basis(const matrix::PID_Matrix<C>& M,
				    std::vector<C>& diagonal);

template<typename C> // find a solution |x| for |A*x==b|
  bool has_solution(const matrix::PID_Matrix<C>& A,
		    matrix::Vector<C> b); // by value
template<typename C> // find a solution |x| for |A*x==b|
  matrix::Vector<C> find_solution(const matrix::PID_Matrix<C>& A,
				  matrix::Vector<C> b); // by value

// inline implementations

template<typename C>
C gcd (matrix::Vector<C> row, matrix::PID_Matrix<C>* col,bool& flip,
       size_t dest)
{ if (col!=nullptr)
    *col = matrix::PID_Matrix<C>(row.size());
  containers::sl_list<size_t> active_entries;
  C min(0); size_t mindex;
  for (size_t j=0; j<row.size(); ++j)
    if (row[j]!=C(0))
    {
      active_entries.push_back(j);
      if (min==C(0) or abs(row[j])<min)
	min=abs(row[mindex=j]);
    }
  if (active_entries.empty())
    return C(0);
  if (row[mindex]<C(0))
  {
    row[mindex]=-row[mindex];
    flip = not flip;
    if (col!=nullptr)
      (*col)(mindex,mindex)=C(-1);
  }

  while (not active_entries.singleton())
  { const size_t cur_col = mindex;
    const C d = row[mindex];
    for (auto it=active_entries.begin(); not active_entries.at_end(it); )
      if (*it==cur_col)
	++it;
      else
      { auto j=*it;
	C q=arithmetic::divide(row[j],d);
	if (col!=nullptr)
	  col->columnOperation(j,cur_col,-q);
	if ((row[j] -= d*q)==C(0))
	  active_entries.erase(it); // this leaves |it| pointing to next
	else
	{ if (row[j]<min)
	    min=row[mindex=j];
	  ++it; // don't forget to increment in this case
	}
      }
    assert(active_entries.singleton() or mindex!=cur_col);
  }

  if (col!=nullptr and mindex!=dest)
  {
    col->swapColumns(dest,mindex);
    flip=not flip;
  }

  return min;
}

/* transform |M| to column echelon form using only PID operations
   postcondition: |M| has unchanged column span, |result.size==M.numColumns|,
   and for |j| and |i=result.n_th(j)|: |M(i,j)>0| and |M(ii,j)==0| for |ii>i|
*/
template<typename C>
  bitmap::BitMap column_echelon(matrix::PID_Matrix<C>& M,
				matrix::PID_Matrix<C>& col,
				bool& flip)
{ using std::abs;
  const size_t n=M.numColumns();
  col=matrix::PID_Matrix<C>(n); // start with identity matrix
  matrix::PID_Matrix<C> ops; // working matrix, accumulates column operations
  flip=false;
  bitmap::BitMap result(M.numRows()); // set of pivot rows found so far
  size_t l=n; // limit of columns yet to consider
  for (size_t i=M.numRows(); i-->0; )
  { int d=gcd(M.partial_row(i,0,l),&ops,flip,l-1);
    assert(ops.numRows()==l and ops.numColumns()==l);
    if (d==0)
      continue; // if partial row was already zero, just skip over current row
    matrix::column_apply(M,ops,0);
    matrix::column_apply(col,ops,0);
    result.insert(i);
    --l; // now we have a pivot in column |l-1|
    assert(M(i,l)==d); // |column_apply| should have achieved this
  } // |for(i)

  while (l-->0)
  { // then erase column |l| from |M|, rotate it in |col| towards the right
    M.eraseColumn(l);
    matrix::Vector<C> cc=col.column(l); // this one too
    for (size_t j=l; j<M.numColumns(); ++j)
      col.set_column(j,col.column(j+1));
    col.set_column(M.numColumns(),cc);
    flip ^= (M.numColumns()-l)%2;
  }

  return result;
} // |column_echelon|

// when |E| is echolon with |pivots|, the following solves by back-substitution
template<typename C>
  matrix::Vector<C> echelon_solve(const matrix::PID_Matrix<C>& E,
				  const bitmap::BitMap& pivots,
				  matrix::Vector<C> b,
				  arithmetic::big_int& f) // needed scale factor
{ assert(b.size()==E.numRows());
  using arithmetic::gcd;
  f=arithmetic::big_int(1);
  matrix::Vector<C> result(E.numColumns());
  size_t j=pivots.size();
  for (size_t i=E.numRows(); i-->0; )
    if (pivots.isMember(i))
    {
      --j;
      assert(E(i,j)>C(0)); // since it is a pivot
      const C d=gcd(b[i],E(i,j)); // we ensure a positive second argument
      assert(d>C(0));
      const C m = b[i]/d; // factor for column |j| in upcoming subtraction
      if (d<E(i,j)) // then division is not exact
      {
	const C q = E(i,j)/d;
	f *= q; // need to scale up |b| by an additional factor |q|
	for (size_t k=0; k<=i; ++k)
	  b[k]*=q;
	for (size_t l=j+1; l<result.size(); ++l)
	  result[l] *= q;
      }
      assert(m==b[i]/E(i,j)); // since |b[i]| now is what was |b[i]*E(i,j)/d|
      result[j] = m;
      for (size_t k=0; k<=i; ++k)
	b[k] -= E(k,j)*m; // subtract off contribution from htis column
      assert(b[i]==C(0)); // that was the point of the subtraction
    }
    else if (b[i]!=C(0))
      throw std::runtime_error("Inconsistent linear system");

  assert(j==0); // every column had its pivot, and |result| is fully defined

  return result;
} // |echelon_solve|


} // |namespace matreduc|
} // |namespace atlas|

#endif
