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
#include "bitmap_fwd.h"

#include <stdexcept>

namespace atlas {

namespace matreduc {

// compute |d=gcd| of all entries of |row|
// also if |col!=nullptr| assign matrix of col operations for |row->[d,0,0,...]|
template<typename C>
  C gcd (matrix::Vector<C> row, matrix::PID_Matrix<C>* col);

// find reduced basis for span of |vectors| and column operations to get there
template<typename C>
  bitmap::BitMap column_echelon(matrix::PID_Matrix<C>& vectors,
				matrix::PID_Matrix<C>& col,
				bool& flip);

template<typename C> // find |row,col| making |row*M*col| (returned) diagonal
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

} // |namespace matreduc|
} // |namespace atlas|

#endif
