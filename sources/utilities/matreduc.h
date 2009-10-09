/*
  matreduc.h

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef   	MATREDUC_H
# define   	MATREDUC_H

#include "matrix_fwd.h"
#include "bitmap_fwd.h"

namespace atlas {

namespace matreduc {

template<typename C>
  bitmap::BitMap column_echelon(matrix::Matrix<C>& vectors);

template<typename C>
matrix::Vector<C> diagonalise(matrix::Matrix<C> M, // by value
			      matrix::Matrix<C>& row,
			      matrix::Matrix<C>& col);
} // |namespace matreduc|
} // |namespace atlas|

#endif 	    /* !ECHELON_H */
