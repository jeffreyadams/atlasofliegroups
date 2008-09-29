/*!
\file
  This is smithnormal.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SMITHNORMAL_H  /* guard against multiple inclusions */
#define SMITHNORMAL_H

#include "matrix_fwd.h"

/******** function declarations **********************************************/

namespace atlas {

namespace smithnormal {

template<typename C>
  void addMultiple(matrix::Vector<C>&, const matrix::Vector<C>&, const C&);

template<typename C>
  void blockReduce(typename std::vector<matrix::Vector<C> >::iterator,
		   matrix::Matrix<C>&);

template<typename C>
  void blockShape(typename std::vector<matrix::Vector<C> >::iterator,
		  matrix::Matrix<C>&);

template<typename C>
  void columnReduce(matrix::Matrix<C>&, size_t);

template <typename C>
  typename matrix::Matrix<C>::index_pair
  findBlockReduction(const matrix::Matrix<C>&);

template <typename C>
  typename matrix::Matrix<C>::index_pair
  findReduction(const matrix::Matrix<C>&);

template<typename C>
  bool hasBlockReduction(const matrix::Matrix<C>&);

template<typename C>
  bool hasReduction(const matrix::Matrix<C>&);

template<typename C>
  void prepareMatrix(typename std::vector<matrix::Vector<C> >::iterator,
		     matrix::Matrix<C>&);

template<typename C>
  void reduce(typename std::vector<matrix::Vector<C> >::iterator,
	      matrix::Matrix<C>&);

template<typename C>
  void rowReduce(typename std::vector<matrix::Vector<C> >::iterator,
		 matrix::Matrix<C>&, size_t);

template<typename C>
  void smithNormal(std::vector<C>&,
		   typename std::vector<matrix::Vector<C> >::iterator,
		   const matrix::Matrix<C>&);

template<typename C>
  void smithNormal(std::vector<C>&,
		   typename std::vector<matrix::Vector<C> >::iterator,
		   const std::vector<matrix::Vector<C> >&);

template<typename C>
  void smithStep(typename std::vector<matrix::Vector<C> >::iterator,
		 matrix::Matrix<C>&);

}

}

#include "smithnormal_def.h"

#endif
