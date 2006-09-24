/*!
\file
  This is smithnormal.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef SMITHNORMAL_H  /* guard against multiple inclusions */
#define SMITHNORMAL_H

#include "matrix_fwd.h"

/******** function declarations **********************************************/

namespace atlas {

namespace smithnormal {

template<typename C>
  void addMultiple(std::vector<C>&, const std::vector<C>&, const C&);

template<typename C>
  void blockReduce(typename std::vector<std::vector<C> >::iterator, 
		   matrix::Matrix<C>&);

template<typename C>
  void blockShape(typename std::vector<std::vector<C> >::iterator, 
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
  void prepareMatrix(typename std::vector<std::vector<C> >::iterator, 
		     matrix::Matrix<C>&);

template<typename C>
  void reduce(typename std::vector<std::vector<C> >::iterator, 
	      matrix::Matrix<C>&);

template<typename C>
  void rowReduce(typename std::vector<std::vector<C> >::iterator, 
		 matrix::Matrix<C>&, size_t);

template<typename C>
  void smithNormal(std::vector<C>&, 
		   typename std::vector<std::vector<C> >::iterator, 
		   const matrix::Matrix<C>&);

template<typename C>
  void smithNormal(std::vector<C>&,
		   typename std::vector<std::vector<C> >::iterator, 
		   const std::vector<std::vector<C> >&);

template<typename C>
  void smithStep(typename std::vector<std::vector<C> >::iterator, 
		 matrix::Matrix<C>&);

}

}

#include "smithnormal_def.h"

#endif
