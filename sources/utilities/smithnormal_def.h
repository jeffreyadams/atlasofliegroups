/*!
\file
  This is smithnormal_def.h.  This file contains a straightforward
  implementation of the smith normal form algorithm for integer-type
  matrices. We have implemented it as a template to get a little bit
  of flexibility in the choice of coefficients.

  Currently the intention is to use it for Cartan matrices and matrices of
  involutions, so the computational burden should be negligible. Of course
  for more serious uses one should look at specialized libraries, or think
  a lot more.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "matrix.h"

/*****************************************************************************

  This file contains a straightforward implementation of the smith normal form
  algorithm for integer-type matrices. We have implemented it as a template
  to get a little bit of flexibility in the choice of coefficients.

  Currently the intention is to use it for Cartan matrices and matrices of
  involutions, so the computational burden should be negligible. Of course
  for more serious uses one should look at specialized libraries, or think
  a lot more.

******************************************************************************/

namespace atlas {

/*****************************************************************************

        Chapter I --- Functions declared in matrix.h

  The following functions are defined :

    - addMultiple(a,b,c) : adds c times the vector b to the vector a;
    - blockReduce(m) : transforms m when it is already in block-shape, but
      there are further reductions in the lower block;
    - blockShape(b,m) : carries out the elementary operations which isolate
      the upper left block;
    - columnReduce(m) : auxiliary to reduce, handles the case of a column
      operation;
    - findReduction(m) : finds the reduction point;
    - hasReduction(m) : tells if there is a reduction point;
    - rowReduce(b,m) : auxiliary to reduce, handles the case of a row
      operation;
    - reduce(b,m) : does the reduction found by hasReduction();
    - smithNormal(invf,b,m) : a naive implementation of Smith normal form;
    - smithNormal(invf,b,b1) : another interface to smithNormal;
    - smithStep(invf,b,m) : carries out one step in the main loop of the
      Smith normal form algorithm;

******************************************************************************/

namespace smithnormal {


/*!
  Adds c times b to a. It is assumed that the size of a is at least the
  size of b.
*/
template<typename C>
void addMultiple(matrix::Vector<C>& a, const matrix::Vector<C>& b, const C& c)
{
  for (size_t j = 0; j < b.size(); ++j)
    a[j] += c*b[j];
}


/*!
  Assuming that hasBlockReduction has returned true, this function adds to
  column zero the column of m containing an element not divisible by m(0,0).
*/
template<typename C> void blockReduce(matrix::Matrix<C>& m)
{
  typedef typename matrix::Matrix<C>::index_pair IP;

  IP k = findBlockReduction(m);
  m.columnOperation(0,k.second,1);
}


/*!
  Finishes off the off-diagonal zeroeing of the first line and the first
  column of m. It is assumed that hasReduction(m) returns false, so that
  m(0,0) divides all entries in the first row and column.
*/
template<typename C>
  void blockShape(typename std::vector<matrix::Vector<C> >::iterator b,
		  matrix::Matrix<C>& m)
{
  // zero first row

  for (size_t j = 1; j < m.rowSize(); ++j) {
    if (m(0,j) == 0)
      continue;
    C c = m(0,j)/m(0,0);
    m.columnOperation(j,0,-c);
  }

  // zero first column

  for (size_t j = 1; j < m.columnSize(); ++j) {
    if (m(j,0) == 0)
      continue;
    C c = m(j,0)/m(0,0);
    m(j,0) = 0; // don't need row operation here
    addMultiple(b[0],b[j],c);
  }
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template<typename C>
  void columnReduce(matrix::Matrix<C>& m, size_t j)
{
  C a = intutils::divide(m(0,j),m(0,0));
  m.columnOperation(j,0,-a);
  m.swapColumns(0,j);
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template <typename C>
  typename matrix::Matrix<C>::index_pair
findBlockReduction(const matrix::Matrix<C>& m)
{
  for (size_t j = 1; j < m.rowSize(); ++j)
    for (size_t i = 1; i < m.columnSize(); ++i)
      if (m(i,j)%m(0,0))
	return std::make_pair(i,j);

  return std::make_pair(0,0); // this should never be reached
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template <typename C>
  typename matrix::Matrix<C>::index_pair
findReduction(const matrix::Matrix<C>& m)
{
  for (size_t j = 1; j < m.rowSize(); ++j)
    if (m(0,j)%m(0,0))
      return std::make_pair(0,j);

  for (size_t i = 1; i < m.columnSize(); ++i)
    if (m(i,0)%m(0,0))
      return std::make_pair(i,0);

  return std::make_pair(0,0); // this should never be reached
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template<typename C>
bool hasBlockReduction(const matrix::Matrix<C>& m)
{
  if (m(0,0) == 1)
    return false;

  for (size_t j = 1; j < m.rowSize(); ++j)
    for (size_t i = 1; i < m.columnSize(); ++i)
      if (m(i,j)%m(0,0))
	return true;

  return false;
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template<typename C>
bool hasReduction(const matrix::Matrix<C>& m)
{
  if (m(0,0) == 1)
    return false;

  for (size_t j = 1; j < m.rowSize(); ++j)
    if (m(0,j)%m(0,0))
      return true;

  for (size_t i = 1; i < m.columnSize(); ++i)
    if (m(i,0)%m(0,0))
      return true;

  return false;
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template<typename C>
  void reduce(typename std::vector<matrix::Vector<C> >::iterator b,
	      matrix::Matrix<C>& m)
{
  typedef typename matrix::Matrix<C>::index_pair IP;

  IP k = findReduction(m);

  if (k.first) { // row operation
    rowReduce(b,m,k.first);
  }
  else { // column operation
    columnReduce(m,k.second);
  }
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template<typename C>
  void prepareMatrix(typename std::vector<matrix::Vector<C> >::iterator b,
		     matrix::Matrix<C>& m)
{
  typedef typename matrix::Matrix<C>::index_pair IP;

  IP k = m.absMinPos();

  // transfer k to (0,0) and change sign if necessary

  if (k.first) { // swap rows
    m.swapRows(0,k.first);
    b[0].swap(b[k.first]);
  }
  if (k.second) { // swap columns
    m.swapColumns(0,k.second);
  }
  if (m(0,0) < 0) { // change sign of column
    m.changeColumnSign(0);
  }
}


/*!
  Does the reduction in the case of a row operation. The reduction consists
  in subtracting from row i the multiple of row 0 which will leave in
  m(i,0) the remainder of the Euclidian division of m(i,0) by m(0,0), and
  then swapping the two rows. One also has to update the basis b accordingly.
*/
template<typename C>
  void rowReduce(typename std::vector<matrix::Vector<C> >::iterator b,
		 matrix::Matrix<C>& m, size_t i)
{
  C a = intutils::divide(m(i,0),m(0,0));
  m.rowOperation(i,0,-a);
  addMultiple(b[0],b[i],a);
  m.swapRows(0,i);
  b[0].swap(b[i]);
}

/*!
  This is a simple implementation of the Smith normal form algorithm, intended
  for use on small matrices such as Cartan matrices.

  It is assumed that the matrix m is given, whose columns express the
  coordinates of a set of vectors relative to the basis b. This relation
  between m and b is preserved throughout the algorithm; the group generated
  by the column combinations does not change. The basis b and the matrix m
  change, until m reaches a diagonal form.

  The invariant factors are appended to invf.
*/
template<typename C>
void smithNormal(std::vector<C>& invf,
		 typename std::vector<matrix::Vector<C> >::iterator b,
		 const matrix::Matrix<C>& m)
{
  // initialize

  matrix::Matrix<C> ml(m);

  // while ml is non-zero, get an additional element in the basis, and an
  // additional invariant factor (this should write down somewhere the new
  // basis in terms of the old one.)

  typename std::vector<matrix::Vector<C> >::iterator bp = b;
  C a = 1;

  while (!ml.isZero()) {
    smithStep(bp,ml);
    C c = ml(0,0);
    ml.eraseRow(0);
    ml.eraseColumn(0);
    ml /= c;
    a *= c;
    invf.push_back(a);
    ++bp;
  }
}

template<typename C>

/*!
  Another interface to the smithNormal procedure, where we are given the
  sublattice through a generating family f, expressed by its coordinates
  in b.
*/
void smithNormal(std::vector<C>& invf,
		 typename std::vector<matrix::Vector<C> >::iterator b,
		 const std::vector<matrix::Vector<C> >& f)
{
  matrix::Matrix<C> m(f);
  smithNormal(invf,b,m);
}


/*!
  This function carries out the main loop in the Smith normal form algorithm.
  Assuming m is non-zero, it finds one additional invariant factor, and a
  corresponding basis vector.
*/
template<typename C>
void smithStep(typename std::vector<matrix::Vector<C> >::iterator b,
	       matrix::Matrix<C>& m)
{

  while(1) {

    // swap rows and columns and change signs if necessary to have upper left
    // element to be > 0 and smallest in absolute value

    prepareMatrix(b,m);

    // check divisibility in row and column

    while (hasReduction(m)) {
      reduce(b,m);
    }

    // zero off-diagonal entries in first row and first column

    blockShape(b,m);

    // reduce remaining block

    if (!hasBlockReduction(m)) // done
      break;

    blockReduce(m);

  }
}

} // namespace smithnormal

} // namespace atlas
