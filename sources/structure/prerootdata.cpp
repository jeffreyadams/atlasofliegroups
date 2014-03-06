/*!
\file
  \brief Constructs a root datum from user interaction: implementation.

  The idea is to construct an abstract root datum (a lattice and a
  subset of "roots," together with the dual lattice and a subset of
  "coroots") specified interactively as a product of simple Lie types,
  then dividing by a specified subgroup of the center of a simply
  connected form.
*/
/*
  This is prerootdata.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include <string> // used implicitly in throwing |std::runtime_error|

#include "prerootdata.h"

#include <cassert>
#include <stdexcept>

#include "matrix.h"
#include "lietype.h"

namespace atlas {

namespace {

  WeightList rootBasis(const LieType& lt, const WeightList&);
  CoweightList corootBasis(const LieType& lt, const CoweightList& lb);

}

/*****************************************************************************

        Chapter I -- The PreRootDatum class

  The PreRootDatum class is here just so that the ingredients for building
  a RootDatum can be gotten conveniently in various ways (for instance,
  interactively). It contains the rank, and bases for the roots and the
  coroots, expressed in the current lattice.

******************************************************************************/

namespace prerootdata {

  // constructor


/*!
\brief Constructs the PreRootDatum whose lattice has basis |b|,
expressed in terms of the simply connected weight lattice basis for |lt|.

More precisely, we build the unique rootdatum whose Cartan matrix is the
standard (Bourbaki) Cartan matrix of the semisimple part of the type |lt|
(i.e., the Cartan matrix of |lt| with any null rows and columns for torus
factors removed), and such that the the vectors |b| express the basis of the
weight lattice $X$ in terms of the fundamental weight basis (dual basis to the
simple coroots mixed with standard basis for the torus factors) of |lt|.

This somewhat convoluted description comes from the way the lattice $X$ may be
chosen (via user interaction) as a sublattice of the lattice for a simply
connected group.

The constructor puts in |d_roots| the list of simple roots expressed in the
basis |b|, and in |d_coroots| the list of simple coroots expressed in the
dual basis.
*/
PreRootDatum::PreRootDatum(const LieType& lt,
			   const WeightList& b)
: d_roots(rootBasis(lt,b)), d_coroots(corootBasis(lt,b)), d_rank(lt.rank())
{}


  // accessors

int_Matrix PreRootDatum::Cartan_matrix() const
{
  int_Matrix Cartan(d_roots.size(),d_coroots.size());

  for (weyl::Generator i = 0; i < d_roots.size(); ++i)
    for (weyl::Generator j = 0; j < d_coroots.size(); ++j)
      Cartan(i,j) = d_roots[i].dot(d_coroots[j]);

  return Cartan;
}

template<typename C>
  void PreRootDatum::simpleReflect(matrix::Vector<C>& v, weyl::Generator s)
  const
{ v -= d_roots[s].scaled(d_coroots[s].dot(v)); }

void PreRootDatum::simple_reflect(weyl::Generator s, LatticeMatrix& M) const
{
  assert(M.numRows()==rank());
  for (unsigned int j=0; j<M.numColumns(); ++j)
  {
    int c=0;
    for (unsigned int i=0; i<rank(); ++i)
      c+= d_coroots[s][i]*M(i,j);
    for (unsigned int i=0; i<rank(); ++i)
      M(i,j) -= d_roots[s][i]*c;
  }
}

void PreRootDatum::simple_reflect(LatticeMatrix& M,weyl::Generator s) const
{
  assert(M.numColumns()==rank());
  for (unsigned int i=0; i<M.numRows(); ++i)
  {
    int c=0;
    for (unsigned int j=0; j<rank(); ++j)
      c+= M(i,j)*d_roots[s][j];
    for (unsigned int j=0; j<rank(); ++j)
      M(i,j) -= c*d_coroots[s][j];
  }
}



  // manipulator

void PreRootDatum::swap(PreRootDatum& other)

{
  d_roots.swap(other.d_roots);
  d_coroots.swap(other.d_coroots);
  std::swap(d_rank,other.d_rank);
}

} // |namespace prerootdata|

/*****************************************************************************

        Chapter II -- Auxiliary functions

  This sections contains the definitions of some auxiliary functions private
  to the present module :

    - WeightList& rootBasis(const CartanMatrix&, const WeightList&):
      writes down the simple roots in the lattice basis;
    - WeightList corootBasis(const LieType& lt,const WeightList& lb): writes
      down the simple coroots for |lt| in the dual lattice basis of |lb|;

******************************************************************************/

namespace {


/*!
 \brief Writes down the simple roots in the lattice basis.

  Given the lattice basis |lb|, expressed in terms of the simple weight basis,
  and a Lie type |lt|, the Cartan matrix of |lt| may be interpreted as giving
  in its _rows_ the coordinates of the simple roots in the simple weight basis
  (this interpretation depends on our definition of the Cartan matrix to
  include in the non-semisimple case null columns at torus factor positions;
  in fact for symmetry null rows are inserted there as well, which we skip).
  This function returns the coordinates of the simple root basis in |lb|, so
  it is simply a base change by left multiplication by the inverse of |lb|.
*/
WeightList rootBasis(const LieType& lt, const WeightList& lb)
{
  assert(lb.size()==lt.rank());
  LatticeCoeff d;
  LatticeMatrix q(lb,lb.size()); q.invert(d);

  if (d==0)
    throw std::runtime_error("Dependent lattice generators");

  // push back simple roots expressed in |lb|

  const size_t rk=lt.rank();
  Weight v(rk); // temporary vector
  WeightList result; result.reserve(lt.semisimple_rank());
  size_t r=0; // row number in Cartan matrix
  for (size_t i=0; i<lt.size(); ++i)
    for (size_t j=0; j<lt[i].rank(); ++j,++r)
      if (lt[i].type()!='T') // skip on torus factors
      {
	for (size_t k=0; k<rk; ++k) // get Cartan entries for root
	  v[k]=lt.Cartan_entry(r,k);
	result.push_back((q*v)/d); // may |throw std::runtime_error|
      }

  return result;
}

/*! \brief Writes down the simple coroots in the dual lattice basis.

  Given the lattice basis |lb|, expressed in terms of the simple weight basis,
  and the Lie type |lt|, this function returns the simple coroots of the
  system. If |q| is the matrix of |lb| in the simple weight basis, its
  transpose is the matrix of the dual basis of the simple weight basis (a
  basis consisting of simple coroots completed with vectors in the radical) on
  the dual basis of the sublattice |lt| (which can be a basis of a _larger_
  lattice than that of the coweights). So all that is necessary is to drop
  from the transposed matrix the columns representing vectors in the radical;
  this is what |lt| is used for. In fact we avoid transposition, copying rows.
*/
CoweightList corootBasis(const LieType& lt, const WeightList& lb)
{
  assert(lb.size()==lt.rank());
  LatticeMatrix q(lb,lb.size()); // square matrix

  CoweightList result; result.reserve(lt.semisimple_rank());
  size_t r=0; // row number in |q|
  for (size_t i=0; i<lt.size(); ++i)
    for (size_t j=0; j<lt[i].rank(); ++j,++r)
      if (lt[i].type()!='T') // skip on torus factors
	result.push_back(q.row(r));

  return result;
}


} // namespace


// template instantiation

namespace prerootdata {

template void PreRootDatum::simpleReflect
  (matrix::Vector<arithmetic::Numer_t>& v, weyl::Generator s) const;

} // |namespace prerootdata|

} // namespace atlas
