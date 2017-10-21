/*
  Constructing a root datum from user interaction: implementation.

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

#include "prerootdata.h"

#include <cassert>
#include <stdexcept>

#include "matrix.h"
#include "lietype.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- The PreRootDatum class

  The PreRootDatum class is here just so that the ingredients for building
  a RootDatum can be gotten conveniently in various ways (for instance,
  interactively). It contains the rank, and bases for the roots and the
  coroots, expressed in the current lattice.

******************************************************************************/

namespace prerootdata {

  // constructor


/*
  Construct the PreRootDatum whose lattice has basis |b|, expressed in terms
  of the simply connected weight lattice basis for |lt|.

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
PreRootDatum::PreRootDatum(const LieType& lt)
  : simple_roots(lt.rank(),lt.semisimple_rank())
  , simple_coroots(simple_roots.numRows(),simple_roots.numColumns(),0)
{
  weyl::Generator s=0; // tracks (co)roots, goes up to semisimple rank
  unsigned int r=0; // |r| indexes rows of |Cartan|, goes up to rank
  for (unsigned int k=0; k<lt.size(); ++k) // run over "simple" factors
    if (lt[k].type()=='T') // only do non-torus factors;
      r+=lt[k].rank();  // skip empty row(s) of Cartan matrix, keep |s|
    else
      for (unsigned int i=lt[k].rank(); i-->0; ++r,++s) // here |s| increases
      { // set simple roots and coroots as for simply connected root datum
	for (unsigned int j=0; j<lt.rank(); ++j)
	  simple_roots(j,s) = lt.Cartan_entry(r,j); // Cartan row to column
	simple_coroots(r,s) = 1; // coroot |s| is canonical basis vector |r|
      }
}


  // accessors

int_Matrix PreRootDatum::Cartan_matrix() const
{
  const auto s = semisimple_rank();
  int_Matrix Cartan(s,s);

  for (weyl::Generator i = 0; i<s; ++i)
    for (weyl::Generator j = 0; j<s; ++j)
      Cartan(i,j) = simple_root(i).dot(simple_coroot(j));

  return Cartan;
}

// replace by root datum for a finite central quotient with weight |sublattice|
PreRootDatum& PreRootDatum::quotient(const LatticeMatrix& sublattice)
{
  const auto r=rank();
  if (sublattice.numRows()!=r or sublattice.numColumns()!=r)
    throw std::runtime_error("Sub-lattice matrix not square of the right size");

  arithmetic::big_int d;
  LatticeMatrix inv=inverse(sublattice,d);

  if (d.is_zero())
    throw std::runtime_error("Dependent lattice generators");

  const auto den = d.int_val(); // denominator; if this throws we're out of luck

  try
  {
    for (unsigned int j=0; j<semisimple_rank(); ++j)
    {
      simple_roots.set_column(j,inv*simple_root(j)/den);
      simple_coroots.set_column(j,sublattice.right_prod(simple_coroot(j)));
    }
    return *this;
  }
  catch (std::runtime_error& e)
  { // relabel |std::runtime_error("Inexact integer division")| from division
    throw std::runtime_error("Sub-lattice does not contain the root lattice");
  }
}

template<typename C>
void PreRootDatum::simple_reflect(weyl::Generator s,matrix::Vector<C>& v)
  const
{ const auto alpha = simple_root(s); // temporarily store this simple root
  v.subtract(alpha.begin(),simple_coroot(s).dot(v));
}

void PreRootDatum::simple_reflect(weyl::Generator s, LatticeMatrix& M) const
{
  assert(M.numRows()==rank());
  for (unsigned int j=0; j<M.numColumns(); ++j)
  {
    int c=0;
    for (unsigned int i=0; i<rank(); ++i)
      c+= simple_coroots(i,s)*M(i,j);
    for (unsigned int i=0; i<rank(); ++i)
      M(i,j) -= simple_roots(i,s)*c;
  }
}

void PreRootDatum::simple_reflect(LatticeMatrix& M,weyl::Generator s) const
{
  assert(M.numColumns()==rank());
  for (unsigned int i=0; i<M.numRows(); ++i)
  {
    int c=0;
    for (unsigned int j=0; j<rank(); ++j)
      c+= M(i,j)*simple_roots(j,s);
    for (unsigned int j=0; j<rank(); ++j)
      M(i,j) -= c*simple_coroots(j,s);
  }
}

} // |namespace prerootdata|

// template instantiation

namespace prerootdata {

template void PreRootDatum::simple_reflect
  (weyl::Generator s,matrix::Vector<arithmetic::Numer_t>& v) const;

} // |namespace prerootdata|

} // |namespace atlas|
