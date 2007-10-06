/*!
\file
\brief Template definitions for the class RootData.
*/
/*
  This is rootdata_def.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "lattice.h"

/******** template definitions from rootdata.h *******************************/

namespace atlas {

namespace rootdata {

/*!
  In this template we assume that I is an InputIterator with value_type
  Weight, and that O is an OutputIterator with the same value_type. We
  assume that the elements of [first,last[ are in the root lattice, and
  we output their expressions in the simple root basis.
*/
template <typename I, typename O>
  void RootDatum::toRootBasis(const I& first, const I& last, O out) const
{
  using namespace latticetypes;

  for (I i = first; i != last; ++i) {
    Weight v(d_semisimpleRank);
    for (size_t j = 0; j < v.size(); ++j) {
      const RatWeight& scw_j = d_coweights[j];
      v[j] = LT::scalarProduct(scw_j,*i);
    }
    *out++ = v;
  }
}


/*!
  In this template we assume that I is an Input_iterator with value_type
  RootNbr, and O an OutputIterator with value_type Weight. We assume that
  |rb| contains a basis of some _sub_ rootsystem of |rd|, and that |I| outputs
  RootNbr's corresponding to roots in that subsystem. Then we write to
  |out| the expression of the root in the basis rb.

  The idea is to use the coroots of the vectors in rb to get the expression
  of the vectors in the simple weight basis; then the Cartan matrix yields
  the expression of rb in that same basis, and it is enough to apply the
  usual baseChange.
*/
template <typename I, typename O>
  void toRootBasis(const I& first, const I& last, O out, const RootList& rb,
		   const RootDatum& rd)
{
  using namespace lattice;
  using namespace latticetypes;

  WeightList wl;
  WeightList wb;

  toSimpleWeights(first,last,back_inserter(wl),rb,rd);
  toSimpleWeights(rb.begin(),rb.end(),back_inserter(wb),rb,rd);

  baseChange(wl.begin(),wl.end(),out,wb.begin(),wb.end());
}

template <typename I, typename O>
  void toSimpleWeights(const I& first, const I& last, O out,
		       const RootList& rb, const RootDatum& rd)

/*!
  In this template we assume that I is an InputIterator with value type
  RootNbr, and that rb is a basis for some _sub_ rootsystem of rd; we
  assume that I outputs values corresponding to roots in the subsystem.
  Then for each v in [first,last[ we output to o the expression of v
  in the simple weight basis of the root system; this is obtained simply
  by taking scalar products with the coroots.
*/

{
  using namespace latticetypes;

  Weight v(rb.size());

  for (I i = first; i != last; ++i) {
    for (unsigned long j = 0; j < rb.size(); ++j)
      v[j] = scalarProduct(rd.root(*i),rd.coroot(rb[j]));
    *out++ = v;
  }

  return;
}

}

}
