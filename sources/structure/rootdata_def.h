/*!
\file
\brief Template definitions for the class RootData.
*/
/*
  This is rootdata_def.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "lattice.h"

/******** template definitions from rootdata.h *******************************/

namespace atlas {

namespace rootdata {

/*!
  In this template we assume that |I| is an |InputIterator| with |value_type|
  |RootNbr|, and that |O| is an |OutputIterator| with |value_type| |Weight|.
  We output via |out| the expressions of the roots numbered in |[first,last[|
  in the simple root basis.
*/
template <typename I, typename O>
  void RootDatum::toRootBasis(const I& first, const I& last, O out) const
{
  for (I i = first; i != last; ++i) {
    *out++ = inSimpleRoots(*i);
  }
}


/*!
  In this template we assume that |I| is an |Input_iterator| with |value_type|
  |RootNbr|, and that |O| is an |OutputIterator| with |value_type| |Weight|.
  We assume that |rb| contains a basis of some _sub_ rootsystem of |rd|, and
  that |I| inputs |RootNbr|s corresponding to roots in that subsystem. Then we
  write to |out| the expression of the root in the basis |rb|.

  The idea is to use the coroots of the vectors in |rb| to get the expression
  of the both the input roots and those from |rb| in the simple weight basis
  for |rb| (this is done by |toSimpleWeights| below); the latter matrix (the
  Cartan matrix of |rb|) being square, we can then apply |baseChange|.

  We could have done this more directly if |d_coweights| had been available
  for arbitrary roots, but the are only stored for simple roots.
*/
template <typename I, typename O>
  void RootDatum::toRootBasis(const I& first, const I& last, O out,
			      const RootList& rb) const
{
  latticetypes::WeightList wl; // roots |[first,last[| in weight basis of |rb|
  latticetypes::WeightList wb; // (square) Cartan matrix of |rb|

  toSimpleWeights(first,last,back_inserter(wl),rb);
  toSimpleWeights(rb.begin(),rb.end(),back_inserter(wb),rb);

  lattice::baseChange(wl.begin(),wl.end(),out,wb.begin(),wb.end());
}


/*!
  In this template we assume that |I| is an |InputIterator| with value type
  |RootNbr|, and that |rb| is a basis for some _sub_ rootsystem of |rd|; we
  assume that |I| inputs values corresponding to roots in the subsystem.
  Then for each |v| in |[first,last[| we output to |o| the expression of |v|
  in the simple weight basis of the root subsystem |rb|; this is obtained
  simply by taking scalar products with the coroots of the roots in |rb|.
*/
template <typename I, typename O>
  void RootDatum::toSimpleWeights(const I& first, const I& last, O out,
				  const RootList& rb) const
{
  latticetypes::Weight v(rb.size());

  for (I i = first; i != last; ++i) {
    for (unsigned long j = 0; j < rb.size(); ++j)
      v[j] = root(*i).scalarProduct(coroot(rb[j]));
    *out++ = v;
  }
}

} // namespace rootdata

} // namespace atlas
