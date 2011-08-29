/*
  This is subdatum.cpp.

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "subdatum.h"

#include "rootdata.h"
#include "tits.h"
#include "complexredgp.h"
#include "realredgp.h"
#include "kgb.h"

#include <cassert>

namespace atlas {

namespace subdatum {

  /*

        class |SubDatum|


  */

SubDatum::SubDatum(RealReductiveGroup& GR,
		   const RatWeight& gamma,
		   KGBElt x)
  : SubSystem(SubSystem::integral(GR.rootDatum(),gamma))
  , base_ww()
  , delta(GR.complexGroup().involutionMatrix(GR.kgb().involution(x)))
  , Tg(static_cast<const SubSystem&>(*this),delta,base_ww)
  , ini_tw()
{ const RootDatum& pd = parent_datum();

  // We'd like to |assert| here that |gamma| is valid, but it's hard to do.

  // we must still modify |delta| by |ww| going up to sub-distinguished invol
  for (size_t i=0; i<base_ww.size(); ++i) // apply in in reverse, to the right
    delta *= pd.root_reflection(parent_nr_simple(base_ww[i]));

  ini_tw = Tg.weylGroup().element(base_ww); // so that |ini_ww*delta=theta|

  // now store in |base_ww| the parent twisted involution for |delta|
  RootNbrList simple_image(pd.semisimpleRank());
  for (weyl::Generator s=0; s<simple_image.size(); ++s)
    simple_image[s] = pd.rootNbr(delta*pd.simpleRoot(s));

  // set the value that will be accessible later as |base_twisted_in_parent()|
  base_ww = rootdata::wrt_distinguished(pd,simple_image);
#ifndef NDEBUG
  const ComplexReductiveGroup& GC=GR.complexGroup();
  // We test that the subdatum shares its most split involution with |GC|
  WeightInvolution M=delta;
  const WeylWord sub_w0 = Weyl_group().word(Weyl_group().longest());
  for (size_t i=0; i<sub_w0.size(); ++i)
    M *= pd.root_reflection(parent_nr_simple(sub_w0[i]));

  assert(M==GC.dualDistinguished().negative_transposed());
#endif
}

WeightInvolution SubDatum::matrix(TwistedInvolution tw)
  const
{
  WeightInvolution theta = delta;
  WeylWord ww = Weyl_group().word(tw);
  for (size_t i=ww.size(); i-->0; )
    theta.leftMult(parent_datum().root_reflection(parent_nr_simple(ww[i])));
  return theta;
}

} // |namespace subdatum|

} // |namespace atlas|
