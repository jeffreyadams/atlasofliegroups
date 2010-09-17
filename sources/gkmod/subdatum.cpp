/*
  This is subdatum.cpp.

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "subdatum.h"
#include "rootdata.h"
#include "realredgp.h"
#include "kgb.h"

#include <cassert>

namespace atlas {

namespace subdatum {

SubSystem::SubSystem(const rootdata::RootDatum& parent,
		     const rootdata::RootList& sub_sys)
  : rootdata::RootSystem(parent.cartanMatrix(sub_sys).transposed())
  , rd(parent)
  , sub_W(rootdata::RootSystem::cartanMatrix())
  , pos_map(numPosRoots(),~0)
  , inv_map(rd.numRoots()+1,~0) // one spare entry for "unfound root in parent"
  , sub_root(sub_sys.size())
{
  for (weyl::Generator i=0; i<sub_sys.size(); ++i)
  {
    rootdata::RootNbr alpha = pos_map[i]=sub_sys[i];
    inv_map[alpha] = simpleRootNbr(i);
    inv_map[rd.rootMinus(alpha)] = numPosRoots()-1-i;

    size_t count=0; weyl::Generator s;
    while (alpha!=rd.simpleRootNbr(s=rd.find_descent(alpha)))
    {
      rd.simple_reflect_root(alpha,s);
      ++count;
    }

    sub_root[i].to_simple.resize(count);
    sub_root[i].reflection.resize(2*count+1);
    sub_root[i].simple=sub_root[i].reflection[count]=s; // set middle letter

    size_t j=count; // redo search loop, now storing values
    for (alpha=sub_sys[i]; j-->0; rd.simple_reflect_root(alpha,s))
    {
      s=rd.find_descent(alpha);
      sub_root[i].to_simple[j]=s;
      sub_root[i].reflection[count+1+j]=sub_root[i].reflection[count-1-j]=s;
    }
    assert(alpha==rd.simpleRootNbr(sub_root[i].simple));
  }

  for (unsigned int i=rank(); i<numPosRoots(); ++i)
  {
    rootdata::RootNbr alpha = posRootNbr(i); // root number in subsystem
    weyl::Generator s = find_descent(alpha);
    simple_reflect_root(alpha,s);
    rootdata::RootNbr beta = pos_map[posRootIndex(alpha)]; // in parent
    pos_map[i] = rd.permuted_root(sub_root[s].reflection,beta);
    inv_map[pos_map[i]] = posRootNbr(i);
    inv_map[rd.rootMinus(pos_map[i])] = numPosRoots()-1-i;
  }
}

SubSystem SubSystem::integral // pseudo contructor for integral system
  (const rootdata::RootDatum& parent, const latticetypes::RatWeight& gamma)
{
  latticetypes::LatticeCoeff n=gamma.denominator();
  const latticetypes::Weight& v=gamma.numerator();
  rootdata::RootSet int_roots(parent.numRoots());
  for (size_t i=0; i<parent.numPosRoots(); ++i)
    if (v.dot(parent.posCoroot(i))%n == 0)
      int_roots.insert(parent.posRootNbr(i));

  // it suffices that simple basis computed below live until end of constructor
  return SubSystem(parent,parent.simpleBasis(int_roots));
}

rootdata::RootNbr SubSystem::parent_nr(rootdata::RootNbr alpha) const
{
  return isPosRoot(alpha) ? pos_map[posRootIndex(alpha)]
    : parent_datum().rootMinus(pos_map[numPosRoots()-1-alpha]) ;
}

prerootdata::PreRootDatum SubSystem::pre_root_datum() const
{
  latticetypes::WeightList roots(rank()),coroots(rank()); // rank of subsystem
  for (weyl::Generator i=0; i<rank(); ++i)
  {
    roots[i]  =rd.root  (parent_nr_simple(i));
    coroots[i]=rd.coroot(parent_nr_simple(i));
  }
  return prerootdata::PreRootDatum(roots,coroots,rd.rank());
}

weyl::Twist SubSystem::twist(const latticetypes::LatticeMatrix& theta,
			     weyl::WeylWord& ww) const
{
  rootdata::RootList Delta(rank());
  for (weyl::Generator i=0; i<rank(); ++i)
  {
    rootdata::RootNbr image =
      inv_map[rd.rootNbr(theta*rd.root(parent_nr_simple(i)))];
    assert(image < numRoots());
    Delta[i] = rootMinus(image); // minus |theta| image of |root(i)|
  }

  weyl::WeylWord wrt = // this is ordered for composition on parent side
    rootdata::wrt_distinguished(*this,Delta); // make |Delta| distinguished

  // twist |Delta| gives |sub|-side fundamental, |theta=parent(wrt).-^Delta|

  weyl::Twist result;
  for (weyl::Generator i=0; i<rank(); ++i)
    result[i] = simpleRootIndex(Delta[i]);

  // we want |ww| such that |-^theta=sub(ww).Delta|; so |ww=twist(wrt^{-1})|
  // the following inversion and twist are neutral if |ww| is twisted involution
  ww.resize(wrt.size()); // we must reverse and twist, to express on our side
  for (size_t i=0; i<wrt.size(); ++i)
    ww[wrt.size()-1-i] = result[wrt[i]];
  return result;
}

// Here we seek twist and |ww| on parent side (dual with respect to |sub|)
weyl::Twist SubSystem::parent_twist(const latticetypes::LatticeMatrix& theta,
				    weyl::WeylWord& ww) const
{
  rootdata::RootList Delta(rank());
  for (weyl::Generator i=0; i<rank(); ++i)
  {
    rootdata::RootNbr image =
      inv_map[rd.rootNbr(theta*rd.root(parent_nr_simple(i)))];
    assert(image < numRoots());
    Delta[i] = image; // minus |theta| image of |root(i)|
  }

  weyl::WeylWord wrt = // this is ordered for composition on parent side
    rootdata::wrt_distinguished(*this,Delta); // make |Delta| distinguished

  weyl::Twist result;
  for (weyl::Generator i=0; i<rank(); ++i)
    result[i] = simpleRootIndex(Delta[i]);

  ww.resize(wrt.size());
  for (size_t i=0; i<wrt.size(); ++i)
    ww[i] = wrt[i];
  return result;
}

latticetypes::LatticeMatrix SubSystem::action_matrix(const weyl::WeylWord& ww)
 const
{
  latticetypes::LatticeMatrix result(rd.rank());
  for (size_t i=0; i<ww.size(); ++i)
    result *= rd.root_reflection(parent_nr_simple(ww[i]));
  return result;
}


gradings::Grading SubSystem::induced(gradings::Grading base_grading) const
{
  base_grading.complement(rd.semisimpleRank()); // flag compact simple roots
  gradings::Grading result;
  for (weyl::Generator s=0; s<rank(); ++s)
  {
    gradings::Grading mask (rd.root_expr(parent_nr_simple(s)));
    result.set(s,mask.dot(base_grading));
  }

  result.complement(rank());
  return result;
}

  /*

        class |SubDatum|


  */

SubDatum::SubDatum(realredgp::RealReductiveGroup& GR,
		   const latticetypes::RatWeight& gamma,
		   kgb::KGBElt x)
  : SubSystem(SubSystem::integral(GR.rootDatum(),gamma))
  , base_ww()
  , delta(GR.complexGroup().involutionMatrix(GR.kgb().involution(x)))
  , Tg(static_cast<const SubSystem&>(*this),delta,base_ww)
  , ini_tw()
{ // we must still modify |delta| by |ww| into the distinguished involution
  for (size_t i=0; i<base_ww.size(); ++i)
    delta *= parent_datum().root_reflection(parent_nr_simple(base_ww[i]));

  ini_tw = Tg.weylGroup().element(base_ww); // store relative word compactly

  const rootdata::RootDatum& pd = parent_datum();
  rootdata::RootList simple_image(pd.semisimpleRank());
  for (weyl::Generator s=0; s<simple_image.size(); ++s)
    simple_image[s] = pd.rootNbr(delta*pd.simpleRoot(s));

  base_ww = rootdata::wrt_distinguished(pd,simple_image);
}

latticetypes::LatticeMatrix SubDatum::involution(weyl::TwistedInvolution tw)
  const
{
  latticetypes::LatticeMatrix theta = delta;
  weyl::WeylWord ww = Weyl_group().word(tw);
  for (size_t i=ww.size(); i-->0; )
    theta.leftMult(parent_datum().root_reflection(parent_nr_simple(ww[i])));
  return theta;
}

} // |namespace subdatum|

} // |namespace atlas|
