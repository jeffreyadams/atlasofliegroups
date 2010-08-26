/*
  This is subdatum.cpp.

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "subdatum.h"

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
      inv_map[rd.rootNbr(theta.apply(rd.root(parent_nr_simple(i))))];
    assert(image < numRoots());
    Delta[i] = rootMinus(image); // minus transposed |theta| image |coroot(i)|
  }

  weyl::WeylWord wrt = // this is ordered for composition on parent side
    rootdata::wrt_distinguished(*this,Delta); // make |Delta| distinguished

  weyl::Twist result;
  for (weyl::Generator i=0; i<rank(); ++i)
    result[i] = simpleRootIndex(Delta[i]);

  // the following inversion and twist are neutral if |ww| is twisted involution
  ww.resize(wrt.size()); // we must reverse and twist, to express on our side
  for (size_t i=0; i<wrt.size(); ++i)
    ww[wrt.size()-1-i] = result[wrt[i]];
  return result;
}

} // |namespace subdatum|

} // |namespace atlas|
