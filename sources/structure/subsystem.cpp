/*
  This is subsystem.cpp.

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "subsystem.h"

#include "ratvec.h"	// acces to infinitesimal character |gamma|
#include "bitmap.h"	// root sets

#include "prerootdata.h"// value returned
#include "rootdata.h"	// |RootSystem| construction and methods
#include "cartanclass.h"// |InvolutionData|
#include "weyl.h"	// subobject

#include <cassert>

namespace atlas {

namespace subsystem {

SubSystem::SubSystem(const RootDatum& parent,
		     const RootNbrList& sub_sys)
  : RootSystem(parent.cartanMatrix(sub_sys).transposed()) // build
  , rd(parent) // share
  , pos_map(numPosRoots(),~0)
  , inv_map(rd.numRoots()+1,~0) // one spare entry for "unfound root in parent"
  , sub_root(numPosRoots())
{
  for (unsigned int i=0; i<numPosRoots(); ++i)
  {
    // first set |pos_map[i]|
    if (i<sub_sys.size()) // simple in subsystem, just use |sub_sys|
      pos_map[i]=sub_sys[i];
    else // use previously stored values in |pos_map|, |sub_root| recursively
    {
      RootNbr sub_alpha = posRootNbr(i);
      weyl::Generator s = find_descent(sub_alpha); // generator for subsystem
      simple_reflect_root(s,sub_alpha);
      RootNbr beta = pos_map[posRootIndex(sub_alpha)]; // |beta| is in parent
      pos_map[i] = rd.permuted_root(sub_root[s].reflection,beta);
    }

    RootNbr alpha = pos_map[i]; // now we use parent numbering
    inv_map[alpha] = posRootNbr(i); // refers simple root |i| in subsystem
    inv_map[rd.rootMinus(alpha)] = rootMinus(inv_map[alpha]); // its negative

    // in the remainder we work in parent datum; must find conjugate to simple
    size_t count=0; weyl::Generator s;
    while (alpha!=rd.simpleRootNbr(s=rd.find_descent(alpha)))
    { // just count the reflections needed to make alpha simple
      rd.simple_reflect_root(s,alpha);
      ++count;
    }

    // now we can dimension our Weyl words, and set |sub_root[i].simple|
    sub_root[i].to_simple.resize(count);
    sub_root[i].reflection.resize(2*count+1);
    sub_root[i].simple=sub_root[i].reflection[count]=s; // set middle letter

    size_t j=count; // redo search loop, storing remaining Weyl word letters
    for (alpha=pos_map[i]; j-->0; rd.simple_reflect_root(s,alpha))
    {
      s=rd.find_descent(alpha);
      sub_root[i].to_simple[j]=s; // write |to_simple| word from right to left
      sub_root[i].reflection[count+1+j]= // and |reflection| from outside in
      sub_root[i].reflection[count-1-j]=s;
    }
    assert(alpha==rd.simpleRootNbr(sub_root[i].simple)); // check |alpha|
  }
}

SubSystem SubSystem::integral // pseudo contructor for integral system
  (const RootDatum& parent, const RatWeight& gamma)
{
  arithmetic::Numer_t n=gamma.denominator(); // signed!
  const Ratvec_Numer_t& v=gamma.numerator();
  RootNbrSet int_roots(parent.numRoots());
  for (size_t i=0; i<parent.numPosRoots(); ++i)
    if (parent.posCoroot(i).dot(v)%n == 0)
      int_roots.insert(parent.posRootNbr(i));

  // it suffices that simpleBasis computed below live until end of constructor
  return SubSystem(parent,parent.simpleBasis(int_roots));
}

RootNbr SubSystem::to_parent(RootNbr alpha) const
{
  RootNbr result = pos_map[rt_abs(alpha)];
  if (is_negroot(alpha))
    result = parent_datum().rootMinus(result);
  return result;
}

PreRootDatum SubSystem::pre_root_datum() const
{
  WeightList roots(rank());  // rank of subsystem
  CoweightList coroots(rank());
  for (weyl::Generator i=0; i<rank(); ++i)
  {
    roots[i]  =rd.root  (parent_nr_simple(i));
    coroots[i]=rd.coroot(parent_nr_simple(i));
  }
  return PreRootDatum(roots,coroots,rd.rank());
}

// compute twist and subsystem twisted involution $ww$ for $-theta^t$
// used in |GlobalTitsGroup| constructor and in |iblock|, |test| cmds
weyl::Twist SubSystem::twist(const WeightInvolution& theta,
			     WeylWord& ww) const
{
  RootNbrList Delta(rank()); // list of subsystem simple images by theta
  for (weyl::Generator i=0; i<rank(); ++i)
  {
    RootNbr image =
      inv_map[rd.root_index(theta*rd.root(parent_nr_simple(i)))];
    assert(image < numRoots());  // |image| is number of image in subsystem
    Delta[i] = rootMinus(image); // |-theta| image of |root(i)|
  }

  WeylWord wrt = // this is ordered according to parent perspective
    rootdata::wrt_distinguished(*this,Delta); // make |Delta| distinguished

  // |Delta| now describes the |sub|-side fundamental involution, a twist

  weyl::Twist result; // the subsystem twist that |Delta| has been reduced to
  for (weyl::Generator i=0; i<rank(); ++i)
    result[i] = RootSystem::simpleRootIndex(Delta[i]);

  // Let |theta_0| be the involution such that |Delta| describes $-theta_0^t$
  // then (for integrality systems) |theta_0| is parent quasi-split involution
  // and $theta=pw.theta_0$ where |pw| is |wrt| in terms of parent generators;
  // equivalently, $-theta^t=Delta.wrt^{-1}$ in terms of the subsystem

  // However, we want |ww| such that $-theta^t=ww.Delta$ (on subsystem side),
  // therefore |ww=twist(wrt^{-1})|. Nonetheless if |theta| is an involution,
  // it is the same to say $-theta^t=ww.Delta$ or $-theta^t=Delta.ww^{-1}$,
  // so following inversion and twist do nothing (|ww| is twisted involution).
  ww.resize(wrt.size()); // we must reverse and twist, to express on our side
  for (size_t i=0; i<wrt.size(); ++i)
    ww[wrt.size()-1-i] = result[wrt[i]];

  return result; // the result returned is the subsystem twist, not |ww|
}

// Here we seek twist and |ww| on parent side (dual with respect to |sub|)
// used in |blocks::param_block::compute_duals| (just for the induced |Twist|)
weyl::Twist SubSystem::parent_twist(const WeightInvolution& theta,
				    WeylWord& ww) const
{
  // beginning is identical to |SubSystem::twist| above
  RootNbrList Delta(rank());
  for (weyl::Generator i=0; i<rank(); ++i)
  {
    RootNbr image =
      inv_map[rd.root_index(theta*rd.root(parent_nr_simple(i)))];
    assert(image < numRoots());
    Delta[i] = image; // PLUS |theta| image of |root(i)|
  }

  // set |ww| to parent-side word, omitting inversion and twist done above
  ww = rootdata::wrt_distinguished(*this,Delta);
  // the above call has made |Delta| sub-distinguished from parent side

  weyl::Twist result;
  for (weyl::Generator i=0; i<rank(); ++i)
    result[i] = RootSystem::simpleRootIndex(Delta[i]);

  return result;
}

// get positive roots by converting the array |pos_map| to a |BitMap|
RootNbrSet SubSystem::positive_roots() const
{ return RootNbrSet(rd.numRoots(),pos_map); }

InvolutionData SubSystem::involution_data(const WeightInvolution& theta) const
{ return InvolutionData(rd,theta,positive_roots()); }

SubSystemWithGroup::SubSystemWithGroup(const RootDatum& parent,
				       const RootNbrList& sub_sys)
  : SubSystem(parent,sub_sys) // build
  , sub_W(RootSystem::cartanMatrix()) // use sub-side Cartan matrix built above
{}

SubSystemWithGroup SubSystemWithGroup::integral // pseudo contructor
  (const RootDatum& parent, const RatWeight& gamma)
{
  arithmetic::Numer_t n=gamma.denominator(); // signed!
  const Ratvec_Numer_t& v=gamma.numerator();
  RootNbrSet int_coroots(parent.numRoots());
  for (size_t i=0; i<parent.numPosRoots(); ++i)
    if (parent.posCoroot(i).dot(v)%n == 0)
      int_coroots.insert(parent.posRootNbr(i));

  // it suffices that simpleBasis computed below live until end of constructor
  return SubSystemWithGroup(parent,parent.simpleBasis(int_coroots));
}

} // |namespace subdatum|

} // |namespace atlas|
