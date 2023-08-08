/*
  This is subsystem.cpp.

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "subsystem.h"

#include "ratvec.h"	// access to infinitesimal character |gamma|
#include "bitmap.h"	// root sets
#include "matreduc.h"   // |diagonalise| used in |codec| constructor

#include "lietype.h"	// returning |ext_gen|
#include "prerootdata.h"// returning |PreRootDatum|
#include "rootdata.h"	// |RootSystem| construction and methods
#include "innerclass.h"	// |integrality_datum_item| construction
#include "cartanclass.h"// |InvolutionData|
#include "weyl.h"	// subobject
#include "repr.h"       // for |repr::codec|

#include <cassert>

namespace atlas {

namespace subsystem {

SubSystem::SubSystem(const RootDatum& parent, const sl_list<RootNbr>& sub_sys)
  : RootSystem(parent.Cartan_matrix(sub_sys.to_vector()), // build new system
	       parent.prefer_coroots()) // no flip, roots will be roots
  , rd(parent) // share
  , which(parent.numPosRoots())
  , pos_map() // will be filled to size |numPosRoots()|
  , inv_map(numPosRoots()) // for |rd.posRootNbr(i)| at |which.position(i)|
  , sub_root(numPosRoots()) // initially empty vector of posroot information
{
  pos_map.reserve(numPosRoots()); // |pos_map| is indexed by \emph{our} posroots
  for (RootNbr alpha : sub_sys)
  {
    pos_map.push_back(alpha); // push simple roots (for subsystem) first
    which.insert(rd.posroot_index(alpha));
  }

  for (unsigned int i=sub_sys.size(); i<numPosRoots(); ++i)
  {
    RootNbr sub_alpha = posRootNbr(i); // a non-simple posroot of subsystem
    weyl::Generator s = find_descent(sub_alpha); // generator for subsystem
    simple_reflect_root(s,sub_alpha); // lower |sub_alpha| in our system
    RootNbr beta = pos_map[posroot_index(sub_alpha)]; // |beta| is in parent
    pos_map.push_back(rd.reflected_root(pos_map[s],beta));
    which.insert(rd.posroot_index(pos_map.back()));
  }

  // now complete setting |inv_map| and |sub_root| arrays
  for (unsigned int i=0; i<numPosRoots(); ++i)
  {
    RootNbr alpha = pos_map[i]; // now we use parent numbering
    assert(rd.is_posroot(alpha)); // conjugating to simple supposes this
    RootNbr alpha_pos = rd.posroot_index(alpha);
    assert(which.isMember(alpha_pos));
    inv_map[which.position(alpha_pos)] = posRootNbr(i);

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
  } // |for(i<numPosRoots())|
}

SubSystem SubSystem::integral // pseudo constructor for integral system
  (const RootDatum& parent, const RatWeight& gamma)
{
  arithmetic::Numer_t n=gamma.denominator(); // signed!
  const Ratvec_Numer_t& v=gamma.numerator();
  RootNbrSet int_poscoroots(parent.numPosRoots());
  for (size_t i=0; i<parent.numPosRoots(); ++i)
    if (parent.posCoroot(i).dot(v)%n == 0)
      int_poscoroots.insert(i);

  // it suffices that |pos_simples| computed below live until end of constructor
  return SubSystem(parent,parent.pos_simples(int_poscoroots));
}

RootNbr SubSystem::to_parent(RootNbr alpha) const
{
  RootNbr result = pos_map[rt_abs(alpha)];
  if (is_negroot(alpha))
    result = parent_datum().rootMinus(result);
  return result;
}

RootNbr SubSystem::from_parent(RootNbr alpha) const
{
  RootNbr alpha_pos = parent_datum().rt_abs(alpha);
  if (alpha_pos<parent_datum().numPosRoots() and which.isMember(alpha_pos))
  {
    RootNbr result = inv_map[which.position(alpha_pos)];
    return parent_datum().is_posroot(alpha) ? result : rootMinus(result);
  }
  return RootNbr(-1);
}

PreRootDatum SubSystem::pre_root_datum() const
{
  auto pr = parent_datum().rank(), sr=rank();
  LatticeMatrix simple_roots(pr,sr);
  LatticeMatrix simple_coroots(pr,sr);

  for (unsigned int j=0; j<sr; ++j)
  {
    simple_roots  .set_column(j,rd.root  (parent_nr_simple(j)));
    simple_coroots.set_column(j,rd.coroot(parent_nr_simple(j)));
  }
  return PreRootDatum(simple_roots,simple_coroots,not prefer_coroots());
}

// get simple roots by converting intial range of |pos_map| to a |BitMap|
RootNbrSet SubSystem::simple_roots() const
{ return RootNbrSet(rd.numRoots(),&pos_map[0],&pos_map[rank()]); }

// get positive roots by converting the array |pos_map| to a |BitMap|
RootNbrSet SubSystem::positive_roots() const
{ return RootNbrSet(rd.numRoots(),pos_map); }

RankFlags SubSystem::singular_generators(const RatWeight& gamma) const
{
  const Ratvec_Numer_t& v=gamma.numerator();
  RankFlags result;
  for (weyl::Generator s=0; s<rank(); ++s)
    result.set(s,simple_coroot(s).dot(v) == 0);

  return result;
}

InvolutionData SubSystem::involution_data(const WeightInvolution& theta) const
{ return InvolutionData(rd,theta,positive_roots()); }

SubSystemWithGroup::SubSystemWithGroup(const RootDatum& parent,
				       const sl_list<RootNbr>& sub_sys)
  : SubSystem(parent,sub_sys) // build
  , sub_W(RootSystem::Cartan_matrix()) // use sub-side Cartan matrix built above
{}

SubSystemWithGroup SubSystemWithGroup::integral // pseudo constructor
  (const RootDatum& parent, const RatWeight& gamma)
{
  arithmetic::Numer_t n=gamma.denominator(); // signed!
  const Ratvec_Numer_t& v=gamma.numerator();
  RootNbrSet int_poscoroots(parent.numPosRoots());
  for (size_t i=0; i<parent.numPosRoots(); ++i)
    if (parent.posCoroot(i).dot(v)%n == 0)
      int_poscoroots.insert(i);

  // it suffices that |pos_simples| computed below live until end of constructor
  return SubSystemWithGroup(parent,parent.pos_simples(int_poscoroots));
}

integral_datum_item::integral_datum_item
  (InnerClass& ic,const RootNbrSet& int_poscoroots)
    : ic(ic)
    , int_sys( ic.root_datum(), ic.root_datum().pos_simples(int_poscoroots) )
    , simple_coroots(int_sys.rank(),ic.rank()) // first is |RootSystem::rank|
{
  for (unsigned i=0; i<simple_coroots.n_rows(); ++i)
    simple_coroots.set_row(i,int_sys.simple_coroot(i));
}


SubSystem integral_datum_item::int_system(const WeylElt& w) const
{ return SubSystem { int_sys.parent_datum(), image_simples(w) }; }

sl_list<RootNbr> integral_datum_item::image_simples(const WeylElt& w) const
{
  WeylWord ww = ic.Weyl_group().word(ic.unfold(w));
  sl_list<RootNbr> result;

  const auto& rd = ic.root_datum();
  for (weyl::Generator s=0; s<int_sys.rank(); ++s)
  {
    RootNbr image = rd.permuted_root(ww,int_sys.parent_nr_simple(s));
    assert(rd.is_posroot(image)); // |ww| must map to integrally dominant
    result.push_back(image);
  }
  result.sort(); // force |integrality_datum| numbering
  return result;
}

int_Matrix integral_datum_item::coroots_matrix(const WeylElt& w) const
{
  auto integral_simples = image_simples(w);

  const auto& rd = ic.root_datum(); // the inner class root datum
  int_Matrix result(integral_simples.size(), rd.rank());
  unsigned i=0;
  for (const auto& alpha : integral_simples)
    result.set_row(i++,rd.coroot(alpha));
  return result;
}

repr::codec integral_datum_item::data (InvolutionNbr inv) const
  { return { ic,inv,simple_coroots }; }

repr::codec integral_datum_item::data (InvolutionNbr inv, const WeylElt& w) const
  { return { ic,inv, coroots_matrix(w) }; }

} // |namespace subdatum|

} // |namespace atlas|
