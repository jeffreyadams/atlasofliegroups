/*
  This is subsystem.cpp.

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "subsystem.h"

#include "ratvec.h"	// acces to infinitesimal character |gamma|
#include "bitmap.h"	// root sets
#include "matreduc.h"   // |diagonalise| used in |codec| constructor

#include "prerootdata.h"// value returned
#include "rootdata.h"	// |RootSystem| construction and methods
#include "innerclass.h"	// |integrality_datum_item| construction
#include "cartanclass.h"// |InvolutionData|
#include "weyl.h"	// subobject

#include <cassert>

namespace atlas {

namespace subsystem {

SubSystem::SubSystem(const RootDatum& parent, const sl_list<RootNbr>& sub_sys)
  : RootSystem(parent.cartanMatrix(sub_sys.to_vector()).transposed(), // build
	       not parent.prefer_coroots()) // flip, since coroots ar now roots
  , rd(parent) // share
  , pos_map()
  , inv_map(rd.numRoots()+1,-1) // one spare entry for "unfound root in parent"
  , sub_root(numPosRoots()) // initially empty vector of posroot information
{
  pos_map.reserve(numPosRoots()); // |pos_map| is indexed by \emph{our} posroots
  for (RootNbr alpha : sub_sys)
    pos_map.push_back(alpha);
  for (unsigned int i=sub_sys.size(); i<numPosRoots(); ++i)
  {
    RootNbr sub_alpha = posRootNbr(i);
    weyl::Generator s = find_descent(sub_alpha); // generator for subsystem
    simple_reflect_root(s,sub_alpha); // lower |sub_alpha| in our system
    RootNbr beta = pos_map[posroot_index(sub_alpha)]; // |beta| is in parent
    pos_map.push_back(rd.reflected_root(pos_map[s],beta));
  }

  // now complete setting |inv_map| and |sub_root| arrays
  for (unsigned int i=0; i<numPosRoots(); ++i)
  {
    RootNbr alpha = pos_map[i]; // now we use parent numbering
    inv_map[alpha] = posRootNbr(i); // refers to posroot |i| in subsystem
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

#if 0 // method is unused
// compute (dual side) twist and subsystem twisted involution |ww| for |theta|
weyl::Twist SubSystem::twist(const WeightInvolution& theta, WeylWord& ww) const
{
  RootNbrList Delta(rank()); // list of subsystem simple images by theta
  for (weyl::Generator i=0; i<rank(); ++i)
  {
    RootNbr image =
      inv_map[rd.root_index(theta*rd.root(parent_nr_simple(i)))];
    assert(image < numRoots());  // |image| is number of image in subsystem
    Delta[i] = rootMinus(image); // |-theta| image of |root(i)|
  }

  WeylWord wrt = // rightmost letter |wrt| applies first to disinguished |Delta|
    rootdata::wrt_distinguished(*this,Delta); // make |Delta| distinguished

  // |Delta| now describes a twist of the subsystem Dynkin diagram

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
  // so following reversal and twist do nothing (|ww| is twisted involution).
  ww.resize(wrt.size()); // we must reverse and twist, to express on our side
  for (size_t i=0; i<wrt.size(); ++i)
    ww[wrt.size()-1-i] = result[wrt[i]];

  return result; // the result returned is the subsystem twist, not |ww|
}
#endif


// get positive roots by converting the array |pos_map| to a |BitMap|
RootNbrSet SubSystem::positive_roots() const
{ return RootNbrSet(rd.numRoots(),pos_map); }

InvolutionData SubSystem::involution_data(const WeightInvolution& theta) const
{ return InvolutionData(rd,theta,positive_roots()); }

SubSystemWithGroup::SubSystemWithGroup(const RootDatum& parent,
				       const sl_list<RootNbr>& sub_sys)
  : SubSystem(parent,sub_sys) // build
  , sub_W(RootSystem::cartanMatrix()) // use sub-side Cartan matrix built above
{}

SubSystemWithGroup SubSystemWithGroup::integral // pseudo contructor
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
  : integral(ic.rootDatum(),ic.rootDatum().pos_simples(int_poscoroots))
  , simple_coroots(integral.rank(),ic.rank())
{
  for (unsigned i=0; i<integral.rank(); ++i)
    simple_coroots.set_row(i,
			   ic.rootDatum().coroot(integral.parent_nr_simple(i)));
}

integral_datum_item::codec::codec
  (const InnerClass& ic,
   unsigned int int_sys_nr, InvolutionNbr inv, const int_Matrix& coroots_mat)
    : coroots_matrix(ic.int_item(int_sys_nr).coroots_matrix())
    , theta_1_image_basis(ic.involution_table().theta_1_image_basis(inv))
    , diagonal(), in(), out()
{
  // get image of $-1$ eigenlattice in int-orth quotient, in coroot coordinates
  int_Matrix A = coroots_mat * theta_1_image_basis, row,col;
  diagonal=matreduc::diagonalise(A,row,col);
  // ensure |diagonal| entries positive, since we shall be reducing modulo them
  if (diagonal.size()>0 and diagonal[0]<0) // only this entry might be negative
  {
    diagonal[0] = -diagonal[0];
    row.rowMultiply(0u,-1); // restablish relation |row*A*col==diagonal|
  }

  auto rank = diagonal.size();
  in  = std::move(row); // keep full coordinate transform
  out = col.block(0,0,col.numRows(),rank); // chop part for final zero entries
}

} // |namespace subdatum|

} // |namespace atlas|
