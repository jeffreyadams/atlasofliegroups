/*!
\file
\brief Implementation of the class RealReductiveGroup.
*/
/*
  This is realredgp.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "realredgp.h"

#include "matreduc.h" // |adapted_basis|

#include "cartanclass.h"  // |Fiber|, and |toMostSplit| function (assertion)
#include "complexredgp.h" // various methods
#include "rootdata.h"     // |refl_prod| function (assertion)
#include "tori.h"         // |tori::RealTorus| used
#include "kgb.h"          // |KGB| constructed
#include "tits.h"         // |TitsGroup| and |TitsElement| access

#include <cassert>

namespace atlas {

/*****************************************************************************

        Chapter I -- The RealReductiveGroup class

******************************************************************************/

namespace realredgp {


/*!
  Synopsis : constructs a real reductive group from the datum of a complex
  reductive group and a real form.
*/
RealReductiveGroup::RealReductiveGroup
  (ComplexReductiveGroup& G_C, RealFormNbr rf)
  : d_complexGroup(G_C)
  , d_realForm(rf)
  , d_connectivity() // wait for most split torus to be constructed below

  , square_class_cocharacter(some_coch(G_C,G_C.xi_square(rf)))
  , torus_part_x0(G_C.x0_torus_part(rf))

  , d_Tg(new // allocate private copy
	 TitsCoset(G_C,grading_of_simples(G_C,square_class_cocharacter)))
  , kgb_ptr(NULL)
  , dual_kgb_ptr(NULL)
  , d_status()
{ construct(); }


RealReductiveGroup::RealReductiveGroup
  (ComplexReductiveGroup& G_C, RealFormNbr rf,
   const RatCoweight& coch, TorusPart x0_torus_part)
  : d_complexGroup(G_C)
  , d_realForm(rf)
  , d_connectivity() // wait for most split torus to be constructed below

  , square_class_cocharacter(coch)
  , torus_part_x0(x0_torus_part)

  , d_Tg(new // allocate private copy
	 TitsCoset(G_C,grading_of_simples(G_C,square_class_cocharacter)))
  , kgb_ptr(NULL)
  , dual_kgb_ptr(NULL)
  , d_status()
{ construct(); }

void RealReductiveGroup::construct()
{
  ComplexReductiveGroup& G_C=d_complexGroup;
  RealFormNbr rf=d_realForm;

  tori::RealTorus msT = G_C.cartan(G_C.mostSplit(rf)).fiber().torus();
  d_connectivity = topology::Connectivity(msT,G_C.rootDatum());

  d_status.set(IsConnected,d_connectivity.component_rank() == 0);
  d_status.set(IsCompact,msT.isCompact());

  d_status.set(IsQuasisplit,rf == G_C.quasisplit());
  d_status.set(IsSplit,msT.isSplit());
  d_status.set(IsSemisimple,G_C.rank() == G_C.semisimpleRank());

#ifndef NDEBUG
  // construct the torus for the most split Cartan
  const Fiber& fundf = G_C.fundamental();
  RootNbrSet so= cartanclass::toMostSplit(fundf,rf,G_C.rootSystem());

  // recompute matrix of most split Cartan
  const RootDatum& rd = G_C.rootDatum();
  tori::RealTorus T1
    (rootdata::refl_prod(so,rd) * G_C.distinguished()); // factors commute

  topology::Connectivity c(T1,rd);
  assert(d_connectivity.component_rank() == c.component_rank());
#endif
}


RealReductiveGroup::~RealReductiveGroup()
{ delete d_Tg; delete kgb_ptr; delete dual_kgb_ptr; }

/******** accessors *********************************************************/



/******** manipulators ******************************************************/


void RealReductiveGroup::swap(RealReductiveGroup& other)
{
  assert(&d_complexGroup==&other.d_complexGroup); // cannot swap references
  std::swap(d_realForm,other.d_realForm);
  d_connectivity.swap(other.d_connectivity);
  std::swap(d_Tg,other.d_Tg);
  std::swap(kgb_ptr,other.kgb_ptr);
  std::swap(dual_kgb_ptr,other.dual_kgb_ptr);
  std::swap(d_status,other.d_status);
}


const RootDatum& RealReductiveGroup::rootDatum() const
  { return d_complexGroup.rootDatum(); }

const TitsGroup& RealReductiveGroup::titsGroup() const
  { return d_Tg->titsGroup(); }

const WeylGroup& RealReductiveGroup::weylGroup() const
  { return d_complexGroup.weylGroup(); }

const TwistedWeylGroup& RealReductiveGroup::twistedWeylGroup() const
  { return d_complexGroup.twistedWeylGroup(); }

BitMap RealReductiveGroup::Cartan_set() const
  { return complexGroup().Cartan_set(d_realForm); }

// Returns Cartan \#cn (assumed to belong to cartanSet()) of the group.
const CartanClass& RealReductiveGroup::cartan(size_t cn) const
  { return d_complexGroup.cartan(cn); }

RatCoweight RealReductiveGroup::g() const
  { return g_rho_check()+rho_check(rootDatum()); }

// the base grading can be computed directly from $g-\rho\check$, as imaginary
// simple roots are simple-imaginary, so dot product flags the compact ones
Grading RealReductiveGroup::base_grading() const
{
  return complexredgp::grading_of_simples(complexGroup(),g_rho_check());
}

size_t RealReductiveGroup::numCartan() const { return Cartan_set().size(); }

size_t RealReductiveGroup::rank() const { return rootDatum().rank(); };

size_t RealReductiveGroup::semisimpleRank() const
  { return rootDatum().semisimpleRank(); }

size_t RealReductiveGroup::numInvolutions()
  { return complexGroup().numInvolutions(Cartan_set()); }

size_t RealReductiveGroup::KGB_size() const
 { return d_complexGroup.KGB_size(d_realForm); }

size_t RealReductiveGroup::mostSplit() const
 { return d_complexGroup.mostSplit(d_realForm); }

/*
  Algorithm: the variable |rset| is first made to flag, among the imaginary
  roots of the fundamental Cartan, those that are noncompact for the chosen
  representative (in the adjoint fiber) of the real form of |G|. The result is
  formed by extracting only the information concerning the presence of the
  \emph{simple} roots in |rset|.
*/
Grading RealReductiveGroup::grading_offset()
{
  RootNbrSet rset= noncompactRoots(); // grading for real form rep
  return cartanclass::restrictGrading(rset,rootDatum().simpleRootList());
}

cartanclass::square_class RealReductiveGroup::square_class() const
  { return d_complexGroup.fundamental().central_square_class(d_realForm); }



const size_t RealReductiveGroup::component_rank() const
  { return d_connectivity.component_rank(); }
const SmallBitVectorList& RealReductiveGroup::dualComponentReps() const
  { return d_connectivity.dualComponentReps(); }

const WeightInvolution& RealReductiveGroup::distinguished() const
  { return d_complexGroup.distinguished(); }

RootNbrSet RealReductiveGroup::noncompactRoots() const
  { return d_complexGroup.noncompactRoots(d_realForm); }


// return stored KGB structure, after generating it if necessary
const KGB& RealReductiveGroup::kgb()
{
  if (kgb_ptr==NULL)
    kgb_ptr = new KGB(*this,Cartan_set(),false); // generate as non-dual
  return *kgb_ptr;
}

// return stored KGB structure, after generating it if necessary
const KGB& RealReductiveGroup::kgb_as_dual()
{
  if (dual_kgb_ptr==NULL)
    dual_kgb_ptr = new KGB(*this,Cartan_set(),true); // generate as dual
  return *dual_kgb_ptr;
}

// return stored Bruhat order of KGB, after generating it if necessary
const BruhatOrder& RealReductiveGroup::Bruhat_KGB()
{
  kgb(); // ensure |kgb_ptr!=NULL|, but we cannot use (|const|) result here
  return kgb_ptr->bruhatOrder(); // get Bruhat order (generate if necessary)
}


//				Functions

TorusPart minimal_torus_part
  (const ComplexReductiveGroup& G, RealFormNbr wrf, RatCoweight coch,
   TwistedInvolution tw, // by value, modified
   const RatCoweight& torus_factor
  )
{
  assert(torus_factor==torus_factor*G.matrix(tw)); // assuming already projected

  RatCoweight diff = (torus_factor-coch).normalize();
  assert (diff.denominator()==1); // since $\exp(2i\pi diff)=1$

  Grading gr = complexredgp::grading_of_simples(G,coch);
  TitsCoset Tc(G,gr);
  const auto& Tg = Tc.titsGroup();
  const auto& W = Tg.weylGroup();
  const auto& i_tab = G.involution_table();
  const TwistedInvolution e; // basis (identity) twisted involution

  TitsElt a (Tg,TorusPart(diff.numerator()),tw);
  while (a.tw()!=e) // move to fundamental fiber
  {
    weyl::Generator s = W.leftDescent(a.tw());
    if (Tg.hasTwistedCommutation(s,a.tw()))
    {
      SmallSubspace mod_space = i_tab.mod_space(i_tab.nr(a.tw()));
      Tc.inverse_Cayley_transform(a,s,mod_space);
    }
    else
      Tc.basedTwistedConjugate(a,s);
  }
  i_tab.reduce(a); // ensure reduction modulo fiber group denominator subgroup

// now find the minimal element in the imaginary Weyl group orbit of |a.t()|
// that induces |G.simple_roots_x0_compact(wrf)| on imaginary simple roots
  const auto base_cpt = compacts_for(G,y_values::exp_pi(coch));
  const auto wrf_cpt = G.simple_roots_x0_compact(wrf);
  const auto afr = G.fundamental().adjointFiberRank();
  const cartanclass::AdjointFiberElt image // which will give required compacts
    ( (base_cpt^wrf_cpt).slice(G.simple_roots_imaginary()), afr);

  const TorusPart t = Tg.left_torus_part(a);
  const auto& fund_fg = G.fundamental().fiberGroup();
  const cartanclass::FiberElt y = fund_fg.toBasis(t);

  const containers::sl_list<TorusPart> candidates =
    complexredgp::preimage(G.fundamental(), G.xi_square(wrf),y,image);

  assert(not candidates.empty());

  // choose minimal torus part among reduced ones giving |wrf_cpt| as compacts
  auto it = candidates.begin();
  auto min = *it;
  while (not (++it).at_end())
    if (t+*it<t+min) // offset |t| does not cancel from this relation!
      min = *it;

  Tg.left_add(min,a);
  assert((i_tab.reduce(a),Tg.left_torus_part(a)==t+min)); // already reduced
  assert(compact_simples(Tc,a,G.simple_roots_imaginary())==wrf_cpt);

  return Tg.left_torus_part(a);
} // |minimal_torus_part|

} // |namespace realredgp|

} // |namespace atlas|
