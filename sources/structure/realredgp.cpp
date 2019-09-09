/*
  This is realredgp.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Implementation of the class RealReductiveGroup.


#include "realredgp.h"

#include "matreduc.h" // |adapted_basis|

#include "cartanclass.h"  // |Fiber|, and |toMostSplit| function (assertion)
#include "innerclass.h"   // various methods
#include "rootdata.h"     // |refl_prod| function (assertion)
#include "topology.h"     // |topology::dual_component_group_basis|
#include "kgb.h"          // |KGB| constructed
#include "tits.h"         // |TitsGroup| and |TitsElement| access

#include <cassert>

namespace atlas {
namespace realredgp {

/*****************************************************************************

        Chapter I -- The RealReductiveGroup class

******************************************************************************/

RealReductiveGroup::RealReductiveGroup (InnerClass& G_C, RealFormNbr rf)
  : d_innerClass(G_C)
  , d_realForm(rf)
  , dual_pi0_gens() // wait for most split torus to be constructed below

  , square_class_cocharacter(some_coch(G_C,G_C.xi_square(rf)))
  , torus_part_x0(G_C.x0_torus_part(rf))

  , d_Tg(new // allocate private copy
	 TitsCoset(G_C,grading_of_simples(G_C,square_class_cocharacter)))
  , kgb_ptr(nullptr)
  , dual_kgb_ptr(nullptr)
  , d_status()
{ construct(); }


RealReductiveGroup::RealReductiveGroup
  (InnerClass& G_C, RealFormNbr rf,
   const RatCoweight& coch, TorusPart x0_torus_part)
  : d_innerClass(G_C)
  , d_realForm(rf)
  , dual_pi0_gens() // wait for most split torus to be constructed below

  , square_class_cocharacter(coch)
  , torus_part_x0(x0_torus_part)

  , d_Tg(new // allocate private copy
	 TitsCoset(G_C,grading_of_simples(G_C,square_class_cocharacter)))
  , kgb_ptr(nullptr)
  , dual_kgb_ptr(nullptr)
  , d_status()
{ construct(); }

void RealReductiveGroup::construct()
{
  const InnerClass& G_C=d_innerClass;
  const RealFormNbr rf=d_realForm;
  const RootDatum& rd = G_C.rootDatum();

  const auto ms_tau = G_C.cartan(G_C.mostSplit(rf)).involution();
  dual_pi0_gens = topology::dual_component_group_basis(ms_tau,rd);

  d_status.set(IsConnected,dual_pi0_gens.empty());
  d_status.set(IsCompact,(ms_tau-1).is_zero());

  d_status.set(IsQuasisplit,rf == G_C.quasisplit());
  d_status.set(IsSplit,(ms_tau+1).is_zero());
  d_status.set(IsSemisimple,G_C.rank() == G_C.semisimpleRank());

#ifndef NDEBUG
  // construct the torus for the most split Cartan
  RootNbrSet noncompacts = G_C.noncompactRoots(rf);
  TwistedInvolution xi; // empty: the distinguished one for the inner class
  RootNbrSet so =
    gradings::max_orth(noncompacts,
		       G_C.involution_data(xi).imaginary_roots(),
		       G_C.rootSystem());

  // recompute matrix of most split Cartan
  auto tau1 = rootdata::refl_prod(so,rd) * G_C.distinguished();
  assert(component_rank() ==
	 topology::dual_component_group_basis(tau1,rd).size());
#endif
}


RealReductiveGroup::~RealReductiveGroup()
{ delete d_Tg; delete kgb_ptr; delete dual_kgb_ptr; }

/******** accessors *********************************************************/



/******** manipulators ******************************************************/


void RealReductiveGroup::swap(RealReductiveGroup& other)
{
  assert(&d_innerClass==&other.d_innerClass); // cannot swap references
  std::swap(d_realForm,other.d_realForm);
  dual_pi0_gens.swap(other.dual_pi0_gens);

  std::swap(square_class_cocharacter,other.square_class_cocharacter);
  std::swap(torus_part_x0,other.torus_part_x0);

  std::swap(d_Tg,other.d_Tg);
  std::swap(kgb_ptr,other.kgb_ptr);
  std::swap(dual_kgb_ptr,other.dual_kgb_ptr);
  std::swap(d_status,other.d_status);
}


const RootDatum& RealReductiveGroup::rootDatum() const
  { return d_innerClass.rootDatum(); }

const TitsGroup& RealReductiveGroup::titsGroup() const
  { return d_Tg->titsGroup(); }

const WeylGroup& RealReductiveGroup::weylGroup() const
  { return d_innerClass.weylGroup(); }

const TwistedWeylGroup& RealReductiveGroup::twistedWeylGroup() const
  { return d_innerClass.twistedWeylGroup(); }

BitMap RealReductiveGroup::Cartan_set() const
  { return innerClass().Cartan_set(d_realForm); }

// Returns Cartan \#cn (assumed to belong to cartanSet()) of the group.
const CartanClass& RealReductiveGroup::cartan(size_t cn) const
  { return d_innerClass.cartan(cn); }

RatCoweight RealReductiveGroup::g() const
  { return g_rho_check()+rho_check(rootDatum()); }

// the base grading can be computed directly from $g-\rho\check$, as imaginary
// simple roots are simple-imaginary, so dot product flags the compact ones
Grading RealReductiveGroup::base_grading() const
{
  return innerclass::grading_of_simples(innerClass(),g_rho_check());
}

size_t RealReductiveGroup::numCartan() const { return Cartan_set().size(); }

size_t RealReductiveGroup::rank() const { return rootDatum().rank(); };

size_t RealReductiveGroup::semisimpleRank() const
  { return rootDatum().semisimpleRank(); }

size_t RealReductiveGroup::numInvolutions()
  { return innerClass().numInvolutions(Cartan_set()); }

size_t RealReductiveGroup::KGB_size() const
 { return d_innerClass.KGB_size(d_realForm); }

size_t RealReductiveGroup::mostSplit() const
 { return d_innerClass.mostSplit(d_realForm); }

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


const size_t RealReductiveGroup::component_rank() const
{ return dual_pi0_gens.size(); }
const SmallBitVectorList& RealReductiveGroup::dualComponentReps() const
  { return dual_pi0_gens; }

const WeightInvolution& RealReductiveGroup::distinguished() const
  { return d_innerClass.distinguished(); }

RootNbrSet RealReductiveGroup::noncompactRoots() const
  { return d_innerClass.noncompactRoots(d_realForm); }


// return stored KGB structure, after generating it if necessary
const KGB& RealReductiveGroup::kgb()
{
  if (kgb_ptr==nullptr)
    kgb_ptr = new KGB(*this,Cartan_set(),false); // generate as non-dual
  return *kgb_ptr;
}

// return stored KGB structure, after generating it if necessary
const KGB& RealReductiveGroup::kgb_as_dual()
{
  if (dual_kgb_ptr==nullptr)
    dual_kgb_ptr = new KGB(*this,Cartan_set(),true); // generate as dual
  return *dual_kgb_ptr;
}

// return stored Bruhat order of KGB, after generating it if necessary
const BruhatOrder& RealReductiveGroup::Bruhat_KGB()
{
  kgb(); // ensure |kgb_ptr!=nullptr|, but we cannot use (|const|) result here
  return kgb_ptr->bruhatOrder(); // get Bruhat order (generate if necessary)
}


//				Functions

TorusPart minimal_torus_part
  (const InnerClass& G, RealFormNbr wrf, const RatCoweight& coch,
   const TwistedInvolution& tw, const RatCoweight& torus_factor
  )
{
  assert(torus_factor==torus_factor*G.matrix(tw)); // assuming already projected

  const auto& rd = G.rootDatum();

  RatCoweight diff = (torus_factor-coch).normalize();
  assert (diff.denominator()==1); // since $\exp(2i\pi diff)=1$

  TorusPart tp(diff.numerator());
  { // move to fundamental fiber
    TitsCoset Tc(G,innerclass::grading_of_simples(G,coch));
    const auto& Tg = Tc.titsGroup();
    TitsElt a (Tg,tp,tw);
    const auto& W = Tg.weylGroup();
    const auto& i_tab = G.involution_table();
    const TwistedInvolution e; // basis (identity) twisted involution

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
    tp = Tg.left_torus_part(a);
  }
  RatCoweight cowt = coch + lift(tp);

  const auto basis = G.fundamental_fiber().simpleImaginary();
  const auto r = basis.size();
  Grading start_grading;
  for (unsigned i=0; i<r; ++i)
    start_grading.set(i,cowt.dot(rd.root(basis[i]))%2==0);

  RankFlags basis_simples{0u},mask;
  { unsigned int i=0;
    while (i<r and rd.is_simple_root(basis[i]))
      basis_simples.set(rd.simpleRootIndex(basis[i++]));
    mask = RankFlags(constants::lMask[i]);
  }

  Grading target // grading of |basis_simples| that we are looking for
    = G.simple_roots_x0_compact(wrf).slice(basis_simples)^mask;

  BitMap seen (1ul<<rd.rank()); seen.insert(tp.data().to_ulong());
  containers::simple_list<TorusPart> candidates;

  { containers::stack<std::pair<TorusPart,Grading> > to_do;
    to_do.push(std::make_pair(tp,start_grading));

    std::vector<TorusPart> m_alpha; m_alpha.reserve(r);
    std::vector<Grading> grading_shift; grading_shift.reserve(r);
    for (unsigned i=0; i<r; ++i)
    {
      const auto& alpha_v = rd.coroot(basis[i]);
      m_alpha.push_back(TorusPart(alpha_v));
      Grading shift;
      for (unsigned j=0; j<r; ++j)
        shift.set(j,alpha_v.dot(rd.root(basis[j]))%2!=0);
      grading_shift.push_back(shift);
    }

    do
    { const auto tp=to_do.top().first; const auto gr=to_do.top().second;
      to_do.pop();
      if ((gr&mask)==target) // record those that are candidate for being x0
	candidates.push_front(tp);
      for(auto it=gr.begin(); it(); ++it)
      { const auto new_tp = tp + m_alpha[*it];
	if (not seen.isMember(new_tp.data().to_ulong()))
	{ seen.insert(new_tp.data().to_ulong());
	  to_do.push(std::make_pair(new_tp,gr^grading_shift[*it]));
	}
      }
    }
  while (not to_do.empty());
  }

  assert(not candidates.empty());

  // choose minimal torus part among reduced ones giving |wrf_cpt| as compacts
  auto it = candidates.begin();
  auto min = *it;
  while (not (++it).at_end())
    if (*it<min)
      min = *it;

  return min;
} // |minimal_torus_part|

} // |namespace realredgp|

} // |namespace atlas|
