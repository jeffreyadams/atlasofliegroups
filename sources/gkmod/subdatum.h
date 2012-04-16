/*!
\file
\brief Class definition for the SubDatum class.
*/
/*
  This is subdatum.h

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SUBDATUM_H  /* guard against multiple inclusions */
#define SUBDATUM_H

#include "atlas_types.h"

#include "weyl.h"
#include "tits.h"
#include "subsystem.h"

namespace atlas {

namespace subdatum {


/* the class |SubDatum| is mathematically much richer than |SubSystem| (also,
   note the latter holds a reference to a parent root\emph{datum}, so the
   terminology is somewhat misleading). While |SubSystem| may be regarded as a
   companion to |RootDatum|, the |SubDatum| class depends on gkmod stuff. Its
   unique constructor requires a root datum involution |theta| (coded by a KGB
   element for a real form) and an infinitesimal character. The latter defines
   a subsystem by integrality, while |theta| defines a subsystem twist and a
   word that expresses it.
 */
class SubDatum : public SubSystemWithGroup
{
  WeylWord base_ww;  // we need this variable mostly in the constructor!
  WeightInvolution delta; // together with twist: what was missing
  TitsGroup Tg; // twisted Weyl group, plus stuff our base class knows
  WeylElt ini_tw; // records involution of initial |x| wrt |SubSystem|

  size_t rank() const; // forbid using this directly
 public:
  SubDatum(RealReductiveGroup& GR,
	   const RatWeight& gamma,
	   KGBElt x);

  const TitsGroup& Tits_group() const {return Tg; }
  TwistedInvolution init_twisted() const { return ini_tw; }
  const WeylWord& base_twisted_in_parent() const { return base_ww; }

  size_t semisimple_rank() const { return SubSystem::rank(); }

  WeightInvolution matrix(TwistedInvolution tw) const;

}; // |class SubDatum|

} // |namespace subdatum|

} // |namespace atlas|

#endif
