/*
  This is bruhat.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Implementation of the class BruhatOrder.
#include "bruhat.h"

namespace atlas {


/*****************************************************************************

        Chapter I -- The BruhatOrder class

******************************************************************************/

namespace bruhat {


// Computes the full poset from the stored Hasse diagram.
void BruhatOrder::fillPoset()
{
  if (d_poset.size()==0)
    poset::Poset(d_hasse).swap(d_poset);
}

} // |namespace bruhat|

} // |namespace atlas|
