/*!
\file
\brief Implementation of the class BruhatOrder.
*/

/*
  This is bruhat.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "bruhat.h"

namespace atlas {


/*****************************************************************************

        Chapter I -- The BruhatOrder class

******************************************************************************/

namespace bruhat {


/*!
  \brief Computes the full poset from the stored Hasse diagram.
*/
void BruhatOrder::fillPoset()
{
  if (d_poset.size()==0)
    poset::Poset(d_hasse).swap(d_poset);
}

} // namespace bruhat

} // namespace atlas
