/*
  This is kltest.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
/*
   Declaration of the functions in namespace |kltest|.

   These functions are designed to test a mathematical assertion used in
   the implementation of the KL algorithm.
*/

#ifndef KLTEST_H  /* guard against multiple inclusions */
#define KLTEST_H

#include "../Atlas.h"

#include "permutations.h"

namespace atlas {

/******** type declarations *************************************************/

namespace kltest {

}

/******** function declarations **********************************************/

namespace kltest {

  bool checkBasePoint(const KGB&);

  void dualityPermutation(Permutation&, const kl::KL_table&);

  bool dualityVerify(const kl::KL_table& kl_tab, const kl::KL_table& dual_kl_tab);
}

/******** type definitions **************************************************/

namespace kltest {

}

}

#endif
