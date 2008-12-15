/*!
\file
\brief Declarations for CartanClass and Fiber classes.
*/
/*
  This is cartanclass_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef CARTANCLASS_FWD_H  /* guard against multiple inclusions */
#define CARTANCLASS_FWD_H

#include <utility>

/******** forward type definitions *******************************************/

namespace atlas {

namespace cartanclass {

  class InvolutionData;
  class Fiber;
  class CartanClass;

  //!\brief Element of the adjoint fiber or adjoint fiber group
  typedef unsigned int AdjointFiberElt;

  //!\brief Element of the fiber group
  typedef unsigned int FiberElt;

  /*!
  \brief Number of a W_imag orbit on the adjoint fiber group.
  */
  typedef unsigned short adjoint_fiber_orbit;
  /*!
  \brief Number of a W_imag orbit on the fiber group.
  */
  typedef unsigned short fiber_orbit;

  /*!\brief
    Identification of a class of real forms \f$f\f$ determined by value
    of \f$f^2\f$
  */
  typedef unsigned short int square_class;

  /*!
  \brief Second number is a possible square in Z(G) of a
  strong real form. First is a W_im orbit on the coset of the
  fiber group corresponding to that square.
  */
  typedef std::pair<fiber_orbit,square_class> StrongRealFormRep;
}

}

#endif
