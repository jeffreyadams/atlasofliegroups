/*!
\file
\brief Class definitions and function declarations for the class InvolutionSet.
*/
/*
  This is involutions.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef INVOLUTIONS_H  /* guard against multiple inclusions */
#define INVOLUTIONS_H

#include "involutions_fwd.h"

#include "complexredgp_fwd.h"
#include "weyl.h"

namespace atlas {

/******** constant declarations *********************************************/

namespace involutions {

  const size_t UndefInvolution = ~0ul;

}

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace involutions {

class InvolutionSet {

 protected:

  size_t d_size;
  size_t d_rank; // semisimple rank
  std::vector<std::vector<size_t> > d_action; // indexed by generator first
  std::vector<size_t> d_cartan;               // indexed by involution
  std::vector<weyl::TwistedInvolution> d_involution;
  std::vector<weyl::TwistedInvolution> d_dualInvolution;

 public:
  // constructors and destructors
  InvolutionSet();

  explicit InvolutionSet(complexredgp::ComplexReductiveGroup&);

  virtual ~InvolutionSet() {}

  // copy, assignment and swap
  void swap(InvolutionSet&);

  // accessors
  size_t action(size_t s, size_t w) const { return d_action[s][w]; }

  const weyl::TwistedInvolution& involution(size_t j) const
    { return d_involution[j]; }

  const weyl::TwistedInvolution& dualInvolution(size_t j) const
    { return d_dualInvolution[j]; }

  size_t involutionNbr(const weyl::TwistedInvolution&,
		       const weyl::TwistedWeylGroup&) const;

  const size_t rank() const { return d_rank; }

  const size_t size() const { return d_size; }

  // manipulators
}; // |class InvolutionSet|

} // |namespace involutions|

} // |namesspace atlas|

#endif
