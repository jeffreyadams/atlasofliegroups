/*
  This is complexredgp_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <iostream>

#include "complexredgp_io.h"

#include "cartan_io.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "ioutils.h"
#include "lietype.h"
#include "matrix.h"
#include "prettyprint.h"
#include "realform.h"
#include "realform_io.h"
#include "rootdata.h"
#include "tags.h"

/*
  This module contains the definition of an i/o interface for a
  ComplexReductiveGroup object. I have now convinced myself that the ownership
  relation should be that the intrface carries an unowned nonconstant pointer
  to the group object. So the group may evolve during its lifetime, being
  passed from one interface to another during its lifetime, and should finally
  be destroyed by its legitimate owner.
*/

namespace atlas {

/*****************************************************************************

        Chapter I --- The Interface class

******************************************************************************/

namespace complexredgp_io {

Interface::Interface(complexredgp::ComplexReductiveGroup& G,
		     const layout::Layout& lo)
  : d_complexGroup(&G)
  , d_layout(lo)
  , d_realFormInterface(G,lo)
  , d_dualRealFormInterface(realform_io::Interface(G,lo,tags::DualTag()))
{}

/******** copy, assignment and swap ******************************************/

/*
  Note: recall that the complex group is not owned; only the pointers are
  swapped.
*/
void Interface::swap(Interface& other)
{
  std::swap(d_complexGroup,other.d_complexGroup);
  d_layout.swap(other.d_layout);
  d_realFormInterface.swap(other.d_realFormInterface);
  d_dualRealFormInterface.swap(other.d_dualRealFormInterface);
}

} // namespace complexredgp_io

/*****************************************************************************

        Chapter II --- Functions declared in complexredgp_io

******************************************************************************/

namespace complexredgp_io {

std::ostream& printBlockSizes(std::ostream& strm, const Interface& CI)
{
  const complexredgp::ComplexReductiveGroup& G = CI.complexGroup();
  const realform_io::Interface rfi = CI.realFormInterface();
  const realform_io::Interface drfi = CI.dualRealFormInterface();

  size_t rf = G.numRealForms();
  size_t drf = G.numDualRealForms();

  matrix::Matrix<unsigned long> block(rf,drf);
  unsigned long maxEntry = 0;

  for (size_t i = 0; i < block.numRows(); ++i)
    for (size_t j = 0; j < block.numColumns(); ++j) {
      block(i,j) = G.block_size(rfi.in(i),drfi.in(j));
      if (block(i,j) > maxEntry)
	maxEntry = block(i,j);
    }

  int width = ioutils::digits(maxEntry,10ul);
  prettyprint::printMatrix(strm,block,width+2);

  return strm;
}



/*
  Synopsis: outputs the gradings corresponding to the various real forms
  defined for Cartan #cn.

  The gradings are output in the same order as the corresponding orbits are
  output in the "cartan" command.
*/
std::ostream& printGradings(std::ostream& strm, size_t cn, const Interface& CI)
{
  const complexredgp::ComplexReductiveGroup& G = CI.complexGroup();
  const cartanclass::CartanClass& cc = G.cartan(cn);

  realform::RealFormList rfl(cc.numRealForms());
  const realform_io::Interface& rfi = CI.realFormInterface();

  for (size_t i = 0; i < rfl.size(); ++i)
    rfl[i] = rfi.out(G.realFormLabels(cn)[i]);

  cartan_io::printGradings(strm,cc.fiber(),rfl,G.rootDatum());

  return strm;
}

}

}
