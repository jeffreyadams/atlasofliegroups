/*
  This is complexredgp_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "complexredgp_io.h"

#include <iostream>

#include "cartan_io.h"	// |printGradings|
#include "cartanclass.h"
#include "complexredgp.h"

#include "ioutils.h"	// |digits|
#include "prettyprint.h" // |printMatrix|
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

Interface::Interface(ComplexReductiveGroup& G,
		     const lietype::Layout& lo)
  : d_complexGroup(&G)
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
  d_realFormInterface.swap(other.d_realFormInterface);
  d_dualRealFormInterface.swap(other.d_dualRealFormInterface);
}

} // namespace complexredgp_io

/*****************************************************************************

        Chapter II --- Functions declared in complexredgp_io

******************************************************************************/

namespace complexredgp_io {

std::ostream& printBlockSizes(std::ostream& strm, Interface& CI)
{
  ComplexReductiveGroup& G = CI.complexGroup();
  const realform_io::Interface rfi = CI.realFormInterface();
  const realform_io::Interface drfi = CI.dualRealFormInterface();

  matrix::Matrix<unsigned long> block(G.numRealForms(),G.numDualRealForms());
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
std::ostream& printGradings(std::ostream& strm, size_t cn, Interface& CI)
{
  ComplexReductiveGroup& G = CI.complexGroup();
  const CartanClass& cc = G.cartan(cn);

  RealFormNbrList rfl(cc.numRealForms());
  const realform_io::Interface& rfi = CI.realFormInterface();

  for (size_t i = 0; i < rfl.size(); ++i)
    rfl[i] = rfi.out(G.realFormLabels(cn)[i]);

  cartan_io::printGradings(strm,cc.fiber(),rfl,G.rootDatum());

  return strm;
}

}

}
