/*
  This is output.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef OUTPUT_H  /* guard against multiple inclusions */
#define OUTPUT_H

#include <iostream>
#include <iosfwd>

#include "../Atlas.h"
#include "tags.h"

namespace atlas {

namespace output {

/******** function declarations *********************************************/

  class FormNumberMap;
  class Interface;

std::ostream& printBlockSizes
  (std::ostream&, InnerClass&, Interface&);


std::ostream& printRealForms(std::ostream& strm, const FormNumberMap& m);

std::ostream&
printCartanClass(std::ostream&,
		 const InnerClass&, size_t, output::Interface&);

std::ostream& printFiber(std::ostream&, const Fiber&,
			 const RealFormNbrList&);

std::ostream& printGradings
  (std::ostream&, InnerClass&, size_t, Interface&);

std::ostream& printGradings(std::ostream&, const Fiber&,
			    const RealFormNbrList&,
			    const RootSystem&);


std::ostream& printBlockStabilizer(std::ostream& strm,
				   RealReductiveGroup& G_R,
				   size_t cn,
				   RealFormNbr dual_rf);

std::ostream& printCartanClasses(std::ostream&,
				 RealReductiveGroup& G_R,
				 output::Interface&);

std::ostream& printCartanOrder(std::ostream&,
			       const RealReductiveGroup&);

std::ostream& printRealWeyl(std::ostream& strm,
			    RealReductiveGroup& G_R,
			    size_t cn);

std::ostream& printStrongReal(std::ostream& strm,
			      InnerClass& G_C,
			      const output::FormNumberMap& rfi,
			      size_t cn);

/******** type definitions **************************************************/

// an object to map an internal form number to an external number and name
class FormNumberMap
{

  RealFormNbrList d_in;
  RealFormNbrList d_out;

  std::vector<std::string> d_name; // real form names indexed by external number

 public:

// constructors and destructors
  FormNumberMap(const InnerClass&, const lietype::Layout&);

  FormNumberMap(const InnerClass&, const lietype::Layout&,
		tags::DualTag);

// accessors
  RealFormNbr in(RealFormNbr external_rf) const {  return d_in[external_rf]; }
  RealFormNbr out(RealFormNbr internal_rf) const { return d_out[internal_rf]; }

  const std::string& type_name(RealFormNbr external_rf) const
  { return d_name[external_rf]; }

  size_t numRealForms() const { return d_in.size(); }

// no manipulators
};

class Interface : public std::pair<FormNumberMap,FormNumberMap>
{
  typedef std::pair<FormNumberMap,FormNumberMap> base;
 public:

// constructors and destructors
  Interface(InnerClass& G, const lietype::Layout& lo)
  : base(FormNumberMap(G,lo),FormNumberMap(G,lo,tags::DualTag()))
  {}


  const output::FormNumberMap& realFormInterface() const
  { return this->first; }

  const output::FormNumberMap& dualRealFormInterface() const
  { return this->second; }

};

} // |namespace output|

} // |namespace atlas|

#endif
