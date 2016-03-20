/*
  This is realform_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef REALFORM_IO_H  /* guard against multiple inclusions */
#define REALFORM_IO_H

#include <iosfwd>
#include <string>

#include "../Atlas.h"
#include "tags.h"

namespace atlas {

/******** type declarations *************************************************/

namespace output {

  class FormNumberMap;

}

/******** functions declarations *********************************************/

namespace output {

  std::ostream& printRealForms(std::ostream&, const FormNumberMap&);

}

/******** type definitions ***************************************************/

namespace output {


// an object to map an internal form number to an external number and name
class FormNumberMap
{

  RealFormNbrList d_in;
  RealFormNbrList d_out;

  std::vector<std::string> d_name; // real form names indexed by external number

 public:

// constructors and destructors
  FormNumberMap(const ComplexReductiveGroup&, const lietype::Layout&);

  FormNumberMap(const ComplexReductiveGroup&, const lietype::Layout&,
	    tags::DualTag);

// copy, assignment and swap
  void swap(FormNumberMap&);

// accessors
  RealFormNbr in(RealFormNbr external_rf) const {  return d_in[external_rf]; }
  RealFormNbr out(RealFormNbr internal_rf) const { return d_out[internal_rf]; }

  const char* typeName(RealFormNbr external_rf) const;

  size_t numRealForms() const { return d_in.size(); }

// no manipulators
};

} // |namespace output|

} // |namespace|

#endif
