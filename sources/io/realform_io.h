/*
  This is realform_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef REALFORM_IO_H  /* guard against multiple inclusions */
#define REALFORM_IO_H

#include <iosfwd>

#include "complexredgp_fwd.h"

#include "lietype.h"
#include "tags.h"

namespace atlas {

/******** type declarations *************************************************/

namespace realform_io {

  class Interface;

}

/******** functions declarations *********************************************/

namespace realform_io {

  std::ostream& printRealForms(std::ostream&, const Interface&);

}

/******** type definitions ***************************************************/

namespace realform_io {

class Interface {

 private:

  RealFormNbrList d_in;
  RealFormNbrList d_out;

  std::vector<std::string> d_name;

 public:

// constructors and destructors
  Interface() {};

  Interface(const ComplexReductiveGroup&, const lietype::Layout&);

  Interface(const ComplexReductiveGroup&, const lietype::Layout&,
	    tags::DualTag);

  ~Interface() {};

// copy, assignment and swap
  void swap(Interface&);

// accessors
  RealFormNbr in(RealFormNbr rf) const {
    return d_in[rf];
  }

  size_t numRealForms() const {
    return d_in.size();
  }

  RealFormNbr out(RealFormNbr rf) const {
    return d_out[rf];
  }

  const char* typeName(RealFormNbr) const;

// manipulators
};

}

}

#endif
