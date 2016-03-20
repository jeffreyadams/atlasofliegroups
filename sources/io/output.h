/*
  This is output.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef OUTPUT_H  /* guard against multiple inclusions */
#define OUTPUT_H



#include "realform_io.h"	// containment of |Interface|s

namespace atlas {

/******** function declarations *********************************************/

namespace output {

std::ostream& printBlockSizes(std::ostream&, Interface&);

std::ostream& printGradings(std::ostream&, size_t, Interface&);

}

/******** type definitions **************************************************/

namespace output {

class Interface {

 private:

  ComplexReductiveGroup* d_complexGroup;

  output::FormNumberMap d_realFormInterface;
  output::FormNumberMap d_dualRealFormInterface;

 public:

// constructors and destructors
  Interface(ComplexReductiveGroup&, const lietype::Layout&);

  ~Interface() {}

// copy, assignment and swap
  void swap(Interface&);

// accessors
  const ComplexReductiveGroup& complexGroup() const {
    return *d_complexGroup;
  }

  const output::FormNumberMap& dualRealFormInterface() const {
    return d_dualRealFormInterface;
  }

  const output::FormNumberMap& realFormInterface() const {
    return d_realFormInterface;
  }

// manipulators
  ComplexReductiveGroup& complexGroup() {
    return *d_complexGroup;
  }

};

}

}

#endif
