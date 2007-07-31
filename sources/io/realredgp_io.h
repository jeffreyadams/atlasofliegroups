/*
  This is realredgp_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef REALREDGP_IO_H  /* guard against multiple inclusions */
#define REALREDGP_IO_H

#include <iosfwd>

#include "realredgp_io_fwd.h"

#include "complexredgp_io_fwd.h"
#include "realform_io.h"
#include "realredgp_fwd.h"

#include "realform.h"

namespace atlas {

/******** function declarations **********************************************/

namespace realredgp_io {

std::ostream& printBlockStabilizer(std::ostream& strm,
				   const realredgp::RealReductiveGroup& G_R,
				   size_t cn,
				   realform::RealForm rf);

std::ostream& printCartanClasses(std::ostream&,
				 const realredgp_io::Interface&);

std::ostream& printCartanOrder(std::ostream&,
			       const realredgp::RealReductiveGroup&);

std::ostream& printRealWeyl(std::ostream& strm,
			    const realredgp::RealReductiveGroup& G_R,
			    size_t cn);

std::ostream& printStrongReal(std::ostream& strm,
			      const realredgp::RealReductiveGroup& G_R,
			      const realform_io::Interface& rfi,
			      size_t cn);

}

/******** type definitions ***************************************************/

namespace realredgp_io {

class Interface {

 private:

  realredgp::RealReductiveGroup* d_realGroup;
  complexredgp_io::Interface* d_complexInterface;

 public:

// constructors and destructors
  Interface():d_realGroup(0),d_complexInterface(0) {}

  Interface(realredgp::RealReductiveGroup& G_R)
    :d_realGroup(&G_R),d_complexInterface(0) {}

  Interface(realredgp::RealReductiveGroup& G_R, complexredgp_io::Interface& CI)
    :d_realGroup(&G_R),d_complexInterface(&CI) {}

// copy, assignment and swap
  void swap(Interface&);

// accessors
  const complexredgp_io::Interface& complexInterface() const {
    return *d_complexInterface;
  }

  const complexredgp::ComplexReductiveGroup& complexGroup() const;

  const realform_io::Interface& realFormInterface() const;

  const realredgp::RealReductiveGroup& realGroup() const {
    return *d_realGroup;
  }

// manipulators
  complexredgp_io::Interface& complexInterface() {
    return *d_complexInterface;
  }

  realredgp::RealReductiveGroup& realGroup() {
    return *d_realGroup;
  }
};

}

}

#endif
