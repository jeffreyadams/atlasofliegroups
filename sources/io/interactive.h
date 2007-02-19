/*
  This is interactive.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef INTERACTIVE_H  /* guard against multiple inclusions */
#define INTERACTIVE_H

#include "interactive_fwd.h"

#include "bitmap_fwd.h"
#include "complexredgp_io_fwd.h"
#include "prerootdata_fwd.h"
#include "realredgp_fwd.h"
#include "rootdata_fwd.h"

#include "error.h"
#include "input.h"
#include "latticetypes.h"
#include "layout.h"
#include "lietype.h"
#include "realform.h"
#include "tags.h"

#include <iostream>

namespace atlas {

/******** type declarations **************************************************/

/* the strange namespace is an historic artifact; the class was moved */

namespace ioutils {

  class OutputFile;

}

/******** function declarations **********************************************/

namespace interactive {

  std::string getFileName(std::string prompt)
    throw(error::InputError);

  void bitMapPrompt(std::string&, const char*, const bitmap::BitMap&);

  void getCartanClass(size_t&, const bitmap::BitMap&,
		      input::InputBuffer&)
    throw(error::InputError);

  void getInnerClass(latticetypes::LatticeMatrix&, layout::Layout&)
    throw(error::InputError);

  void getInteractive(lietype::LieType&) throw(error::InputError);

  void getInteractive(lietype::InnerClassType&, const lietype::LieType&)
    throw(error::InputError);

  void getInteractive(prerootdata::PreRootDatum&, latticetypes::WeightList&,
		      const lietype::LieType&) throw(error::InputError);

  void getInteractive(realform::RealForm&, const complexredgp_io::Interface&)
    throw(error::InputError);

  void getInteractive(realform::RealForm&, const complexredgp_io::Interface&,
		      tags::DualTag)
    throw(error::InputError);

  void getInteractive(realform::RealForm&, const complexredgp_io::Interface&,
		      const realform::RealFormList&, tags::DualTag)
    throw(error::InputError);

  void getInteractive(realredgp::RealReductiveGroup&,
		      complexredgp_io::Interface&)
    throw(error::InputError);

  void getInteractive(complexredgp_io::Interface&)
    throw(error::InputError);

  void getInteractive(unsigned long&, const char*, unsigned long)
    throw(error::InputError);

  void getInteractive(unsigned long&, const char*, const bitmap::BitMap&,
		      input::InputBuffer* linep = 0)
    throw(error::InputError);

  input::InputBuffer& inputLine();
}

/******** type definitions ***************************************************/

namespace ioutils {

class OutputFile {
 private:
  std::ostream* d_stream;
  bool d_foutput;
 public:
  OutputFile() throw(error::InputError);
  ~OutputFile();
  template<typename T> std::ostream& operator<< (const T& arg)
    {return *d_stream << arg;}
  operator std::ostream& () {return *d_stream;}
};


class InputFile {
 private:
  std::ifstream* d_stream;
 public:
  InputFile(std::string prompt,
            std::ios_base::openmode mode
	      =std::ios_base::in | std::ios_base::binary)
    throw(error::InputError);
  ~InputFile();
  operator std::ifstream& () {return *d_stream;}
};

}

}

#endif
