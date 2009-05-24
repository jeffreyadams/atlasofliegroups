/*
  This is interactive.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef INTERACTIVE_H  /* guard against multiple inclusions */
#define INTERACTIVE_H

#include "interactive_fwd.h"

#include "bitmap_fwd.h"
#include "complexredgp_fwd.h"
#include "complexredgp_io_fwd.h"
#include "prerootdata_fwd.h"
#include "realredgp_fwd.h"
#include "rootdata_fwd.h"

#include "error.h"
#include "input.h"
#include "latticetypes.h"
#include "lietype.h"
#include "realform.h"
#include "standardrepk.h"
#include "tags.h"

namespace atlas {

/******** type declarations **************************************************/

/* the strange namespace is an historic artifact; the class was moved */

namespace ioutils {

  class OutputFile;

}

/******** function declarations **********************************************/

namespace interactive {

std::string getFileName(const std::string& prompt)
    throw(error::InputError);

bool open_binary_file(std::ofstream& block_out,const std::string& prompt);

 void bitMapPrompt(std::string&, const char*, const bitmap::BitMap&);

 size_t get_Cartan_class(const bitmap::BitMap& cs) throw(error::InputError);

  void getInteractive(lietype::LieType&) throw(error::InputError);

  void getInteractive(lietype::InnerClassType&, const lietype::LieType&)
    throw(error::InputError);

  void getInteractive(prerootdata::PreRootDatum&,
		      latticetypes::WeightList&,
		      const lietype::LieType&) throw(error::InputError);

  void getInteractive(realform::RealForm&, const complexredgp_io::Interface&)
    throw(error::InputError);

  void getInteractive(realform::RealForm&, const complexredgp_io::Interface&,
		      tags::DualTag)
    throw(error::InputError);

  void getInteractive(realform::RealForm&, const complexredgp_io::Interface&,
		      const realform::RealFormList&, tags::DualTag)
    throw(error::InputError);

  realredgp::RealReductiveGroup getRealGroup(complexredgp_io::Interface&)
    throw(error::InputError);

  void getInteractive(complexredgp::ComplexReductiveGroup*&,
		      complexredgp_io::Interface*&)
    throw(error::InputError);

  unsigned long get_bounded_int(input::InputBuffer& ib,
				const char* prompt,
				unsigned long limit)
    throw(error::InputError);

  void getInteractive(unsigned long& v, const char* prompt,
		      const bitmap::BitMap& choices,
		      input::InputBuffer* linep = 0)
    throw(error::InputError);

  latticetypes::Weight get_weight(input::InputBuffer& ib,
				  const char* prompt,
				  size_t rank)
    throw(error::InputError);

  standardrepk::StandardRepK get_standardrep
    (const standardrepk::KhatContext& khc)
    throw(error::InputError);

  input::InputBuffer& common_input();
  input::InputBuffer& sr_input();

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
