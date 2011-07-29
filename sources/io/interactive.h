/*
  This is interactive.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006--2011 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef INTERACTIVE_H  /* guard against multiple inclusions */
#define INTERACTIVE_H

#include "atlas_types.h"

#include <string>
#include <ios>

#include "error.h"
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

  void bitMapPrompt(std::string&, const char*, const BitMap&);

  size_t get_Cartan_class(const BitMap& cs) throw(error::InputError);

  void getInteractive(LieType&) throw(error::InputError);

  void getInteractive(InnerClassType&, const LieType&)
    throw(error::InputError);

  void getInteractive(PreRootDatum&,
		      WeightList&,
		      const LieType&) throw(error::InputError);

  void getInteractive(RealFormNbr&, const complexredgp_io::Interface&)
    throw(error::InputError);

  void getInteractive(RealFormNbr&, const complexredgp_io::Interface&,
		      tags::DualTag)
    throw(error::InputError);

  void getInteractive(RealFormNbr&, const complexredgp_io::Interface&,
		      const RealFormNbrList&, tags::DualTag)
    throw(error::InputError);

  RealReductiveGroup getRealGroup(complexredgp_io::Interface&)
    throw(error::InputError);

  void getInteractive(ComplexReductiveGroup*&,
		      complexredgp_io::Interface*&)
    throw(error::InputError);

  unsigned long get_bounded_int(input::InputBuffer& ib,
				const char* prompt,
				unsigned long limit)
    throw(error::InputError);

  unsigned long get_int_in_set(const char* prompt,
			       const BitMap& choices,
			       input::InputBuffer* linep = 0)
    throw(error::InputError);

  Weight get_weight(input::InputBuffer& ib,
				  const char* prompt,
				  size_t rank)
    throw(error::InputError);

  RatWeight get_ratweight(input::InputBuffer& ib,
					const char* prompt,
					size_t rank)
    throw(error::InputError);

  standardrepk::StandardRepK get_standardrep
    (const standardrepk::SRK_context& c)
    throw(error::InputError);

  repr::StandardRepr get_repr(const repr::Rep_context& c)
    throw(error::InputError);

  input::InputBuffer& common_input();
  input::InputBuffer& sr_input();

  void get_parameter(RealReductiveGroup& GR,
		     KGBElt& x,
		     Weight& lambda_rho,
		     RatWeight& gamma);

  void getInteractive(atlas::Parabolic &psg, size_t rank) throw(error::InputError);

}

/******** type definitions ***************************************************/

namespace ioutils {

class OutputFile
{
  std::ostream* d_stream;
  bool d_foutput;
 public:
  OutputFile() throw(error::InputError);
  ~OutputFile();
  template<typename T> std::ostream& operator<< (const T& arg)
    {return *d_stream << arg;}
  operator std::ostream& () {return *d_stream;}
}; // |class OutputFile|


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
