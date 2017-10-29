/*
  This is interactive.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006--2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef INTERACTIVE_H  /* guard against multiple inclusions */
#define INTERACTIVE_H

#include "../Atlas.h"

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

  // most of these functions may |throw error::InputError|, however exception
  // specifications are not given (we don't claim absence of other exceptions)
  std::string getFileName(const std::string& prompt);

  bool open_binary_file(std::ofstream& block_out,const std::string& prompt);

  bool get_yes_or_no(const char* prompt);

  void bitMapPrompt(std::string&, const char*, const BitMap&);

  size_t get_Cartan_class(const BitMap& cs);

  void get_group_type(InnerClass*&, output::Interface*&,
     lietype::Layout& layout, LatticeMatrix& basis);

  void getInteractive(LieType&);

  PreRootDatum get_pre_root_datum(LatticeMatrix& basis, const LieType& lt);
  WeightInvolution
    getInnerClass(lietype::Layout& lo, const LatticeMatrix& basis);
  void getInteractive(InnerClassType&, const LieType&);

  RealFormNbr get_real_form(output::Interface&);

  RealFormNbr get_dual_real_form(output::Interface&,
				 const InnerClass& G,
				 RealFormNbr rf);


  unsigned long get_bounded_int(input::InputBuffer& ib,
				const char* prompt,
				unsigned long limit);

  unsigned long get_int_in_set(const char* prompt,
			       const BitMap& choices,
			       input::InputBuffer* linep = 0);

  Weight get_weight(input::InputBuffer& ib, const char* prompt, size_t rank);

  RatWeight get_ratweight
    (input::InputBuffer& ib, const char* prompt, size_t rank);

  standardrepk::StandardRepK get_standardrep
    (const standardrepk::SRK_context& c);

  StandardRepr get_repr(const Rep_context& c);

  input::InputBuffer& common_input();
  input::InputBuffer& sr_input();

  SubSystemWithGroup get_parameter(RealReductiveGroup& GR,
				   KGBElt& x,
				   Weight& lambda_rho,
				   RatWeight& gamma);

  void getInteractive(atlas::Parabolic &psg, size_t rank);

// get second distinguished involution that commutes with the inner class one
  WeightInvolution get_commuting_involution
    (const lietype::Layout& lo, const LatticeMatrix& basis);

}

/******** type definitions ***************************************************/

namespace ioutils {

// the next class emulates an (impossible) derivation from std::ostream&
class OutputFile
{
  std::ostream* d_stream;
  bool d_foutput;
 public:
  OutputFile();
  ~OutputFile();
  template<typename T> std::ostream& operator<< (const T& arg)
    {return *d_stream << arg;}
  operator std::ostream& () {return *d_stream;}
  bool is_std_cout () const { return not d_foutput; }
}; // |class OutputFile|


class InputFile {
 private:
  std::ifstream* d_stream;
 public:
  InputFile(std::string prompt,
            std::ios_base::openmode mode
	      =std::ios_base::in | std::ios_base::binary);
  ~InputFile();
  operator std::ifstream& () {return *d_stream;}
};

}

}

#endif
