/*
  This is interactive.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

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

namespace atlas {

/******** function declarations **********************************************/

namespace interactive {

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

}

#endif
