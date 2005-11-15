/*
  This is basic_io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef BASIC_IO_H  /* guard against multiple inclusions */
#define BASIC_IO_H

#include <iosfwd>
#include <vector>

#include "abelian_fwd.h"
#include "bitset_fwd.h"
#include "latticetypes_fwd.h"
#include "rootdata_fwd.h"
#include "weyl_fwd.h"

#include "lietype.h"

/******** function declarations *********************************************/

namespace atlas {

namespace basic_io {

// types from abelian
std::ostream& operator<< (std::ostream&, const abelian::GrpArr&);

// types from bitset
template<size_t d>
  std::ostream& operator<< (std::ostream&, const bitset::BitSet<d>&);

// types from bitvector
template<size_t dim>
  std::ostream& operator<< (std::ostream&, const bitvector::BitVector<dim>&);

// types from latticetypes
std::ostream& operator<< (std::ostream&, const latticetypes::LatticeElt&);

// types from lietype
std::ostream& operator<< (std::ostream& strm, 
			  const lietype::SimpleLieType& slt);

std::ostream& operator<< (std::ostream&, const lietype::LieType&);

// types from weyl
std::ostream& operator<< (std::ostream&, const weyl::WeylWord&);

// other functions
template<typename I>
std::ostream& seqPrint(std::ostream&, const I&, const I&, 
		       const char* sep = ",", const char* pre = "", 
		       const char* post = "");

}

}

#include "basic_io_def.h"

#endif
