/*
  This is dynkin.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef DYNKIN_H  /* guard against multiple inclusions */
#define DYNKIN_H

#include <map>

#include "latticetypes_fwd.h"

#include "bitset.h"
#include "lietype.h"
#include "setutils.h"

/******** type declarations *************************************************/

namespace atlas {

namespace dynkin {

typedef std::pair<size_t, size_t> Edge;
typedef unsigned Multiplicity;

class DynkinDiagram;

}

/******** function declarations *********************************************/

namespace dynkin {

  void bourbaki(setutils::Permutation&, const DynkinDiagram&);

  void components(bitset::RankFlagsList&, const DynkinDiagram&);

  void components(bitset::RankFlagsList&, const DynkinDiagram&,
		  const bitset::RankFlags&); // to be implemented

  void lieType(lietype::LieType&, const latticetypes::LatticeMatrix&);

  void normalize(setutils::Permutation&, const DynkinDiagram&);

}

/******** type definitions **************************************************/

namespace dynkin {

class DynkinDiagram {

 private:

  std::vector<bitset::RankFlags> d_star;
  std::map<Edge,Multiplicity> d_label;

 public:

// constructors and destructors

  DynkinDiagram() {}

  explicit DynkinDiagram(const latticetypes::LatticeMatrix&);

  DynkinDiagram(const bitset::RankFlags&, const DynkinDiagram&);

  ~DynkinDiagram() {}

// accessors

  bitset::RankFlags component(size_t) const;

  bitset::RankFlags extremities() const;

  bool isConnected() const;

  bool isSimplyLaced() const;

  bool isString() const;

  Edge labelEdge() const;

  Multiplicity maxMultiplicity() const;

  size_t node() const;

  size_t rank() const {
    return d_star.size();
  }

  const bitset::RankFlagsList& star() const {
    return d_star;
  }

  const bitset::RankFlags& star(size_t j) const {
    return d_star[j];
  }

// manipulators
};

}

}

#endif
