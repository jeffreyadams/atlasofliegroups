/*!
\file
\brief Class definitions and function declarations for DynkinDiagram
*/
/*
  This is dynkin.h
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

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

  /*!
  Since a Dynkin diagram will be a bitset (a subset of RANK_MAX
  elements, the possible vertices 0, 1,...,RANK_MAX-1), an Edge is a
  pair of numbers (between 0 and RANK_MAX-1).
  */
typedef std::pair<size_t, size_t> Edge;

  /*!
  The Multiplicity of an Edge should be 1, 2, or 3.
  */
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
  /*!
  \brief A Dynkin diagram of at most RANK_MAX vertices.

  The collection of vertices is represented as a BitSet d_star, regarded as a
  subset of the integers 0,1,...,RANK_MAX-1.  These are the vertices.
  An Edge is a pair of vertices.  The map d_label attaches to each of
  certain Edges an unsigned integer Multiplicity.
  */
class DynkinDiagram {

 private:

  /*!
  Subset of the integers 0,1,...,RANK_MAX-1 corresponding to the
  vertices of the Dynkin diagram.
  */
  std::vector<bitset::RankFlags> d_star;

  /*!
  Map from certain pairs of vertices (the edges of the Dynkin diagram)
  to unsigned integers (their multiplicities).
  */
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
