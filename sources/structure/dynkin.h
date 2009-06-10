/*!
\file
\brief Class definitions and function declarations for DynkinDiagram
*/
/*
  This is dynkin.h
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef DYNKIN_H  /* guard against multiple inclusions */
#define DYNKIN_H

#include <utility> // for |std::pair|
#include <vector>

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

} // |namespace dynkin|

/******** function declarations *********************************************/

namespace dynkin {

  // find (semisimple) Lie type given by Cartan matrix
  lietype::LieType Lie_type(const latticetypes::LatticeMatrix& cm);

  // same, also set |pi| to permutation "straightening" each diagram component
  lietype::LieType Lie_type(const latticetypes::LatticeMatrix& cm,
			    bool Bourbaki, bool check,
			    setutils::Permutation& pi);

  bitset::RankFlagsList components(const DynkinDiagram& d);

  setutils::Permutation normalize(const DynkinDiagram& d);

  setutils::Permutation bourbaki(const DynkinDiagram& d);

}

/******** type definitions **************************************************/

namespace dynkin {
  /*!
  \brief A Dynkin diagram of at most RANK_MAX vertices.

  The collection of vertices is represented as vector, each of whose elements
  is a BitSet d_star[i], which regarded as a subset of the integers
  0,1,...,RANK_MAX-1 gives the neighbours of vertex i.
  An Edge is a pair of indices designating vertices. The map d_label attaches
  to certain Edges an unsigned integer Multiplicity.
  */
class DynkinDiagram {

 private:

  /*!
  Adjacency relations: d_star[i].test(j) tells whether i and j are neighbors
  */
  std::vector<bitset::RankFlags> d_star;

  /*!
  \brief List of edges from longer to shorter node, with edge label (2 or 3)
  */
  std::vector<std::pair<Edge,Multiplicity> > d_downedge;

 public:

// constructors and destructors

  DynkinDiagram() {}

  explicit DynkinDiagram(const latticetypes::LatticeMatrix&);

  DynkinDiagram(const bitset::RankFlags&, const DynkinDiagram&);

  ~DynkinDiagram() {}

// accessors

  int cartanEntry(size_t i,size_t j) const;

  bitset::RankFlags component(size_t) const;

  Multiplicity edgeMultiplicity(size_t i,size_t j) const
    { return cartanEntry(i,j)*cartanEntry(j,i); }

  bitset::RankFlags extremities() const;

  bool isConnected() const;

  bool isSimplyLaced() const;

  bool isString() const;

  size_t num_labels() const { return d_downedge.size(); }
  Edge labelEdge(size_t i=0) const ;

  Multiplicity maxMultiplicity() const;

  size_t node() const;

  size_t rank() const {
    return d_star.size();
  }

  const bitset::RankFlagsList& star() const { return d_star; }

  bitset::RankFlags star(size_t j) const { return d_star[j]; }

// manipulators
};

}

}

#endif
