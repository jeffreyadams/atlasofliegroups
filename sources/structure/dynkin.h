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

#include "atlas_types.h"
#include "lietype_fwd.h"

#include "bitset.h"
#include "permutations.h"

/******** type declarations *************************************************/

namespace atlas {

namespace dynkin {

  /*!
  Since a Dynkin diagram will be a bitset (a subset of RANK_MAX
  elements, the possible vertices 0, 1,...,RANK_MAX-1), an Edge is a
  pair of numbers (between 0 and RANK_MAX-1).
  */
typedef std::pair<unsigned char, unsigned char> Edge;

  /*!
  The Multiplicity of an Edge should be 1, 2, or 3.
  */
typedef unsigned char Multiplicity;

class DynkinDiagram;

} // |namespace dynkin|

/******** function declarations *********************************************/

namespace dynkin {

  // find (semisimple) Lie type given by Cartan matrix
  LieType Lie_type(const int_Matrix& cm);

  // same, also set |pi| to permutation "straightening" each diagram component
  LieType Lie_type(const int_Matrix& cm,
		   bool Bourbaki, bool check,
		   permutations::Permutation& pi);

  RankFlagsList components(const DynkinDiagram& d);

  permutations::Permutation normalize(const DynkinDiagram& d);

  permutations::Permutation bourbaki(const DynkinDiagram& d);

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
  std::vector<RankFlags> d_star;

  /*!
  \brief List of edges from longer to shorter node, with edge label (2 or 3)
  */
  std::vector<std::pair<Edge,Multiplicity> > d_downedge;

 public:

// constructors and destructors

  DynkinDiagram() {}

  explicit DynkinDiagram(const int_Matrix& Cartan);

  DynkinDiagram(const RankFlags& selection, const DynkinDiagram& d);

  ~DynkinDiagram() {}

// accessors

  int cartanEntry(size_t i,size_t j) const;

  RankFlags component(size_t) const;

  Multiplicity edgeMultiplicity(size_t i,size_t j) const
    { return cartanEntry(i,j)*cartanEntry(j,i); }

  RankFlags extremities() const;

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

  const RankFlagsList& star() const { return d_star; }

  RankFlags star(size_t j) const { return d_star[j]; }

// manipulators
};

}

}

#endif
