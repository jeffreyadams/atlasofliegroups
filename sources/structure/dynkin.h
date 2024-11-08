/*
  This is dynkin.h
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2014,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Class definitions and function declarations for DynkinDiagram

#ifndef DYNKIN_H  /* guard against multiple inclusions */
#define DYNKIN_H

#include <utility> // for |std::pair|
#include <vector>

#include "Atlas.h"

#include "bitset.h"
#include "lietype.h"
#include "permutations.h"
#include "sl_list.h"

namespace atlas {
namespace dynkin {

/******** type declarations *************************************************/

/*
  Since a Dynkin diagram will be a bitset (a subset of RANK_MAX elements, the
  possible vertices 0, 1,...,RANK_MAX-1), an Edge is a pair of numbers
  (between 0 and RANK_MAX-1), which numbers will fit in |unsigned char|
 */
typedef std::pair<unsigned char, unsigned char> Edge;

class DynkinDiagram;

/******** function declarations *********************************************/

  // the following two functions are main entry points to this module
  // find (semisimple) Lie type given by Cartan matrix
  LieType Lie_type(const int_Matrix& cm);

  // same, also set |pi| to permutation "straightening" each diagram component
  LieType Lie_type(const int_Matrix& cm, Permutation& pi);

  void permute(const Permutation& pi, DynkinDiagram& D);

/******** type definitions **************************************************/

/*
  A Dynkin diagram of at most RANK_MAX vertices.

  The collection of vertices is represented as vector, each of whose elements
  is a BitSet d_star[i], which regarded as a subset of the integers
  0,1,...,RANK_MAX-1 gives the neighbours of vertex i. An Edge is a pair of
  indices designating vertices. The map d_label attaches to certain Edges an
  unsigned integer Multiplicity.
*/
class DynkinDiagram
{
 public:

  struct comp_info
  {
    lietype::TypeLetter type; unsigned char offset;
    RankFlags support;
    std::vector<unsigned char> position;

    // only initialise |support| to singleton; |DinkinDiagram| ctor does rest
  comp_info(unsigned i) : type('T'), offset(-1), support(1u<<i), position() {}
    unsigned rank() const { return position.size(); }
  };

 private:
  // Adjacency relations: |d_star[i].test(j)|: whether |i| and |j| are neighbors
  std::vector<RankFlags> d_star;

  // List of edges from longer to shorter node, with edge label (2 or 3)
  sl_list<std::pair<Edge,short int> > down_edges;

  sl_list<comp_info> comps;

 public:

// constructors and destructors

  DynkinDiagram() {} // allow classes to contain a diagram, and assign it late
  explicit DynkinDiagram(const int_Matrix& Cartan);

// accessors

  unsigned int rank() const { return d_star.size(); }

  LieType type() const;
  Permutation perm() const;

  bool is_simply_laced() const { return down_edges.empty(); }

  int Cartan_entry(unsigned int i,unsigned int j) const;
  int edge_multiplicity(unsigned int i,unsigned int j) const
    { return Cartan_entry(i,j)*Cartan_entry(j,i); }
  bool are_adjacent(unsigned int i, unsigned int j) const
    { return star(i).test(j); }

  const containers::sl_list<comp_info>& components() const { return comps; }

  // these are mostly internal methods, but cannot be private in current setup
  RankFlags star(unsigned int j) const { return d_star[j]; } // copies it

  // transformer returns related Dynkin diagram for delta-fixed weight space
  // (only Coxeter diagram matters for us, but this maps A5->B3, A6->B3, D5->C4)
  DynkinDiagram folded(const ext_gens& orbits) const;

// manipulators:

private: // functions that are only meaningful for connected diagrams

  void classify(const int_Matrix& Cartan);
  void classify(const int_Matrix& Cartan,comp_info& ci) const;
}; // |class DynkinDiagram|

} // |namespace dynkin|
} // |namespace atlas|

#endif
