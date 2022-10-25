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

#include "../Atlas.h"

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

// The Multiplicity of an Edge should be 1, 2, or 3.
typedef unsigned char Multiplicity;

class DynkinDiagram;

/******** function declarations *********************************************/

  // the following two functions are main entry points to this module
  // find (semisimple) Lie type given by Cartan matrix
  LieType Lie_type(const int_Matrix& cm);

  // same, also set |pi| to permutation "straightening" each diagram component
  LieType Lie_type(const int_Matrix& cm,
		   bool Bourbaki, // whether Bourbaki order is requested
		   bool check, // if true, be prepared for arbitrary |cm|
		   Permutation& pi);


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
  std::vector<std::pair<Edge,Multiplicity> > d_downedge;

  containers::sl_list<comp_info> comps;

 public:

// constructors and destructors

  DynkinDiagram() {} // allow classes to contain a diagram, and assign it late
  explicit DynkinDiagram(const int_Matrix& Cartan);

// accessors

  unsigned int rank() const { return d_star.size(); }

  // the next two attributes are only valid after calling |classisfy()|
  LieType type() const;
  Permutation perm() const;

  bool is_simply_laced() const { return edge_label()==1; }

  int Cartan_entry(unsigned int i,unsigned int j) const;
  Multiplicity edge_multiplicity(unsigned int i,unsigned int j) const
    { return Cartan_entry(i,j)*Cartan_entry(j,i); }
  bool are_adjacent(unsigned int i, unsigned int j) const
    { return star(i).test(j); }

  const containers::sl_list<comp_info>& components() const { return comps; }

  // these are mostly internal methods, but cannot be private in current setup
  RankFlags extremities() const; // find nodes of degree 1 or 0
  Edge labelled_edge() const; // find edge with label (not 1), assumed present
  unsigned int fork_node() const; // find node of degree>2, or rank() if none
  RankFlags star(unsigned int j) const { return d_star[j]; } // copies it

  // transformers (return related Dynkin diagrams)
  DynkinDiagram subdiagram(RankFlags selection) const;
  DynkinDiagram folded(const ext_gens& orbits) const;

  LieType classify_semisimple(Permutation& pi) const
  { pi=perm(); return type(); }

// manipulators:

  friend LieType Lie_type(const int_Matrix&); // uses |component_kind()|
private: // functions that are only meaningful for connected diagrams

  Multiplicity edge_label() const; // (maximal) label of edge, or 1 if none
  lietype::TypeLetter component_kind() const;
  SimpleLieType classify_simple(Permutation& pi) const;
  DynkinDiagram& classify(const int_Matrix& Cartan);
  void classify(const int_Matrix& Cartan,comp_info& ci) const;
}; // |class DynkinDiagram|

} // |namespace dynkin|
} // |namespace atlas|

#endif
