/*
  This is preorder.h

  Copyright (C) 2025 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef PREORDER_H  /* guard against multiple inclusions */
#define PREORDER_H

/*
  This unit implements a simple class in the spirit of |Poset|, to represent a
  preorder on an initial segment of the natural numbers, or more generally a
  transitive relation (where neither reflexivity nor irreflexivity is imposed).

  As in |poset| we use bitmaps for every element to record the elements
  reachable from it. Being transitive is not a class invariant: one can
  introduce an arbitrary relation, and then construct the transitive closure for
  it, which is in fact the main operation implemented.

  The dense implementation and relative inefficiency of the closure operation
  (inevitable due to the sheer amount of pairs that need consideration) make
  this class mostly suitable for relatively small base sets; the
  |graph::OrientedGraph| provides algorithms to get information about the
  preorder defined by such a graph, such as the cells of mutual reachability,
  without generating that preorder explitily.
*/

#include <vector>
#include "bitmap.h"
#include "sl_list.h"
#include "graph_fwd.h"

namespace atlas {

namespace preorder {


class Preorder
{
public:
  using Elt = unsigned int; using List = containers::sl_list<Elt>;

private: // data
  struct entry
  {
    bitmap::BitMap out; // elements reachable from here
    Elt ref; // if not |nil| then an equivalent element
    entry () : out(), ref(-1) {}
    entry (Elt i, bitmap::BitMap&& out) : out(std::move(out)),ref(i) {}
  };
  std::vector<entry> data; //
  bool ready;

public:
  explicit Preorder(Elt n=0);
  explicit Preorder(const graph::OrientedGraph& g);

  // Construct a |Preorder| from a generating relation
  template<typename I> // input iterator to some container of |Elt|
  explicit Preorder(I begin, I end)
    : data (std::distance(begin,end)), ready(false)
  {
    Elt i=0;
    for (auto it=begin; it!=end; ++it,++i)
      data[i] = entry(i,bitmap::BitMap(size(),it->begin(),it->end()));
  }

  // accessors
  Elt size() const { return data.size(); }
  Elt root(Elt x) const { while (data[x].ref!=x) x = data[x].ref; return x; }

  bool relates(Elt x, Elt y) const
  { return x<data.size() and data[root(x)].out.isMember(y); }

  containers::sl_list<List> cliques() const; // in transitive closure graph

  Preorder reversed () const;
  // manipulators
  Preorder& closure (); // compute transitive closure and return |*this|
}; // |class Preorder|

} // |namespace preorder|

} // |namespace atlas|

#endif
