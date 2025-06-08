/*
  This is preorder.cpp.

  This file contains a simple-minded implementation of a preorder structure,
  where one explicitly keeps the bit-matrix of the order relation.

  Copyright (C) 2025 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include <cassert>
#include "preorder.h"
#include "graph.h"

namespace atlas {

namespace preorder {

Preorder::Preorder(Elt n)
  : data(n)
  , ready(true) // empty relation is vacuously transitive
{ for (Elt i=0; i<n; ++i)
    data[i] = entry(i,bitmap::BitMap(n));
}

Preorder::Preorder(const graph::OrientedGraph& g)
  : data(g.size()), ready(false)
{
  for (Elt i=0; i<size(); ++i)
  {
    const auto& L=g.edgeList(i);
    data[i]=entry(i,bitmap::BitMap(size(),L.begin(),L.end()));
  }
}

Preorder Preorder::reversed() const
{
  Preorder result(size());
  for (Elt i=0; i<size(); ++i)
    for (Elt x : data[i].out)
      result.data[x].out.insert(i);
  result.ready = ready;
  if (ready) // then equivalence classe are preserved
    for (Elt i=0; i<size(); ++i)
      result.data[i].ref = data[i].ref;
  return result;
}

/*
  Computing the closure can simply exhaust the operations of |m |= out[x]| where
  |m| is a bitmap entry of |out| and |m.isMember(x)|. However this can take a
  long time to settle large cliques, since we must repeatedly loop over all
  their |out| entries until each of them stabilises, which involves very
  repetitive and not very useful work. But we can detect loops that lead to such
  cliques early, and we know that all |out| bitmaps will eventually be equal, so
  we use the |ref| auxiliary field to record this, and store only one bitmap for
  each group of currently equivalenced entries. Indeed, the set |ref| links
  define a rooted tree for each such group, the references pointing to the root
  where the bitmap for the group is held. The tree structure evolves to a star
  becuase every time we traverse the tree to the root, we redirect all |ref|
  entries used to point directly to the root.
 */
Preorder& Preorder::closure ()
{
  if (ready)
    return *this;
  do
  {
    ready=true; // start optimistically
    for (Elt i=0; i<size(); ++i)
    {
      entry& d = data[i];
      if (d.ref!=i) // only isolated or group-root entries done
	continue;
      for (Elt x : d.out) // |d.out| changes in loop; only speeds up things
      {
	{
	  Elt y = data[x].ref;
	  if (y==i)
	    continue; // skip this element of the group of |i|
	  while (data[y].ref!=y)
	    y = data[y].ref;
	  while (x!=y) // short-circuit all |ref| values visited to |y|
	  { Elt& link = data[x].ref;
	    x = link; link=y;
	  }
	  assert(data[x].out.capacity()==size()); // bitmap at |x| still present
	} // now |x| has caught up with |y| we can forget the latter
	if (x==i) // it was indirectly an element of the group of |i|
	  continue;
	if (data[x].out.isMember(i)) // then we found a loop linking |i| and |x|
	{
	  ready = ready and
	    data[x].out==data[i].out; // whether no new relations are added
	  data[x].out |= data[i].out; // transfer |i| to |x|
	  data[i].ref = x;
	  data[i].out = bitmap::BitMap(0); // no longer use this one
	  break;
	}
	if (not data[i].out.contains(data[x].out))
	{
	  ready=false;
	  data[i].out |= data[x].out;
	}
      } // |for (x)|
    } // |for (i)|
  } // |do|
  while (not ready);
  return *this;
}

auto Preorder::cliques() const -> containers::sl_list<List>
{
  containers::sl_list<List> result;
  assert(ready); // so we can be a |const| method
  for (Elt x=0; x<size(); ++x)
    if (data[x].ref==x) // only consider root elements of cliques
    {
      if (not data[x].out.isMember(x))
	continue; // without a self-loop, there is no actual clique
      auto& L=result.push_back(List{}); // no need to push |x|, is in |out[x]|
      for (auto y : data[x].out)
	if (root(y)==x)
	  L.push_back(y);
    }
  return result;
}


} // |namespace preorder|

} // |namespace atlas|
