/*
  This is graph.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "graph.h"

#include <set>
#include <stack>

#include "bitmap.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- The OrientedGraph class

  ... explain here when it is stable ...

******************************************************************************/

namespace {

using namespace graph;
  
void getClass(const OrientedGraph&, Vertex, bitmap::BitMap&,
	      partition::Partition&, OrientedGraph* p);

}

namespace graph {

/******** accessors **********************************************************/
void OrientedGraph::cells(partition::Partition& pi, OrientedGraph* p) const

/*
  Synopsis: puts in pi the partition of the vertex set into left cells, and,
  if p != 0, puts in p the graph of the order relation induced on the quotient.

  Explanation: define a preorder relation on the vertices by setting x <= y 
  iff there is an oriented path from x to y. This function puts in pi the 
  partition function corresponding to the equivalence classes of this preorder.
  We use the Tarjan algorithm, explained in one of the Knuth books, but which 
  I learned from Bill Casselman.

  The vertices for which the partition function is already defined will
  be said to be dealt with. These are marked off in an auxiliary bitmap.
  The algorithm goes as follows. Start with the first vertex that has not
  been dealt with, say x0. We will have to deal (potentially) with the set
  of all vertices visible from x0 (i.e >= x0). There will be a stack of
  vertices, corresponding to a path originating from x0, such that at each
  point in time, each vertex y >= x0 which is not dealt with will be
  equivalent to an element of the stack; the least such (in the stack
  ordering) will be called y_min, so that for instance x0_min = 0. We
  record these values in an additional table, initialized to some inaccessible
  value (here we use the size of the poset)

  Now let x be at the top of the stack. Look at the edges originating from
  x. Ignore the ones that go to vertices which are dealt with. If there
  are no edges left, x is minimal and is a class by itself; we can take
  it off, mark it as dealt with, and continue. Otherwise, run through the
  edges one by one. Let the edge be x->y. If y_min is defined, this means
  that y has already been examined, and is not dealt with. But each such
  element is equivalent to an element in the active stack, so y_min should
  be one of the elements in the active stack, hence x is visible from y_min:
  in other words, x and y are equivalent, and we set x_min = y_min if
  y_min < x_min. Otherwise, y is seen for the first time; then we just put
  it on the stack. When we are done with the edges of x, we have now the
  value of x_min which is the inf over the edges originating from x of
  the y_min. If this value is equal to the stack-position of x, we see that x 
  is minimal in its class, and we get a new class by taking all the successors
  of x not already dealt with. We then move to the parent of x and continue
  the process there.
*/

{
  using namespace bitmap;

  BitMap b(size());
  std::vector<Vertex> min(size(),size());

  pi.resize(size());

  std::vector<Vertex> v(1);
  std::vector<const EdgeList*> elist(1);
  std::vector<size_t> ecount(1);

  for (Vertex x = 0; x < size(); ++x) {

    if (b.isMember(x)) /* x is dealt with */
      continue;

    v[0] = x;
    elist[0] = &edgeList(x);
    ecount[0] = 0;
    min[x] = 0;
    size_t t = 1; // stack position

    while(t) {

      Vertex y = v[t-1];
      Vertex z;
      const EdgeList& e = *elist[t-1];

      for (; ecount[t-1] < e.size(); ++ecount[t-1]) {
	z = e[ecount[t-1]];
	if (b.isMember(z))
	  continue;
	if (min[z] == size()) // z is new
	  goto add_path;
	if (min[y] > min[z])
	  min[y] = min[z];
      }

      // at this point we have exhausted the edges of y
      if (min[y] == t-1) { // take off class
	getClass(*this,y,b,pi,p);
      }
      else if (min[y] < min[v[t-2]]) // if t=1, previous case holds
	min[v[t-2]] = min[y];
      --t;
      continue;

    add_path:
      v.resize(t+1);
      elist.resize(t+1);
      ecount.resize(t+1);
      v[t] = z;
      elist[t] = &edgeList(z);
      ecount[t] = 0;
      min[z] = t;
      ++t;
    }
  }

  return;
}

}

/*****************************************************************************

        Chapter II -- Functions local to this module

  ... explain here when it is stable ...

******************************************************************************/

namespace {

void getClass(const OrientedGraph& g, Vertex y, bitmap::BitMap& b,
	      partition::Partition& pi, OrientedGraph* p)

/*
  Synopsis: marks off the class of y in b, and, if p !=0, extends the poset
  p to accomodate the new vertex

  Precondition: y is the minimal element in its class, in terms of the
  numbering of the vertices;

  After the element y has been identified as minimal among the elements not
  already marked in b, this function marks off the equivalence class of y;
  these are just the elements visible from y and not already marked in b.
  The class is also written as a new class in pi.

  NOTE: the new class is always maximal in the poset-so-far, so only one
  new row of edges needs to be added.
*/

{
  std::stack<Vertex> toDo;

  toDo.push(y);
  b.insert(y);
  pi.newClass(y);
  unsigned long c = pi(y);

  std::set<Vertex> a;

  while (not toDo.empty()) {

    Vertex x = toDo.top();
    toDo.pop();
    const EdgeList& e = g.edgeList(x);

    for (size_t j = 0; j < e.size(); ++j) {
      Vertex z = e[j];
      if (b.isMember(z)) {
	if (p and (pi(z) < c)) /* add a new edge to p */
	  a.insert(pi(z));
	continue;
      } else {
	toDo.push(z);
	b.insert(z);
	pi.addToClass(c,z);
      }
    }
  }

  if (p) { // write new edges
    p->resize(c+1);
    EdgeList& f = p->edgeList(c);
    std::copy(a.begin(),a.end(),back_inserter(f));
  }

  return;
}

}

}
