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
  been dealt with, say x0. Starting from x0, we are going to examine all
  vertices y thaat are visible from x0 (i.e., >= x0). There will be a stack of
  vertices, corresponding to a path originating from x0, such that at each
  point in time, each vertex y >= x0 which has been examined will be
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

  NOTE: I believe that in this process it is not guaranteed that x_min is
  the actual minimal element of the stack equivalent to x (this is also why
  we cannot take off a class simply by looking at x_min). The real 
  interpretation of x_min is: the smallest element of the stack that is visible
  from y using edges already processed. The important issue here is whether
  x_min < x or not, and for this we only need to examine whether there is
  a strict succesor y of x with y_min < x.
*/

{
  using namespace bitmap;

  BitMap b(size());
  std::vector<Vertex> min(size(),size());

  pi.resize(size());

  std::stack<Vertex> v;
  std::stack<size_t> ecount;

  for (Vertex x0 = 0; x0 < size(); ++x0) {

    if (b.isMember(x0)) /* x0 is dealt with */
      continue;

    v.push(x0);
    ecount.push(0);
    min[x0] = 0;
    size_t t = 1; // stack position

    while(not v.empty()) {

      Vertex x = v.top();
      Vertex y;
      const EdgeList& e = edgeList(x);

      for (; ecount.top() < e.size(); ++ecount.top()) {
	y = e[ecount.top()];
	if (b.isMember(y))
	  continue;
	if (min[y] == size()) // y is new
	  goto add_to_stack;
	if (min[x] > min[y])
	  min[x] = min[y];
      }

      // at this point we have exhausted the edges of x
      v.pop();
      ecount.pop();
      if (min[x] == t-1) { // take off class
	getClass(*this,x,b,pi,p);
      }
      --t;
      continue;

    add_to_stack:
      v.push(y);
      ecount.push(0);
      min[y] = t;
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
