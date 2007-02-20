/*!
\file
  This is graph.cpp
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "graph.h"

#include <set>
#include <stack>
#include <cassert>


#include "bitmap.h"
#include "tags.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- The OrientedGraph class

******************************************************************************/


namespace graph {

/******** accessors **********************************************************/

/*!
  OrientedGraph::cells puts in pi the partition of the vertex set into strong
  components for the graph, and if gr!=NULL also puts in *gr the graph of the
  partial order relation induced by our graph on the quotient.

  Explanation: define a preorder relation on the vertices by setting x <= y
  iff there is an oriented path from x to y. This function puts in pi the
  partition function corresponding to the equivalence classes of this
  preorder, which are called strong components of the directed graph (two
  vertices are strongly connected if there are oriented paths between them in
  both directions). We use R. E. Tarjan's algorithm, explained in Knuths's
  book The Stanford GraphBase, pp. 512-519 (and also in his future volume 4 of
  The Art of Computer Programming). Fokko explained he had learned this
  algorithm from Bill Casselman, but unfortunately his rendition was not
  correct (even though it gave the right results in the vast majority of the
  cases) so I had to adapt it, MvL.

  The strong components are going to be discovered in a highest-first
  topological order ("topological order" means some total order compatible
  with a given partial order): when a strong component is discovered all
  outgoing edges from it lead to components that were previously discovered.
  Since the elements of the already discovered components are marked, and
  edges leading to them are henceforth ignored, the problem at each moment is
  to try to isolate a strong component with no outgoing edges. Eventually we
  prefer numbering in the opposite ordering to have the highest strong
  component last, whence we perform a reversal at the end.

  The vertices will be traversed in the order of depth-first search. When they
  are first seen we say the become "active", and some data associated to them
  is entered on a local vector structure |active| that behaves somewhat like a
  stack (but we need more liberal access than |std::stack| would provide).
  When all its descendants have been traverse we shall call a vertex "mature",
  but in this case it remains in the |active structure. Finally when a strong
  component is discovered, its elements become "settled"; the information in
  |active| about them is then discarded. For each vertex |x| we keep a value
  |rank[x]| which is |0| until |x| becomes active, at which point it becomes
  the sequence number of |x| during the traversal (for once this numbering
  starts from 1), and when |x| gets settled we set |rank[x]=infinity|. Nothing
  particular happens when |x| becomes mature.

  Extra information is kept for all vertices in |active| (this infomation is
  accessible for the current vertex during traversal, but is not (easily)
  accessible given just some Vertex value). We record the Vertex value |v| of
  the current vertex, the location of its |parent| in |active| to be able to
  continue the traversal when the vertex matures, and the index |next_edge| up
  to which its edges have already been considered (for immature vertices).
  Finally and most importantly, we record for |v| the minimal sequence number
  |min| of any active vertex that is either |v| itself or can be reached in
  one step from a mature descendant |d| of |v| (to be exact, this information
  is incorporated into the |min| value of |v| once the branch of the
  depth-first search tree from |v| containing |d| matures). When a vertex |v|
  matures, its unsettled descendants are precisely the active vertices not
  older than |v|, and |v| is the head of a strong componenent if and only if
  its |min| value equals its own sequence number |rank[v]|; this means that
  none of its descendants (including |v| itself) can reach an active vertex
  older than |v|. When this happens, the unsettled descendants of |v| form the
  final segment of |active| starting at |v|.

  During the depth-first traversal, there are (apart from the settled vertices
  that are ignored) two kind of vertex encounters: the initial ones for which
  the edge leading to it becomes part of the search tree, and the encounters
  of a vertex that is already active ("cross edges"). In the former case we
  make the new vertex active, add its information, and make it current. In the
  latter case we just take note, in the value of |min| for the current vertex,
  of the fact that the other vertex can be reached . Finally when no edges for
  the current vertex |x| remain, we must test whether it heads a strong
  component, and if not update the |min| information of its parent; after that
  the parent becomes current (again). It remains to prove that the criterion
  |min==rank[v]| used for decicding whether |v| heads a strong component is
  always correct; we shall do that below.
*/

namespace {
typedef size_t seqno; // sequence number in depth-first traversal
typedef size_t work_addr; // reference to an active vertex by location

struct info
{
  Vertex v;    // identification of vertex in the graph
  work_addr parent; // location of parent in depth-first traversal
  size_t next_edge; // index in edge list of next edge to consider
  seqno min;        // minimal rank of vertex reachable as indicated above

  // constructor
  info (const Vertex x, work_addr p, seqno rank)
    : v(x), parent(p), next_edge(0), min(rank) { }
};

} // namespace

void OrientedGraph::cells(partition::Partition& pi, OrientedGraph* gr) const
{
  std::vector<seqno> rank(size(),0);

  const work_addr nil= size(); // impossible index into |active|
  const seqno infinity= size()+1; // impossible sequence number

  std::vector<info> active;

  pi.resize(size());           // |pi| will be a partition of the vertex set
  if (gr!=NULL) gr->resize(0); // start induced graph with a clean slate

  for (Vertex x0 = 0; x0 <size(); ++x0) // find all CONNECTED graph components
    if (rank[x0]<infinity)              // each time this holds gives a new one
    {
      seqno count=1;
      rank[x0]=count++;
      active.push_back(info(x0,nil,rank[x0])); // x0 has no parent
      work_addr cur_pos=0; // point current position to x0

      while(cur_pos!=nil)
	{
	next_x:
	  info& x = active[cur_pos];
	  const EdgeList& edges = edgeList(x.v);

	  while (x.next_edge < edges.size())
	    {
	      Vertex y = edges[x.next_edge++];
	      if (rank[y]==0) // y is a fresh vertex
		{
		  work_addr y_pos=active.size();
		  rank[y]=count++;
		  active.push_back(info(y,cur_pos,rank[y]));
		  cur_pos=y_pos;
		  goto next_x;
		}
	      else // |y| was seen before (cross edge), or |y| is settled
		{  // if |y| is settled nothing will happen
		  if (rank[y] < x.min) // then record that we can reach y
		    x.min = rank[y];
		}
	    } // while (x.next_edge < edges.size())

	  // at this point we have exhausted the edges of |x.v|, it matures
	  work_addr new_pos=x.parent; // we will back-up to parent of |x| next

	  if (x.min == rank[x.v]) // no older vertex reachable from |x|
	    { // split off strong component
	      pi.newClass(x.v);  // x will be added again in loop, no harm done
	      unsigned long c= pi(x.v);
	      std::vector<const EdgeList*> out; // to gather outgoing edges
	      out.reserve(active.size()-cur_pos);
	      for (work_addr i=cur_pos; i<active.size(); ++i)
		{
		  Vertex y=active[i].v; // the first time |y==x.v|
		  pi.addToClass(c,y);
		  rank[y]=infinity;     // |y| is now settled
		  out.push_back(&edgeList(y));
		}
	      // now remove |x| and its descendance from |active|
	      active.erase(active.begin()+cur_pos,active.end());

	      if (gr!=NULL) gr->addLinks(out,pi);
	    }
	  else // |x| matures but does not head a new strong component
	    {  // note that |x| cannot be |x0|, so active[x.parent] exists
	      if (x.min < active[x.parent].min) // then update parent info
		active[x.parent].min=x.min; // what x sees, its parent sees
	    }
	  cur_pos=new_pos;
	}  // while(cur_pos!=nil)

      assert(active.empty());

    } //for (x0) if (rank[x0]<infinity)

  // now reverse the numbering of the classes in the partition
  if (false)
  {
    unsigned long last=pi.classCount()-1;
    std::vector<unsigned long> f(pi.size());
    for (size_t i=0; i<f.size(); ++i) f[i]=last-pi(i); // opposite renumbering

    // replace |pi| by its reversely numbered equivalent
    partition::Partition(f,tags::UnnormalizedTag()).swap(pi);

    if (gr!=NULL) gr->reverseNumbering();
  }
}

/*
  We must still show that when |x.min==rank[x.v]| above holds, then |x.v| and
  all of its unsettled descendants form a strong component. The propagation of
  the |min| value to parents shows that |x.min==rank[x.v]| implies that none
  of the descendants of |x.v| gives access to any active vertices older than
  |x.v| (i.e., with a smaller sequence number). So it remains to show that
  from any such proper active descendant |y| one _can_ reach |x.v| itself.
  Since |y| did not head a strong component (otherwise it would have become
  settled) it gives access to an active vertex strictly older than itself,
  which must also be a descendant of |x.v|. Iterating, we finally reach |x.v|.

  We must also show that removing the strong component does not change any of
  the remaining |min| values (of course they stay the same, but their
  interpretation as minimal sequence number bla bla should also remain valid).
  But this is simply true because none of the removed vertices gives access to
  any vertices that remain active, and any removed vertex is newer than any
  remaining vertex.
*/


void OrientedGraph::reverseEdges()
{
  std::vector<EdgeList> new_edges(size(),EdgeList());
  for (Vertex x=0; x<size(); ++x)
    {
      const EdgeList& el = edgeList(x);
      for (size_t j=0; j<el.size(); ++j)
	new_edges[el[j]].push_back(x);
    }

  d_edges=new_edges; // assignment removes spare capacity in lists, swap won't
}

void OrientedGraph::reverseNumbering()
{
  Vertex last=size()-1;
  std::vector<EdgeList> new_edges(size(),EdgeList());

  for (Vertex x=0; x<size(); ++x)
    {
      const EdgeList& el = edgeList(x);
      EdgeList& new_list= new_edges[last-x]; new_list.reserve(el.size());
      for (size_t j=el.size(); j-->0; )
	new_list.push_back(last-el[j]);
    }
  d_edges.swap(new_edges); // swap preferred here: lists capacities are tight
}

/*
  The auxiliary method |addLinks| is called above to extend the induced graph
  on the quotient structure (note how in the public accessor |cells| we made
  use of our right to call private methods to call the manipulator |addLinks|
  on a _different_ OrientedGraph!). The argument |out| gives a vector of
  EdgeList pointers, whose destination vertices must be translated through the
  parition function |pi| to obtain a set of existing vertices of the induced
  graph, to which a freshly created vertex should be linked. There will almost
  certainly be internal links in the component added, which means the new
  vertex created will also appear among the images by |pi|, but it will be
  ignored (no loop in the induced graph is created).

  Precondition: all values y=(*out[i])[j] (which are vertices of some _other_
  graph) have their value |pi(y)| already defined, and |pi(y)<=c| where
  |c=size()| is the current size of the (induced) graph |*this| (it is
  therefore also the number of the new vertex that |addLinks will create).
*/

void OrientedGraph::addLinks
  (const std::vector<const EdgeList*>& out, const partition::Partition& pi)
{
  Vertex c=newVertex(); // extend graph by vertex; it now has |c+1| vertices
  bitmap::BitMap seen(c+1);
  for (size_t i=0; i<out.size(); ++i)
    for (size_t j = 0; j < (*out[i]).size(); ++j)
      seen.insert(pi((*out[i])[j]));
  seen.remove(c); // exclude any edge from class |c| to itself
  EdgeList& e = edgeList(c); // the new edge list to define
  e.reserve(seen.size());

  // by bitmap iterator's semantics we can thus add edges to set bit positions:
  std::copy(seen.begin(),seen.end(),back_inserter(e));

}

} // namespace graph

} // namespace atlas
