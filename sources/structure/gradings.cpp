/*
  This is gradings.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "gradings.h"

#include <functional>

#include "comparison.h"
#include "partition.h"

#include "bitvector.h"
#include "rootdata.h"

/*
  A grading is always a $\Z/2Z$ grading of a (sub)system of imaginary roots.
  The grading value 0 is called "compact", while 1 is called "noncompact".

  Gradings are stored as BitSet<RANK_MAX> values, and in depends on the
  context what their actual size is, and which set of imaginary roots is
  graded (either the imaginary subsystem at a given involution, or the span of
  the delta-fixed simple roots.

  Gradings used to be extensively employed, but they have been supplanted by
  the finer notion of strong involutions; these define gradings but in a
  non-injective way, so while gradings are their most visible aspect, they do
  not in general supply sufficient information for computations for instance
  for modeling the set K\G/B. And since |Grading| is a typedef rather than a
  class, there is little left in this module to be defined.
*/

namespace atlas {
namespace gradings {

/*****************************************************************************

        Chapter I -- Functions declared in gradings.h

******************************************************************************/


/* find maximal long orthogonal set of roots in |subsys| so that, with
   |non_compact| as initial set of noncompact roots, each root is non-compact,
   and passing through these roots only compact roots remain in the orthogonal
 */
RootNbrSet max_orth(const RootNbrSet& non_compact,
		    const RootNbrSet& subsys,
		    const RootSystem& rs)
{
  RootNbrSet ncr = non_compact; // working variables
  RootNbrSet sub = subsys;

  ncr &= rs.posRootSet(); // restricting to positive roots is safe
  sub &= rs.posRootSet(); // and improves performance

  RootNbrSet result(rs.numRoots());
  while (not ncr.empty())
  {
    RootNbr alpha = ncr.front();
    result.insert(alpha);
    for (RootNbrSet::iterator it = sub.begin(); it(); ++it)
      if (not rs.is_orthogonal(alpha,*it))
	sub.remove(*it); // remove non-orthogonal roots as candidates
      else if (rs.sum_is_root(alpha,*it))
	ncr.flip(*it); // and flip compactness status of these short roots

    ncr &= sub; // ensure removed roots, and |alpha|, are no longer considered
  }

  return rs.long_orthogonalize(result); // ensure \emph{strong} orthogonality
}




/*****************************************************************************

        Chapter II -- Auxiliary classes

******************************************************************************/

/*****************************************************************************

        Chapter III -- The GradingCompare class

  This is a function object that compares gradings in a size-first fashion.

******************************************************************************/

/*
  Compare |lhs| and |rhs|, returns whether |lhs<rhs| as gradings

  By definition this holds iff either |lhs| has fewer set bits, or they have
  the number, and |lhs<rhs| holds in the usual (bitwise lexicographic) sense.
*/
bool GradingCompare::operator() (const Grading& lhs, const Grading& rhs)
{
  size_t lhc=lhs.count(), rhc=rhs.count();
  return lhc!=rhc ? lhc<rhc : lhs<rhs;
}


/****************************************************************************/


#if 0 // the remaining functions are never used, but are retained for reference

/*
  Whether |v| is noncompact w.r.t. the grading |g|.

  NOTE : it is essential that |v| is expressed in the root basis in which
  |g| is also given!

  Given this, the condition is just that the sum of the coordinates of |v|
  that are flagged by set bits of |g| is odd.
*/
bool isNonCompact(const rootdata::Root& v, const Grading& g)
{
  unsigned s=0;

  for (auto it=g.begin(); it(); ++it)
    s += v[*it];

  return s%2!=0;
}

/* Here is an example of how one could use |isNonCompact|.

  This function puts in |cr| the list of all roots in |rd| that are
  noncompact for the grading |g| of the full root system.

  An important point is _not_ to call |lattice::baseChange| (as was originally
  done in this example) to convert roots to the the "basis" of simple roots,
  since that set need not have full rank (in the non-semisimple case), and
  |baseChange| cannot handle that. However |RootSystem::root_expr| is safe.
*/

RootNbrSet noncompact_roots(const Grading& g, const RootSystem& rs)
{
  RootNbrSet result(rs.numRoots());
  for (RootNbr i=0; i<rs.numRoots(); ++i)
    result.set_to(i,isNonCompact(rs.root_expr(i),g));
  return result;
}



/*
  A GradingAction |a| is a function object, to be used for the determination
  of $W$-orbits on gradings. The result of calling |a(s,g)|, where |s| is a
  generator of $W$ and |g| a grading of the root datum used to construct |a|,
  is the the grading |sg| obtained from the one encoded by |g| by the action
  of the simple reflection number |s|, converted (back) to unsigned long. This
  grading is defined by $sg(i)=g(si)$ for all simple roots |i|, where |si| is
  the (in general non-simple) root obtained from applying |s| to |i|. To
  compute |g(si)| one needs to express |si| in the simple root basis.

  Note that since we are working modulo 2, $g(si)$ is either $g(i)$ or
  $g(i)+g(s)$, where $s$ is now interpreted as the simple roots corresponding
  to the reflection; the second case applies if and only if the Cartan matrix
  has an odd entry at (i,s), since
  $s(\alpha_i)=\alpha_i-\left<\alpha_i,\alpha_s^\vee\right>\alpha_s$. So if
  |g(s)==0| there is nothing to change for any |i|, while if |g(s)==1| we must
  add (modulo 2) column |s| of the Cartan matrix to |g|. For this reason we
  store upon construction the transpos Cartan matrix modulo 2.

  NOTE : the conversion to/from |unsigned long| depends on the fact that
  |RANK_MAX| does not exceed the number of bits in an |unsigned long|. It
  would have been more natural to to give argument and result as |Grading|.
  The restriction however is already built into |partition::makeOrbit|, which
  can only make orbits on sets whose points are represented by |unsigned long|.
*/
namespace {

class GradingAction
{
  GradingList d_cartan; // transpose Cartan matrix reduced modulo 2

  bool cartan(unsigned long i, unsigned long j) const
  { return d_cartan[j].test(i); }
public:
  GradingAction(const WeightList& wts)
    : d_cartan(wts.begin(),wts.end()) {}

  // apply simple generator |s| to grading encoded by |g|
  unsigned long operator()(size_t s, unsigned long g) const
  {
    Grading gr(g); // convert from |unsigned long|
    if (gr.test(s))
      gr ^= d_cartan[s];
    return gr.to_ulong(); // convert (back) to |unsigned long|
  }

}; // |class GradingAction|

} // |namespace|

/*
  This function returns a set of representatives of $W$-conjugacy classes
  of $Z/2Z$-gradings on the root system |rs|, which arise from a $Z$-grading.

  Such a grading is determined freely by its value at each simple root;
  therefore a grading is represented by RankFlags (=BitSet<RANK_MAX>).

  We return the first representative in each class with the least possible
  number of set bits. This code appears to be nowhere used.
*/
GradingList grading_classes(const RootSystem& rs)
{
  Partition pi = partition::orbits
    (GradingAction(rs.Cartan_matrix().columns()),rs.rank(),1ul<<rs.rank());

  // sort orbits
  GradingList result; result.reserve(pi.classCount());
  for (Partition::iterator i(pi); i(); ++i)
    result.push_back
      (Grading(*std::min_element
	       (i->first,
		i->second,
		comparison::compare(std::ptr_fun(&bits::bitCount)))
	       ));

  return result;
}

/*
  This function returns the equations saying that the grading should
  be noncompact (have value 1) at each element of |nc|.
*/
BinaryEquationList noncompact_eqns(const WeightList& nc)
{
  BinaryEquationList result; result.reserve(nc.size());
  for (size_t i=0; i<nc.size(); ++i)
  {
    BinaryEquation e(nc[i]); // reduce mod 2
    e.pushBack(true); // add a bit 1, saying that value should be noncompact
    result.push_back(e);
  }
  return result;
}

/*
  Transform the grading |gr| of |rl| according to |so|.

  Precondition: |gr| is a grading of the roots in |rl|; |so| is a set of
  long-orthogonal roots, orthogonal to the roots in |rl|.

  Assume that all roots in |rl| and |so| are imaginary for some involution,
  and that a grading is given that matches |gr| on |rl|, and for which all
  roots of |so| are noncompact (grading 1). Then by Cayley-transforming though
  the roots of |so| (in any order, since they are long-orthogonal) one
  gets an involution for which the roots in |rl| are still imaginary (those in
  |so| have become real), and a grading of those roots that possibly differs
  from |gr|; this function transforms the grading |gr| into that new grading.

  Formula: the rule for Cayley-transforming a grading through an imaginary
  noncompact root $\alpha$ in |so| is that the grading of $\beta$ in |rl| is
  flipped if and only if $\alpha+\beta$ is a root. So in all the grading of
  $\beta$ changes iff this condition is met for an odd number of $\alpha$s.
*/
void transform_grading(Grading& gr,
	const RootNbrList& rl, const RootNbrSet& so, const RootSystem& rs)
{
  for (size_t i = 0; i < rl.size(); ++i)
    for (RootNbrSet::iterator jt=so.begin(); jt(); ++jt)
      if (rs.sumIsRoot(rl[i],*jt))
	gr.flip(i);
}

/*
  In this function we assume that |o| is an orthogonal set of roots in the
  root system |rs|. Let |o_orth| be the set of roots of |rs| orthogonal to
  |o|. Then we add to |eqn| the equations which express the fact that after
  descent through |o|, the roots in |o_orth| have to be all compact (i.e., the
  grading on |rs| becomes the trivial grading on |o_orth|.)

  Recall that when descending through a short noncompact root |a|, if |b| is a
  short root in |rs| such that |a +- b| is a root (in other words, |a| and |b|
  are short roots generating a subsystem of type $B_2$), then the descent
  through |a| _changes_ the grading of |b|. So for every root $\alpha$ of
  |o_orth|, the grading must take the value given by the parity of the number
  of roots $\beta$ in |o| that with $\alpha$ so span a type $B_2$
  system.

  So the approach here is that we take a basis |oob| of the root system in
  |o_orth|, and compute the degree of the elements in that basis through this
  rule. This then gives us the equations, to the number of rank(o_orth), that
  we have to add to |eqn|.
*/
  void add_compact_eqns(BinaryEquationList& eqn,
	const RootNbrSet& o, const RootNbrSet& subsys, const RootSystem& rs)
{
  // make orthogonal root system

  RootNbrSet o_orth = rootdata::makeOrthogonal(o,subsys,rs);

  // find a basis, since having the desired grading on a basis will suffice
  RootNbrList oob = rs.simpleBasis(o_orth);

  // find signs

  Grading g; // compact grading on basis |oob|, "length" is |oob.size()|
  transform_grading(g,oob,o,rs);

  // express oob in subsys's root basis
  RootNbrList rsb = rs.simpleBasis(subsys); // (caller knew this value)

  int_VectorList woob; // "weight" form of |oob|, on basis |rsb|

  rs.toRootBasis(oob.begin(),oob.end(),back_inserter(woob),rsb);

  // write equations


  for (unsigned long i = 0; i < woob.size(); ++i)
  {
    BinaryEquation e(woob[i]); // reduce left hand side modulo 2
    e.pushBack(g[i]); // value grading should take on root |oob[i]| (rhs)
    eqn.push_back(e);
  }
}


/*
  Given a RootSet |subsys|, which is assumed to hold a root subsystem in |rs|,
  and an orthogonal subsystem |o| inside |subsys|, this function returns a
  grading |g| of the root system |subsys| (expressed in terms of its canonical
  basis) for which |o| is a maximal orthogonal set of noncompact roots (or
  more precisely, for which (1) all roots in |o| are noncompact for |gr|, and
  (2) the descent of |gr| through the roots of |o| to the subsystem of
  |subsys| orthogonal to |o| leads to a compact grading of that subsystem.)
  The answer is given as the |RootSet| which flags the noncompact roots.

  The condition given amounts to solving a system of linear equations modulo
  2. In addition to the values 1 of |gr| on |o| given by (1), condition (2)
  says that for any root $\alpha$ of |subsys| orthogonal to |o|, the value
  $gr(\alpha)$ should be equal modulo 2 to the number of roots $\beta$
  of |o| such that $\alpha+\beta$ is a root (of |rs|). Thus the values of
  |gr| are fixed on |o| and on a basis of the root subsystem of |subsys|
  orthogonal to |o|.

  NOTE: the solution is not unique in general. However the grading that is
  found should be unique up to conjugacy. If the subspace orthogonal to |o| is
  spanned by the roots it contains, then we have as many equations as
  unknowns; but even then the system can be degenerate. This happens already
  in case B2, for the split Cartan.
*/
RootNbrSet grading_for_orthset
  (const RootNbrSet& o, const RootNbrSet& subsys, const RootSystem& rs)
{
  // find basis |rb| for |subsys|
  RootNbrList rb = rs.simpleBasis(subsys);

  // express roots of |o| in basis |rb| (of |subsys|) just found
  int_VectorList wo; // will hold "weights" of length |rb.size()|
  rs.toRootBasis(o.begin(),o.end(),back_inserter(wo),rb);

  // now do the same for \emph{all} roots in |subsys|
  int_VectorList wrs; // will hold "weights" of length |rb.size()|
  rs.toRootBasis(subsys.begin(),subsys.end(),back_inserter(wrs),rb);

  // write equations for the roots in o
  BinaryEquationList eqn = noncompact_eqns(wo);
  // now |eqn| describes gradings with value 1 on elements of |o|.

  // add equations for the maximality
  add_compact_eqns(eqn,o,subsys,rs);
  // now |eqn| also requires specific values on roots orthogonal to |o|

  // solve system

  BinaryEquation gc(rb.size());
  firstSolution(gc,eqn);

  Grading g = gc.data();

  // return the grading

  RootNbrSet result;

  for (unsigned long i = 0; i < wrs.size(); ++i)
    if (isNonCompact(wrs[i],g))
      result.insert(subsys.n_th(i));

  return result;
}

/*
  Returns in |gt| a maximal set of strongly orthogonal noncompact roots for
  the grading |g|.

  The algorithm is simple: pick any noncompact root; find the orthogonal
  system; write down the grading on the orthogonal; and repeat. When done,
  transform the system (if it isn't already) into a strongly orthogonal one.
*/
RootNbrSet grading_type(const Grading& g, const RootSystem& rs)
{ // start with set of all roots
  RootNbrSet all_roots(rs.numRoots()); all_roots.fill();

  return max_orth(noncompact_roots(g,rs),all_roots,rs);
}
#endif

} // |namespace gradings|
} // |namespace atlas|
