/*!
\file
\brief Implementation for the type Grading.
*/
/*
  This is gradings.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "gradings.h"

#include <functional>

#include "latticetypes_fwd.h"

#include "comparison.h"
#include "partition.h"
#include "rootdata.h"

/*
   A grading is always a $\Z/2Z$ grading of a (sub)system of imaginary roots.
   The grading value 0 is called "compact", while 1 is called "noncompact"
*/

namespace atlas {
namespace gradings {
namespace {


  /* The class |GradingAction| serves to provide an action function argument
     to |partition::makeOrbits|.
  */
  class GradingAction
  {
    GradingList d_cartan; // transpose Cartan matrix reduced modulo 2

    bool cartan(unsigned long i, unsigned long j) const
    { return d_cartan[j].test(i); }
  public:
    GradingAction(const latticetypes::WeightList& wts)
      : d_cartan(wts.begin(),wts.end()) {}
    // apply simple generator |s| to grading encoded by |g|
    unsigned long operator()(size_t s, unsigned long g) const;

  };

} // |namespace|
} // |namespace gradings|

/*****************************************************************************

        Chapter I -- Functions declared in gradings.h

******************************************************************************/

namespace gradings {

/* find maximal long orthogonal set of roots in |subsys| so that, with
   |non_compact| as initial set of noncompact roots, each root is non-compact,
   and passing through these roots only compact roots remain in the orthogonal
 */
rootdata::RootSet max_orth(const rootdata::RootSet& non_compact,
			   const rootdata::RootSet& subsys,
			   const rootdata::RootSystem& rs)
{
  rootdata::RootSet ncr = non_compact; // working variables
  rootdata::RootSet sub = subsys;

  ncr &= rs.posRootSet();    // restricting to positive roots is safe
  sub &= rs.posRootSet(); // and improves performance

  rootdata::RootSet result(rs.numRoots());
  while (not ncr.empty())
  {
    rootdata::RootNbr alpha = ncr.front();
    result.insert(alpha);
    for (rootdata::RootSet::iterator it = sub.begin(); it(); ++it)
      if (not rs.isOrthogonal(alpha,*it))
	sub.remove(*it);
      else if (rs.sumIsRoot(alpha,*it))
	ncr.flip(*it);

    ncr &= sub; // this removes in particular |alpha|
  }

  return rs.long_orthogonalize(result);
}



#if 0 // the following functions are never used

/*!
  Synopsis: tells whether |v| is noncompact w.r.t. the grading |g|.

  NOTE : it is essential that |v| is expressed in the root basis in which
  |g| is also given!

  Given this, the condition is just that the sum of the coordinates of |v|
  that are flagged by set bits of |g| is odd.
*/
bool isNonCompact(const rootdata::Root& v, const Grading& g)
{
  bool b = false;

  for (unsigned long i = 0; i < v.size(); ++i)
    b ^= g[i] and v[i]%2!=0;

  return b;
}

/* Here is an example of how one could use |isNonCompact|.

  This function puts in |cr| the list of all roots in |rd| that are
  noncompact for the grading |g| of the full root system.

  An important point is _not_ to call |lattice::baseChange| (as was originally
  done in this example) to convert roots to the the "basis" of simple roots,
  since that set need not have full rank (in the non-semisimple case), and
  |baseChange| cannot handle that. However |RootSystem::root_expr| is safe.
*/

rootdata::RootSet
noncompact_roots(const Grading& g, const rootdata::RootSystem& rs)
{
  rootdata::RootSet result(rs.numRoots());
  for (rootdata::RootNbr i=0; i<rs.numRoots(); ++i)
    result.set_to(i,isNonCompact(rs.root_expr(i),g));
  return result;
}


/*!
  This function returns a set of representatives of $W$-conjugacy classes
  of $Z/2Z$-gradings on the root system |rs|, which arise from a $Z$-grading.

  Such a grading is determined freely by its value at each simple root;
  therefore a grading is represented by bitset::RankFlags (=BitSet<RANK_MAX>).

  We return the first representative in each class with the least possible
  number of set bits. This code appears to be nowher used.
*/
GradingList grading_classes(const rootdata::RootSystem& rs)
{
  partition::Partition pi = partition::orbits
    (GradingAction(rs.cartanMatrix().columns()),rs.rank(),1ul<<rs.rank());

  // sort orbits
  GradingList result; result.reserve(pi.classCount());
  for (partition::PartitionIterator i(pi); i(); ++i)
    result.push_back
      (Grading(*std::min_element
	       (i->first,
		i->second,
		comparison::compare(std::ptr_fun(&bits::bitCount)))
	       ));

  return result;
}

/*!
  This function returns the equations saying that the grading should
  be noncompact (have value 1) at each element of |nc|.
*/
bitvector::BinaryEquationList
noncompact_eqns(const latticetypes::WeightList& nc)
{
  bitvector::BinaryEquationList result; result.reserve(nc.size());
  for (size_t i=0; i<nc.size(); ++i)
  {
    bitvector::BinaryEquation e(nc[i]); // reduce mod 2
    e.pushBack(true); // add a bit 1, saying that value should be noncompact
    result.push_back(e);
  }
  return result;
}

/*!
  \brief Transforms the grading |gr| of |rl| according to |so|.

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
  noncompact root \f$\alpha\f$ in |so| is that the grading of \f$\beta\f$ in
  |rl| is flipped if and only if \f$\alpha+\beta\f$ is a root. So in all the
  grading of \f$\beta\f$ changes iff this condition is met for an odd number
  of \f$\alpha\f$s.
*/
void transform_grading(gradings::Grading& gr,
		      const rootdata::RootList& rl,
		      const rootdata::RootSet& so,
		      const rootdata::RootSystem& rs)
{
  for (size_t i = 0; i < rl.size(); ++i)
    for (rootdata::RootSet::iterator jt=so.begin(); jt(); ++jt)
      if (rs.sumIsRoot(rl[i],*jt))
	gr.flip(i);
}

/*!
  In this function we assume that |o| is an orthogonal set of roots in the
  root system |rs|. Let |o_orth| be the set of roots of |rs| orthogonal to
  |o|. Then we add to |eqn| the equations which express the fact that after
  descent through |o|, the roots in |o_orth| have to be all compact (i.e., the
  grading on |rs| becomes the trivial grading on |o_orth|.)

  Recall that when descending through a short noncompact root |a|, if |b| is a
  short root in |rs| such that |a +- b| is a root (in other words, |a| and |b|
  are short roots generating a subsystem of type $B_2$), then the descent
  through |a| _changes_ the grading of |b|. So for every root \f$\alpha\f$ of
  |o_orth|, the grading must take the value given by the parity of the number
  of roots \f$\beta\f$ in |o| that with \f$\alpha\f$ so span a type $B_2$
  system.

  So the approach here is that we take a basis |oob| of the root system in
  |o_orth|, and compute the degree of the elements in that basis through this
  rule. This then gives us the equations, to the number of rank(o_orth), that
  we have to add to |eqn|.
*/
  void add_compact_eqns(bitvector::BinaryEquationList& eqn,
			const rootdata::RootSet& o,
			const rootdata::RootSet& subsys,
			const rootdata::RootSystem& rs)
{
  // make orthogonal root system

  rootdata::RootSet o_orth = rootdata::makeOrthogonal(o,subsys,rs);

  // find a basis, since having the desired grading on a basis will suffice
  rootdata::RootList oob = rs.simpleBasis(o_orth);

  // find signs

  Grading g; // compact grading on basis |oob|, "length" is |oob.size()|
  transform_grading(g,oob,o,rs);

  // express oob in subsys's root basis
  rootdata::RootList rsb = rs.simpleBasis(subsys); // (caller knew this value)

  latticetypes::WeightList woob; // "weight" form of |oob|, on basis |rsb|

  rs.toRootBasis(oob.begin(),oob.end(),back_inserter(woob),rsb);

  // write equations


  for (unsigned long i = 0; i < woob.size(); ++i)
  {
    bitvector::BinaryEquation e(woob[i]); // reduce left hand side modulo 2
    e.pushBack(g[i]); // value grading should take on root |oob[i]| (rhs)
    eqn.push_back(e);
  }
}


/*!
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
  says that for any root \f$\alpha\f$ of |subsys| orthogonal to |o|, the value
  \f$gr(\alpha)\f$ should be equal modulo 2 to the number of roots \f$\beta\f$
  of |o| such that \f$\alpha+\beta\f$ is a root (of |rs|). Thus the values of
  |gr| are fixed on |o| and on a basis of the root subsystem of |subsys|
  orthogonal to |o|.

  NOTE: the solution is not unique in general. However the grading that is
  found should be unique up to conjugacy. If the subspace orthogonal to |o| is
  spanned by the roots it contains, then we have as many equations as
  unknowns; but even then the system can be degenerate. This happens already
  in case B2, for the split Cartan.
*/
rootdata::RootSet grading_for_orthset
  (const rootdata::RootSet& o,
   const rootdata::RootSet& subsys,
   const rootdata::RootSystem& rs)
{
  // find basis |rb| for |subsys|
  rootdata::RootList rb = rs.simpleBasis(subsys);

  // express roots of |o| in basis |rb| (of |subsys|) just found
  latticetypes::WeightList wo; // will hold "weights" of length |rb.size()|
  rs.toRootBasis(o.begin(),o.end(),back_inserter(wo),rb);

  // now do the same for \emph{all} roots in |subsys|
  latticetypes::WeightList wrs; // will hold "weights" of length |rb.size()|
  rs.toRootBasis(subsys.begin(),subsys.end(),back_inserter(wrs),rb);

  // write equations for the roots in o
  bitvector::BinaryEquationList eqn = noncompact_eqns(wo);
  // now |eqn| describes gradings with value 1 on elements of |o|.

  // add equations for the maximality
  add_compact_eqns(eqn,o,subsys,rs);
  // now |eqn| also requires specific values on roots orthogonal to |o|

  // solve system

  bitvector::BinaryEquation gc(rb.size());
  firstSolution(gc,eqn);

  Grading g = gc.data();

  // return the grading

  rootdata::RootSet result;

  for (unsigned long i = 0; i < wrs.size(); ++i)
    if (isNonCompact(wrs[i],g))
      result.insert(subsys.n_th(i));

  return result;
}

/*!
  Returns in |gt| a maximal set of strongly orthogonal noncompact roots for
  the grading |g|. This function is currently unused [MvL 14 oct 2007]

  The algorithm is simple: pick any noncompact root; find the orthogonal
  system; write down the grading on the orthogonal; and repeat. When done,
  transform the system (if it isn't already) into a strongly orthogonal one.
*/
rootdata::RootSet grading_type(const Grading& g,
			       const rootdata::RootSystem& rs)
{ // start with set of all roots
  rootdata::RootSet all_roots(rs.numRoots()); all_roots.fill();

  return max_orth(noncompact_roots(g,rs),all_roots,rs);
}
#endif

} // |namespace gradings|

/*****************************************************************************

        Chapter II -- Auxiliary classes

******************************************************************************/

/*!
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
  \f$s(\alpha_i)=\alpha_i-\left<\alpha_i,\alpha_s^\vee\right>\alpha_s\f$. So if
  |g(s)==0| there is nothing to change for any |i|, while if |g(s)==1| we must
  add (modulo 2) column |s| of the Cartan matrix to |g|. For this reason we
  store upon construction the transpos Cartan matrix modulo 2.

  NOTE : the conversion to/from |unsigned long| depends on the fact that
  |RANK_MAX| does not exceed the number of bits in an |unsigned long|. It
  would have been more natural to to give argument and result as |Grading|.
  The restriction however is already built into |partition::makeOrbit|, which
  can only make orbits on sets whose points are represeted by |unsigned long|.
*/
namespace gradings {
namespace {

unsigned long GradingAction::operator() (size_t s, unsigned long g) const
{
  Grading gr(g); // convert from |unsigned long|
  if (gr.test(s))
    gr ^= d_cartan[s];

  return gr.to_ulong(); // convert (back) to |unsigned long|
}


} // |namespace|
} // |namespace gradings|

/*****************************************************************************

        Chapter III -- The GradingCompare class

  This is a function object that compares gradings in a size-first fashion.

******************************************************************************/

namespace gradings {


/*!
  Synopsis: compares |lhs| and |rhs|, returns whether |lhs<rhs| as gradings

  By definition this holds iff either |lhs| has fewer set bits, or they have
  the number, and |lhs<rhs| holds in the usual (bitwise lexicographic) sense.
*/
bool GradingCompare::operator() (const Grading& lhs, const Grading& rhs)
{
  size_t lhc=lhs.count(), rhc=rhs.count();
  return lhc!=rhc ? lhc<rhc : lhs<rhs;
}

} // |namespace gradings|
} // |namespace atlas|
