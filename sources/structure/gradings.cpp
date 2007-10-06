/*!
\file
\brief Implementation for the type Grading.
*/
/*
  This is gradings.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "gradings.h"

#include <functional>

#include "comparison.h"
#include "latticetypes.h"
#include "lattice.h"
#include "partition.h"
#include "rootdata.h"

/*
   A grading is always a $\Z/2Z$ grading of a (sub)system of imaginary roots.
   The grading value 0 is called "compact", while 1 is called "noncompact"
*/

namespace atlas {

namespace {

  using namespace constants;
  using namespace gradings;
  using namespace latticetypes;
  using namespace rootdata;

  // local function declarations

  void noncompactEquations(BinaryEquationList&, const WeightList&);
  void compactEquations(BinaryEquationList&, const RootList&, const RootList&,
			const RootDatum&);
  // local class declarations: |BitCount| and |GradingAction|

  /* The class |BitCount| makes an adaptable unary function object (for use by
     comparison::Compare) from |bits::bitCount|. The |argument_type| actually
     represents a |Grading|, but at the point of use, the gradings in question
     are actually encoded as |unsigned long| values inside a |Partition|.
  */
  struct BitCount
    : public std::unary_function<unsigned long, unsigned long>
  {
    result_type operator() (argument_type a) const
    { return bits::bitCount(a); }
  };


  /* The class |GradingAction| serves to provide an action functino argument
     to |partition::makeOrbits|.
  */
  class GradingAction
    : public std::binary_function<unsigned long, unsigned long, unsigned long>
  // derivation is not really necessary here
  {
    GradingList d_cartan; // transpose Cartan matrix reduced modulo 2

    bool cartan(unsigned long i, unsigned long j) const {
      return d_cartan[j].test(i);
    }
  public:
    GradingAction(const RootDatum& rd);
    // apply simple generator |s| to grading encoded by |g|
    unsigned long operator()(unsigned long s, unsigned long g) const;

  };

}

/*****************************************************************************

        Chapter I -- Functions declared in gradings.h

******************************************************************************/

namespace gradings {

bool isNonCompact(const rootdata::Root& v, const Grading& g)

/*!
  Synopsis: tells whether |v| is noncompact w.r.t. the grading |g|.

  NOTE : it is essential that |v| is expressed in the root basis in which
  |g| is also given! Use |baseChange| (declared in lattice.h) if necessary.

  Given this, the condition is just that the sum of the coordinates of |v|
  that are flagged by set bits of |gr| is odd.
*/

{
  bool b = false;

  for (unsigned long j = 0; j < v.size(); ++j)
    if (v[j]%2!=0)
      b ^= g[j];

  return b;
}

  /* Here is an example of how to use |isNonCompact| */

void noncompactRoots(rootdata::RootList& ncr, const Grading& g,
		     const rootdata::RootDatum& rd)

{
  using namespace lattice;

  // need to have roots in root basis

  WeightList rl;
  baseChange(rd.beginRoot(),rd.endRoot(),back_inserter(rl),
	     rd.beginSimpleRoot(),rd.endSimpleRoot());

  // make list of compact roots

  for (unsigned long j = 0; j < rl.size(); ++j)
    if (isNonCompact(rl[j],g))
      ncr.push_back(j);
}

void compactRoots(rootdata::RootList& cr, const Grading& g,
		  const rootdata::RootDatum& rd)

/*!
  This function puts in |cr| the list of roots in |rd| that are compact for
  the grading |g|.
*/

{
  using namespace lattice;

  // need to have roots in root basis

  WeightList rl;
  baseChange(rd.beginRoot(),rd.endRoot(),back_inserter(rl),
	     rd.beginSimpleRoot(),rd.endSimpleRoot());

  // make list of compact roots

  for (unsigned long j = 0; j < rl.size(); ++j)
    if (not isNonCompact(rl[j],g))
      cr.push_back(j);
}

void makeGradings(GradingList& gl, const rootdata::RootDatum& rd)

/*!
  This function puts in |gl| a set of representatives of $W$-conjugacy classes
  of $Z/2Z$-gradings on the root system of |rd|, which arise from a
  $Z$-grading.

  Such a grading is determined freely by its value at each simple root;
  therefore a grading is represented by bitset::RankFlags (=BitSet<RANK_MAX>).

  We return the first representative in each class with the least possible
  number of set bits
*/

{
  partition::Partition pi;
  partition::makeOrbits
    (pi,GradingAction(rd),rd.semisimpleRank(),(1UL) << rd.semisimpleRank());

  // sort orbits

  for (partition::PartitionIterator i(pi); i(); ++i) {
    gl.push_back(Grading(*std::min_element
			 (i->first,i->second,comparison::compare(BitCount()))
			 ));

  }
}

void findGrading(RootSet& ncr, const RootList& o, const RootList& rs,
		 const RootDatum& rd)

/*!
  Given a RootList |rs|, which is assumed to hold a root subsystem in |rd|,
  and an orthogonal subsystem |o| inside |rs|, this function returns a grading
  |g| of the root system |rs| (expressed in terms of its canonical basis) for
  which |o| is a maximal orthogonal set of noncompact roots (or more
  precisely, for which (1) all roots in |o| are noncompact for |gr|, and (2)
  the descent of |gr| through the roots of |o| to the subsystem of |rs|
  orthogonal to |o| leads to a compact grading of that subsystem.) The answer
  is written directly in the RootSet |ncr|, which flags the noncompact roots.

  The condition given amounts to solving a system of linear equations modulo
  2. In addition to the values 1 of |gr| on |o| given by (1), condition (2)
  says that for any root \f$\alpha\f$ of |rs| orthogonal to |o|, the value
  \f$gr(\alpha)\f$ should be equal modulo 2 to the number of roots \f$\beta\f$ of |o|
  such that \f$\alpha+\beta\f$ is a root (of |rs|). Thus the values of |gr| are
  fixed on |o| and on a basis of the root subsystem of |rs| orthogonal to |o|.

  NOTE: the solution is not unique in general. However the grading that is
  found should be unique up to conjugacy. If the subspace orthogonal to |o| is
  spanned by the roots it contains, then we have as many equations as
  unknowns; but even then the system can be degenerate. This happens already
  in case B2, for the split Cartan.
*/

{
  using namespace lattice;
  using namespace rootdata;

  // find basis |rb| for |rs|

  RootList rb; rootdata::rootBasis(rb,rs,rd);

  // express roots of |o| in basis |rb| (of subsystem |rs| of |rd|) just found
  WeightList wo; // will hold "weights" of length |rb.size()| (actually roots)
  toRootBasis(o.begin(),o.end(),back_inserter(wo),rb,rd);

  // do the same for all roots in |rs|
  WeightList wrs; // will hold "weights" of length |rb.size()| (actually roots)
  toRootBasis(rs.begin(),rs.end(),back_inserter(wrs),rb,rd);


  /* As "equation" is a $Z/2Z$-linear one for gradings of the root lattice
     modulo 2 of the root subsystem |rs| with basis |rb|. Each equation
     describes an element of that modular lattice, and a value in $Z/2Z$ that
     the grading should take on that element. Equations are represented by the
     type |latticetypes::BinaryEquation|, which allows for |RANK_MAX+1| bits,
     of which |rb.size()+1| bits are used. The first |rb.size()| bits give the
     (mod 2) coefficients of root vector the left hand side of the equation,
     the final bit gives the right hand side.
   */

  // write equations for the roots in o

  BinaryEquationList eqn; noncompactEquations(eqn,wo);
  // now |eqn| describes gradings with value 1 on elements of |o|.

  // add equations for the maximality
  compactEquations(eqn,o,rs,rd);
  // now |eqn| also requires specific values on roots orthogonal to |o|

  // solve system

  BinaryEquation gc(rb.size());
  firstSolution(gc,eqn);

  Grading g = gc.data();

  // write the grading in ncr

  for (unsigned long j = 0; j < wrs.size(); ++j)
    if (isNonCompact(wrs[j],g))
      ncr.insert(rs[j]);

  // |ncr| now holds a solution
}

void gradingType(rootdata::RootList& gt, const Grading& g,
		 const rootdata::RootDatum& rd)

/*!
  Returns in gt a maximal set of strongly orthogonal noncompact roots for the
  grading.

  The algorithm is simple: pick any noncompact root; find the orthogonal
  system; write down the grading on the orthogonal; when done, transform the
  system into a strongly orthogonal one if needed.
*/

{
  using namespace lattice;
  using namespace bitmap;

  typedef BitMap::iterator BI;

  // need to have roots in root basis

  WeightList wl;
  baseChange(rd.beginRoot(),rd.endRoot(),back_inserter(wl),
	     rd.beginSimpleRoot(),rd.endSimpleRoot());

  // find grading of all roots

  BitMap nc(wl.size());

  for (unsigned long j = 0; j < wl.size(); ++j)
    if (isNonCompact(wl[j],g))
      nc.insert(j);

  // find maximal orthogonal set

  RootList rl;

  for (unsigned long j = 0; j < wl.size(); ++j)
    rl.push_back(j);

  while (!nc.empty()) {

    // insert new wlwment in gt

    RootNbr n = nc.front();
    gt.push_back(n);

    // put in rl the orthogonal to rl.front()

    RootList orth;
    for (unsigned long j = 0; j < rl.size(); ++j) {
      if (rd.isOrthogonal(n,rl[j]))
	orth.push_back(rl[j]);
      else
	nc.remove(rl[j]);
    }
    rl.swap(orth);

    // find grading

    for (unsigned long j = 0; j < rl.size(); ++j)
      if (sumIsRoot(rd.root(n),rd.root(rl[j]),rd)) { // change parity
	if (nc.isMember(rl[j]))
	  nc.remove(rl[j]);
	else
	  nc.insert(rl[j]);
      }
  }

  // correct gt to be made of strongly orthogonal roots

  for (unsigned long i = 1; i < gt.size(); ++i)
    for (unsigned long j = 0; j < i; ++j)
      if (sumIsRoot(rd.root(gt[i]),rd.root(gt[j]),rd)) {
	// replace gt[i] and gt[j] by their sum and difference
	// note that they generate a little B2 system
	Weight a(rd.root(gt[i]));
	a += rd.root(gt[j]);
	RootNbr n = rd.rootNbr(a);
	a = rd.root(gt[i]);
	a -= rd.root(gt[j]);
	gt[i] = n;
	gt[j] = rd.rootNbr(a);
      }

  return;
}

} // namespace gradings

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
  \f$s(\alpha_i)=\alpha_i-\<alpha_i,\alpha_s^\vee>\alpha_s\f$. So if |g(s)==0|
  there is nothing to change for any |i|, while if |g(s)==1| we must add
  (modulo 2) column |s| of the Cartan matrix to |g|. For this reason we store
  upon construction the transpos Cartan matrix modulo 2.

  NOTE : the conversion to/from |unsigned long| depends on the fact that
  |RANK_MAX| does not exceed the number of bits in an |unsigned long|. It
  would have been more natural to to give argument and result as |Grading|.
  The restriction however is already built into |partition::makeOrbit|, which
  can only make orbits on sets whose points are represeted by |unsigned long|.
*/

namespace {

GradingAction::GradingAction(const RootDatum& rd)
  : d_cartan(rd.semisimpleRank())
{
  for (unsigned long i = 0; i < rd.semisimpleRank(); ++i)
    for (unsigned long j = 0; j < rd.semisimpleRank(); ++j)
      d_cartan[j].set(i,rd.cartan(i,j)%2!=0); // reduce modulo 2, transposing
}

unsigned long GradingAction::operator() (unsigned long s, unsigned long g)
  const
{
  Grading gr(g); // convert from |unsigned long|
  if (gr.test(s))
    gr ^= d_cartan[s];

  return gr.to_ulong(); // convert (back) to |unsigned long|
}


}

/*****************************************************************************

        Chapter III -- The GradingCompare class

  This is a function object that compares gradings in a size-first fashion.

******************************************************************************/

namespace gradings {

bool GradingCompare::operator() (const Grading& lhs, const Grading& rhs)

/*!
  Synopsis: compares |lhs| and |rhs|, returns whether |lhs<rhs| as gradings

  By definition this holds iff either |lhs| has fewer set bits, or they have
  the number, and |lhs<rhs| holds in the usual (bitwise lexicographic) sense.
*/

{
  size_t lhc=lhs.count(), rhc=rhs.count();
  return lhc!=rhc ? lhc<rhc : lhs<rhs;
}

}

/*****************************************************************************

        Chapter III -- Auxiliary functions

******************************************************************************/

namespace {

void noncompactEquations(BinaryEquationList& eqn, const WeightList& nc)

/*!
  This function adds to eqn the equations saying that the grading should
  be noncompact (have value 1) at each element of |nc|.
*/

{
  for (unsigned long j = 0; j < nc.size(); ++j) {
    BinaryEquation e(nc[j]); // reduce mod 2
    e.pushBack(true); // add a bit 1, saying that value should be noncompact
    eqn.push_back(e);
  }
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
  of roots \f$\beta\f$ in |o| that with \f$\alpha\f$ so span a type $B_2$ system.

  So the approach here is that we take a basis |oob| of the root system in
  |o_orth|, and compute the degree of the elements in that basis through this
  rule. This then gives us the equations, to the number of rank(o_orth), that
  we have to add to |eqn|.
*/
void compactEquations(BinaryEquationList& eqn, const RootList& o,
		      const RootList& rs, const RootDatum& rd)
{
  using namespace lattice;

  // make orthogonal root system

  RootList o_orth;
  makeOrthogonal(o_orth,o,rs,rd);

  // find basis

  RootList oob; // basis of |o_orth|
  rootdata::rootBasis(oob,o_orth,rd);

  // find signs

  Grading g; // grading on basis |oob|, therefore of length |oob.size()|

  for (unsigned long i = 0; i < oob.size(); ++i)
    for (unsigned long j = 0; j < o.size(); ++j)
      if (rd.sumIsRoot(oob[i],o[j]))
	g.flip(i); // use |flip|, not |set|, as same |i| can have more flips

  // express oob in rs's root basis
  RootList rsb; // simple basis of |rs| (could have been transmitted as arg)
  rootBasis(rsb,rs,rd);

  WeightList woob; // |oob| in form of "weights", expressed on basis |rsb|

  toRootBasis(oob.begin(),oob.end(),back_inserter(woob),rsb,rd);

  // write equations


  for (unsigned long j = 0; j < woob.size(); ++j) {
    BinaryEquation e(woob[j]); // reduce left hand side modulo 2
    e.pushBack(g[j]); // value grading should take on root |oob[j]| (rhs)
    eqn.push_back(e);
  }
}

}

}
