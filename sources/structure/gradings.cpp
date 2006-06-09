/*!
\file
\brief Implementation for the type Grading.
*/
/*
  This is gradings.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "gradings.h"

#include <functional>

#include "comparison.h"
#include "latticetypes.h"
#include "lattice.h"
#include "partition.h"
#include "rootdata.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace {

  using namespace constants;
  using namespace gradings;
  using namespace latticetypes;
  using namespace rootdata;

  // local function declarations

  void noncompactEquations(LongComponentList&, const WeightList&);
  void compactEquations(LongComponentList&, const RootList&, const RootList&, 
			const RootDatum&);
  void toGrading(Grading&, const LongComponent&);

  // local class declarations

  class BitCount 
    :public std::unary_function<unsigned long, unsigned long> {

  public:
  // accessors
    result_type operator() (argument_type) const;
  };

  class GradingAction {
  private:
    GradingList d_cartan;
    size_t d_rank;
  public:
  // constructors and destructors
    GradingAction() {}
    GradingAction(const RootDatum& rd);
    ~GradingAction() {}
  // accessors
    unsigned long operator()(unsigned long, unsigned long) const;
    bool cartan(unsigned long i, unsigned long j) const {
      return d_cartan[i].test(j);
    }
  };

}

/*****************************************************************************

        Chapter I -- Functions declared in gradings.h

  ... explain here when it is stable ...

******************************************************************************/

namespace gradings {

void findGrading(RootSet& ncr, const RootList& o, const RootList& rs,
		 const RootDatum& rd)

/*!
  Given a RootList rs, which is assumed to hold a root system, and an
  orthogonal subsystem o in rs, this function returns a grading of the
  root system rs (expressed in terms of its canonical basis) for which
  o is a maximal orthogonal set of noncompact roots (or more precisely,
  for which all roots in o are imaginary, and the descent leads to a
  compact root system in the orthogonal of o.) The answer is written
  directly in the RootSet ncr, which flags the noncompact roots.

  The first condition amounts to solving a system of linear equations modulo 2.
  The second condition in fact just adds some new conditions : indeed, it
  amounts to saying that after applying the descent rules to the root
  system orthogonal to o, we should get a compact root system; therefore
  we can determine the value of the grading on a basis of the orthogonal
  root system.

  NOTE : this grading should be unique up to conjugacy! However the
  solution is not unique in general. If the orthogonal to o is spanned by the
  roots it contains, we have as many equations as unknowns; but even then
  the system can be degenerate. This happens already in case B2, for the
  split Cartan.
*/

{
  using namespace lattice;
  using namespace rootdata;

  // find basis for rs

  RootList rb;
  rootBasis(rb,rs,rd);

  // write roots in rootbasis

  WeightList wo;
  toRootBasis(o.begin(),o.end(),back_inserter(wo),rb,rd);

  WeightList wrs;
  toRootBasis(rs.begin(),rs.end(),back_inserter(wrs),rb,rd);

  // write equations for the roots in o

  LongComponentList eqn;

  noncompactEquations(eqn,wo);

  // add equations for the maximality

  compactEquations(eqn,o,rs,rd);

  // solve system

  LongComponent gc(rb.size());
  firstSolution(gc,eqn);

  Grading g = gc.data();

  // write the grading in ncr

  for (unsigned long j = 0; j < wrs.size(); ++j)
    if (isNonCompact(wrs[j],g))
      ncr.insert(rs[j]);

  return;

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

void compactRoots(rootdata::RootList& cr, const Grading& g, 
		  const rootdata::RootDatum& rd)

/*!
  This function puts in cr the list of roots in rd that are compact for
  the grading g.
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

  return;
}

bool isNonCompact(const rootdata::Root& v, const Grading& g)

/*!
  Synopsis: tells whether v is noncompact w.r.t. the grading g.

  NOTE : it is essential that v is expressed in the root basis in which
  g is also given! Use baseChange (declared in lattice.h) if necessary.
*/

{
  bool b = false;

  for (unsigned long j = 0; j < v.size(); ++j)
    if (v[j]%2)
      b ^= g[j];

  return b;
}

void makeGradings(GradingList& gl, const rootdata::RootDatum& rd)

/*!
  This function puts in gl a set of representatives of W-conjugacy classes
  of Z_2-gradings on the root system of rd, which arise from a Z-grading.
  Such a grading is defined by the sign of each simple root; therefore
  the datum for a grading fits inside a BitSet<RANK_MAX>.

  We return a representative in each class with the least possible number
  of set bits.
*/

{
  using namespace comparison;
  using namespace partition;

  GradingAction a(rd);
  Partition pi;
  makeOrbits(pi,a,rd.semisimpleRank(),(1UL) << rd.semisimpleRank());

  // sort orbits

  for (PartitionIterator i(pi); i(); ++i) {
    std::vector<unsigned long> o(i->first,i->second);
    std::stable_sort(o.begin(),o.end(),compare(BitCount()));
    gl.push_back(Grading(o.front()));
  }

  return;
}

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

  return;
}

}

/*****************************************************************************

        Chapter II -- Auxiliary classes

  ... explain here when it is stable ...

******************************************************************************/

namespace {

BitCount::result_type BitCount::operator() (argument_type g) const

/*!
  Counts the number of set bits in g.
*/

{
  Grading gr(g);
  return gr.count();
}

GradingAction::GradingAction(const RootDatum& rd)
  :d_cartan(rd.semisimpleRank()),d_rank(rd.semisimpleRank())

/*!
  The idea is that d_cartan holds the reduced Cartan matrix.
*/

{
  for (unsigned long i = 0; i < rd.semisimpleRank(); ++i)
    for (unsigned long j = 0; j < rd.semisimpleRank(); ++j)
      if (rd.cartan(i,j)%2)
	d_cartan[i].set(j);
}

unsigned long GradingAction::operator() (unsigned long s, unsigned long g)
  const

/*!
  GradingAction is a function object, to be used for the determination of
  W-orbits on gradings. The result of opereator() on (s,g) is the action
  of the simple reflection #s on the grading j. This is obtained by letting
  (sg)(i) = g(si), where si is the in general non-simple root obtained
  from applying s to i; from the expression of si in the simple root basis
  we immediately get its grading.

  Note that since we are working modulo 2, si is either i or i+i_s, where
  i_s is the simple reflection corresponding to s; to tell which we only
  have to look at the Cartan matrix mod 2, which is stored in d_cartan.

  NOTE : this depends on the fact that RANK_MAX does not exceed the number
  of bits in an ulong.
*/

{
  Grading result;
  Grading gr(g);

  for (unsigned long j = 0; j < d_rank; ++j) { // compute bit #j of result
    bool b = gr.test(j);
    b ^= cartan(j,s) & gr.test(s);
    if (b)
      result.set(j);
  }

  return result.to_ulong();
}

void toGrading(Grading& g, const LongComponent&  gc)

/*!
  This is just a little function which extracts the first RANK_MAX bits
  from gc.
*/

{}

}

/*****************************************************************************

        Chapter III -- The GradingCompare class

  This is a function object that compares gradings in a size-first fashion.

******************************************************************************/

namespace gradings {

bool GradingCompare::operator() (const Grading& lhs, const Grading& rhs)

/*!
  Synopsis: compares lhs and rhs.

  We have lhs < rhs iff either lhs has fewer set members, or they have the
  number, and lhs < rhs in the usual sense.
*/

{
  if (lhs.count() < rhs.count())
    return true;

  if (lhs.count() > rhs.count())
    return false;

  return lhs < rhs;
}

}

/*****************************************************************************

        Chapter III -- Auxiliary functions

  ... explain here when it is stable ...

******************************************************************************/

namespace {

void noncompactEquations(LongComponentList& eqn, const WeightList& nc)

/*!
  This function adds to eqn the equations saying that the grading should
  be noncompact for each element of nc.
*/

{
  using namespace lattice;

  LongComponent e;

  for (unsigned long j = 0; j < nc.size(); ++j) {
    mod2(e,nc[j]);
    e.pushBack(true);
    eqn.push_back(e);
  }

  return;
}

void compactEquations(LongComponentList& eqn, const RootList& o, 
		      const RootList& rs, const RootDatum& rd)

/*!
  In this function we assume that o is an orthogonal set of roots in the
  root system rs. Let o_orth be its orthogonal. Then we add to eqn the 
  equations which express the fact that after descent through o, the roots in 
  o_orth have to be all compact (i.e., the grading on rs becomes the trivial 
  grading on o_orth.)

  Recall that when we descend through a short noncompact root a, and if
  b is a short root in rs s.t. a +- b is a root (in other words, a and b
  generate a system of type B2), then the descent through a _changes_ the
  grading of b.

  So the approach here is that we take a basis of the rootsystem in o_orth,
  and compute the degree of the elements in that basis through this rule. This
  then gives us the equations, to the number of rank(o_orth), that we have
  to add to eqn.
*/

{
  using namespace lattice;

  // make orthogonal root system

  RootList o_orth;
  makeOrthogonal(o_orth,o,rs,rd);

  // find basis

  RootList oob;
  rootBasis(oob,o_orth,rd);

  // find signs

  Grading g;

  for (unsigned long i = 0; i < oob.size(); ++i)
    for (unsigned long j = 0; j < o.size(); ++j)
      if (sumIsRoot(oob[i],o[j],rd))
	g.flip(i);

  // express oob in rs's rootbasis

  RootList rsb;
  rootBasis(rsb,rs,rd);

  WeightList woob;

  toRootBasis(oob.begin(),oob.end(),back_inserter(woob),rsb,rd);

  // write equations

  LongComponent e;

  for (unsigned long j = 0; j < woob.size(); ++j) {
    mod2(e,woob[j]);
    e.pushBack(g[j]);
    eqn.push_back(e);
  }

  return;
}

}

}
