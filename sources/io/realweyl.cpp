/*
  This is realweyl.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include "realweyl.h"

#include "bitvector.h"
#include "cartanclass.h"
#include "latticetypes.h"
#include "rootdata.h"
#include "weyl.h"
#include "weylsize.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace {
void orthogonalMAlpha(rootdata::RootList&, unsigned long, 
		      const cartanclass::Fiber&, const rootdata::RootDatum&);
void rGenerators(latticetypes::ComponentList&, const rootdata::RootList&,
		 const cartanclass::Fiber&, const rootdata::RootDatum&);
}

/*****************************************************************************

        Chapter I -- The RealWeyl class

  ... explain here when it is stable ...

******************************************************************************/

namespace realweyl {

RealWeyl::RealWeyl(const cartanclass::CartanClass& cc, 
		   unsigned long x, unsigned long y,
		   const rootdata::RootDatum& rd, const weyl::WeylGroup& W)
  :d_group(&W)

{
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tags;

  d_imaginary = cc.simpleImaginary();
  d_real = cc.simpleReal();
  d_complex = cc.simpleComplex();

  RootDatum drd(rd,DualTag());

  lieType(d_complexType,d_complex,rd);
  lieType(d_imaginaryType,d_imaginary,rd);
  lieType(d_realType,d_real,drd);

  const Fiber& f = cc.fiber();
  const Fiber& df = cc.dualFiber();

  // construct basis of imaginary compact roots
  RootSet rs;

  f.noncompactRootSet(rs,x);
  ~rs;
  rs &= f.imaginaryRootSet();
  rootBasis(d_imaginaryCompact,rs,rd);
  lieType(d_imaginaryCompactType,d_imaginaryCompact,rd);

  // construct imaginary R-group
  orthogonalMAlpha(d_imaginaryOrth,x,f,rd);
  rGenerators(d_imaginaryR,d_imaginaryOrth,f,rd);

  // construct basis of real compact roots
  df.noncompactRootSet(rs,y);
  ~rs;
  rs &= df.imaginaryRootSet();
  rootBasis(d_realCompact,rs,drd);
  lieType(d_realCompactType,d_realCompact,drd);

  // construct real R-group
  orthogonalMAlpha(d_realOrth,y,df,drd);
  rGenerators(d_realR,d_realOrth,df,drd);
}

}

/*****************************************************************************        Chapter II -- The RealWeylGenerators class

  ... explain here when it is stable ...

******************************************************************************/

namespace realweyl {

RealWeylGenerators::RealWeylGenerators(const RealWeyl& rw,
				       const cartanclass::CartanClass& cc, 
				       const rootdata::RootDatum& rd)
  :d_group(&rw.weylGroup())
  

/*
  Synopsis: extracts the various generator lists from the data contained in rw.
*/

{
  using namespace rootdata;
  using namespace weyl;

  const WeylGroup& W = *d_group;
  WeylElt e; // the identity element

  size_t nic = rw.numImaginaryCompact();
  d_imaginaryCompact.assign(nic,e);

  for (size_t j = 0; j < nic; ++j) {
    WeylWord ww;
    toWeylWord(ww,rw.imaginaryCompact(j),rd);
    W.prod(d_imaginaryCompact[j],ww);
  }

  size_t nir = rw.numImaginaryR();
  d_imaginaryR.assign(nir,e);

  for (size_t j = 0; j < nir; ++j) {
    const latticetypes::Component& c = rw.imaginaryR(j);
    for (size_t i = 0; i < c.size(); ++i)
      if (c.test(i)) {
	WeylWord ww;
	toWeylWord(ww,rw.imaginaryOrth(i),rd);
	W.prod(d_imaginaryR[j],ww);
      }
  }

  size_t ni = rw.numImaginary();
  d_imaginary.assign(ni,e);

  for (size_t j = 0; j < ni; ++j) {
    WeylWord ww;
    toWeylWord(ww,rw.imaginary(j),rd);
    W.prod(d_imaginary[j],ww);
  }

  size_t nrc = rw.numRealCompact();
  d_realCompact.assign(nrc,e);

  for (size_t j = 0; j < nrc; ++j) {
    WeylWord ww;
    toWeylWord(ww,rw.realCompact(j),rd);
    W.prod(d_realCompact[j],ww);
  }

  size_t nrr = rw.numRealR();
  d_realR.assign(nrr,e);

  for (size_t j = 0; j < nrr; ++j) {
    const latticetypes::Component& c = rw.realR(j);
    for (size_t i = 0; i < c.size(); ++i)
      if (c.test(i)) {
	WeylWord ww;
	toWeylWord(ww,rw.realOrth(i),rd);
	W.prod(d_realR[j],ww);
      }
  }

  size_t nr = rw.numReal();
  d_real.assign(nr,e);

  for (size_t j = 0; j < nr; ++j) {
    WeylWord ww;
    toWeylWord(ww,rw.real(j),rd);
    W.prod(d_real[j],ww);
  }

  size_t nc = rw.numComplex();
  d_complex.assign(nc,e);

  for (size_t j = 0; j < nc; ++j) {
    WeylWord ww;
    rootdata::RootNbr rn = rw.complex(j);
    toWeylWord(ww,rn,rd);
    W.prod(d_complex[j],ww);
    rn = cc.rootInvolution(rn);
    toWeylWord(ww,rn,rd);
    W.prod(d_complex[j],ww);
  }

  return;
}
}

/*****************************************************************************

        Chapter III -- Functions declared in realweyl.h

  ... explain here when it is stable ...

******************************************************************************/

namespace realweyl {

void blockStabilizerSize(size::Size& c, const RealWeyl& rw)

/*
  Synopsis: puts in c the size of the block stabilizer for rw.
*/

{
  using namespace size;
  using namespace weylsize;

  c.reset();
  Size cp;

  weylSize(cp,rw.complexType());
  c *= cp;
  weylSize(cp,rw.imaginaryCompactType());
  c *= cp;
  weylSize(cp,rw.realCompactType());
  c *= cp;

  size_t twoPower = rw.numImaginaryR();
  twoPower += rw.numRealR();
  c.twoShift(twoPower);

  return;
}

void dualRealWeylSize(size::Size& c, const RealWeyl& rw)

/*
  Synopsis: puts in c the size of the dual real weyl group for rw.
*/

{
  using namespace size;
  using namespace weylsize;

  c.reset();
  Size cp;

  weylSize(cp,rw.complexType());
  c *= cp;
  weylSize(cp,rw.imaginaryType());
  c *= cp;
  weylSize(cp,rw.realCompactType());
  c *= cp;

  size_t twoPower = rw.numRealR();
  c.twoShift(twoPower);

  return;
}

void realWeylSize(size::Size& c, const RealWeyl& rw)

/*
  Synopsis: puts in c the size of the real weyl group for rw.
*/

{
  using namespace size;
  using namespace weylsize;

  c.reset();
  Size cp;

  weylSize(cp,rw.complexType());
  c *= cp;
  weylSize(cp,rw.imaginaryCompactType());
  c *= cp;
  weylSize(cp,rw.realType());
  c *= cp;

  size_t twoPower = rw.numImaginaryR();
  c.twoShift(twoPower);

  return;
}

}

/*****************************************************************************

        Chapter IV -- Local functions

  ... explain here when it is stable ...

******************************************************************************/

namespace {

void orthogonalMAlpha(rootdata::RootList& rl, unsigned long x, 
		      const cartanclass::Fiber& f, 
		      const rootdata::RootDatum& rd)

/*
  Synopsis: puts in rl the list of positive imaginary roots that are orthogonal
  to the sum of positive imaginary compact roots for the fiber.

  Explanation: the imaginary roots orthogonal to 2 rho_ic constitute a set
  of strongly orthogonal roots (they must all be non-compact, and therefore
  no sum of two of them can be a root.) So their root system is an A_1^n.

  NOTE: they are even super-orthogonal, see IC4, prop 3.20.
*/

{
  using namespace latticetypes;
  using namespace rootdata;

  Weight tworho_ic;
  compactTwoRho(tworho_ic,x,f,rd);

  RootSet rs;

  // put positive imaginary noncompact roots in rs
  f.noncompactRootSet(rs,x);
  rs &= rd.posRootSet();

  // keep the ones that are orthogonal to tworho_ic
  RootSet::iterator rs_end = rs.end();

  for (RootSet::iterator i = rs.begin(); i != rs_end; ++i)
    if (rd.isOrthogonal(tworho_ic,*i))
      rl.push_back(*i);

  return;
}

void rGenerators(latticetypes::ComponentList& cl, const rootdata::RootList& rl,
		 const cartanclass::Fiber& f, const rootdata::RootDatum& rd)

/*
  Synopsis: puts in cl a basis of the R-group for this torus and this
  real form.

  Precondition: rl contains the list of positive imaginary roots orthogonal
  to the sum of positive imaginary compact roots (those are noncompact,
  and pairwise strongly orthogonal);

  Explanation: let n be the number of elements in rl. We want to look at
  the map (Z_2)^n -> elts. of order two in f's fiber group, that takes each 
  element in rl to the corresponding m_alpha; the R-group is generated by the 
  products of simple reflections in rl corresponding to the combinations that 
  lie in the _kernel_ of this map.
*/

{
  using namespace bitvector;
  using namespace constants;
  using namespace latticetypes;

  size_t rln = rl.size();
  BitMatrix<RANK_MAX> m(f.fiberRank(),rln);

  for (size_t j = 0; j < rln; ++j) {
    Component v(f.fiberRank());
    f.mAlpha(v,rd.coroot(rl[j]));
    for (size_t i = 0; i < f.fiberRank(); ++i)
      if (v[i])
	m.set(i,j);
  }

  m.kernel(cl);

  return;
}

}

}
