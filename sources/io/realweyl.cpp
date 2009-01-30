/*
  This is realweyl.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "realweyl.h"

#include "bitvector.h"
#include "cartanclass.h"
#include "latticetypes.h"
#include "rootdata.h"
#include "weyl.h"
#include "weylsize.h"


namespace atlas {

namespace {
void orthogonalMAlpha(rootdata::RootList&,
		      unsigned long,
		      const cartanclass::Fiber&,
		      const rootdata::RootDatum&);
void rGenerators(latticetypes::SmallBitVectorList&,
		 const rootdata::RootList&,
		 const cartanclass::Fiber&,
		 const rootdata::RootDatum&);
}

/*****************************************************************************

        Chapter I -- The RealWeyl class

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
  rootdata::rootBasis(d_imaginaryCompact,f.compactRoots(x),rd);
  lieType(d_imaginaryCompactType,d_imaginaryCompact,rd);

  // construct imaginary R-group
  orthogonalMAlpha(d_imaginaryOrth,x,f,rd);
  rGenerators(d_imaginaryR,d_imaginaryOrth,f,rd);

  // construct basis of real compact roots
  rootBasis(d_realCompact,df.compactRoots(y),drd);
  lieType(d_realCompactType,d_realCompact,drd);

  // construct real R-group
  orthogonalMAlpha(d_realOrth,y,df,drd);
  rGenerators(d_realR,d_realOrth,df,drd);
}

}

/*****************************************************************************        Chapter II -- The RealWeylGenerators class

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
    W.mult(d_imaginaryCompact[j],rd.reflectionWord(rw.imaginaryCompact(j)));
  }

  size_t nir = rw.numImaginaryR();
  d_imaginaryR.assign(nir,e);

  for (size_t j = 0; j < nir; ++j) {
    const latticetypes::SmallBitVector& c = rw.imaginaryR(j);
    for (size_t i = 0; i < c.size(); ++i)
      if (c[i]) {
	W.mult(d_imaginaryR[j],rd.reflectionWord(rw.imaginaryOrth(i)));
      }
  }

  size_t ni = rw.numImaginary();
  d_imaginary.assign(ni,e);

  for (size_t j = 0; j < ni; ++j) {
    W.mult(d_imaginary[j],rd.reflectionWord(rw.imaginary(j)));
  }

  size_t nrc = rw.numRealCompact();
  d_realCompact.assign(nrc,e);

  for (size_t j = 0; j < nrc; ++j) {
    W.mult(d_realCompact[j],rd.reflectionWord(rw.realCompact(j)));
  }

  size_t nrr = rw.numRealR();
  d_realR.assign(nrr,e);

  for (size_t j = 0; j < nrr; ++j) {
    const latticetypes::SmallBitVector& c = rw.realR(j);
    for (size_t i = 0; i < c.size(); ++i)
      if (c[i]) {
	W.mult(d_realR[j],rd.reflectionWord(rw.realOrth(i)));
      }
  }

  size_t nr = rw.numReal();
  d_real.assign(nr,e);

  for (size_t j = 0; j < nr; ++j) {
    W.mult(d_real[j],rd.reflectionWord(rw.real(j)));
  }

  size_t nc = rw.numComplex();
  d_complex.assign(nc,e);

  for (size_t j = 0; j < nc; ++j) {
    rootdata::RootNbr rn = rw.complex(j);
    W.mult(d_complex[j],rd.reflectionWord(rn));
    rn = cc.involution_image_of_root(rn);
    W.mult(d_complex[j],rd.reflectionWord(rn));
  }

  return;
}
}

/*****************************************************************************

        Chapter III -- Functions declared in realweyl.h

******************************************************************************/

namespace realweyl {


/*
  Synopsis: puts in c the size of the block stabilizer for rw.
*/
void blockStabilizerSize(size::Size& c, const RealWeyl& rw)
{
  c =  weylsize::weylSize(rw.complexType());
  c *= weylsize::weylSize(rw.imaginaryCompactType());
  c *= weylsize::weylSize(rw.realCompactType());

  size_t twoPower = rw.numImaginaryR();
  twoPower += rw.numRealR();
  c.twoShift(twoPower);

}


/*
  Synopsis: puts in c the size of the real weyl group for rw.
*/
void realWeylSize(size::Size& c, const RealWeyl& rw)
{
  c =  weylsize::weylSize(rw.complexType());
  c *= weylsize::weylSize(rw.imaginaryCompactType());
  c *= weylsize::weylSize(rw.realType());

  size_t twoPower = rw.numImaginaryR();
  c.twoShift(twoPower);
}


/*
  Synopsis: puts in c the size of the dual real weyl group for rw.
*/
void dualRealWeylSize(size::Size& c, const RealWeyl& rw)
{
  c =  weylsize::weylSize(rw.complexType());
  c *= weylsize::weylSize(rw.imaginaryType());
  c *= weylsize::weylSize(rw.realCompactType());

  size_t twoPower = rw.numRealR();
  c.twoShift(twoPower);
}

} // |namespace realweyl|

/*****************************************************************************

        Chapter IV -- Local functions

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

  Weight tworho_ic = compactTwoRho(x,f,rd);

  // put positive imaginary noncompact roots in rs
  RootSet rs = f.noncompactRoots(x);
  rs &= rd.posRootSet();

  // keep the ones that are orthogonal to tworho_ic
  RootSet::iterator rs_end = rs.end();

  for (RootSet::iterator i = rs.begin(); i != rs_end; ++i)
    if (rd.isOrthogonal(tworho_ic,*i))
      rl.push_back(*i);
}

void rGenerators(latticetypes::SmallBitVectorList& cl,
		 const rootdata::RootList& rl,
		 const cartanclass::Fiber& f,
		 const rootdata::RootDatum& rd)

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
    SmallBitVector v = f.mAlpha(rd.coroot(rl[j])); // |v.size()=f.fiberRank()|
    for (size_t i = 0; i < v.size(); ++i)
      m.set(i,j,v[i]);
  }

  m.kernel(cl);
}

}

}
