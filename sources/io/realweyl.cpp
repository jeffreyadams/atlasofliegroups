/*
  This is realweyl.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "realweyl.h"

#include "bitvector.h"
#include "cartanclass.h"
#include "rootdata.h"
#include "weyl.h"
#include "weylsize.h"


namespace atlas {

namespace {
void orthogonalMAlpha(RootNbrList&,
		      unsigned long,
		      const Fiber&,
		      const RootDatum&);
void rGenerators(SmallBitVectorList&,
		 const RootNbrList&,
		 const Fiber&,
		 const RootDatum&);
}

/*****************************************************************************

        Chapter I -- The RealWeyl class

******************************************************************************/

namespace realweyl {

RealWeyl::RealWeyl(const CartanClass& cc,
		   unsigned long x, unsigned long y,
		   const RootDatum& rd, const WeylGroup& W)
  :d_group(&W)

{
  using namespace cartanclass;
  using namespace rootdata;
  using namespace tags;

  d_imaginary = cc.simpleImaginary();
  d_real = cc.simpleReal();
  d_complex = cc.simpleComplex();

  RootDatum drd(rd,DualTag());

  d_complexType = rd.Lie_type(d_complex);
  d_imaginaryType = rd.Lie_type(d_imaginary);
  d_realType = drd.Lie_type(d_real);

  const Fiber& f = cc.fiber();
  const Fiber& df = cc.dualFiber();

  // construct basis of imaginary compact roots
  d_imaginaryCompact= rd.simpleBasis(f.compactRoots(x));
  d_imaginaryCompactType = rd.Lie_type(d_imaginaryCompact);

  // construct imaginary R-group
  orthogonalMAlpha(d_imaginaryOrth,x,f,rd);
  rGenerators(d_imaginaryR,d_imaginaryOrth,f,rd);

  // construct basis of real compact roots
  d_realCompact = drd.simpleBasis(df.compactRoots(y));
  d_realCompactType = drd.Lie_type(d_realCompact);

  // construct real R-group
  orthogonalMAlpha(d_realOrth,y,df,drd);
  rGenerators(d_realR,d_realOrth,df,drd);
}

}

/*****************************************************************************        Chapter II -- The RealWeylGenerators class

******************************************************************************/

namespace realweyl {

RealWeylGenerators::RealWeylGenerators(const RealWeyl& rw,
				       const CartanClass& cc,
				       const RootDatum& rd)
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
    const SmallBitVector& c = rw.imaginaryR(j);
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
    const SmallBitVector& c = rw.realR(j);
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
    RootNbr rn = rw.complex(j);
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

/*
  Synopsis: puts in rl the list of positive imaginary roots that are orthogonal
  to the sum of positive imaginary compact roots for the fiber.

  Explanation: the imaginary roots orthogonal to 2 rho_ic constitute a set
  of strongly orthogonal roots (they must all be non-compact, and therefore
  no sum of two of them can be a root.) So their root system is an A_1^n.

  NOTE: they are even super-orthogonal, see IC4, prop 3.20.
*/
void orthogonalMAlpha(RootNbrList& rl, unsigned long x,
		      const Fiber& f,
		      const RootDatum& rd)
{
  Weight tworho_ic = compactTwoRho(x,f,rd);

  // put positive imaginary noncompact roots in rs
  RootNbrSet rs = f.noncompactRoots(x);
  rs &= rd.posRootSet();

  // keep the ones that are orthogonal to tworho_ic
  for (RootNbrSet::iterator it=rs.begin(); it(); ++it)
    if (rd.isOrthogonal(tworho_ic,*it))
      rl.push_back(*it);
}

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
void rGenerators(SmallBitVectorList& cl,
		 const RootNbrList& rl,
		 const Fiber& f,
		 const RootDatum& rd)
{
  size_t rln = rl.size();
  BinaryMap m(f.fiberRank(),rln);

  for (size_t j=0; j<rln; ++j)
  {
    SmallBitVector v = // |v.size()=f.fiberRank()|
      f.mAlpha(rd.coroot(rl[j]));
    for (size_t i = 0; i < v.size(); ++i)
      m.set(i,j,v[i]);
  }

  m.kernel(cl);
}

}

}
