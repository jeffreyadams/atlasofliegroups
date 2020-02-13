/*
  This is output.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2016,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "output.h"

#include <iostream>
#include <sstream>
#include <cassert>


#include "arithmetic.h"	// |remainder|
#include "partition.h"
#include "poset.h"
#include "tags.h"


#include "lietype.h"
#include "dynkin.h"
#include "gradings.h"
#include "innerclass.h"
#include "cartanclass.h"
#include "tori.h"
#include "realredgp.h"
#include "realweyl.h"	// |RealWeyl| class
#include "realweyl_io.h" // |printBlockStabilizer|

#include "basic_io.h"	// |seqPrint|
#include "ioutils.h"	// |foldLine|, |digits|
#include "prettyprint.h" // |printMatrix|
#include "graph.h"	// |OrientedGraph|

namespace atlas {

namespace output {

//			 Functions declared in output.h


  std::ostream& printBlockSizes(std::ostream& strm,
				const InnerClass& G,Interface& CI)
{
  const FormNumberMap rfi = CI.realFormInterface();
  const FormNumberMap drfi = CI.dualRealFormInterface();

  matrix::Matrix<arithmetic::big_int>
    block(G.numRealForms(),G.numDualRealForms());
  arithmetic::big_int maxEntry(0);

  for (size_t i = 0; i < block.numRows(); ++i)
    for (size_t j = 0; j < block.numColumns(); ++j)
    {
      block(i,j) = G.block_size(rfi.in(i),drfi.in(j));
      if (block(i,j) > maxEntry)
	maxEntry = block(i,j);
    }

  int width = ioutils::digits(maxEntry,10ul);
  prettyprint::printMatrix(strm,block,width+2);

  return strm;
}

namespace {
  // a class to help sorting real forms, based on "compactness" and grading
  struct RealFormData
  {
    RealFormNbr number;
    Grading grading;
    unsigned int depth;

    RealFormData(RealFormNbr rf, Grading gr, const RootNbrSet& so)
    : number(rf),grading(gr),depth(so.size()) {}

    bool operator< (const RealFormData& other) const
    { return depth != other.depth ? depth<other.depth : grading<other.grading; }

  };


  std::ostream& printType(std::ostream& out,
			  const Grading& gr,
			  const lietype::Layout& lo);

} // |namespace|


/*
  Build the standard real form interface for G and lo.

  Explanation: d_in[rf] is the inner numbering (the one used by G) of the
  weak real form rf; d_out[rf] is the outer numbering (the one used by the
  interface) of rf --- so d_in and d_out are inverses of each other. The
  names are names for the corresponding Lie algebras.
*/
FormNumberMap::FormNumberMap(const InnerClass& G,
			     const lietype::Layout& lo)
: d_in(G.numRealForms()), d_out(G.numRealForms()), d_name(G.numRealForms())
{
  const size_t nrf = G.numRealForms();
  const auto& rs = G.rootSystem();
  const auto& rd = G.rootDatum();

  std::vector<RealFormData> rf_data; rf_data.reserve(nrf);

  for (RealFormNbr rf = 0; rf<nrf; ++rf) // rf is internal real form number
  {
    RootNbrSet so = gradings::max_orth
      (G.noncompactRoots(rf),
       InvolutionData(rd,G.distinguished()).imaginary_roots(), rs);
    Grading gr = cartanclass::specialGrading
      (G.weak_real_partition(),rf,G.simple_roots_imaginary());
    rf_data.push_back(RealFormData(rf,gr,so));
  }

  std::sort(rf_data.begin(),rf_data.end());

  for (size_t i = 0; i<nrf; ++i) // now |i| is external number for |rf_data[i]|
  {
    d_in[i] = rf_data[i].number;
    d_out[d_in[i]] = i;
  }

  // write names
  std::ostringstream os;

  for (size_t i = 0; i<nrf; ++i)
  {
    os.str("");
    printType(os,rf_data[i].grading,lo);
    d_name[i] = os.str();
  }
} // |FormNumberMap::FormNumberMap|


// like the previous constructor, but for the _dual_ real forms.
FormNumberMap::FormNumberMap(const InnerClass& G,
			     const lietype::Layout& lo, tags::DualTag)
: d_in(G.numDualRealForms())
, d_out(G.numDualRealForms())
, d_name(G.numDualRealForms())
{
  const size_t ndrf = G.numDualRealForms();
  const auto& drd = G.dualRootDatum();
  const auto& drs = G.dualRootSystem();
  const lietype::Layout dlo = dual(lo);

  std::vector<RealFormData> rf_data; rf_data.reserve(ndrf);

  for (RealFormNbr drf = 0; drf<ndrf; ++drf)
  {
    RootNbrSet so = gradings::max_orth
      (G.parity_coroots(drf),
       InvolutionData(drd,G.dualDistinguished()).imaginary_roots(), drs);
    Grading gr = cartanclass::specialGrading
      (G.dual_weak_real_partition(),drf,G.simple_roots_real());
    rf_data.push_back(RealFormData(drf,gr,so));
  }

  std::sort(rf_data.begin(),rf_data.end());

  for (size_t i=0; i<ndrf; ++i)
  {
    d_in[i] = rf_data[i].number;
    d_out[d_in[i]] = i;
  }

  // write names
  std::ostringstream os;
  for (size_t i=0; i<ndrf; ++i)
  {
    os.str("");
    printType(os,rf_data[i].grading,dlo);
    d_name[i] = os.str();
  }
} // |FormNumberMap::FormNumberMap|, dual version


std::ostream& printRealForms(std::ostream& strm, const FormNumberMap& I)
{
  for (size_t i = 0; i < I.numRealForms(); ++i)
    std::cout << i << ": " << I.type_name(i) << std::endl;

  return strm;
}

/*
  Output the gradings corresponding to the various real forms
  defined for Cartan #cn.

  The gradings are output in the same order as the corresponding orbits are
  output in Fokko's "cartan" command.
*/
std::ostream& printGradings(std::ostream& strm,
			    const InnerClass& G, size_t cn, Interface& CI)
{
  const CartanClass& cc = G.cartan(cn);

  RealFormNbrList rfl(cc.numRealForms());
  const output::FormNumberMap& rfi = CI.realFormInterface();

  for (cartanclass::adjoint_fiber_orbit i = 0; i < rfl.size(); ++i)
    rfl[i] = rfi.out(G.realFormLabels(cn)[i]);

  printGradings(strm,cc.fiber(),rfl,G.rootDatum());

  return strm;
}


// Print information about the Cartan class #cn.
std::ostream& printCartanClass(std::ostream& strm,
			       const InnerClass& G,
			       size_t cn,
			       output::Interface& CI)
{
  const RootSystem& rs = G.rootDatum();

  const CartanClass& cc = G.cartan(cn);
  const auto& tau = cc.involution();
  const Fiber& f = cc.fiber();

  prettyprint::printTorusType(strm,tau) << std::endl;

  {
    std::ostringstream os;
    os << "canonical twisted involution: ";
    prettyprint::printWeylElt(os,G.involution_of_Cartan(cn),G.weylGroup());
    ioutils::foldLine(strm,os.str(),"",",") << std::endl;
  }

  size_t orbit_size=cc.orbitSize();
  strm << "twisted involution orbit size: " << orbit_size
       << "; fiber size: " << f.fiberSize()
       << "; strong inv: " << orbit_size*f.fiberSize()
       <<std::endl;

  // print type of imaginary root system
  LieType ilt = rs.subsystem_type(cc.simpleImaginary());

  if (ilt.size() == 0)
    strm << "imaginary root system is empty" << std::endl;
  else
    strm << "imaginary root system: " << ilt << std::endl;

  // print type of real root system
  LieType rlt = rs.subsystem_type(cc.simpleReal());

  if (rlt.size() == 0)
    strm << "real root system is empty" << std::endl;
  else
    strm << "real root system: " << rlt << std::endl;

  // print type of complex root system
  LieType clt = rs.subsystem_type(cc.simpleComplex());

  if (clt.size() == 0)
    strm << "complex factor is empty" << std::endl;
  else
    strm << "complex factor: " << clt << std::endl;

  RealFormNbrList rfl(cc.numRealForms());
  const output::FormNumberMap& rfi = CI.realFormInterface();

  for (cartanclass::adjoint_fiber_orbit i = 0; i < rfl.size(); ++i)
    rfl[i] = rfi.out(G.realFormLabels(cn)[i]);

  printFiber(strm,f,rfl);

  return strm;
}


// Print the fiber data.
std::ostream& printFiber(std::ostream& strm, const Fiber& f,
			 const RealFormNbrList& rfl)
{
  const Partition& pi = f.weakReal();
  unsigned long c = 0;

  for (Partition::iterator i(pi); i(); ++i,++c)
  {
    std::ostringstream os;
    os << "real form #";
    os << rfl[c] << ": ";
    basic_io::seqPrint(os,i->first,i->second,",","[","]");
    os << " (" << i->second - i->first << ")" << std::endl;
    ioutils::foldLine(strm,os.str(),"",",");
  }

  return strm;
}

/*
  Print the gradings of the simple-imaginary roots corresponding
  to the various real forms defined for |f|.

  Precondition: |rfl| contains the outer numbering of the real forms;

  The gradings are output in the same order as the orbit corresponding to the
  real form is output in the "cartan" command.
*/
std::ostream& printGradings(std::ostream& strm, const Fiber& f,
			    const RealFormNbrList& rfl,
			    const RootSystem& rs)
{
  typedef std::vector<unsigned long>::const_iterator VI;
  const auto afr = f.adjointFiberRank();

  const RootNbrList& si = f.simpleImaginary();
  int_Matrix cm = rs.cartanMatrix(si);
  dynkin::DynkinDiagram d(cm);
  Permutation a = dynkin::bourbaki(d);
  a.inv_conjugate(cm);

  strm << "cartan matrix of imaginary root system is:" << std::endl;
  prettyprint::printMatrix(strm,cm);

  const Partition& pi = f.weakReal();
  unsigned long c = 0;

  for (Partition::iterator it(pi); it(); ++it)
  {
    std::ostringstream os;
    os << "real form #";
    os << rfl[c] << ": ";

    const VI i_first = it->first;
    const VI i_last = it->second;

    os << '[';

    for (VI j = i_first; j != i_last; ++j)
    {
      if (j != i_first)
	os << ',';
      cartanclass::AdjointFiberElt afe(RankFlags(*j),afr);
      Grading gr= a.pull_back(f.grading(afe));
      prettyprint::prettyPrint(os,gr,f.simpleImaginary().size());
    }

    os << ']' << std::endl;

    std::string line = os.str();
    ioutils::foldLine(strm,line,"",",");
    ++c;
  }

  return strm;
}



std::ostream& printBlockStabilizer(std::ostream& strm,
				   RealReductiveGroup& G_R,
				   size_t cn, RealFormNbr drf)

{
  const InnerClass& G_C = G_R.innerClass();
  const RootDatum& rd = G_R.rootDatum();
  const WeylGroup& W = G_R.weylGroup();

  RealFormNbr rf = G_R.realForm();

  cartanclass::AdjointFiberElt x = G_C.representative(rf,cn);
  cartanclass::AdjointFiberElt y = G_C.dualRepresentative(drf,cn);
  const CartanClass& cc = G_C.cartan(cn);

  realweyl::RealWeyl rw(cc,x,y,rd,W);
  realweyl::RealWeylGenerators rwg(rw,cc,rd);

  realweyl_io::printBlockStabilizer(strm,rw,rwg);

  // check if the size is correct
  size::Size c;
  realweyl::blockStabilizerSize(c,rw);
  c *= G_C.fiberSize(rf,cn);
  c *= G_C.dualFiberSize(drf,cn);
  c *= cc.orbitSize();
  assert(c == W.order());

  return strm;
}

// Print information about all the Cartan classes for |G|.
std::ostream& printCartanClasses(std::ostream& strm,
				 RealReductiveGroup& G,
				 output::Interface& G_CI)
{
  const BitMap& b = G.innerClass().Cartan_set(G.realForm());
  bool first = true;

  for (BitMap::iterator i = b.begin(); i(); ++i)
  {
    if (first)
      first = false;
    else
      strm << std::endl << std::endl;
    strm << "Cartan #" << *i << ":" << std::endl;
    output::printCartanClass(strm,G.innerClass(),*i,G_CI);
  }

  return strm;
}

/* the function below is much like |kgb_io::printBruhatOrder|, but cannot
   call that function, which needs a different type, nameley |BruhatOrder|.
 */

// Print the Hasse diagram of the Cartan ordering of G_R.
std::ostream& printCartanOrder(std::ostream& strm,
			       const RealReductiveGroup& G_R)
{
  const poset::Poset& p = G_R.innerClass().Cartan_ordering();

  // get Hasse diagram of the Cartan classes for |G_R|
  graph::OrientedGraph g = p.hasseDiagram(G_R.mostSplit());

  strm << "0:" << std::endl; // this is the only minimal element

  // all covering relations in |g| are grouped by lower (covered) element
  for (size_t j = 1; j < g.size(); ++j)
  {
    const graph::EdgeList& e = g.edgeList(j);
    if (not e.empty()) // suppress other non-covered elements: not for |G_R|
      basic_io::seqPrint(strm << j << ": ",e.begin(),e.end()) << std::endl;
  }

  return strm;
}


/*
  Print the real Weyl group corresponding to Cartan #cn.

  Precondition: cartan #cn is defined for this real form.
*/
std::ostream& printRealWeyl(std::ostream& strm,
			    RealReductiveGroup& G_R, // modifiable for |cartan|
			    size_t cn)
{
  const InnerClass& G_C = G_R.innerClass();

  RealFormNbr rf = G_R.realForm();

  const RootDatum& rd = G_C.rootDatum();
  const WeylGroup& W = G_C.weylGroup();
  const CartanClass& cc = G_C.cartan(cn);
  const auto dafr = cc.dualFiber().adjointFiberRank();

  cartanclass::AdjointFiberElt x = G_C.representative(rf,cn);
  cartanclass::AdjointFiberElt y (RankFlags(0),dafr); // dual quasisplit form

  realweyl::RealWeyl rw(cc,x,y,rd,W);
  realweyl::RealWeylGenerators rwg(rw,cc,rd);

  realweyl_io::printRealWeyl(strm,rw,rwg);

  // check if the size is correct
  size::Size c;
  realweyl::realWeylSize(c,rw);
  c *= G_C.fiberSize(rf,cn);
  c *= cc.orbitSize();
  assert(c == W.order());

  return strm;
}

/*
  Synopsis: outputs information about the strong real forms of G.

  Explanation: the inverse image in \X of a class of weak real forms of G is
  of the form Z.\X(z), where z is an admissible value for x^2 for a strong real
  form of that class; so there is associated to the class of weak real forms
  a coset (1+\delta)(Z).z in Z^\delta. All the various \X(z) for all possible
  choices of z are isomorphic by Z-translation. Therefore it is enough to
  describe the combinatorial structure of one of them, for each class of weak
  real forms. This is a finite problem in all cases.

  We output the orbits of W_im in X(z), which correspond to the various strong
  real forms; we label them with the corresponding weak real form.
*/
std::ostream& printStrongReal(std::ostream& strm,
			      const InnerClass& G_C,
			      const output::FormNumberMap& rfi,
			      size_t cn)
{
  const CartanClass& cc = G_C.cartan(cn);
  const RealFormNbrList& rfl = G_C.realFormLabels(cn);

  size_t n = cc.numRealFormClasses();

  if (n>1)
    strm << "there are " << n << " real form classes:\n" << std::endl;

  for (cartanclass::square_class csc=0; csc<n; ++csc)
  {
    // print information about the square of real forms, in center
    {
      RealFormNbr wrf=cc.fiber().realFormPartition().classRep(csc);
      RealFormNbr fund_wrf= rfl[wrf]; // lift weak real form to |fund|
      cartanclass::square_class f_csc=G_C.xi_square(fund_wrf);

      // having the square class number of the fundamental fiber, get grading
      RatCoweight coch = some_coch(G_C,f_csc);
      Grading base_grading = grading_of_simples(G_C,coch);

      RatWeight z (G_C.rank());
      for (Grading::iterator it=base_grading.begin(); it(); ++it)
	z += G_C.rootDatum().fundamental_coweight(*it);

      Ratvec_Numer_t& zn = z.numerator();
      for (size_t i=0; i<z.size(); ++i)
        zn[i]=arithmetic::remainder(zn[i],z.denominator());
      strm << "class #" << f_csc
	   << ", possible square: exp(2i\\pi(" << z << "))" << std::endl;
    }

    const Partition& pi = cc.fiber_partition(csc);

    unsigned long c = 0;

    for (Partition::iterator i(pi); i(); ++i,++c)
    {
      std::ostringstream os;
      RealFormNbr rf = rfl[cc.toWeakReal(c,csc)];
      os << "real form #" << rfi.out(rf) << ": ";
      basic_io::seqPrint(os,i->first,i->second,",","[","]")
	<< " (" << i->second - i->first << ")" << std::endl;
      ioutils::foldLine(strm,os.str(),"",",");
    }

    if (n>1)
      strm << std::endl;
  }

  return strm;
}



//				Local functions

namespace {



// Print the complex form of |slt| (actually of two such factors).
std::ostream& printComplexType(std::ostream& strm,
			       const SimpleLieType& slt)
{
  size_t rk = slt.rank();
  switch (slt.type())
  {
  case 'A':
    strm << "sl(" << rk+1 << ",C)";
    break;
  case 'B':
    strm << "so(" << 2*rk+1 << ",C)";
    break;
  case 'C':
    strm << "sp(" << 2*rk << ",C)";
    break;
  case 'D':
    strm << "so(" << 2*rk << ",C)";
    break;
  case 'E':
    strm << "e" << rk << "(C)";
    break;
  case 'F':
  case 'f':
    strm << "f4(C)";
    break;
  case 'G':
  case 'g':
    strm << "g2(C)";
    break;
  case 'T':
    strm << "gl(1,C)";
    if (rk > 1) // this can no longer occur
    {
      strm << "^" << rk;
    }
    break;
  default:
    assert(false && "unexpected type in printComplexType");
    break;
  }

  return strm;
}

// we often need to print either '(n)' or '(n-m,m)' in |print_real_form_name|
// and when we do, we want them weakly decreasing, and suppressing a final 0
inline std::ostream& split(std::ostream& s,size_t n,size_t m)
{
  s << '(';
  if (m==0)
    s << n;
  else
  {
    if (2*m>n)
      m=n-m; // make sure larger number is listed first
    s << n-m << ',' << m;
  }
  return s << ')';
}

/*
  Print the name of the real form of simple type |slt| represented by |gr|,
  which is a special grading at the fundamental fiber of the simple roots (only
  imaginary ones give the actual grading, but complex ones retain a 0 bit).

  (Complex real forms do not involve a single simple type, do not come here.)

  As explained in |cartanclass::specialGrading|, we can find for any real form a
  corresponding grading of the simple roots that marks at most one root for each
  simple factor of the Lie type, and from that information one can deduce an
  identfying name for the real form (among other ones in its inner class).

  This function assumes it is called with such a grading, and that the simple
  roots of |lt| are in the standard order. We cannot however rely in all cases
  on the noncompact root being a specific one for the real form (some of them
  allow for more than one position), since (1) no effort was done in
  |specialGrading| to ensure such a choice, and (2) simple roots may have been
  in an unusual order at the point of selecting |gr|, so it would have been hard
  to do correctly at that time anyway. If that order was indeed unusual, the
  bits of |gr| have been since permuted to correct for it, so interpretation of
  bit positions becomes straightforward here.

  There is a subtle point though: the straightening of the Dynkin diagram was
  done without regard to the inner class involution. This means that we may have
  a usual numbering of |slt|, but in an unequal rank case with (effectively) an
  unusual automorphism determining which simple roots are imaginary (and
  therefore subject to actual grading). The problem arises only in unequal rank
  D4 and is easy to deal with: any nonzero grading defines so(5,3) there rather
  than so(7,1).
*/
std::ostream& print_real_form_name(std::ostream& strm, const Grading& gr,
				   const SimpleLieType& slt,
				   const lietype::TypeLetter ic)
{
  size_t rk = slt.rank();
  bool gr_trivial = gr.count()==0;
  size_t m = gr_trivial ? 0 : gr.firstBit()+1; // usually relevant for printing

  switch (slt.type())
  {
  case 'A':
    if (rk == 1)
      strm << (gr_trivial ? "su(2)" : "sl(2,R)");
    else
    {
      size_t n=rk+1;
      if (ic=='c')
	split(strm << "su",n,m);
      else // unequal rank $A_n$ case
	if (rk%2!=0 and gr_trivial) // both parts of this condition needed
	  strm << "sl(" << n/2 << ",H)";
	else
	  strm << "sl(" << n << ",R)";
    }
    break;
  case 'B':
    split(strm << "so",2*rk+1,2*m); // surprisingly easy to state
    break;
  case 'C':
    if (m == rk)
      strm << "sp(" << 2*rk << ",R)";
    else
      split (strm << "sp",rk,m);
    break;
  case 'D':
    {
      size_t n=2*rk;
      if (ic=='c' or (rk%2==0 and ic=='s')) // equal rank case
	if (m<rk-1) // so(2p,2q) form
	  split (strm << "so",n,2*m);
	else if (rk%2!=0)
	  strm << "so*(" << n << ")";
	else // so* type with label depending on m and parity of |rk/2|
	  strm << "so*(" << n << ((rk%4==0)==(m==rk) ? ")[1,0]" : ")[0,1]");
      else // unequal rank case
	if (rk>4)
	  split(strm << "so",n,2*m+1);
        else // unequal rank D4, we could have $m>=2$
	  strm << (gr_trivial ? "so(7,1)" : "so(5,3)");
    }
    break;
  case 'E':
    strm << 'e' << rk; // it always starts like this
    if (gr_trivial and (ic=='c' or rk>6)) // for |rk>6|, |'s'| means |'c'|
      break; // from |switch|: compact reak form
    strm << '(';
    if (rk==6) // E6, noncompact forms
      if (ic=='c')
	strm << (m==1 or m==6 ? "so(10).u(1)" : "su(6).su(2)");
      else // unequal rank
	strm << (gr_trivial ? "f4" : "R");
    else if (rk==7) // E7, noncompact forms
      strm << (m==7 ? "e6.u(1)" : m==2 or m==5 ? "R" : "so(12).su(2)");
    else // E8, noncompact forms: e8(e7.su(2)) positions 2,3,6,7, e8(R): 0,1,4,5
      strm <<(gr.any(RankFlags(0xCC)) ? "e7.su(2)" : "R");
    strm << ')';
    break;
  case 'f':
    m=5-m; // reverse interpretation of |m| for type 'f'
  case 'F':
    strm << (gr_trivial ? "f4" : m>=3 ? "f4(so(9))" : "f4(R)");
   break;
  case 'g': // for type G root order is unused
  case 'G':
    strm << (gr_trivial ? "g2" : "g2(R)");
    break;
  case 'T':
    strm << (ic=='c' ? "u(1)" : "gl(1,R)");
    if (rk > 1) // this can no longer occur
      strm << "^" << rk;
    break;
  default: // cannot happen
    assert(false && "unexpected type in printSimpleType");
    break;
  }

  return strm;
}


/*
  Print descriptive name of the real form represented by |d_gr|.

  Precondition: |d_gr| contains a grading of the simple roots (effectively only
  grading those that are imaginary), which represents the given real form, and
  which is special in that it contains a mimimal number of roots marked
  noncompact; this implies there is percisely one such simple root for each
  non-quasicompact non-complex non-torus simple factor of the Lie type.

  Torus factors $T_1$ may occur in |lo.d_type|, which have no corresponding bits
  in the grading |d_gr|, but such factors should match letters in |lo.d_inner|
  (either |'c'|, |'C'| or |'s'|, where |'C'| accounts for two factors $T_1$),
  which will be transmitted to |print_real_form_name| so that the correct Lie
  algebra name for the torus factors can be printed.
*/
std::ostream& printType(std::ostream& strm,
			const Grading& d_gr,
			const lietype::Layout& lo)
{
  const LieType& lt=lo.d_type;
  const InnerClassType& ict=lo.d_inner;

  // adapt grading to standard ordering of type |lt|
  Grading gr = lo.d_perm.pull_back(d_gr);

  size_t i = 0; // index of simple factor in |lt|

  for (size_t j = 0; j<ict.size(); gr>>=lt[i].semisimple_rank(),++i,++j)
  {
    SimpleLieType slt = lt[i];
    if (ict[j] == 'C')
    {
      printComplexType(strm,slt);
      gr >>= slt.semisimple_rank();  ++i; // skip one more simple Lie type
    }
    else
    {
      Grading grs = gr;
      grs.truncate(slt.semisimple_rank());
      print_real_form_name(strm,grs,slt,ict[j]);
    }
    if (j < ict.size()-1)
      strm << ".";
  }

  return strm;
}


} // |namespace|

} // |namespace output|

} // |namespace atlas|
