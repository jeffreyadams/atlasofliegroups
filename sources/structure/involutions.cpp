/*
  This is involutions.cpp.

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include <iostream>

#include "involutions.h"

#include "arithmetic.h"
#include "matreduc.h"
#include "rootdata.h"
#include "weyl.h"
#include "y_values.h"
#include "tits.h"
#include "innerclass.h"
#include "lattice.h"

#include "kgb.h" // for |KGB_elt_entry|

namespace atlas {

namespace involutions {

// ------------------------------ InvolutionData ------------------------------

InvolutionData::InvolutionData(const RootDatum& rd,
			       const WeightInvolution& theta)
  : root_perm(rd.rootPermutation(theta))
  , d_imaginary(rd.numRoots()) // we can only dimension root sets for now
  , d_real(rd.numRoots())
  , d_complex(rd.numRoots())
  , d_simpleImaginary()        // here even dimensioning is pointless
  , d_simpleReal()
{ classify_roots(rd); }

void InvolutionData::classify_roots(const RootSystem& rs)
{
  for (RootNbr alpha = 0; alpha<rs.numRoots(); ++alpha)
    if (root_perm[alpha] == alpha)
      d_imaginary.insert(alpha);
    else if (root_perm[alpha] == rs.rootMinus(alpha))
      d_real.insert(alpha);
    else
      d_complex.insert(alpha);

  // find simple-imaginary roots
  d_simpleImaginary=rs.simpleBasis(imaginary_roots());
  d_simpleReal=rs.simpleBasis(real_roots());
}

// the follwing constructor intersects all root sets with those of |sub|, but
// uses the full parent datum for numbering roots, and in |root_involution()|
InvolutionData::InvolutionData(const RootDatum& rd,
			       const WeightInvolution& theta,
			       const RootNbrSet& positive_subsystem)
: root_perm(rd.rootPermutation(theta))
, d_imaginary(rd.numRoots())
, d_real(d_imaginary.capacity())
, d_complex(d_imaginary.capacity())
, d_simpleImaginary()        // here even dimensioning is pointless
, d_simpleReal()
{
  for (RootNbrSet::const_iterator it=positive_subsystem.begin(); it(); ++it)
  {
    RootNbr alpha = *it;
    RootNbr minus_alpha = rd.rootMinus(alpha);
    if (root_perm[alpha] == alpha)
    {
      d_imaginary.insert(alpha);
      d_imaginary.insert(minus_alpha);
    }
    else if (root_perm[alpha] == minus_alpha)
    {
      d_real.insert(alpha);
      d_real.insert(minus_alpha);
    }
    else
    {
      d_complex.insert(alpha);
      d_complex.insert(minus_alpha);
    }
  }
  // find simple-imaginary roots
  d_simpleImaginary=rd.simpleBasis(imaginary_roots());
  d_simpleReal=rd.simpleBasis(real_roots());
}

InvolutionData::InvolutionData(const RootSystem& rs,
			       const RootNbrList& s_image)
  : root_perm(rs.extend_to_roots(s_image))
  , d_imaginary(rs.numRoots()) // we can only dimension root sets for now
  , d_real(rs.numRoots())
  , d_complex(rs.numRoots())
  , d_simpleImaginary()
  , d_simpleReal()
{ classify_roots(rs); }

InvolutionData InvolutionData::build
  (const RootSystem& rs,
   const TwistedWeylGroup& W,
   const TwistedInvolution& tw)
{ return InvolutionData(rs,W.simple_images(rs,tw)); }

void InvolutionData::swap(InvolutionData& other)
{
  root_perm.swap(other.root_perm);
  d_imaginary.swap(other.d_imaginary);
  d_real.swap(other.d_real);
  d_complex.swap(other.d_complex);
  d_simpleImaginary.swap(other.d_simpleImaginary);
  d_simpleReal.swap(other.d_simpleReal);
}


void InvolutionData::cross_act(const Permutation& root_reflection)
{
  root_reflection.renumber(root_perm); // first half of conjugation
  // then do |root_reflection.permute(root_perm)|, exploiting involution
  RootNbr beta;
  for (RootNbr alpha=0; alpha<root_perm.size(); ++alpha)
    if (alpha<(beta=root_reflection[alpha])) // do each 2-cycle just once
      std::swap(root_perm[alpha],root_perm[beta]);

  d_imaginary = root_reflection.renumbering(d_imaginary);
  d_real = root_reflection.renumbering(d_real);
  d_complex = root_reflection.renumbering(d_complex);
  d_simpleImaginary = root_reflection.renumbering(d_simpleImaginary);
  d_simpleReal = root_reflection.renumbering(d_simpleReal);
}



// ------------------------------ InvolutionTable -----------------------------



// describe involution by Cartan class 0 involution |cross|, and Cayley roots
RootNbrList Cayley_roots(const TwistedInvolution& tw,
			 const RootSystem& rs,
			 const TwistedWeylGroup& W,
			 TwistedInvolution& cross) // gets twisted-conjugated
{
  weyl::InvolutionWord inv_ex = W.involution_expr(tw);

  RootNbrList result; // roots whose reflections give Cayley part
  result.reserve(rs.rank()); // an upper bound to size needed

  for (size_t i=inv_ex.size(); i-->0; )
    if (inv_ex[i]>=0) // Cayley transform by simple root
      result.push_back(rs.simpleRootNbr(inv_ex[i]));
    else // cross action by simple root
    {
      const weyl::Generator s=~inv_ex[i];
      W.twistedConjugate(cross,s); // record cross action
      for (size_t i=0; i<result.size(); ++i)
	rs.simple_reflect_root(s,result[i]); // |s|-reflect the |result| roots
    }

  return result;
}


InvolutionTable::InvolutionTable
  (const RootDatum& r, const WeightInvolution& d,  const TwistedWeylGroup& t)
: rd(r), delta(d), tW(t)
, pool(), hash(pool), data()
, torus_simple_reflection()
{
  torus_simple_reflection.reserve(tW.rank());

  for (weyl::Generator s=0; s<tW.rank(); ++s)
  {
    WeightInvolution M(rd.rank()); // identity
    rd.simple_reflect(s,M);        // now it is the reflection matrix on $X^*$
    torus_simple_reflection.push_back(BinaryMap(M.transposed())); // |^M%2|
  }

}

// a method used to create the first element of a Cartan orbit of involutions
InvolutionNbr InvolutionTable::add_involution
  (const TwistedInvolution& canonical)
{
  InvolutionNbr result=hash.match(weyl::TI_Entry(canonical));
  if (result<data.size()) return result; // skip if |canonical| already known

  // compute the involution matrix |theta| using |Cayley_roots|
  const WeylGroup& W= tW.weylGroup();
  TwistedInvolution cross; // records $w\in W$ conjugating |delta| to |theta|
  RootNbrList Cayleys = Cayley_roots(canonical,rd,tW,cross);

  WeightInvolution theta = delta;
  W.act(rd,cross,theta); // taking |cross| as |WeylElt| multiplies theta<-delta
  for (auto& alpha : Cayleys)
    rd.reflect(alpha,theta);

  int_Matrix A=theta; // will contain |id-theta|, later row-saturated
  A.negate() += 1;

  int_Matrix col; bool flip;
  matreduc::column_echelon(A,col,flip);  // now |A| holds basis for image
  int_Matrix M_real = col.inverse().block(0,0,A.numColumns(),A.numRows());

  unsigned int W_length=W.length(canonical);
  unsigned int length = (W_length+Cayleys.size())/2;
  data.push_back(record(theta,InvolutionData(rd,theta),
			M_real, A, // |lift_mat|
			length,W_length,
			tits::fiber_denom(theta))); // |mod_space| for |x|
  assert(data.size()==hash.size());

  return result;
}

// extend from known element |n| to conjugate by |s|
InvolutionNbr InvolutionTable::add_cross(weyl::Generator s, InvolutionNbr n)
{
  weyl::TI_Entry tw = involution(n);
  int d =tW.twistedConjugate(tw,s); // modify |tw|, |d| is W-length difference
  InvolutionNbr result=hash.match(tw);
  if (result<data.size()) return result;

  data.push_back(data[n]); // start out with a copy of the data
  record& me=data.back();

  rd.simple_reflect(s,me.theta);
  rd.simple_reflect(me.theta,s); // not |twisted(s)|: |delta| is incorporated
  rd.simple_reflect(me.M_real,s); // apply $s$ before |M_real|
  rd.simple_reflect(s,me.lift_mat); // and apply it after |lift_mat|
  me.id.cross_act(rd.simple_root_permutation(s));
  me.length   += d/2;
  me.W_length += d;
  me.mod_space.apply(torus_simple_reflection[s]);

  assert(data.size()==hash.size());

  return result;
}

bool
InvolutionTable::is_complex_simple(InvolutionNbr n,weyl::Generator s) const
{ return complex_roots(n).isMember(rd.simpleRootNbr(s)); }

bool
InvolutionTable::is_imaginary_simple(InvolutionNbr n,weyl::Generator s) const
{ return imaginary_roots(n).isMember(rd.simpleRootNbr(s)); }

bool
InvolutionTable::is_real_simple(InvolutionNbr n,weyl::Generator s) const
{ return real_roots(n).isMember(rd.simpleRootNbr(s)); }

// whether a complex root |alpha| is a descent at involution |n|
bool InvolutionTable::complex_is_descent(InvolutionNbr n,RootNbr alpha) const
{ // whether involution inverses positivity when applied to |alpha|
  return rd.is_posroot(alpha)==rd.is_negroot(root_involution(n,alpha));
}

// this method makes involution table usable in X command, even if inefficient
KGB_elt_entry InvolutionTable::x_pack(const GlobalTitsElement& x) const
{
  const TwistedInvolution& tw= x.tw();
  InvolutionNbr i = nr(tw);
  assert(i<hash.size());
  RatWeight wt = x.torus_part().log_2pi();
  // we need projector modulo kernel of |theta^tr+1|, cf. constructor
  int_Matrix A = matrix(i);
  A.transpose() += 1;
  int_Matrix projector = lattice::row_saturate(A);
  Ratvec_Numer_t p = projector * wt.numerator();

  // reduce modulo integers and return
  for (size_t j=0; j<p.size(); ++j)
    p[j]= arithmetic::remainder(p[j],wt.denominator());
  return KGB_elt_entry(RatWeight(p,wt.denominator()).normalize(),x);
}

bool
InvolutionTable::x_equiv(const GlobalTitsElement& x0,
			 const GlobalTitsElement& x1) const
{
  if (x0.tw()!=x1.tw())
    return false;

  InvolutionNbr i = nr(x0.tw());
  assert(i<hash.size());
  RatWeight wt = x0.torus_part().log_2pi()-x1.torus_part().log_2pi();

  // we need projector modulo kernel of |theta^tr+1|, cf. constructor
  int_Matrix A = matrix(i);
  A.transpose() += 1;
  int_Matrix projector = lattice::row_saturate(A);
  Ratvec_Numer_t p = projector * wt.numerator();

  for (size_t i=0; i<p.size(); ++i)
    if (p[i]%wt.denominator()!=0)
      return false;

  return true;
}


// choose unique representative for real projection class of a rational weight
void InvolutionTable::real_unique(InvolutionNbr inv, RatWeight& y) const
{
  const record& rec=data[inv];
  Ratvec_Numer_t v = rec.M_real * y.numerator();
  // reduce $v=(1-theta)y$ expressed on image $(1-theta)X^*$-basis modulo 2
  for (unsigned i=0; i<v.size(); ++i)
    v[i]= arithmetic::remainder(v[i],2*y.denominator());

  y.numerator()= rec.lift_mat * v; // original |y| now "mapped to" |(1-theta)*y|
  (y/=2).normalize(); // and this gets us back to the class of the original |y|
}

Weight InvolutionTable::y_lift(InvolutionNbr inv, TorusPart y_part) const
{
  const record& rec=data[inv];
  Weight result(rec.lift_mat.numRows(),0);
  // set |result=lift_mat*lift| where |lift| is the integer lift of |y_part|
  for (unsigned j=0; j<y_part.size(); ++j)
    if (y_part[j])
      for (unsigned i=0; i<result.size(); ++i)
	result[i] += rec.lift_mat(i,j);
  return result;
}

// ------------------------------ Cartan_orbit --------------------------------



Cartan_orbit::Cartan_orbit(InvolutionTable& i_tab,
			   const InnerClass& G,
			   CartanNbr cn)
  : Cartan_class_nbr(cn)
  , start(i_tab.size())
  , size(G.cartan(cn).orbitSize())
{
  const TwistedInvolution& canonical = G.involution_of_Cartan(cn);
  assert(i_tab.unseen(canonical)); // should only generate orbit once

  i_tab.reserve(start+size);
  for (InvolutionNbr i=i_tab.add_involution(canonical); i<i_tab.size(); ++i)
    for (weyl::Generator s=0; s<i_tab.semisimple_rank(); ++s)
      i_tab.add_cross(s,i);

  assert(i_tab.size()==end()); // |size| new elements were found

} // |Cartan_orbit::Cartan_orbit|


// ------------------------------ Cartan_orbits -------------------------------


const CartanNbr undefined = ~0;

void Cartan_orbits::set_size(CartanNbr n_Cartans)
{ Cartan_index.resize(n_Cartans,undefined); orbit.reserve(n_Cartans); }

void Cartan_orbits::add(const InnerClass& G, CartanNbr cn)
{
  assert(cn<Cartan_index.size());
  if (Cartan_index[cn]!=undefined)
    return; // class was already added before, so nothing to do
  Cartan_index[cn]=orbit.size(); // if not, it will be added at this position

  // now actually generate the involutions associated to this Cartan class
  orbit.push_back(Cartan_orbit(static_cast<InvolutionTable&>(*this),G,cn));
}

unsigned int Cartan_orbits::locate(InvolutionNbr i) const
{
  unsigned int low=0, high=orbit.size();
  assert(high>0 and i<orbit.back().end());

  while (high-low>1)
  {
    CartanNbr mid=(low+high)/2;
    if (i<orbit[mid].start)
      high=mid;
    else
      low=mid;
  }
  assert(low<orbit.size());
  assert(orbit[low].contains(i));
  return low;
}

// whether |i<j|
bool Cartan_orbits::comparer::operator()
  (InvolutionNbr i, InvolutionNbr j) const
{
  int d = t.length(i)-t.length(j);
  if (d!=0 or (d = t.Weyl_length(i)-t.Weyl_length(j))!=0)
    return d<0;
  else
    return t.involution(i)<t.involution(j);
}

void InvolutionTable::reduce(TitsElt& a) const
{
  InvolutionNbr i=nr(a.tw());
  a.reduce(mod_space(i));
}


size_t Cartan_orbits::total_size(const BitMap& Cartan_classes) const
{
  size_t s=0;
  for (CartanNbr cn=0; cn<orbit.size(); ++cn)
    s+=operator[](cn).size;
  return s;
}

} // |namsepace|

} // |namespace atlas|
