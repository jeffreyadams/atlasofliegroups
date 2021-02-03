/*
  This is involutions.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef INVOLUTIONS_H  /* guard against multiple inclusions */
#define INVOLUTIONS_H

#include <vector>
#include <cassert>
#include <functional>

#include "../Atlas.h"

#include "hashtable.h"   // containment
#include "permutations.h"// containment root permutation in |InvolutionData|
#include "bitmap.h"      // containment root sets in |InvolutionData|

#include "weyl.h"        // containment of |WI_Entry|
#include "subquotient.h" // containment of |SmallSubspace|

#include "matrix.h"      // inlined matrix arithmetic

/* The purpose of this module is to provide a central registry of (twisted)
   involutions, in the form of a hash table to encode them by numbers, and
   supplementary information in the form of a table indexed by those numbers.
   This information, which includes root classification and the (somewhat
   voluminous) involution matrix, is generated as soon a an involution is
   registered here, which happens for whole twisted conjugacy classes at a
   time, through a call to |Cartan_orbits::add| defined below. If user code
   should need additional information associated involutions, it might define
   an array indexed by |InvolutionNbr|; but currently this happens nowhere.
 */

namespace atlas {

namespace involutions {


// this class gathers some information associated to a root datum involution
// the data depends only on the permutation of the |RootSystem| by |theta|
class InvolutionData
{
  Permutation root_perm; // permutation of all roots
  RootNbrSet d_imaginary, d_real, d_complex;
  RootNbrList d_simpleImaginary; // imaginary roots simple wrt subsystem
  RootNbrList d_simpleReal; // real roots simple wrt subsystem
 public:
  InvolutionData(const RootDatum& rd,
		 const WeightInvolution& theta);
  InvolutionData(const RootDatum& rd,
		 const WeightInvolution& theta,
		 const RootNbrSet& positive_subsystem);
  InvolutionData(const RootSystem& rs,
		 const RootNbrList& simple_images);
  static InvolutionData build(const RootSystem& rs,
			      const TwistedWeylGroup& W,
			      const TwistedInvolution& tw);
  void swap(InvolutionData&);
  // manipulators
private:
  void classify_roots(const RootSystem& rs); // workhorse for contructors
public:
  void cross_act(const Permutation& r_perm); // change (cheaply) to conjugate

  //accessors
  const Permutation& root_involution() const { return root_perm; }
  RootNbr root_involution(RootNbr alpha) const  { return root_perm[alpha]; }
  const RootNbrSet& imaginary_roots() const  { return d_imaginary; }
  const RootNbrSet& real_roots() const       { return d_real; }
  const RootNbrSet& complex_roots() const    { return d_complex; }
  size_t imaginary_rank() const { return d_simpleImaginary.size(); }
  const RootNbrList& imaginary_basis() const
    { return d_simpleImaginary; }
  RootNbr imaginary_basis(size_t i) const
    { return d_simpleImaginary[i]; }
  size_t real_rank() const { return d_simpleReal.size(); }
  const RootNbrList& real_basis() const { return d_simpleReal; }
  RootNbr real_basis(size_t i) const { return d_simpleReal[i]; }

}; // |class InvolutionData|

typedef unsigned int InvolutionNbr;

class InvolutionTable
{
 public: // there is no danger to expose these constant references
  const RootDatum& rd;
  const WeightInvolution& delta;
  const TwistedWeylGroup& tW;

 private:
  weyl::TI_Entry::Pooltype pool;
  HashTable<weyl::TI_Entry, InvolutionNbr> hash;

  struct record
  {
    InvolutionData id; // stuff that does not involve weight coordinates
    WeightInvolution theta;
    int_Matrix M_real; // $1-\theta$; then expression in a basis of its image
    int_Matrix lift_mat; // image basis $1-\theta$: |lift_mat*M_real==1-theta|
    unsigned int length;
    unsigned int W_length;
    SmallSubspace mod_space; // for |x|

  record(const WeightInvolution& inv,
	 const InvolutionData& inv_d,
	 const int_Matrix& Mre, const int_Matrix& lm,
	 unsigned int l,
	 unsigned int Wl,
	 const SmallSubspace& ms)
  : id(inv_d), theta(inv)
  , M_real(Mre), lift_mat(lm)
  , length(l), W_length(Wl), mod_space(ms) {}
  }; // |struct record|

  std::vector<record> data;
  std::vector<BinaryMap> torus_simple_reflection;

 public:
  InvolutionTable // constructor; starts without any involutions
    (const RootDatum& , const WeightInvolution&,  const TwistedWeylGroup&);

  //accessors
  size_t size() const { return pool.size(); }
  bool unseen(const TwistedInvolution& tw) const
  { return hash.find(weyl::TI_Entry(tw))==hash.empty; }
  InvolutionNbr nr(const TwistedInvolution& tw) const
  { return hash.find(weyl::TI_Entry(tw)); }

  unsigned int semisimple_rank() const { return tW.rank(); }

  const weyl::TI_Entry& involution(InvolutionNbr n) const
  { assert(n<size()); return pool[n]; }

  const WeightInvolution& matrix(InvolutionNbr n) const
  { assert(n<size()); return data[n].theta; }
  const WeightInvolution& matrix(const TwistedInvolution& tw) const
  { return matrix(nr(tw)); }

  unsigned int length(InvolutionNbr n) const
  { assert(n<size()); return data[n].length; }
  unsigned int Weyl_length(InvolutionNbr n) const
  { assert(n<size()); return data[n].W_length; }

  unsigned int length(const TwistedInvolution& tw) const
  { return length(nr(tw)); }
  unsigned int Weyl_length(const TwistedInvolution& tw) const
  { return Weyl_length(nr(tw)); }

  const Permutation& root_involution(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.root_involution(); }
  RootNbr root_involution(InvolutionNbr n,RootNbr alpha) const
  { assert(n<size()); return data[n].id.root_involution(alpha); }
  const RootNbrSet& imaginary_roots(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.imaginary_roots(); }
  const RootNbrSet& real_roots(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.real_roots(); }
  const RootNbrSet& complex_roots(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.complex_roots(); }

  size_t imaginary_rank(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.imaginary_rank(); }
  const RootNbrList& imaginary_basis(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.imaginary_basis(); }
  RootNbr imaginary_basis(InvolutionNbr n,weyl::Generator i) const
  { assert(n<size()); return data[n].id.imaginary_basis(i); }
  size_t real_rank(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.real_rank(); }
  const RootNbrList& real_basis(InvolutionNbr n) const
  { assert(n<size()); return data[n].id.real_basis(); }
  RootNbr real_basis(InvolutionNbr n,weyl::Generator i) const
  { assert(n<size()); return data[n].id.real_basis(i); }

  bool is_complex_simple(InvolutionNbr n,weyl::Generator s) const;
  bool is_imaginary_simple(InvolutionNbr n,weyl::Generator s) const;
  bool is_real_simple(InvolutionNbr n,weyl::Generator s) const;

  bool complex_is_descent(InvolutionNbr n,RootNbr alpha) const;

  void reduce(TitsElt& a) const;

  const SmallSubspace& mod_space(InvolutionNbr n) const
  { assert(n<size()); return data[n].mod_space; }

  KGB_elt_entry x_pack(const GlobalTitsElement& x) const; // for X only; slow
  bool x_equiv(const GlobalTitsElement& x0,const GlobalTitsElement& x1) const;

  // functionality for |y| values, represented as |TorusPart| (a small bitset)

  // uniquely chosen representative modulo $(1-\theta)X^*$ of $(1-\theta)/2*y$
  void real_unique(InvolutionNbr inv, RatWeight& y) const;

  const int_Matrix& to_1_theta_image_coordinates(InvolutionNbr inv) const
  { return data[inv].M_real; }
  const int_Matrix& theta_1_image_basis(InvolutionNbr inv) const
  { return data[inv].lift_mat; }

  // pack $(1-\theta)\lambda_rho$ into |TorusPart|: |lift_mat| coordinates mod 2
  TorusPart y_pack(InvolutionNbr inv, const Weight& lambda_rho) const
  { return TorusPart(data[inv].M_real * lambda_rho); }

  // find |(1-theta)*lam_rho| for some |lam_rho| with |ypack(i,lam_rho)=y_part|
  Weight y_lift(InvolutionNbr i, TorusPart y_part) const;

  // apply |delta| to |y_part| at |i0|, the result being at |i1==delta*i0*delta|
  // this is used to twist a parameter by |delta|, which affects its involution
  TorusPart y_act(InvolutionNbr i0, InvolutionNbr i1, // source, destination
		  TorusPart y_part, const WeightInvolution& delta) const
  { return y_unlift(i1,delta*y_lift(i0,y_part)); }


  // the following produces a light-weight function object calling |involution|
  class mapper
  // : public std::unary_function<InvolutionNbr,const weyl::TI_Entry& >
  {
    const InvolutionTable& t;
  public:
  mapper(const InvolutionTable* tab) : t(*tab) {}
    const weyl::TI_Entry& operator() (InvolutionNbr n) const
    { return t.involution(n); }
  };
  mapper as_map() const { return mapper(this); }

  // manipulators

  // these methods construct/propagate information at individual involutions
  InvolutionNbr add_involution(const TwistedInvolution& tw);
  InvolutionNbr add_cross(weyl::Generator s, InvolutionNbr n);

  void reserve(size_t s) { pool.reserve(s); }

 private: // auxiliary for |y_act| above
  // like |y_pack|, but direct |lift| left-inverse: |y_unlift(i,y_lift(i,y))==y|
  // effectively do |y_pack(i,lifted/2)|, but avoid half-integer coordinates
  TorusPart y_unlift(InvolutionNbr inv, const Weight& lifted) const
  { return TorusPart // contructor reduces coordinates modulo 2
      ((data[inv].M_real*lifted) / 2); // division precedes that reduction
  }



}; // |class InvolutionTable|

struct Cartan_orbit
{
  unsigned int Cartan_class_nbr;
  InvolutionNbr start,size;

  Cartan_orbit(InvolutionTable& i_tab,const InnerClass& G, CartanNbr cn);

  bool contains(InvolutionNbr i) const { return i-start<size; }
  InvolutionNbr end() const { return start+size; }

}; // |struct Cartan_orbit|


// we organize everything by Cartan classes of involutions

class Cartan_orbits : public InvolutionTable
{
  std::vector<Cartan_orbit> orbit; // orbits of Cartan classes that were added
  std::vector<unsigned int> Cartan_index; // maps Cartan number to its position

  unsigned int locate(InvolutionNbr i) const;
public:
  Cartan_orbits (const RootDatum& rd, const WeightInvolution& theta,
		 const TwistedWeylGroup& tW)
  : InvolutionTable(rd,theta,tW), orbit(), Cartan_index() { }

// manipulators

  void set_size(CartanNbr n_Cartans); // resize once number of Cartans is known
  void add(const InnerClass& G, CartanNbr cn);

// accessors

  // this method allows mostly finding the range of a given Cartan class
  const Cartan_orbit& operator[](CartanNbr cn) const
  { assert(Cartan_index[cn]!=CartanNbr(~0u)); return orbit[Cartan_index[cn]]; }

  CartanNbr Cartan_class(InvolutionNbr i) const
  { return orbit[locate(i)].Cartan_class_nbr; }
  CartanNbr Cartan_class(const TwistedInvolution& tw) const
  { return Cartan_class(nr(tw)); }

  size_t total_size(const BitMap& Cartan_classes) const;

 // make class useful as comparison object, for |std::sort|
  class comparer
  {
    const Cartan_orbits& t;
  public:
  comparer(const Cartan_orbits* o) : t(*o) {}
    // whether involution |i| less than |j| by length; Weyl length; |i<j|
    bool operator() (InvolutionNbr i, InvolutionNbr j) const;
  };
  comparer less() const { return comparer(this); }

}; // |class Cartan_orbits|




} // |namespace involutions|

} // |namespace atlas|

#endif
