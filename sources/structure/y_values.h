/*
  This is y_values.h, for representing |y| components of module parameters

  Copyright (C) 2011-2014 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef Y_VALUES_H  /* guard against multiple inclusions */
#define Y_VALUES_H

#include <cassert>

#include "../Atlas.h"

#include "bitvector.h" // use of |SmallBitVector| constructor
#include "ratvec.h"   // containment

namespace atlas {

namespace y_values {


inline TorusElement exp_pi(const RatWeight& r);
inline TorusElement exp_2pi(const RatWeight& r);


/* The following classes are newer than |TitsElement|, |TitsGroup| and
   |TitsCoset|, and less efficient, but provide a framework in which te sets
   modelled by the latter classes can be fitted. They are actually used mostly
   for handling |y| values (whence the name of this module), in which case one
   has elements of the torus of the dual group. The terminology is adapted to
   this point of view: we store a rational weight (rather than coweight).
 */


/* An element of finite order in $H$. This is like a rational vector in the
   coordinates, taken modulo integers. But we store twice that value (among
   other things this facilitates adding elements of $H(2)$), so out bijection
   is |exp_pi| rather than |exp_2pi|, and coordinates are taken modulo $2\Z$.
 */
class TorusElement
{
  RatWeight repr; // represents $\exp(\pi i repr)$, no factor 2!

  TorusElement(const RatWeight& r,tags::UnnormalizedTag)
   : repr(r) {} // raw constructor, like exp(pi i r), but no modular reduction
 public:
  explicit TorusElement(size_t rank) : repr(rank) {} // identity torus element

  // $\exp(\pi\ii r )$ or $\exp(2\pi\ii r )$
  TorusElement(const RatWeight& r,bool two);

  // accessors
  size_t rank () const { return repr.size(); }

  RatWeight log_pi(bool normalize) const; // return the stored rational vector
  RatWeight log_2pi() const; // value halved: to be interpreted "mod Z^rank"

  // more often it will be practical to have access to that "mod 2Z^rank" form
  const RatWeight& as_Qmod2Z() const { return repr; }

  bool operator== (const TorusElement& a) const { return repr==a.repr; }
  bool operator!= (const TorusElement& a) const { return repr!=a.repr; }
  bool operator<  (const TorusElement& a) const { return repr<a.repr; }

  TorusElement operator +(const TorusElement& t) const;
  TorusElement operator -(const TorusElement& t) const;

  TorusElement& operator +=(const TorusElement& t)
  { return *this= operator +(t); }
  TorusElement& operator -=(const TorusElement& t)
  { return *this= operator -(t); }

  // this method is to be used only at weights |alpha| taking value +1 or -1
  bool negative_at(const Coweight& alpha) const
  {
    // the following asserts internally that |evaluate_at(alpha)| is integer
    return repr.dot(alpha)%2!=0; // true if evaluates to odd integer
  }

  // evaluation giving rational number modulo 2, represented in interval [0,2)
  Rational evaluate_at(const SmallBitVector& alpha) const;

  // evaluation giving rational number modulo 2, represented in interval [0,2)
  Rational evaluate_at(const Coweight& alpha) const;

// a method for rapidly doing imaginary cross action (for completing fiber)
// requires simple-imaginary |alpha| that is integral on |t| (however unless
// all roots are integral on |t|, these conditions might be contradictory!)
  TorusElement simple_imaginary_cross
    (const RootDatum& dual_rd, // dual for pragmatic reasons
     RootNbr alpha) const; // any simple-imaginary root

  // manipulators

  TorusElement& operator+=(TorusPart v); // arg by value since it is small
  TorusElement& reduce(); // reduce entires mod $2\Z$

  void simple_reflect(const PreRootDatum& prd, weyl::Generator s);
  void reflect(const RootDatum& rd, RootNbr alpha);
  void act_by(const  WeightInvolution& delta);

  TorusElement& left_symmetrise(const WeightInvolution& delta)
  { symmetrise(delta,repr); return this->reduce(); }
  TorusElement& right_symmetrise(const WeightInvolution& delta)
  { symmetrise(repr,delta); return this->reduce(); }
}; // |class TorusElement|

inline TorusElement exp_pi(const RatWeight& r) { return TorusElement(r,false); }
inline TorusElement exp_2pi(const RatWeight& r) { return TorusElement(r,true); }

// a |y| value can be hashed as (fingerprint,InvolutionNbr) pair
struct y_entry
{
  TorusElement t_rep; // a representative torus element, ignored in test
  InvolutionNbr nr;
  RatWeight fingerprint; // charcterizes the torus element

  // obligatory fields for hashable entry
  typedef std::vector<y_entry> Pooltype;
  size_t hashCode(size_t modulus) const; // hash function
  bool operator !=(const y_entry& y) const; // tests |nr| and |fingerprint|

  y_entry (const RatWeight& f, InvolutionNbr i, const TorusElement& t)
  : t_rep(t), nr(i), fingerprint(f) {}
  const TorusElement& repr() const { return t_rep; }

}; //  |struct y_entry|

// maps from $\gamma-\lambda$ shift to integral-only $y$ bits and back
// coding modulo sum of $\Q$-orthogonal of integral coroots and $(1-\theta)X^*$
class y_codec
{ // below $r$ is the rank (of $X^*$) and $d$ the number of relevant $y$ bits
  int_Matrix coder; // size $d\times r$
  int_Matrix decoder; // size $r\times d$, lift $y$ to $\gamma-\lambda$ shift

 public:
  y_codec(const InnerClass& ic,
	  InvolutionNbr inv,
	  const CoweightList& integrally_simple_coroots);

  // methods to inherit coing information to a diffent involution
  y_codec cross(weyl::Generator s, const y_codec& yc);
  y_codec Cayley_descent(weyl::Generator s, const y_codec& yc);

  SmallBitVector encode(const RatWeight& gamma_lambda_shift)
  { const RatWeight prod = (coder*gamma_lambda_shift).normalize();
    assert(prod.denominator()==1);
    return SmallBitVector { prod.numerator() };
  }
  Weight decode (SmallBitVector y_bits)
  { Weight result(decoder.numRows(),0);
    for (auto it = y_bits.data().begin(); it(); ++it)
      result += decoder.column(*it);
    return result;
  }
}; // |y_codec|


//				*** Functions ***

// whether each simple root has integral evaluation on a torus element
bool is_central(const LatticeMatrix& simple_roots, const TorusElement& e);

// For a $\xi$-stable torus element, find $\xi$-stable pre-image by $\exp_1$
RatCoweight stable_log(const TorusElement& t, CoweightInvolution xi);

} // |namespace y_values|

} // |namespace atlas|

#endif
