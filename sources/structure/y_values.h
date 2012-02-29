/*
  This is involurions.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef Y_VALUES_H  /* guard against multiple inclusions */
#define Y_VALUES_H

#include "atlas_types.h"

#include "ratvec.h"   // containment

namespace atlas {

namespace y_values {


inline TorusElement exp_pi(const RatWeight& r);
inline TorusElement exp_2pi(const RatWeight& r);


/* The following classes are newer than |TitsElement|, |TitsGroup| and
   |TitsCoset|, and less efficient, but provide a framework in which te sets
   modelled by the latter classes can be fitted.
*/

/* An element of finite order in $H$. Externally this is like a rational vector
   modulo integers in the coordinates. Internally we store twice the numerator
   of the rational vector, which facilitates adding halves to its coordinates.
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

  RatWeight log_pi(bool normalize) const; // return the stored rational vector
  RatWeight log_2pi() const; // value halved: to be interpreted "mod Z^rank"

  // more often it will be practical to have acces to that "mod 2Z^rank" form
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
    { return repr.scalarProduct(alpha)%2!=0; }

  // evaluation giving rational number modulo 2 (|negative_at| iff equals 1)
  Rational evaluate_at(const Coweight& alpha) const;

// a method for rapidly doing imaginary cross action (for completing fiber)
// requires simple-imaginary |alpha| that is integral on |t| (however unless
// all roots are integral on |t|, these conditions might be contradictory!)
  TorusElement simple_imaginary_cross
    (const RootDatum& dual_rd, // dual for pragmatic reasons
     RootNbr alpha) const; // any simple-imaginary root

  // manipulators

  TorusElement& operator+=(TorusPart v); // arg by value since it is small

  // the following method assumes |prd| is on dual side with respect to torus
  void simple_reflect(const PreRootDatum& prd, weyl::Generator s);
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
  TorusElement repr() const { return t_rep; }

}; //  |struct y_entry|


} // |namespace y_values|

} // |namespace atlas|

#endif
