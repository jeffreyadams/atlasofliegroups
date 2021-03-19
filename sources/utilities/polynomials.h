/*
  This is polynomials.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2016,2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/* Class definitions and function declarations for the class Polynomial. */

#ifndef POLYNOMIALS_H  /* guard against multiple inclusions */
#define POLYNOMIALS_H

#include "polynomials_fwd.h"

#include <limits>
#include <vector>
#include <iostream>

namespace atlas {

namespace polynomials {

/******** constant declarations *********************************************/


  const Degree MinusOne = ~ Degree(0); // -1 as unsigned, the degree of Zero



/******** function declarations *********************************************/


template<typename C>
  bool compare(const Polynomial<C>&, const Polynomial<C>&); // lex ordering

// the following operations will throw on overflow, but not on other errors
// underflow and division by zero trigger failing assertions while debugging
// the rationale is that only overflow should be plausible in production code

// basic operations which test for overflow/underflow, assuming |C| is unsigned
template<typename C>
  void safe_add(C&, C); // version of |operator+=| testing for overflow

template<typename C>
void safe_subtract(C&, C);// version of |operator-=| asserting no underflow

template<typename C>
  void safe_multiply(C&, C); // version of |operator*=| testing for overflow

template<typename C>
  void safe_divide(C&, C); // version of |operator/=| asserting validity


/******** type definitions **************************************************/


/*
  Polynomials with coefficients in |C|.

  The coefficient type |C| could be a signed or unsigned integral type, or any
  type providing ring operations (operator+, operator*= etc.), for instance
  modular integers. Moreover comparisons (<, ==, !=) should be defined for STL
  use, although < need not have any particular mathematical sense
*/
template<typename C> class Polynomial
{ // to some extent we do as if the type is derived from |std::vector<C>|
  using data_type = std::vector<C>;
  data_type d_data; // as this unique data member

 public:

// constructors and destructors
  Polynomial() : d_data() {} // zero polynomial
  explicit Polynomial(C c);  // constant polynomial, |c==0| handled correctly
  Polynomial(Degree d, C c); // initialised to $cX^d$ (with |c!=0|)
  Polynomial(Degree d,const Polynomial& Q); // initialised to $X^d Q$

// copy and move constructors are implicitly declared, their definition defaulted

template <typename U>
  Polynomial(const Polynomial<U>& src) : d_data(src.begin(),src.end()) { }

// accessors
  const C& operator[] (Degree i) const; // get coefficient $X^i$
  const C coef (Degree i) const { return i>=d_data.size() ? C(0) : d_data[i]; }

  bool operator== (const Polynomial& q) const { return d_data == q.d_data; }
  bool operator!= (const Polynomial& q) const { return d_data != q.d_data; }

/*
  Operator < is the default from the standard library < on vector.

  The comparison operation below is only defined in order to allow ordered
  data types containing polynomials, such as |std::set<Polynomial<int> >|.
  Currently no such types are used in the Atlas (but initially they were).
*/
  bool operator< (const Polynomial& q) const { return d_data < q.d_data; }

  // iterator methods |begin| and |end| are as if we derived from |d_data|
  typename data_type::const_iterator begin() const
  { return d_data.begin();}
  typename data_type::const_iterator end() const { return d_data.end();}

  Degree degree() const { return d_data.size()-1; }
  Degree size() const { return d_data.size(); }

  // because of reduction (automatic removal leading zeros), this is easy:
  bool isZero() const { return size() == 0; }

  // |degree()| being unsigned but possibly -1, compare using this safe method:
  bool degree_less_than (Degree d) const { return size()<=d; }

  bool multi_term () const; // whether more than one term is nonzero (printing)

  // write polynomial as $(1+cX)Q+rX^d$, and return $r$
  C up_remainder(C c, Degree d) const; // assumes $d\geq degree()$

// manipulators
  C& operator[] (Degree j); // non-const version of above
  typename data_type::iterator begin()  { return d_data.begin();}
  typename data_type::iterator end()    { return d_data.end();}
  data_type&& data() && { return std::move(d_data); } // pilfer contents

  Polynomial& operator+= (const Polynomial& q);

  Polynomial& operator-= (const Polynomial& q);

  Polynomial& subtract_from (const Polynomial& p); // *this = p - *this

  Polynomial& operator*= (C);
  Polynomial& operator/= (C);
  Polynomial operator* (C c) const { return Polynomial (*this)*=c; }

  Polynomial operator* (const Polynomial& q) const;
  Polynomial& operator*= (const Polynomial& q)
  { (*this)=operator*(q); return *this; } // move-assign from |operator*| value
  Polynomial operator+ (const Polynomial& q) const
    { // avoid requiring reallocation during addition
      return d_data.size()>=q.d_data.size()
	? Polynomial(*this)+=q : Polynomial(q)+=*this;
    }
  Polynomial operator- (const Polynomial& q) const
    { // avoid requiring reallocation during addition
      return d_data.size()>=q.d_data.size()
	? Polynomial(*this)-=q : (Polynomial(q)*=C(-1))+=*this;
    }
  Polynomial operator- () const { return Polynomial(*this)*= C(-1); }

  // as |up_remainder| above, but also change polynomial into quotient $Q$
  // the coefficient in degree |d| may leave a remainder which is returned
  C factor_by_1_plus_q(Degree d); // with precondition |this->degree()<=d|

  // similarly upward divide by $1+q^k$, return remainder in degree |d|
  C factor_by_1_plus_q_to_the(Degree k,Degree d);

  std::ostream& print(std::ostream& strm, const char* x) const;

protected:
  void resize (Degree d) { d_data.resize(d,C(0)); }
  void adjustSize(); // shrink |d_data| to make leading coefficient nonzero

 }; // |template<typename C> class Polynomial|

/*
  The following class template assumes |C| is an _unsigned_ integer type,
  for which |std::numeric_limits<C>| is defined; it then provides safe
  versions of shifted additive operations
 */
template <typename C>
  class Safe_Poly : public Polynomial<C>
{
  typedef Polynomial<C> base;
 public:
  Safe_Poly() : base() {} // zero polynomial
  explicit Safe_Poly(C c) : base(c) {}
  Safe_Poly(Degree d, C c) : base(d,c) {}

  // unlike |operator+| etc., the following test for negative coefficients
  void safeAdd(const Safe_Poly& p, Degree d, C c); // *this += c*q^d*p
  void safeAdd(const Safe_Poly& p, Degree d = 0);  // *this += q^d*p
  void safeDivide(C c);  // *this = *this/c
  void safe_quotient_by_1_plus_q(Degree delta);  // *this = (*this + mq^d)/(q+1)

  void safeSubtract(const Safe_Poly& p, Degree d, C c);
  void safeSubtract(const Safe_Poly& p, Degree d = 0 );

}; // |template <typename C> class Safe_Poly|

/* For somewhat unclear reasons, both |PolEntry| and |SafePolEntry| templates
   use as |Pooltype| a |std::vector| of a base type rather than a vector of
   themselves. Since the base type has no |hashCode| member, we define hash
   functions for them as free functions |hash_code|, that will override the more
   generally templated one provided in hastable.h
*/

template <typename U>
  size_t hash_code(const Safe_Poly<U>& P, size_t modulus)
{ if (P.isZero()) return 0;
  Degree i = P.degree();
  size_t h = P[i];
  while (i-->0)
    h= (h<<21)+(h<<13)+(h<<8)+(h<<5)+h+P[i];
  return h & (modulus-1);
}

template <typename C>
  size_t hash_code(const Polynomial<C>& P, size_t modulus)
{ if (P.isZero()) return 0;
  Degree i = P.degree();
  size_t h = P[i]; // start with leading coefficient, converted to unsigned
  while (i-->0)
    h= (h<<21)+(h<<13)+(h<<8)+(h<<5)+h+P[i];
  return h & (modulus-1);
}

template <typename C> class PolEntry : public Polynomial<C>
{
  using Pol = Polynomial<C>;
public:
  // constructors
  PolEntry() : Pol() {} // default constructor builds zero polynomial
  PolEntry(const Pol& p) : Pol(p) {} // lift polynomial to this class

  // members required for an Entry parameter to the HashTable template
  using Pooltype = std::vector<Pol>;  // associated storage type
  size_t hashCode(size_t modulus) const // hash function; use the one above
  { return hash_code(static_cast<const Pol&>(*this),modulus); }

  // compare polynomial with one from storage
  bool operator!=(typename Pooltype::const_reference e) const
  { const Pol& P=*this;
    if (P.degree()!=e.degree()) return true;
    if (P.isZero()) return false; // since degrees match
    for (Degree i=0; i<=P.degree(); ++i)
      if (P[i]!=e[i]) return true;
    return false; // no difference found
  }

}; // |template <typename C> class PolEntry|

template <typename U> class SafePolEntry : public Safe_Poly<U>
{
  using SafePol = Safe_Poly<U>;
public:
  // constructors
  SafePolEntry() : SafePol() {} // default constructor builds zero polynomial
  SafePolEntry(const SafePol& p) : SafePol(p) {} // lift polynomial to this class

  SafePolEntry(SafePolEntry&&) = default;

  // members required for an Entry parameter to the HashTable template
  using Pooltype = std::vector<SafePol>;  // associated storage type
  size_t hashCode(size_t modulus) const // hash function; use the one above
  { return hash_code(static_cast<const SafePol&>(*this),modulus); }

  // compare polynomial with one from storage
  bool operator!=(typename Pooltype::const_reference e) const
  { const SafePol& P=*this;
    if (P.degree()!=e.degree()) return true;
    if (P.isZero()) return false; // since degrees match
    for (Degree i=0; i<=P.degree(); ++i)
      if (P[i]!=e[i]) return true;
    return false; // no difference found
  }

}; // |class KLPolEntry|

} // |namespace polynomials|

} // |namespace atlas|

#include "polynomials_def.h"

#endif
