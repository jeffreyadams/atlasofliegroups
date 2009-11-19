/*!
\file
\brief Class definitions and function declarations for the class Polynomial.
*/
/*
  This is polynomials.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef POLYNOMIALS_H  /* guard against multiple inclusions */
#define POLYNOMIALS_H

#include "polynomials_fwd.h"

#include <limits>
#include <vector>

namespace atlas {

namespace polynomials {

/******** constant declarations *********************************************/


  const Degree MinusOne = ~ Degree(0); // -1 as unsigned, the degree of Zero



/******** function declarations *********************************************/


template<typename C>
  bool compare(const Polynomial<C>&, const Polynomial<C>&); // lex ordering

// basic operations which test for overflow/underflow, assuming |C| is unsigned
template<typename C>
  void safeAdd(C&, C); // version of |operator+=| testing for overflow

template<typename C>
  void safeProd(C&, C); // version of |operator*=| testing for overflow

template<typename C>
void safeSubtract(C&, C);// version of |operator-=| testing for underflow


/******** type definitions **************************************************/


/*! \brief Polynomials with coefficients in |C|

  The coefficient type |C| could be a signed or unsigned integral type, or any
  type providing ring operations (operator+, operator*= etc.), for instance
  modular integers. Moreover comparisons (<, ==, !=) should be defined for STL
  use, although < need not have any particular mathematical sense
*/
template<typename C> class Polynomial
{
  std::vector<C> d_data;

 public:

// constructors and destructors
  Polynomial() : d_data() {} // zero polynomial
  explicit Polynomial(C c);
  Polynomial(Degree d, C c); // initialised to $cX^d$ (with |c!=0|)

// copy, assignment (both defaults will do) and swap
  void swap(Polynomial& other) { d_data.swap(other.d_data); }

// accessors
  C operator[] (Degree i) const { return d_data[i]; } // get coefficient $X^i$

  bool operator== (const Polynomial& q) const { return d_data == q.d_data; }
  bool operator!= (const Polynomial& q) const { return d_data != q.d_data; }

/*!
\brief Operator < is the default from the standard library < on vector.

  The comparison operation below is only defined in order to allow ordered
  data types containing polynomials, such as |std::set<Polynomial<int> >|.
  Currently no such types are used in the Atlas (but initially they were).
*/
  bool operator< (const Polynomial& q) const { return d_data < q.d_data; }

  Degree degree() const { return d_data.size()-1; }
  Degree size() const { return d_data.size(); }

  bool isZero() const { return size() == 0; } // because of reduction

// manipulators
  C& operator[] (Degree j) { return d_data[j]; } // non-const version of above

  Polynomial& operator+= (const Polynomial& q);

  Polynomial& operator-= (const Polynomial& q);

  Polynomial& subtract_from (const Polynomial& p); // *this = p - *this

  Polynomial& operator*= (C);
  Polynomial operator* (C c) const { return Polynomial (*this)*=c; }

  Polynomial operator* (const Polynomial& q) const;
  Polynomial operator+ (const Polynomial& q) const
    { return Polynomial(*this)+=q; }
  Polynomial operator- (const Polynomial& q) const
    { return Polynomial(*this)-=q; }

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
  explicit Safe_Poly(Degree d, C c) : base(d,c) {}

  // unlike |operator+| etc., the following test for negative coefficients
  void safeAdd(const Safe_Poly& p, Degree d, C c); // *this += c*q^d*p
  void safeAdd(const Safe_Poly& p, Degree d = 0);  // *this += q^d*p

  void safeSubtract(const Safe_Poly& p, Degree d, C c);
  void safeSubtract(const Safe_Poly& p, Degree d = 0 );

}; // |template <typename C> class Safe_Poly|

} // |namespace polynomials|

} // |namespace atlas|

#include "polynomials_def.h"

#endif
