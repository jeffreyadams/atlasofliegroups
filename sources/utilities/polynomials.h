/*!
\file
\brief Class definitions and function declarations for the class Polynomial.
*/
/*
  This is polynomials.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef POLYNOMIALS_H  /* guard against multiple inclusions */
#define POLYNOMIALS_H

#include "polynomials_fwd.h"

#include <limits>
#include <vector>

namespace atlas {

/******** constant declarations *********************************************/

namespace polynomials {

  const Degree MinusOne = ~ Degree(0); // -1 as unsigned, the degree of Zero

}

/******** function declarations *********************************************/

namespace polynomials {

template<typename C>
  bool compare(const Polynomial<C>&, const Polynomial<C>&);

template<typename C>
void safeAdd(C&, C);

template<typename C>
void safeProd(C&, C);

template<typename C>
void safeSubtract(C&, C);

}

/******** type definitions **************************************************/

namespace polynomials {

  /*!  \brief Polynomials with coefficients in C, which is expected to be a
       standard unsigned type.

       The coefficient type C must support addition (+, +=), multiplication
       (*, *=), subtraction (unary and binary -, -=), and
       std::numeric_limits<C> (used to test for overflow in the safeAdd
       operation); moreover comparisons (<, ==, !=) should be defined,
       although < need not have any particular mathematical sense
  */
template<typename C> class Polynomial {

 private:

  std::vector<C> d_data;

  void adjustSize();

 public:

// constructors and destructors
  Polynomial() {} // zero polynomial Zero

  explicit Polynomial(Degree d);

// copy, assignment and swap
  void swap(Polynomial& other) {
    d_data.swap(other.d_data);
  }

// accessors
  C operator[] (Degree j) const {
    return d_data[j];
  }

  bool operator== (const Polynomial& q) const {
    return d_data == q.d_data;
  }

  bool operator != (const Polynomial& q) const {
    return d_data != q.d_data;
  }


  /*!
\brief Operator < is the default from the standard library < on vector.

  The comparison operation below is only defined in order to allow ordered
  data types containing polynomials, such as |std::set<Polynomial<int> >|.
  Currently no such types are used in the Atlas (but initially they were).
  */
  bool operator< (const Polynomial& q) const { return d_data < q.d_data; }

  Degree degree() const { return d_data.size()-1; }

  bool isZero() const { return d_data.size() == 0; } // because of reduction

// manipulators
  C& operator[] (Degree j) { return d_data[j]; } // non-const version of above

  Polynomial& operator+= (const Polynomial& q);

  Polynomial& operator-= (const Polynomial& q);

  Polynomial& subtract_from (const Polynomial& p); // *this = p - *this

  Polynomial& operator*= (C);

  Polynomial operator* (const Polynomial& q) const;

  void safeAdd(const Polynomial& p, Degree d, C c); // *this += c*q^d*p
  void safeAdd(const Polynomial& p, Degree d = 0);  // *this += q^d*p

  void safeSubtract(const Polynomial& p, Degree d, C c);
  void safeSubtract(const Polynomial& p, Degree d = 0 );
};

}

}

#include "polynomials_def.h"

#endif
