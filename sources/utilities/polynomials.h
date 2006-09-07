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

  const Degree UndefDegree = ~0ul;

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

  /*!  \brief Polynomials with coefficients in C, which must be a
       standard unsigned type.

       The coefficient type C must support addition, multiplication,
       subtraction, and std::numeric_limits<C> (used to test for
       overflow in the safeAdd operation).
  */
template<typename C> class Polynomial {

 private:

  std::vector<C> d_data;

  void adjustSize();

 public:

// constructors and destructors
  Polynomial() {}

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

  /*!
\brief Operator < is the default from the standard library < on vector.

The ordering in the boolean compare(P,Q) is a much more useful and
  natural one. According to Fokko, the reason for using this unnatural
  comparison operator is probably that profiling suggested that the
  compare function was using a lot of time in Kazhdan-Lusztig
  computations.  (Each new KL polynomial must be inserted into the set
  of 
  */
bool operator< (const Polynomial& q) const { return d_data < q.d_data;
  }

  Degree degree() const {
    return d_data.size()-1;
  }

  bool isZero() const {
    return d_data.size() == 0;
  }

// manipulators
  C& operator[] (Degree j) {
    return d_data[j];
  }

  Polynomial& operator+= (const Polynomial&);

  Polynomial& operator-= (const Polynomial&);

  Polynomial& operator*= (C);

  Polynomial& safeAdd(const Polynomial&, Degree d = 0, C c = 1);

  Polynomial& safeSubtract(const Polynomial&, Degree d = 0, C c = 1);
};

}

}

#include "polynomials_def.h"

#endif
