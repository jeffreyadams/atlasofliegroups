/*!
\file
  This is polynomials.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

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

  bool operator< (const Polynomial& q) const {
    return d_data < q.d_data;
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
