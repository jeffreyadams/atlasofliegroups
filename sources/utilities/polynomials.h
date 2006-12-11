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

  const Degree MinusOne = ~0; // -1 as unsigned, the degree of Zero

}

/******** function declarations *********************************************/

namespace polynomials {

template<typename C>
  bool compare(const Polynomial<C>&, const Polynomial<C>&);

template<typename C>
void safeAdd(C&, C);

template<typename C>
void safeSubtract(C&, C);

}

/******** type definitions **************************************************/

namespace polynomials {

  /*!  \brief Polynomials with coefficients in C, which must be a
       modular integral type

       The coefficient type C must support addition (+,+=),
       multiplication (*=) and subtraction (unary and binary -),
       and comparison (<, ==, !=) (where < need not have mathematical sense)
  */
template<typename C> class Polynomial {

 private:

  std::vector<C> d_data;

  void adjustSize();

 public:

// constructors and destructors
  Polynomial() {} // zero polynomial Zero

  explicit Polynomial(Degree d); // build X^d. If d==~0 (==-1) gives Zero

// copy, assignment and swap
  void swap(Polynomial& other) {
    d_data.swap(other.d_data);
  }

// accessors

/* N.B. result was changed from 'C' to 'const C&', in order to be able to set
   const C* p_ptr= &p[0] when p is a constant polynomial */
  const C& operator[] (Degree j) const {  return d_data[j]; }

  bool operator== (const Polynomial& q) const { return d_data == q.d_data; }

  /*!
\brief Operator < is the default from the standard library < on vector.

  The ordering in the boolean compare(P,Q) is a much more useful and
  natural one. According to Fokko, the reason for using this unnatural
  comparison operator is probably that profiling suggested that the
  compare function was using a lot of time in Kazhdan-Lusztig
  computations.  (Each new KL polynomial must be inserted into the set
  of all polynomials, using a binary search like procedure, which uses '<')
  */
  bool operator< (const Polynomial& q) const { return d_data < q.d_data; }

  Degree degree() const { return d_data.size()-1; }

  bool isZero() const { return d_data.size() == 0; } // because of reduction

// manipulators
  C& operator[] (Degree j) { return d_data[j]; } // non-const version of above

  /* the names of the following methods are misnomers now, since there is no
     safety concern, but we prefer not to introduce useless changes to kl.cpp
  */

  void safeAdd(const Polynomial&, Degree d = 0, C c = C(1));

  void safeSubtract(const Polynomial& q, Degree d = 0, C c = C(1))
    { safeAdd(q,d,-c); }

  // and some variants using PolRef<C> instead of const Polynomial&

  void safeAdd(PolRef<C>, Degree d = 0, C c = C(1));

  void safeSubtract(PolRef<C> q, Degree d = 0, C c = C(1))
    { safeAdd(q,d,-c); }

};

/* the class PolRef<C> simulates 'const Polynomial<C>&' when no actual such
   polynomial is present in memory, a sequence containing all non-zero
   coefficients is
*/

template <typename C>
class PolRef
  {
    const C* p; // pointer to (virtual) constant term coefficient
    unsigned short int d_begin; // offset of first stored coefficient
    unsigned short int d_end; // degree+1, offset beyond last coefficient

    explicit PolRef(const Polynomial<C>& q) // easy: a polynomial is present
      : p(&q[0]), d_begin(0), d_end(q.degree()+1) {}
  public:
    PolRef(const C* q, unsigned short int begin, unsigned short int end)
      : p(q), d_begin(begin), d_end(end) {}

    // accessors
    bool isZero()  const { return d_begin>=d_end; } // no coefficients present

    Degree degree() const { return isZero() ? MinusOne : d_end-1; }

    Degree begin() const { return d_begin; } // lowest nonzero coefficient
    Degree end()   const { return d_end; }   // degree+1

    // selection supposes test i<end() has already been done
    C operator[] (size_t i) const { return i<d_begin ? C(0) : p[i]; }


    // when we really need a real Polynomial<C>, we shall construct one

    //    operator Polynomial<C> () const ; temporarily removed
    Polynomial<C> freeze ()const
      {
	// for efficiency reasons we wish to always 'return result'
	// when degree()==MinusOne, we shall have result=Zero
	Polynomial<C> result(degree());

        // now coeffs except last are C(0), no need to clear until d_begin

	for (Degree i=d_begin; i<d_end; ++i) // loop skipped if isZero()
	  result[i]=p[i];
	return result; // unique return allows optimised allocation of result
      }
  };

} // namespace polynomials

} // namespace atlas

#include "polynomials_def.h"

#endif
