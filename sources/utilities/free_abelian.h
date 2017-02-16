/*
  This is free_abelian.h

  Copyright (C) 2008 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include <map>

#ifndef FREE_ABELIAN_H  /* guard against multiple inclusions */
#define FREE_ABELIAN_H

#include "free_abelian_fwd.h"

namespace atlas {

namespace free_abelian {


/******** type declarations **************************************************/


/******** function declarations **********************************************/

/******** type definitions ***************************************************/

// A class that represents an element of the free abelian group on the set T
// or, with |C=polynomials::polynomial<int>| counting with $q$-multiplicities
// |T| value type, |C| coefficient type, |Compare| equivalence test for |T|
template<typename T, typename C, typename Compare>
struct Free_Abelian : public std::map<T,C,Compare>
{
  typedef C coef_t;
  typedef std::map<T,coef_t,Compare> base; // the (base) reresentation type
  typedef typename base::iterator iterator;
  typedef typename base::const_iterator const_iterator;

  Free_Abelian() : base(Compare()) {} // default |Compare| value for base

  Free_Abelian(Compare c) : base(c) {} // here a specific |Compare| is used

  explicit Free_Abelian(const base& m) : base(m) {} // promote base to derived

  explicit Free_Abelian(const T& p, Compare c=Compare()) // create a monomial
  : base(c)
  { base::insert(std::make_pair(p,coef_t(1L))); }

  Free_Abelian(const T& p,C m, Compare c=Compare()) // mononomial (single term)
  : base(c)
  { if (m!=C(0))
      base::insert(std::make_pair(p,m));
  }

  // convert other agregate of (monomial,coefficient) pairs to |Free_Abelian|
  // warning: current implementation allows coefficients |C(0)| to slip through
  template<typename InputIterator> // iterator over (T,coef_t) pairs
  Free_Abelian(InputIterator first, InputIterator last, Compare c=Compare())
    : base(first,last,c) {}

  Free_Abelian& add_term(const T& p, C m);
  Free_Abelian& operator+=(const T& p) { return add_term(p,C(1)); }
  Free_Abelian& operator-=(const T& p) { return add_term(p,C(-1)); }

  Free_Abelian& add_multiple(const Free_Abelian& p, C m);

  Free_Abelian& operator+=(const Free_Abelian& p)
  { if (base::empty())
      return *this =(p); // assign, avoiding work on initial addition to empty
    return add_multiple(p,C(1));
  }

  Free_Abelian& operator-=(const Free_Abelian& p)
  { return add_multiple(p,C(-1)); }

  C operator[] (const T& t) const // find coefficient of |t| in |*this|
  {
    typename base::const_iterator p=base::find(t);
    return p==base::end() ? C(0) : p->second;
  }

}; // |class Free_Abelian|

/* When we also want a multiplication, |T| must have operator+=. Rather than
   defining operator*= in Free_Abelian (possibly unused), we prefer to derive.
*/
template<typename T, typename C=long int, typename Compare=std::less<T> >
  struct Monoid_Ring : public Free_Abelian<T,C,Compare>
{
  typedef C coef_t;
  typedef Free_Abelian<T,C,Compare> base;
  typedef typename base::iterator iterator;
  typedef typename base::const_iterator const_iterator;

  Monoid_Ring() : base(Compare()) {}
  Monoid_Ring(Compare c) : base(c) {}

  explicit Monoid_Ring(const typename base::base& m) : base(m) {}
  explicit Monoid_Ring(T p, Compare c=Compare()) : base(p,c) {}
  Monoid_Ring(T p,C m, Compare c=Compare()) : base(p,m,c) {}

  template<typename InputIterator>
  Monoid_Ring(InputIterator first, InputIterator last, Compare c=Compare())
    : base(first,last,c) {}

  Monoid_Ring operator*(const Monoid_Ring& p);

  // add to |*this| a multiple of |p| by the mono-nomial $cX^{\\{expon}}$
  Monoid_Ring& add_multiple(const Monoid_Ring& p, C m,const T& expon);

}; // |class Monoid_Ring|


} // |namespace free_abelian|

} // |namespace atlas|

#include "free_abelian_def.h"

#endif
