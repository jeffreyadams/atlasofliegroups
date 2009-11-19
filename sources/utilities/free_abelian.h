/*!
\file
  This is free_abelian.h
*/
/*
  Copyright (C) 2008 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <map>

#ifndef FREE_ABELIAN_H  /* guard against multiple inclusions */
#define FREE_ABELIAN_H

namespace atlas {


/******** type declarations **************************************************/


/******** function declarations **********************************************/

namespace free_abelian {

}

/******** type definitions ***************************************************/

namespace free_abelian {

// A class that represents an element of the free abelian group on the set T
// in fact, taking |C=polynomials::polynomial<int>| gives $q$-multiplicities
template<typename T, typename C=long int, typename Compare=std::less<T> >
struct Free_Abelian : public std::map<T,C,Compare>
{
  typedef C coef_t;
  typedef std::map<T,coef_t,Compare> base;
  typedef typename base::iterator iterator;
  typedef typename base::const_iterator const_iterator;

  Free_Abelian() : base(Compare()) {}

  Free_Abelian(Compare c) : base(c) {}

  explicit Free_Abelian(const base& m) : base(m) {}

  explicit Free_Abelian(T p, Compare c=Compare()) : base(c)
  { base::insert(std::make_pair(p,1L)); }

  Free_Abelian(T p,C m, Compare c=Compare()) : base(c)
  { base::insert(std::make_pair(p,m)); }

  template<typename InputIterator>
  Free_Abelian(InputIterator first, InputIterator last, Compare c=Compare())
    : base(first,last,c) {}

  Free_Abelian& add_multiple(const Free_Abelian& p, C m);

  Free_Abelian& operator+=(const Free_Abelian& p)
  { return add_multiple(p,C(1)); }

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

  Monoid_Ring& add_multiple(const Monoid_Ring& p, C m,const T& expon);

}; // |class Monoid_Ring|


} // namespace free_abelian

} // namespace atlas

#include "free_abelian_def.h"

#endif
