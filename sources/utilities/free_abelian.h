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

// A class that represent an element of the free abelian group on the set T
template<typename T, typename Compare=std::less<T> >
struct Free_Abelian : public std::map<T,long int,Compare>
{
  typedef long int coef_t;
  typedef std::map<T,coef_t,Compare> base;

  Free_Abelian() : base(Compare()) {}

  Free_Abelian(Compare c) : base(c) {}

  explicit Free_Abelian(const base& m) : base(m) {}

  explicit Free_Abelian(T p, Compare c=Compare()) : base(c)
  { base::insert(std::make_pair(p,1L)); }

  Free_Abelian(T p,coef_t m, Compare c=Compare()) : base(c)
  { base::insert(std::make_pair(p,m)); }

  template<typename InputIterator>
  Free_Abelian(InputIterator first, InputIterator last, Compare c=Compare())
    : base(first,last,c) {}

  Free_Abelian& add_multiple(const Free_Abelian& p, coef_t m);

  Free_Abelian& operator+=(const Free_Abelian& p)
  { return add_multiple(p,1); }

  Free_Abelian& operator-=(const Free_Abelian& p)
  { return add_multiple(p,-1); }

  coef_t operator[] (const T& t) const // find coefficient of |t| in |*this|
  {
    typename base::const_iterator p=base::find(t);
    return p==base::end() ? coef_t(0) : p->second;
  }

}; // Free_Abelian

} // namespace free_abelian

} // namespace atlas

#include "free_abelian_def.h"

#endif
