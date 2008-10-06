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
template<typename T>
struct Free_Abelian : public std::map<T,long int>
{
  typedef long int coef_t;
  typedef std::map<T,coef_t> base;

  Free_Abelian() : base() {}

  explicit Free_Abelian(const base& m) : base(m) {}

  explicit Free_Abelian(T p) : base()
  { base::insert(std::make_pair(p,1L)); }

  Free_Abelian(T p,coef_t m) : base()
  { base::insert(std::make_pair(p,m)); }

  Free_Abelian& add_multiple(const Free_Abelian& p, coef_t m);

  Free_Abelian& operator+=(const Free_Abelian& p)
  { return add_multiple(p,1); }

  Free_Abelian& operator-=(const Free_Abelian& p)
  { return add_multiple(p,-1); }

}; // Free_Abelian

} // namespace free_abelian

} // namespace atlas

#include "free_abelian_def.h"

#endif
