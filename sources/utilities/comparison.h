/*
  This is comparison.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef COMPARISON_H  /* guard against multiple inclusions */
#define COMPARISON_H

/******** type declarations **************************************************/

namespace atlas {

namespace comparison {

  template<typename F> class Compare;

}

namespace comparison {

template<typename F>
  inline Compare<F> compare(F f) {return Compare<F>(f);}

template<typename F>
class Compare {
 private:
  F* d_f;
 public:
  explicit Compare(F f):d_f(&f) {}

  bool operator() (typename F::argument_type x, 
		   typename F::argument_type y) const {
    typename F::result_type lhs = (*d_f)(x);
    typename F::result_type rhs = (*d_f)(y);
    return  lhs < rhs;
  }
};

}

}

#endif
