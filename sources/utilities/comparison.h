/*!
\file
  This is comparison.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef COMPARISON_H  /* guard against multiple inclusions */
#define COMPARISON_H

/******** type declarations **************************************************/

namespace atlas {

namespace comparison {

  template<typename F> class Compare;

}

namespace comparison {

/* The template Compare turns a unary function object into a comparison class
   The comparison is done by comparing the images under the function object
   by means of operator<.

   So if |f| is a unary function object, one can set

   Compare<F> comp=Compare<F>(f); // or Compare<F> comp(f);

   after which |comp(x,y)| is equivalent to |f(x)<f(y)|. Rather than being
   stored in a variable, |Compare<F>(f)| is usually passed as an argument (of
   type |const Compare<F>&|) to a (STL) function needing a comparison object.
 */

// template function allowing |compare(f)| to abbreviate |Compare<F>(f)|
template<typename F>
  inline Compare<F> compare(const F& f) { return Compare<F>(f);}

template<typename F> // F is the type of a unary function object
class Compare {
 private:
  const F& d_f;
 public:
  explicit Compare(const F& f):d_f(f) {}

  bool operator() (typename F::argument_type x,
		   typename F::argument_type y) const
  {
    return  d_f(x) < d_f(y); // F::result_type should define operator<
  }
};

}

}

#endif
