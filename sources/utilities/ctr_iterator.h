/*!
\file
  This is ctr_iterator.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef CTR_ITERATOR_H  /* guard against multiple inclusions */
#define CTR_ITERATOR_H

#include <iterator>

namespace atlas {

/*! The purpose of this type is to be able to apply STL algorithms like
  |std::find|, in situations where the values compared are really integers but
  not actually all stored in memory. This is achieved by having the iterators
  themselves containing the values, so that |operator*| just delivers the
  value stored. To be of any interest, the comparison function should have
  some special context dependent meaning, probably involving table lookup. It
  is a major inconvenience that, since we are comparing integers, the values
  compared against must be accessible in this way, in other words stored in
  the same table that we are (implicitly) searching. For this reason these
  iterators are currently unsused by atlas, and might well remain so forever.

  It is assumed that U is a (typically unsigned) integral type.
*/

/******** type declarations *************************************************/

namespace ctr_iterator {

  template<typename U> class CounterIterator;

}

/******** function declarations **********************************************/

/******** type definitions ***************************************************/

namespace ctr_iterator {

template<typename U> class CounterIterator
  // derive so that associated types are defined (no pointer and ref allowed)
: public std::iterator<std::random_access_iterator_tag,U,ptrdiff_t,void,void>
{

 private:

  U d_val;

 public:

// constructors and destructors
  explicit CounterIterator(U n) : d_val(n) {}

  ~CounterIterator() {}

// assignment
  CounterIterator& operator= (const CounterIterator& i) {
    d_val = i.d_val;
    return *this;
  }

// accessors
  U operator* () const {
    return d_val;
  }

  bool operator== (const CounterIterator& i) const {
    return d_val == i.d_val;
  }

  bool operator!= (const CounterIterator& i) const {
    return d_val != i.d_val;
  }

  // manipulators
  CounterIterator& operator++ () {
    ++d_val;
    return *this;
  }

  CounterIterator operator++ (int) {
     return CounterIterator(d_val++);
  }

  CounterIterator& operator-- () {
    --d_val;
    return *this;
  }

  CounterIterator operator-- (int) {
     return CounterIterator(d_val--);
  }

  CounterIterator& operator+= (difference_type n) {
    d_val += n;
    return *this;
  }

  CounterIterator operator+ (difference_type n) {
    return CounterIterator(d_val + n);
  }

  CounterIterator& operator-= (difference_type n) {
    d_val -= n;
    return *this;
  }

  CounterIterator operator- (difference_type n) {
    return CounterIterator(d_val - n);
  }

  difference_type operator- (const CounterIterator& i) const {
    return d_val - i.d_val;
  }

  U operator[] (difference_type n) {
    return d_val + n;
  }

};

template<typename U>
inline CounterIterator<U> operator+
  (typename CounterIterator<U>::difference_type n, CounterIterator<U> i) {
  return i+n;
}

}

}

#endif
