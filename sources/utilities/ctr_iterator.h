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

namespace atlas {

/*!
  The purpose of this type is to be able to apply stl algorithms like find,
  where the iterator is really an integer, ranging over a non-allocated range.
  In other words, the iterator itself holds its value. It is assumed that
  U is a (typically unsigned) integral type.
*/

/******** type declarations *************************************************/

namespace ctr_iterator {

  template<typename U> class CounterIterator;

}

/******** function declarations **********************************************/

/******** type definitions ***************************************************/

namespace ctr_iterator {

template<typename U> class CounterIterator {

 private:

  U d_val;

 public:

// associated types
  typedef std::random_access_iterator_tag iterator_category;
  typedef U value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  explicit CounterIterator(U n):d_val(n) {}

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
