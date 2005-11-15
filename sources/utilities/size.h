/*
  This is size.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef SIZE_H  /* guard against multiple inclusions */
#define SIZE_H

#include "constants.h"

/******** type declarations **************************************************/

namespace atlas {

namespace size {

  template<typename C> class SizeType;
  // this may have to be modified if RANK_MAX increases
  typedef signed char BaseType;
  typedef unsigned char UnsignedBaseType;
  typedef SizeType<BaseType> Size;

}

/******** constant declarations **********************************************/

namespace size {

  // useless value
  template<unsigned long n> class PrimesMax {
  public:
    static const unsigned long value = 0ul;
  };


  // predefined values for likely instances of RANK_MAX
  template<> class PrimesMax<8> {
  public:
    static const unsigned long value = 4ul;
  };


  template<> class PrimesMax<16> {
  public:
    static const unsigned long value = 7ul;
  };


  template<> class PrimesMax<32> {
  public:
    static const unsigned long value = 11ul;
  };


  template<> class PrimesMax<64> {
  public:
    static const unsigned long value = 18ul;
  };


  const size_t PRIMES_MAX = PrimesMax<constants::RANK_MAX>::value;
}

/******** function declarations **********************************************/

namespace size {

  template<typename C> void factorial (SizeType<C>&, unsigned long);
  unsigned long prime(size_t);

}

/******** type definitions ***************************************************/

namespace size {

template<typename C> class SizeType {

 private:
  C d_data[PRIMES_MAX];

 public:
// constructors and destructors
  SizeType() {
    memset(d_data,0,PRIMES_MAX);
  }

  explicit SizeType(unsigned long);

  ~SizeType()
    {}

// copy and assignment
  SizeType(const SizeType& a) {
    memcpy(d_data,a.d_data,PRIMES_MAX);
  }

  SizeType& operator=(const SizeType& a) {
    memcpy(d_data,a.d_data,PRIMES_MAX); return *this;
  }

// accessors
  C operator[] (size_t j) const {
    return d_data[j];
  }

  bool operator== (const SizeType& c) const {
    return !memcmp(d_data,c.d_data,PRIMES_MAX);
  }

  bool operator!= (const SizeType& c) const {
    return memcmp(d_data,c.d_data,PRIMES_MAX);
  }

  bool hasOverflow() const;

  bool hasOverflow(size_t) const;

  unsigned long piece(size_t) const;

  unsigned long toUlong() const;

// manipulators
  C& operator[] (size_t j) {
    return d_data[j];
  }

  SizeType& operator*= (const SizeType&);

  SizeType& operator*= (unsigned long);

  SizeType& operator/= (const SizeType&);

  void reset() {
    memset(d_data,0,PRIMES_MAX);
  }

  void twoShift(C n) { // multiplication by 2^n; prime #0 is 2
    d_data[0] += n;
  }
};

}

}

#include "size_def.h"

#endif
