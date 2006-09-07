/*
  This is typenumber.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef TYPENUMBER_HPP
#define TYPENUMBER_HPP

#include <cstddef>

namespace atlas {

/******** constant declarations **********************************************/

namespace typenumber {

enum Types { Int, Ulong, WeylElt, Unknown };

}

/******** type declarations **************************************************/

namespace typenumber {

template <typename T> struct TypeNumber;
class TypeData;

}

/******** functions declarations *********************************************/

namespace typenumber {

const char* name(size_t);

}

/******** type definitions ***************************************************/

namespace typenumber {

static size_t UnknownNumber = Unknown;

template<typename T> struct TypeNumber {
  static size_t number() {
    static size_t nbr = UnknownNumber++;
    return nbr;
  };
};

template<> struct TypeNumber<int> {
  static size_t number() {
    return Int;
  };
};

template<> struct TypeNumber<unsigned long> {
  static size_t number() {
    return Ulong;
  };
};

class TypeData {

 private:

  size_t d_size;
  size_t d_number;

  TypeData() {}

 public:

  template<typename T> TypeData(T t)
    :d_size(sizeof(T)),d_number(TypeNumber<T>::number())
    {}

// accessors
  size_t size() const {
    return d_size;
  }

  size_t number() const {
    return d_number;
  }
};

}

}

#endif
