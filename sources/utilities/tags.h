/*!
\file
\brief Definition of dummy argument tags used for constructor overloading.
*/
/*
  This is tags.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef TAGS_H  /* guard against multiple inclusions */
#define TAGS_H

/******** type declarations *************************************************/

namespace atlas {

namespace tags {

/******** general tags for structure theory **********************************/

  /*! Dummy argument to distinguish constructors, etc. using iterators to pass
multiple arguments form those with non-iterator arguments that would otherwise
match template

Typical use is the Matrix constructors defined in
sources/utilities/matrix.h
  */
  struct IteratorTag {};

  /*!
\brief Dummy argument to distinguish two constructors for Partition.
  */
  struct UnnormalizedTag {};

/******** tags for structure theory ******************************************/

  // To distinguish constructor for an adjoint group from copy-constructor
  struct AdjointTag {};

  // To distinguish constructor for derived group from copy-constructor.
  struct DerivedTag {};

  // To distinguish constructor for derived group from copy-constructor.
  struct SimplyConnectedTag {};

  // To distinguish constructor for dual object from copy-constructor
  struct DualTag {};

  // To distinguish an expermental function (constructor) from one it is 
  // ultimately destined to replace.
  struct NewTag {};

}

}

#endif
