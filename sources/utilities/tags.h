/*
   Definition of dummy argument tags used for constructor overloading.
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

  // Dummy argument to distinguish two constructors for |Partition|.
  struct UnnormalizedTag {};

/******** tags for structure theory ******************************************/

  // To distinguish constructor for derived group
  struct DerivedTag {};

  // To distinguish constructor for group modulo the central torus
  struct CoderivedTag {};

  // To distinguish constructor for dual object from copy-constructor
  struct DualTag {};

}

}

#endif
