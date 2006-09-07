/*!
\file
\brief Definition of dummy argument tags used for constructor overloading.
*/
/*
  This is tags.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef TAGS_H  /* guard against multiple inclusions */
#define TAGS_H

/******** type declarations *************************************************/

namespace atlas {

namespace tags {

/******** general tags for structure theory **********************************/

  /*!
\brief Dummy argument to distinguish constructors, etc. using
iterators to pass multiple arguments.

Typical use is the Matrix constructors defined in
sources/utilities/matrix_def.h. 
  */
  struct IteratorTag {};

  /*!
\brief Dummy argument to distinguish two constructors for Partition.
  */
  struct UnnormalizedTag {};

/******** tags for structure theory ******************************************/

  /*!
\brief Dummy argument to distinguish constructors, etc. for an adjoint
group.

[Not used 7/24/06. DV.]
  */
  struct AdjointTag {};

  /*!
\brief Dummy argument to distinguish constructors, etc. for derived
group.

[Not used 7/24/06. DV.]
  */
  struct DerivedTag {};

  /*!
\brief Dummy argument to distinguish constructors, etc. for dual group objects.

The presence of the argument DualTag says that a constructor should
make an object for the dual group (as for example in the constructor
RealTorus declared in the file sources/structure/tori.h).
  */
  struct DualTag {};

  /*!
\brief Dummy argument to distinguish constructors, etc. for a quasisplit
group.

[Not used 7/24/06. DV.]
  */
  struct QuasisplitTag {};

  /*!
\brief Dummy argument to distinguish constructors, etc. for a simply connected
group.

[Not used 7/24/06. DV.]
  */
  struct SimplyConnectedTag {};

}

}

#endif
