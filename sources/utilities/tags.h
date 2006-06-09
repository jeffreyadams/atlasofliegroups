/*!
\file
  This is tags.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef TAGS_H  /* guard against multiple inclusions */
#define TAGS_H

/******** type declarations *************************************************/

namespace atlas {

namespace tags {

/******** general tags for structure theory **********************************/
  struct IteratorTag {};
  struct UnnormalizedTag {};

/******** tags for structure theory ******************************************/
  struct AdjointTag {};
  struct DerivedTag {};
  struct DualTag {};
  struct QuasisplitTag {};
  struct SimplyConnectedTag {};

}

}

#endif
