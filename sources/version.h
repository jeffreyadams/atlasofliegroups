/*
  This is version.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef VERSION_H  /* guard against multiple inclusions */
#define VERSION_H

namespace atlas {

namespace version {

const char* const NAME = "the Atlas of Lie Groups and Representations";
const char* const VERSION = "0.6"; // last advanced January, 2016
const char* const COMPILEDATE = __DATE__ " at " __TIME__;

}

}

#endif
