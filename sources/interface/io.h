/*
  This is io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef IO_H  /* guard against multiple inclusions */
#define IO_H

#include <iosfwd>

/******** constants ********************************************************/

namespace atlas {

namespace io {

  extern const char* MESSAGE_DIR; // a modifyable pointer to const char

}

/******** function declarations ********************************************/

namespace io {

  std::ostream& printFile(std::ostream& strm, const char* text_file_name);
  std::ostream& printFile(std::ostream& strm,
			  const char* text_file_name, const char* dir_name);

}

}

#endif
