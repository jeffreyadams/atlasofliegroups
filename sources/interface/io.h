/*
  This is io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For copyright and license information see the LICENSE file
*/

#ifndef IO_H  /* guard against multiple inclusions */
#define IO_H

/* the following macro should be set to an absolute path from the atlas
   makefile, but if not we use a relative path from the atlas directory */
#ifndef MESSAGE_DIR_MACRO
#define MESSAGE_DIR_MACRO "messages/"
#endif

#include <iosfwd>

/******** constants ********************************************************/

namespace atlas {

namespace io {

  const char* const MESSAGE_DIR = MESSAGE_DIR_MACRO;

}

/******** function declarations ********************************************/

namespace io {

  std::ostream& printFile(std::ostream& strm, const char* text_file_name);
  std::ostream& printFile(std::ostream& strm,
			  const char* text_file_name, const char* dir_name);

}

}

#endif
