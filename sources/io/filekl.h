
#ifndef FILEKL_H
#define FILEKL_H


#include <iosfwd>

#include "../Atlas.h"

namespace atlas {
  namespace filekl {
    
    const BlockElt no_good_ascent = UndefBlock-1;
     // value flagging that no good ascent exists
    const unsigned int magic_code=0x06ABdCF0; 

    
    void write_block_file(const Block& block, std::ostream& out);
    
    void write_matrix_file(const kl::KL_table& klc, std::ostream& out);
    
    void write_KL_store(const kl::KLStore& store, std::ostream& out);

  }
}
#endif

