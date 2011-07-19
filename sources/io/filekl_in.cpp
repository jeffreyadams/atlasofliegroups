
#include "filekl_in.h"
#include <stdexcept>

#include "basic_io.h"

#include <fstream>

#include "blocks.h"
#include "bitset.h"
namespace atlas {
  namespace filekl {

    using basic_io::read_bytes;

    
    const BlockElt UndefBlock= ~BlockElt(0);
    const BlockElt noGoodAscent= UndefBlock-1;
    
    const unsigned int magic_code=0x06ABdCF0; 


    
    BlockElt
    block_info::primitivize(BlockElt x, BlockElt y) const
    {
      RankFlags d=descent_set[y];
    start:
      if (x>=y) return x; // possibly with |x==UndefBlock|
      const ascent_vector& ax=ascents[x];
      for (size_t s=0; s<rank; ++s)
        if (d[s] and ax[s]!=noGoodAscent)
        { x=ax[s]; goto start; } // this should raise $x$, now try another step
      return x; // no raising possible, stop here
    }
    
    bool
    block_info::is_primitive(BlockElt x, const RankFlags d) const
    {
      const ascent_vector& ax=ascents[x];
      for (size_t s=0; s<ax.size(); ++s)
        if (d[s] and ax[s]!=noGoodAscent)
          return false;
      return true;
      // now |d[s]| implies |ascents[s]==noGoodAscent| for all simple roots |s|
    }
    
    const prim_list& block_info::prims_for_descents_of(BlockElt y)
    { RankFlags d=descent_set[y];
      unsigned long s=d.to_ulong();
      prim_list& result=primitives_list[s];
      if (result.empty())
      { for (BlockElt x=0; x<size; ++x)
          if (is_primitive(x,d))
    	result.push_back(x);
        prim_list(result).swap(result); // reallocate to fit snugly
      }
      return result;
    }
    
    block_info::block_info(std::ifstream& in)
      : rank(), size(), max_length(), start_length()
      , descent_set(), ascents(), primitives_list() // don't initialize yet
    {
      in.seekg(0,std::ios_base::beg); // ensure we are reading from the start
      size=read_bytes<4>(in);
      rank=read_bytes<1>(in);
      max_length=read_bytes<1>(in);
    
      // read intervals of block elements for each length
      start_length.resize(max_length+2);
      start_length[0]=BlockElt(0);
      for (size_t i=1; i<=max_length; ++i) start_length[i]=read_bytes<4>(in);
      start_length[max_length+1]=size;
    
      // read descent sets
      descent_set.reserve(size);
      for (BlockElt y=0; y<size; ++y)
        descent_set.push_back(RankFlags(read_bytes<4>(in)));
    
      // read ascent table
      ascents.reserve(size);
      for (BlockElt x=0; x<size; ++x)
        {
          ascents.push_back(ascent_vector());
          ascent_vector& a=ascents.back();
          a.reserve(rank);
          for (size_t s=0; s<rank; ++s)
    	a.push_back(read_bytes<4>(in));
        }
    
      // lists of primitives lazily computed, on demand by |prims_for_descents_of|
      primitives_list.resize(1ul<<rank); // create $2^{rank}$ empty vectors
    }
    
    size_t matrix_info::length (BlockElt y) const
    { return
        std::upper_bound(block.start_length.begin(),block.start_length.end(),y)
        -block.start_length.begin() // index of first element of length |l(y)+1|
        -1; // now we have just |l(y)|
    }
    
    RankFlags matrix_info::descent_set (BlockElt y) const
      { return block.descent_set[y]; }
    
    void matrix_info::set_y(BlockElt y)
    {
      if (y==cur_y)
        { matrix_file.seekg(cur_row_entries); return; }
      cur_y=y;
      const prim_list& weak_prims = block.prims_for_descents_of(y);
      cur_strong_prims.resize(0);
      // restart building from scratch, but don't deallocate storage
    
      matrix_file.seekg(row_pos[y]);
      size_t n_prim=read_bytes<4>(matrix_file);
    
    
      for (size_t i=0; i<n_prim; i+=32)
        {
          unsigned int chunk=read_bytes<4>(matrix_file);
          for (size_t j=0; chunk!=0; ++j,chunk>>=1) // and certainly |j<32|
            if ((chunk&1)!=0) cur_strong_prims.push_back(weak_prims[i+j]);
        }
    
      {
    #ifndef NDEBUG
        const BlockElt* i=
          // point after first block element of length of |y|
          std::upper_bound(&block.start_length[0]
      		    ,&block.start_length[block.max_length+1]
      		    ,y);
        const BlockElt* first=
          // point to first of |weak_prims| of that length
          std::lower_bound(&weak_prims[0],&weak_prims[weak_prims.size()],*(i-1));
        if (cur_strong_prims.back()!=*first)
          {
    	std::cerr << "For y=" << y << ", first weak prim found " << *first
    		  << " does not match computed value "
    		  << cur_strong_prims.back() << ".\n";
    	throw std::runtime_error("Panic");
          }
    #endif
        cur_strong_prims.back()=y; // replace by |y|
      }
    
      cur_row_entries=matrix_file.tellg();
    }
    
    KLIndex matrix_info::find_pol_nr(BlockElt x,BlockElt y)
    {
      set_y(y);
      x_prim=block.primitivize(x,y);
      if (x_prim>=y)
        return KLIndex(x_prim==y ? 1 : 0); // primitivisation copped out
      strong_prim_list::const_iterator it=
        std::lower_bound(cur_strong_prims.begin(),cur_strong_prims.end(),x_prim);
      if (it==cur_strong_prims.end() or *it!=x_prim)
        return KLIndex(0); // not strong
    
      matrix_file.seekg(4*size_t(it-cur_strong_prims.begin()),std::ios_base::cur);
      return KLIndex(read_bytes<4>(matrix_file));
    }
    
    BlockElt matrix_info::prim_nr(unsigned int i,BlockElt y)
    { const prim_list& weak_prims = block.prims_for_descents_of(y);
      const BlockElt* it=
        // point after first block element of length of |y|
        std::upper_bound(&block.start_length[0]
    		    ,&block.start_length[block.max_length+1]
    		    ,y);
      const BlockElt* first=
        // point to first of |weak_prims| of that length
        std::lower_bound(&weak_prims[0],&weak_prims[weak_prims.size()],*(it-1));
      size_t limit=first-&weak_prims[0]; // limiting value for |i|
      if (i<limit) return weak_prims[i];
      else if (i==limit) return y;
      std::cerr << "Limit is " << limit << ".\n";
      throw std::runtime_error("Weakly primitive element index too large");
    }
    
    matrix_info::matrix_info
      (std::ifstream& block_file,std::ifstream& m_file)
    : matrix_file(m_file) // store reference to the matrix file
      , block(block_file) // read in block information
      , row_pos(block.size) // dimension these vectors
      , cur_y(UndefBlock), cur_strong_prims(), cur_row_entries(0)
    {
      matrix_file.seekg(0,std::ios_base::beg);
      if (read_bytes<4>(matrix_file)==magic_code)
      { matrix_file.seekg(-4*std::streamoff(block_size()),std::ios_base::end);
        std::streamoff cumul=0;
        for (BlockElt y=0; y<block_size(); ++y)
        { cumul+= 4*std::streamoff(read_bytes<4>(matrix_file));
          row_pos[y] = cumul;
        }
      }
    
      else // read old file format, scanning the whole file
      {
        size_t l=0; // length of y
        matrix_file.seekg(0,std::ios_base::beg);
        for (BlockElt y=0; y<block.size; ++y)
        {
          while (y>=block.start_length[l+1]) ++l;
    
          if (read_bytes<4>(matrix_file)!=y and y!=0)
          { std::cerr << y << std::endl;
            throw std::runtime_error ("Alignment problem");
          }
          row_pos[y]= matrix_file.tellg(); // record position after row number
          size_t n_prim=read_bytes<4>(matrix_file);
    
          { size_t n_strong_prim=0;
            { // compute number of entries to skip, while reading bitmap
    	  static const unsigned int ulsize=sizeof(unsigned long int);
    	  size_t count=n_prim/(8*ulsize); // number of |unsigned long|s to read
    
    	  while (count-->0)
    	    n_strong_prim+=bits::bitCount(read_bytes<ulsize>(matrix_file));
    
    	  count=(n_prim%(8*ulsize)+31)/32; // maybe some tetra-byte(s) left
    	  while (count-->0)
    	    n_strong_prim+=bits::bitCount(read_bytes<4>(matrix_file));
    	}
    
    	// now skip over matrix entries
    	matrix_file.seekg(4*std::streamoff(n_strong_prim)
    			 ,std::ios_base::cur);
          }
          if (not matrix_file.good())
          { std::cerr << y << std::endl;
    	throw std::runtime_error ("Premature end of file");
          }
        } // for (BlockElt y...)
      } // |if (...==magic_code)|
    
      block_file.close(); // success, we no longer need the block file
    }
    
    polynomial_info::polynomial_info (std::ifstream& coefficient_file)
    : file(coefficient_file), n_pols(read_bytes<4>(file))
    , coef_size(), index_begin(file.tellg()), coefficients_begin()
    { file.seekg(10,std::ios_base::cur); // skip initial 2 indices
      coef_size=read_bytes<5>(file); // size of the |One|
      file.seekg(index_begin+5*n_pols,std::ios_base::beg);
      n_coef=read_bytes<5>(file)/coef_size;
      coefficients_begin=file.tellg();
    }
    
    polynomial_info::~polynomial_info() { file.close(); }
    
    size_t polynomial_info::degree(KLIndex i) const
    { if (i<2) return i-1; // quit exit for Zero and One
      file.seekg(index_begin+5*i,std::ios_base::beg);
      ullong index=read_bytes<5>(file);
      ullong next_index=read_bytes<5>(file);
      size_t length=(next_index-index)/coef_size;
      return length-1;
    }
    
    std::vector<size_t> polynomial_info::coefficients(KLIndex i) const
    { file.seekg(index_begin+5*i,std::ios_base::beg);
      ullong index=read_bytes<5>(file);
      ullong next_index=read_bytes<5>(file);
      size_t length=(next_index-index)/coef_size;
    
      std::vector<size_t> result(length);
      file.seekg(coefficients_begin+index,std::ios_base::beg);
    
      for (size_t i=0; i<length; ++i)
        result[i]=basic_io::read_var_bytes(coef_size,file);
    
      return result;
    }
    
    size_t polynomial_info::leading_coeff(KLIndex i) const
    { if (i<2) return i; // this makes "leading coefficient" of Zero return 0
      file.seekg(index_begin+5*(i+1),std::ios_base::beg);
      ullong next_index=read_bytes<5>(file);
      file.seekg(coefficients_begin+next_index-coef_size,std::ios_base::beg);
      return basic_io::read_var_bytes(coef_size,file);
    }
    
    ullong polynomial_info::coeff_start(KLIndex i) const
    { file.seekg(index_begin+5*i,std::ios_base::beg);
      ullong index=read_bytes<5>(file);
      return index/coef_size;
    }
    
    cached_pol_info::cached_pol_info(std::ifstream& coefficient_file)
      : polynomial_info(coefficient_file)
      , cache(n_polynomials()-2)
    {
      for (KLIndex i=2; i<n_polynomials(); ++i)
      {
    #ifdef VERBOSE
        if ((i&0xFFF)==0) std::cerr << i << '\r';
    #endif
        size_t d=polynomial_info::degree(i);
        if ((d&~degree_mask)!=0)
          throw std::runtime_error("Degree found too large (>=32)");
        cache[i-2]=d;
      }
    #ifdef VERBOSE
        std::cerr << n_polynomials()-1 << '\n';
    #endif
    }
    
    size_t cached_pol_info::degree (KLIndex i) const
    {
      return i<2 ? i-1 : cache[i-2]&degree_mask ;
    }
    
    size_t cached_pol_info::leading_coeff (KLIndex i) const
    {
      if (i<2) return i;
      if ((cache[i-2]&~degree_mask)!=0) return cache[i-2]/(degree_mask+1);
      size_t lc=polynomial_info::leading_coeff(i); // look up in file
      if (lc<=255/(degree_mask+1)) cache[i-2] |= lc*(degree_mask+1);
      return lc;
    }
    
    progress_info::progress_info(std::ifstream& file)
    : first_pol()
    { file.seekg(0,std::ios_base::end); // measure |file|
      if (file.tellg()%12!=0)
        throw std::runtime_error("Row file size not a multiple of 12");
      BlockElt size= file.tellg()/12;
      first_pol.reserve(size+1); first_pol.push_back(0);
      file.seekg(0,std::ios_base::beg); // rewind
      for (BlockElt y=0; y<size; ++y)
      { file.seekg(8,std::ios_base::cur); // skip ahead
        first_pol.push_back(first_pol.back()+read_bytes<4>(file));
      }
      file.close();
    }
    
    BlockElt progress_info::first_row_for_pol(KLIndex i) const
    { const KLIndex* p=std::upper_bound(&first_pol[1],&*first_pol.end(),i);
      return p-&first_pol[1];
    }

  }
}

