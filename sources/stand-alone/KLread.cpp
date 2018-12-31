#include  <string>
#include  <vector>
#include  <iostream>
#include  <stdexcept>
#include  <fstream>
#include  <bitset>
#include  <cassert>
#include  <algorithm>
#include  <sstream>
#include  <cctype>
#include  <memory>


typedef unsigned long long int ullong; 

typedef ullong KLIndex;
 // must be |long long| to avoid problems when multiplied

class polynomial_info
{ std::ifstream& file; // non-owned reference

  KLIndex n_pols;  // must be 64-bits, for computing |5*n_pols|
  unsigned int coef_size;
  ullong n_coef;
  std::streamoff index_begin, coefficients_begin;

  const size_t max_length;
public:
  polynomial_info(std::ifstream& coefficient_file, size_t deg_limit);

  KLIndex n_polynomials() const { return n_pols; }
  unsigned int coefficient_size() const { return coef_size; }
  ullong n_coefficients() const { return n_coef; }

  std::vector<ullong> coefficients(KLIndex i) const;
};

typedef unsigned int BlockElt;
typedef std::bitset<32> RankFlags; // we can go up to rank 32
typedef std::vector<RankFlags> descent_set_vector; // indexed by block element

typedef std::vector<BlockElt> ascent_vector;       // indexed by simple root
typedef std::vector<ascent_vector> ascent_table;   // indexed by block element

typedef std::vector<BlockElt> prim_list;           // list of weak primitives
typedef std::vector<prim_list> prim_table;
 

struct block_info
{
  unsigned int rank;
  size_t size;
  unsigned int max_length;
  std::vector<BlockElt> start_length;
   // size |max_length+2|; defines intervals for each length

  descent_set_vector descent_set;     // $453060*4    =  1812240$ bytes:  2~MB

private:
  ascent_table ascents;               // $453060*8*4  = 14497920$ bytes: 13~MB
  prim_table primitives_list;         // $256*15000*4 = 15360000$ bytes: 14~MB

public:
  block_info(std::ifstream& in); // constructor reads file

  BlockElt primitivise(BlockElt x, const RankFlags d) const;
  const prim_list& prims_for_descents_of(BlockElt y);
private:
  bool is_primitive(BlockElt x, const RankFlags d) const;
  void compute_prim_table();
};

typedef prim_list strong_prim_list;
 

class matrix_info
{
  std::fstream& matrix_file;

  block_info block;

  std::vector<std::streampos> row_pos;     // $453060*8= 3624480$ bytes: 3 MB

// data for currently selected row~|y|
  BlockElt cur_y;		       // $4$ bytes
  strong_prim_list cur_strong_prims;   // $6400*4 = 25600$ bytes (but max 2MB)
  std::streampos cur_row_entries;
   // 8 bytes; points to start of indices for |cur_strong_prims|
//private methods
  matrix_info(const matrix_info&); // copying forbidden
  void set_y(BlockElt y);          // install |cur_y| and dependent data

public:
  BlockElt x_prim; // public variable that is set by |find_pol_nr|
  enum mode { old, revised, transform };
  matrix_info(std::ifstream* block_file,std::fstream* m_file,bool new_format);
  ~matrix_info() { delete &matrix_file; }

  BlockElt block_size() const { return block.size; }
  KLIndex find_pol_nr(BlockElt x,BlockElt y);
  // sets |y|, so not |const|
  std::streamoff row_offset(BlockElt y) const { return row_pos[y]; }
  BlockElt prim_nr(unsigned int i,BlockElt y);
    // find primitive element by its index~|i|
};

class progress_info
{
  std::vector<KLIndex> first_pol;
    // count distinct polynomials in rows before this one
public:
  progress_info(std::ifstream& progress_file);

  BlockElt block_size() const { return first_pol.size()-1; }
  KLIndex first_new_in_row(BlockElt y) const // |y==block_size()| is allowed
  { return first_pol[y]; }
  BlockElt first_row_for_pol(KLIndex i) const;
   // $i={}$\# polynomials allowed
};

const BlockElt UndefBlock= ~BlockElt(0);
const BlockElt noGoodAscent= UndefBlock-1;

const std::ios_base::openmode binary_in=
			    std::ios_base::in
			  | std::ios_base::binary;

const std::ios_base::openmode binary_in_out=
			    std::ios_base::in
			  | std::ios_base::out
			  | std::ios_base::binary;

const unsigned int magic_code=0x06ABdCF0;
const unsigned int work_in_progress=0x76543210;

bool verbose=true;

unsigned int add_bits(unsigned long int x);

template <unsigned int n>
inline ullong read_bytes(std::istream& in)
{
  return static_cast<unsigned char>(in.get())+(read_bytes<n-1>(in)<<8);
}

template<>
inline ullong read_bytes<1>(std::istream& in)
{
  return static_cast<unsigned char>(in.get());
}

ullong read_var_bytes(unsigned int n,std::istream& in)
{ switch(n)
  { case 1: return read_bytes<1>(in);
    case 2: return read_bytes<2>(in);
    case 3: return read_bytes<3>(in);
    case 4: return read_bytes<4>(in);
    case 5: return read_bytes<5>(in);
    case 6: return read_bytes<6>(in);
    case 7: return read_bytes<7>(in);
    case 8: return read_bytes<8>(in);
    default: throw std::runtime_error("Illegal read_var_bytes");
  }
}

void write_int(unsigned int n, std::ostream& out)
{ out.put(char(n&0xFF)); n>>=8;
out.put(char(n&0xFF)); n>>=8;
out.put(char(n&0xFF)); n>>=8;
out.put(char(n));
}

polynomial_info::polynomial_info
  (std::ifstream& coefficient_file, size_t deg_limit)
: file(coefficient_file), n_pols(read_bytes<4>(file))
, coef_size(), index_begin(file.tellg()), coefficients_begin()
, max_length(deg_limit)
{ file.seekg(10,std::ios_base::cur); // skip initial 2 indices
  coef_size=read_bytes<5>(file); // size of the |One|
  file.seekg(index_begin+5*n_pols,std::ios_base::beg);
  n_coef=read_bytes<5>(file)/coef_size;
  coefficients_begin=file.tellg();
}

std::vector<ullong> polynomial_info::coefficients(KLIndex i) const
{ file.seekg(index_begin+5*i,std::ios_base::beg);
  ullong index=read_bytes<5>(file);
  ullong next_index=read_bytes<5>(file);
  size_t length=(next_index-index)/coef_size;
  if (length>max_length)
    {
      std::cerr << "Found degree " << length-1 
		<< ", for polynomial #" << i << ".\n";
      throw std::runtime_error("Degree exceeds limit given");
    }

  std::vector<ullong> result(length);
  file.seekg(coefficients_begin+index,std::ios_base::beg);

  for (size_t i=0; i<length; ++i)
    result[i]=read_var_bytes(coef_size,file);

  if (not file.good())
  {
      file.clear();
      throw std::runtime_error("Input error reading coefficients");
    }

  return result;
}

BlockElt block_info::primitivise(BlockElt x, const RankFlags d) const
{ start:
  const ascent_vector& ax=ascents[x];
  for (size_t s=0; s<rank; ++s)
    if (d[s] and ax[s]!=noGoodAscent)
    { if (ax[s]==UndefBlock) return UndefBlock;
      x=ax[s]; goto start; // this should raise x, now try another step
    }
  return x; // no raising possible, stop here
}

bool  block_info::is_primitive(BlockElt x, const RankFlags d) const
{
  const ascent_vector& ax=ascents[x];
  for (size_t s=0; s<ax.size(); ++s)
    if (d[s] and ax[s]!=noGoodAscent) return false;
  return true;
  // now |d[s]| implies |ascents[s]==noGoodAscent| for all simple roots |s|
}

const prim_list& block_info::prims_for_descents_of(BlockElt y)
{ RankFlags d=descent_set[y];
  unsigned long s=d.to_ulong();
  prim_list& result=primitives_list[s];
  if (result.empty())
  { for (BlockElt x=0; x<size; ++x)
      if (is_primitive(x,d)) result.push_back(x);
    prim_list(result).swap(result); // reallocate to fit snugly
  }
  return result;
}

block_info::block_info(std::ifstream& in)
  : rank(), size(), max_length(), start_length()
  , descent_set(), ascents(), primitives_list() // don't initialise yet
{
  in.seekg(0,std::ios_base::beg); // ensure we are reading from the start
  size=read_bytes<4>(in);
  rank=read_bytes<1>(in);
  max_length=read_bytes<1>(in);

  if (verbose)
    std::cout << "rank=" << rank << ", block size=" << size 
	      << ", maximal length=" << max_length << ".\n";

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

  primitives_list.resize(1ul<<rank); // create $2^{rank}$ empty vectors
}

void matrix_info::set_y(BlockElt y)
{
  if (y==cur_y) { matrix_file.seekg(cur_row_entries); return; }
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
    const BlockElt* i= // point after first block element of length of |y|
      std::upper_bound(&block.start_length[0]
  		    ,&block.start_length[block.max_length+1]
  		    ,y);
    const BlockElt* first= // point to first of |weak_prims| of that length
      std::lower_bound(&weak_prims[0],&weak_prims[weak_prims.size()],*(i-1));
    assert(cur_strong_prims.back()==*first);
    cur_strong_prims.back()=y; // replace by |y|
  }

  cur_row_entries=matrix_file.tellg();
}

KLIndex matrix_info::find_pol_nr(BlockElt x,BlockElt y)
{
  set_y(y);
  RankFlags d=block.descent_set[y];
  x_prim=block.primitivise(x,d);
  if (x_prim==UndefBlock) return KLIndex(0); // primitivisation copped out
  strong_prim_list::const_iterator it=
    std::lower_bound(cur_strong_prims.begin(),cur_strong_prims.end(),x_prim);
  if (it==cur_strong_prims.end() or *it!=x_prim)
    return KLIndex(0); // not strong

  matrix_file.seekg(4*size_t(it-cur_strong_prims.begin()),std::ios_base::cur);
  return KLIndex(read_bytes<4>(matrix_file));
}

BlockElt matrix_info::prim_nr(unsigned int i,BlockElt y)
{ const prim_list& weak_prims = block.prims_for_descents_of(y);
  const BlockElt* it= // point after first block element of length of |y|
    std::upper_bound(&block.start_length[0]
		    ,&block.start_length[block.max_length+1]
		    ,y);
  const BlockElt* first= // point to first of |weak_prims| of that length
    std::lower_bound(&weak_prims[0],&weak_prims[weak_prims.size()],*(it-1));
  size_t limit=first-&weak_prims[0]; // limiting value for |i|
  if (i<limit) return weak_prims[i];
  else if (i==limit) return y;
  std::cerr << "Limit is " << limit << ".\n";
  throw std::runtime_error("Weakly primitive element index too large");
}

matrix_info::matrix_info
  (std::ifstream* block_file,std::fstream* m_file, bool new_format)
: matrix_file(*m_file) // store reference to the matrix file
  , block(*block_file) // read in block information
  , row_pos(block.size) // dimension these vectors
  , cur_y(UndefBlock), cur_strong_prims(), cur_row_entries(0)
{
  if (new_format)
    
    { matrix_file.seekg(-4*std::streamoff(block_size()),std::ios_base::end);
      if (verbose) std::cerr << "Reading indices into the matrix file... ";
      std::streamoff cumul=0;
      for (BlockElt y=0; y<block_size(); ++y)
      { cumul+= 4*std::streamoff(read_bytes<4>(matrix_file));
        row_pos[y] = cumul;
      }
      if (verbose) std::cerr << "done.\n";
    }
  else
    
    { if (verbose)
        std::cerr << "Starting to scan matrix file by 'rows'.\n";
    
      size_t l=0; // length of y
      matrix_file.seekg(0,std::ios_base::beg);
      for (BlockElt y=0; y<block.size; ++y)
        { 
          while (y>=block.start_length[l+1])
          { ++l;
            if (verbose)
              std::cerr << "length " << l << " starts at y=" << y << std::endl;
          }
          if (verbose) std::cerr << y << '\r';
    
          if (read_bytes<4>(matrix_file)!=y and y!=0)
          { std::cerr << y << std::endl;
            throw std::runtime_error ("Alignment problem");
          }
          row_pos[y]= matrix_file.tellg(); // record position after row number
          size_t n_prim=read_bytes<4>(matrix_file);
    
    #if 0
          
          { const prim_list& weak_prims = block.prims_for_descents_of(y);
            prim_list::const_iterator i= // find limit of |length<l| values
                std::lower_bound(weak_prims.begin(),weak_prims.end()
          		      ,block.start_length[l]);
          
          
            if (n_prim!=size_t(i-weak_prims.begin())+1)
              { std::cerr << y << std::endl;
                throw std::runtime_error ("Primitive count problem");
              }
          }
    #endif
          
          { size_t n_strong_prim=0;
            
            { static const unsigned int ulsize=sizeof(unsigned long int);
            
              size_t count=n_prim/(8*ulsize);
              // number of |unsigned long int|s to read
            
              while (count-->0)
                n_strong_prim+=add_bits(read_bytes<ulsize>(matrix_file));
            
              count=(n_prim%(8*ulsize)+31)/32;
               // maybe some tetra-byte(s) left to read
              while (count-->0)
                n_strong_prim+=add_bits(read_bytes<4>(matrix_file));
            }
          
            matrix_file.seekg(4*std::streamoff(n_strong_prim),std::ios_base::cur);
          }
          if (not matrix_file.good())
          { std::cerr << y << std::endl;
    	  throw std::runtime_error ("Premature end of file");
    	}
        }
    
        if (verbose) std::cerr << "\nDone scanning matrix file.\n";
    }
  delete block_file; // success, we no longer need the block file
}

unsigned int add_bits(unsigned long int x)
{ static const unsigned long int b0= ~(~0ul/3);
  // |0xAAAA|\dots ; flags odd bit positions
  static const unsigned long int b1= ~(~0ul/5);
  // |0xCCCC|\dots ; flags positions $\cong2,3 \pmod4$
  static const unsigned long int b2= ~(~0ul/17);
  // |0xF0F0|\dots ; flags positions $\cong4$--$7\pmod8$
  static const unsigned long int ones= ~0ul/255;
  // |0x0101|\dots ; flags the low bit of each octet
  static const unsigned int high_byte_shift=8*(sizeof(unsigned long int)-1);

  x-=(x&b0)>>1;          // replace pairs of bits $10\to01$ and $11\to10$
  x=(x&~b1)+((x&b1)>>2);
   // sideways add 2 groups of pairs of bits to 4-tuples of bits
  x += x>>4;
   // the sums of octets (bytes) are now in lower 4-tuples of those octets
  return (x&~b2)*ones >> high_byte_shift;
   // add lower 4-tuples of bytes in high octet, and extract
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
  if (p==&*first_pol.end())
    throw std::runtime_error("Polynomial index out of range");
  return p-&first_pol[1];
}

BlockElt locate_KL_polynomial (KLIndex i,matrix_info& mi,BlockElt y)
{
  for (BlockElt x=0; x<=y; ++x)
    if (mi.find_pol_nr(x,y)==i) return x;
  if (i==0 and y+1<mi.block_size())
    return mi.x_prim=y+1; // an entry 0 above the diagonal
  throw std::runtime_error("Polynomial could not be located");
}

std::pair<BlockElt,BlockElt> locate_KL_polynomial
  (KLIndex i,matrix_info& mi, const progress_info& pi)
{ BlockElt y=pi.first_row_for_pol(i);
  return std::make_pair(locate_KL_polynomial(i,mi,y),y);
}

void usage()
{ std::cerr <<
  "Usage: KLread [-q] [-l degree] [blockfile matrixfile [rowfile]] polfile\n";
  exit(1);
}

int main(int argc,char** argv)
{ std::string program_name(*argv);
  --argc; ++argv; // read and skip program name

  if (argc==0 or std::string(*argv)=="-help") usage();
  if (std::string(*argv)=="-q") { verbose=false; --argc; ++argv;}

  size_t degree_limit=32; //default value is OK for split~$E_8$
  if (argc>=2 and std::string(*argv)=="-l")
    // then override default |degree_limit|
  { std::istringstream arg_text(argv[1]);
    arg_text >> degree_limit;
    if (arg_text.good()) { argc-=2; argv+=2; }
    else { std::cerr << "Non-numeric argument following -l\n"; exit(1); }
  }

  std::unique_ptr<matrix_info> mi; // unique pointer guarantees clean-up at end
  std::unique_ptr<progress_info> row_info;
  std::ifstream coef_file;

  
  if (argc>=3)
  { argc-=2; // we will consume at least two arguments
    
    {
      std::ifstream* block_file=new std::ifstream;
    std::fstream* matrix_file=new std::fstream;
    block_file->open(*argv++,binary_in);
      matrix_file->open(*argv++,binary_in);
      if (block_file->is_open() and matrix_file->is_open())
      { matrix_info::mode format; 
                                  { unsigned int code=read_bytes<4>(*matrix_file);
                                    format= code==magic_code ? matrix_info::revised :  matrix_info::old;
                                    if (format== matrix_info::old and program_name=="KLwrite")
                                      format= matrix_info::transform;
                                    if (code==work_in_progress)
                                    { if (program_name=="KLwrite")
                                        std::cout << "Reattempting conversion of matrix file.\n";
                                      else
                                      { std::cout << "Broken matrix file, retry the conversion.\n";
                                       exit(1);
                                      }
                                    }
                                    std::cout << "Matrix file format: "
                                      << ( format== matrix_info::old ? "old"
                                         : format== matrix_info::revised ? "new"
                                         : "updating to new" )
                                      << ".\n";
                                  }
        if (format==matrix_info::transform)
          
          { matrix_file->close();
            matrix_file->open(argv[-1],binary_in_out);
            if (not matrix_file->is_open())
            {  std::cerr << "failed to open file '"
          		<< argv[-1] << "' for writing.\n";
                delete block_file; delete matrix_file; exit(1);
            }
          }
        mi=std::unique_ptr<matrix_info> 
           (new matrix_info(block_file,matrix_file,format==matrix_info::revised));
        if (format==matrix_info::transform)
          
          {
            std::streamoff here=matrix_file->tellg();
            matrix_file->seekg(0,std::ios_base::beg); // reset to beginning of file
            write_int(work_in_progress,*matrix_file); // to set ``work in progress''
            matrix_file->seekg(0,std::ios_base::end); // then append to end of file
            if (here!=std::streamoff(matrix_file->tellg()))
              // but we should have been there already
            { std::cerr << "Not at end of file after reading all parts: " << here 
              << "!=" << std::streamoff(matrix_file->tellg()) << ".\n";
              exit(1);
            }
            std::cout << "Appending information to matrix file.\n";
          
            write_int(mi->row_offset(0)/4,*matrix_file);
              // first value written is not a difference
            for (BlockElt y=1; y<mi->block_size(); ++y)
            { if (mi->row_offset(y)%4!=0)
                throw std::runtime_error("Matrix row alignment error");
              write_int((mi->row_offset(y)-mi->row_offset(y-1))/4,*matrix_file);
            }
          
            matrix_file->seekg(0,std::ios_base::beg); // reset to beginning of file
            write_int(magic_code,*matrix_file);
            std::cout << "Conversion of matrix file successfully completed.\n";
          }
      }
      else
        { std::cerr << "failed to open file '"
    		<< argv[block_file->is_open()? -1 : -2 ] << "'.\n";
          delete block_file; delete matrix_file; exit(1);
        }
    }
    if (argc==2)
    { --argc; // we consume a third argument
      
      { std::ifstream row_file(*argv++,binary_in);
        if (row_file.is_open())
        {
          row_info=std::unique_ptr<progress_info>(new progress_info(row_file));
          if (row_info->block_size()!=mi->block_size())
            throw std::runtime_error("Block size mismatch for row file");
        }
        else
          std::cerr << "failed to open file '" << argv[-1] << "', doing without.\n";
      }
    }
  }
  if (argc!=1) usage();
  coef_file.open(argv[0],binary_in);
  if (not coef_file.is_open())
    { std::cerr << "Open failed"; exit(1); }
  
  polynomial_info pol(coef_file,degree_limit);
  
  std::cout << "Coefficient size " << pol.coefficient_size() << ".\n" 
            << pol.n_polynomials() << " polynomials, "
  	  << pol.n_coefficients() << " coefficients.\n";


  
  do
  {
    try
    { if (mi.get()==NULL)
        // if auto-pointer |mi| is unset, no matrix information is present
        std::cout << "index: ";
      else
        std::cout << "give block elements x,y, or polynomial index i as #i: ";
  
      while(isspace(std::cin.peek())) std::cin.get();
  
      KLIndex i;
      if (std::cin.peek()=='#' ? std::cin.get(), true : mi.get()==NULL )
      { 
        {
          std::cin >> i;
          if (not std::cin.good())
            
            { if (std::cin.eof()) break;
              std::string s; std::cin.clear(); std::cin>>s;
              if (s=="quit") break;
              throw std::runtime_error("Non-numeric input");
            }
          if (i>=pol.n_polynomials())
          { std::cerr << "Limit is " << pol.n_polynomials()-1 << ".\n";
            throw std::runtime_error("Index too large");
          }
        }
        if ((mi.get()!=NULL and std::cin.peek()==':') or row_info.get()!=NULL)
        { static const BlockElt UndefBlock = ~0u;
          BlockElt x=UndefBlock,y;
  	if (std::cin.peek()==':' or std::cin.peek()=='>')
          { bool once=std::cin.peek()==':';
  	  
  	  {
  	    std::cin.get(); std::cin >> y;
  	    if (not std::cin.good())
  	    { std::cin.clear(); throw std::runtime_error("Non-numeric input"); }
  	    if (y>=mi->block_size())
  	      throw std::runtime_error("Row number too large");
  	  }
  	  if (once)
  	    x=locate_KL_polynomial(i,*mi,y);
  	  else
  	  { while(++y<mi->block_size())
  	    { try { x=locate_KL_polynomial(i,*mi,y); break;}
  	      catch(std::runtime_error) { std::cerr << y+1 << '\r'; }
  	    }
              if (x==UndefBlock)
                throw std::runtime_error("Not found");
            }
  	}
          else
          { std::pair<BlockElt,BlockElt> p=
              locate_KL_polynomial(i,*mi,*row_info);
          x=p.first; y=p.second;
  	}
          std::cout << "P_{" << x << ',' << y <<  "}=P_{"
  	          << mi->x_prim << ',' << y << "}:\n";
        }
      }
      else
        
        {
          int c; BlockElt x,y;
          while(ispunct(c=std::cin.peek()) or isspace(c)) std::cin.get();
          // skip spaces/punctuation
          std::cin >> x;
          if (not std::cin.good()) // non-numeric input
            
            { if (std::cin.eof()) break;
              std::string s; std::cin.clear(); std::cin>>s;
              if (s=="quit") break;
              throw std::runtime_error("Non-numeric input");
            }
          
          bool convert=std::cin.peek()=='#';
          if (convert) std::cin.get(); 
        
          while(ispunct(c=std::cin.peek()) or isspace(c)) std::cin.get();
          // skip spaces/punctuation
          std::cin >> y;
          if (not std::cin.good())
            throw std::runtime_error("Failure reading y");
        
          
          if (convert) x=mi->prim_nr(x,y);
          if (x>=mi->block_size())
            throw std::runtime_error("First parameter too large");
          if (y>=mi->block_size())
            throw std::runtime_error("Second parameter too large");
        
          if (x>y)
            throw std::runtime_error
              ("Result null by triangularity."); // not really an error
        
          i=mi->find_pol_nr(x,y);
        
          if (mi->x_prim==UndefBlock) // then |primitivise| hit a real non-parity case
            throw std::runtime_error
              ("Result is null because raising the first argument " 
                  "reaches a real non-parity case."); // not really an error
          else
            std::cout << "P_{" << x << ',' << y << "}=P_{" << mi->x_prim << ',' << y 
        	      << "}=polynomial #" << i << ':' << std::endl;
        }
      
      {
        std::vector<ullong> coefficients(pol.coefficients(i));
      
        bool first=true;
        for (size_t i=coefficients.size(); i-->0;)
          if (coefficients[i]!=0)
          { if (first) first=false; else std::cout << " + ";
            if (coefficients[i]!=1 or i==0) std::cout << coefficients[i];
            std::cout << (i==0? "" : "q");
            if (i>1)
            {	if (i<10) std::cout << '^' << i;
      	else std::cout << "^{" << i << '}';
            }
          }
        if (coefficients.size()==0) std::cout << 0;
          // print something for the null polynomial
        else if (coefficients.size()>1)
        // for non-constant polynomials show value at $q=1$
        { ullong sum=0; // sum of coefficients
          for (size_t i=0; i<coefficients.size(); ++i) sum+=coefficients[i];
          std::cout << "; value at q=1: " << sum;
        }
        std::cout << '.' << std::endl;
      }
    }
    catch (std::runtime_error& e) // try again after runtime errors
    { std::cerr << e.what() << std::endl; }
    std::cin.clear(); // clear any previous error condition
    while(std::cin.peek()!=EOF and std::cin.get()!='\n') {}
    // skip to the end of the line
  }
  while (true); 
}


