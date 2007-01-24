#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <bitset>
#include <stdexcept>
#include <cassert>
#include <cctype>

typedef unsigned long long int ullong; // sufficiently long unsigned type


typedef unsigned int BlockElt;
typedef std::bitset<32> RankFlags; // we can go up to rank 32
typedef std::vector<RankFlags> descent_set_vector; // indexed by block element
typedef std::vector<BlockElt> ascent_vector;       // indexed by simple root
typedef std::vector<ascent_vector> ascent_table;   // indexed by block element
typedef std::vector<BlockElt> prim_list;           // list of weak primitives
typedef std::vector<prim_list> prim_table; // effectively indexed by RankFlags

typedef prim_list strong_prim_list; // strongly primitive elements for fixed y

typedef std::vector<bool> prim_bitmap;     // has compact implementation

typedef ullong KLIndex; // must be long long to avoid problems when multiplied
const KLIndex UndefKLIndex=~KLIndex(0);

const std::ios_base::openmode binary_out=
			    std::ios_base::out
			  | std::ios_base::trunc
			  | std::ios_base::binary;

const std::ios_base::openmode binary_in=
			    std::ios_base::in
			  | std::ios_base::binary;

ullong read_bytes(unsigned int n, std::istream& in)
{
  if (n==0) return 0;
  else
    {
      char c; in.get(c); unsigned char low=c;
      return low+(read_bytes(n-1,in)<<8);
    }
}

/* After each variable estimates of size for split E8, on a 32-bit machine.
   Estimates are based on average of 12800 (weakly) primitive elements per
   descent set and half that (6400) of strongly primitive elements per block
   element y
*/
struct block_info
{
  unsigned int rank;
  size_t size;
  unsigned int max_length;
  std::vector<BlockElt> start_length; // of size max_length+2;

  descent_set_vector descent_set;      // 453060*4     =  1812240 bytes :  2 MB
  ascent_table ascents;                // 453060*8*4   = 14497920 bytes : 13 MB

//std::vector<prim_bitmap> prim_map;   // 453060*1600  =724896000 bytes :691 MB
  prim_table primitives_list;          // 256*12800*4  = 13107200 bytes : 12 MB

  block_info(std::istream& in); // constructor reads file
};


const BlockElt UndefBlock= ~BlockElt(0);
const BlockElt noGoodAscent= UndefBlock-1;

BlockElt primitivise(BlockElt x, const RankFlags d, const ascent_table& t)
{
  unsigned int rank=t[0].size(); // any index gives the same value
 start:
  const ascent_vector& tx=t[x];
  for (size_t s=0; s<rank; ++s)
    if (d[s] and tx[s]!=noGoodAscent)
      {
	if (tx[s]==UndefBlock) return UndefBlock;
	x=tx[s]; goto start; // this should raise x, now try another step
      }
  return x; // no raising possible, stop here
}

bool is_primitive(BlockElt x, const RankFlags d, const ascent_table& t)
{
  const ascent_vector& tx=t[x];
  for (size_t s=0; s<tx.size(); ++s)
    if (d[s] and tx[s]!=noGoodAscent) return false;
  return true; // now d[s] implies t[s]==noGoodAscent for all simple roots s
}

void make_prim_table
  (const ascent_table& t, unsigned int rank, prim_table& result)
{
  size_t powerset_rank= 1ul<<rank;
  size_t block_size=t.size();
  result.resize(powerset_rank); // fills result with empty vectors
  for (size_t s=0; s<powerset_rank; ++s)
    {
#ifdef VERBOSE
      std::cerr << "Constructing weakly primitive elements for descent set "
		<< s << '\r';
#endif
      RankFlags d(s);
      for (BlockElt x=0; x<block_size; ++x)
	if (is_primitive(x,d,t)) result[s].push_back(x);
      prim_list(result[s]).swap(result[s]); // reallocate to fit snugly
    }
#ifdef VERBOSE
      std::cerr << "\nDone constructing weakly primitive elements.\n";
#endif
}

block_info::block_info(std::istream& in)
  : rank(), size(), max_length(), start_length()
  , descent_set(), ascents(), primitives_list() // don't initialise yet
{
  in.seekg(0,std::ios_base::beg); // ensure we are reading from the start
  size=read_bytes(4,in);
  rank=read_bytes(1,in);
  max_length=read_bytes(1,in);

#ifdef VERBOSE
  std::cout << "rank=" << rank << ", block size=" << size
	    << ", maximal length=" << max_length << ".\n";
#endif

  start_length.resize(max_length+2);
  start_length[0]=BlockElt(0);
  for (size_t i=1; i<=max_length; ++i) start_length[i]=read_bytes(4,in);
  start_length[max_length+1]=size;

  descent_set.reserve(size);
  for (BlockElt y=0; y<size; ++y)
    descent_set.push_back(RankFlags(read_bytes(4,in)));

  ascents.reserve(size);
  for (BlockElt x=0; x<size; ++x)
    {
      ascents.push_back(ascent_vector());
      ascent_vector& a=ascents.back();
      a.reserve(rank);
      for (size_t s=0; s<rank; ++s)
	a.push_back(read_bytes(4,in));
    }
  make_prim_table(ascents,rank,primitives_list);
}

// after each variable estimates size for split E8, on 32-bit machine
// estimates based on average of 12800 primitive elements per descent set
// and half that (6400) of strongly primitive elements per block element y

class matrix_info
{
  std::ifstream& matrix_file;

  block_info block;

  std::vector<std::streampos> row_pos; // 453060*8     =  3624480 bytes :  3 MB

  BlockElt cur_y;		       // 4 bytes
  strong_prim_list cur_strong_prims;   // 6400*4 = 25600 bytes  (but max 2MB)
  std::streampos cur_row_entries;     // 8 bytes
//prim_bitmap cur_prim;		       // 12800/8 =   1600 bytes (max 56633)
//std::vector<size_t> cur_offset;      // 12800*4  =102400 bytes (but max 2MB)

  matrix_info(const matrix_info&); // copying forbidden
  void set_y(BlockElt y);          // install |cur_y| and dependent data

public:
  matrix_info(std::ifstream* block_file,std::ifstream* m_file);
  ~matrix_info() { delete &matrix_file; }

  BlockElt block_size() const { return block.size; }
  KLIndex find_pol_nr(BlockElt x,BlockElt y,BlockElt& xx);
  // sets y, so not const
};

void matrix_info::set_y(BlockElt y)
{
  if (y==cur_y) { matrix_file.seekg(cur_row_entries); return; }
  cur_y=y;
  const prim_list& weak_prims
    = block.primitives_list[block.descent_set[y].to_ulong()];
  cur_strong_prims.resize(0); // restart fom scratch, but don't destroy space

  matrix_file.seekg(row_pos[y]);
  unsigned int n_prim=read_bytes(4,matrix_file);
  for (size_t i=0; i<n_prim; i+=32)
    {
      unsigned int chunk=read_bytes(4,matrix_file);
      for (size_t j=0; chunk!=0; ++j,chunk>>=1) // and certainly j<32
	if ((chunk&1)!=0) cur_strong_prims.push_back(weak_prims[i+j]);
    }

  // the last bit 1 was for |y|, but we pushed the first element of its length
  // in |weak_prims|; now check that statement, and replace the element by |y|
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

KLIndex matrix_info::find_pol_nr(BlockElt x,BlockElt y,BlockElt& xx)
{
  set_y(y);
  RankFlags d=block.descent_set[y];
  xx=primitivise(x,d,block.ascents);
  if (xx==UndefBlock) return KLIndex(0); // primitivisation copped out
  strong_prim_list::const_iterator it=
    std::lower_bound(cur_strong_prims.begin(),cur_strong_prims.end(),xx);
  if (it==cur_strong_prims.end() or *it!=xx) return KLIndex(0); // not strong

  matrix_file.seekg(4*(it-cur_strong_prims.begin()),std::ios_base::cur);
  return KLIndex(read_bytes(4,matrix_file));
}

/* To prepare for the constructor of matrix_info, we need to sum all the bits
   in a word (sideways addition of bits). The following procedure is
   attributed by Knuth (pre-fascicle 1a to volume 4 of TAOCP, p 11) to D. B.
   Gillies and J. C. P. Miller. The idea is to successively replace groups of
   2, 4, and 8 bits by their sideways sum, and then to use a judicious
   multiplication and extraction to add all the bytes.
*/

const unsigned long int b0= ~(~0ul/3);  // 0xAAAA... ; flags odd bit positions
const unsigned long int b1= ~(~0ul/5);  // 0xCCCC... ; flags bit pos 2,3 (mod4)
const unsigned long int b2= ~(~0ul/17); // 0xF0F0... ; flags bit pos 4-7 (mod8)
const unsigned long int ones= ~0ul/255; // 0x0101... ; flags bytes' low-bits
const unsigned int high_byte_shift=8*(sizeof(unsigned long int)-1);

unsigned int add_bits(unsigned long int x)
{
  x-=(x&b0)>>1;          // replace pairs of bits 10->01 and 11->10
  x=(x&~b1)+((x&b1)>>2); // sideways add 2 groups of pairs of bits to 4-tuples
  x += x>>4;             // sums of 8-tuples are now in lower 4-tuples of bytes
  return (x&~b2)*ones >> high_byte_shift; // add bytes in high byte and extract
}

matrix_info::matrix_info(std::ifstream* block_file,std::ifstream* m_file)
  : matrix_file(*m_file)
  , block(*block_file) // read in block information
  , row_pos(block.size)
  , cur_y(UndefBlock), cur_strong_prims(), cur_row_entries(0)
{
  static const unsigned int ulsize=sizeof(unsigned long int);

#ifdef VERBOSE
      std::cerr << "Starting to scan matrix file by 'rows'.\n";
#endif

  size_t l=0; // length of y
  matrix_file.seekg(0,std::ios_base::beg);
  for (BlockElt y=0; y<block.size; ++y)
    {
#ifdef VERBOSE
      std::cerr << y << '\r';
#endif

      while (y>=block.start_length[l+1])
	{
	  ++l;
#ifdef VERBOSE
	  std::cerr << "length " << l << " starts at y=" << y << std::endl;
#endif
	}
      // now block.start_length[l] <= y < block.start_length[l+1]

      if (read_bytes(4,matrix_file)!=y)
	{ std::cerr << y << std::endl;
	  throw std::runtime_error ("Alignment problem");
	}

      row_pos[y]=matrix_file.tellg(); // record position where row info starts
      unsigned int n_prim=read_bytes(4,matrix_file);

      // test that this number is the number of weak primitive elements for y
      {
	const prim_list& weak_prims
	  = block.primitives_list[block.descent_set[y].to_ulong()];
	prim_list::const_iterator i= // find limit of length<l values
	  std::lower_bound(weak_prims.begin(),weak_prims.end()
			  ,block.start_length[l]);
	// test that number is that of weak_prims of length<l, plus 1 for y
	if (n_prim!=size_t(i-weak_prims.begin())+1)
	{
	  std::cerr << y << ": " << n_prim << "!="
		    << (i-weak_prims.begin())+1
		    << std::endl;
	  throw std::runtime_error ("Primitive count problem");
	}
      }

      size_t count=n_prim/(8*ulsize); // number of unsigned longs to read
      std::streamoff n_strong_prim=0;

      while (count-->0)
	n_strong_prim+=add_bits(read_bytes(ulsize,matrix_file));

      count=(n_prim%(8*ulsize)+31)/32; // maybe some tetrabyte(s) left to read
      while (count-->0)
	n_strong_prim+=add_bits(read_bytes(4,matrix_file));

      matrix_file.seekg(4*n_strong_prim,std::ios_base::cur);
      if (not matrix_file.good())
	{ std::cerr << y << std::endl;
	  throw std::runtime_error ("Premature end of file");
	}
    } // for y

#ifdef VERBOSE
      std::cerr << "\nDone scanning matrix file.\n";
#endif

  delete block_file; // success, we no longer need the block file
}


int main(int argc,char** argv)
{
  --argc; ++argv;

  std::auto_ptr<matrix_info> mi; // guarantees clean-up at end

  if (argc>=3)
    {
      argc-=2; // we will consume two arguments
      std::ifstream* block_file=new std::ifstream;
      std::ifstream* matrix_file=new std::ifstream;
      block_file->open(*argv++,binary_in);
      matrix_file->open(*argv++,binary_in);
      if (block_file->is_open() and matrix_file->is_open())
	mi=std::auto_ptr<matrix_info>(new matrix_info(block_file,matrix_file));
      else
	{
	  std::cerr << "failed to open file '"
		    << argv[block_file->is_open()? -1 : -2 ] << "'.\n";
	  delete block_file; delete matrix_file; exit(1);
	}
    }

  std::ifstream coef_file;
  if (argc==1) coef_file.open(argv[0],binary_in);
  else
    {
      std::string file_name;
      std::cout << "File name: " ;
      std::cin >> file_name;
      coef_file.open(file_name.c_str(),binary_in);
    }
  if (not coef_file.is_open())
    { std::cerr << "Open failed"; exit(1); }

  KLIndex n_polynomials=read_bytes(4,coef_file); // must be 64-bits, for 5*...

  std::streamoff index_begin=coef_file.tellg();
  read_bytes(10,coef_file); // skip initial 2 indices
  unsigned int coefficient_size=read_bytes(5,coef_file); // size of One
  std::cout << "Coefficient size " << coefficient_size << ".\n";

  coef_file.seekg(index_begin+5*n_polynomials,std::ios_base::beg);
  ullong n_coefficients =read_bytes(5,coef_file)/coefficient_size;
  std::streamoff coefficients_begin=coef_file.tellg();

  std::cout << n_polynomials << " polynomials, "
	    << n_coefficients << " coefficients.\n";


  while (true)
    {
      if (mi.get()==NULL)
	std::cout << "index: ";
      else
	std::cout << "give block elements x,y, or polynomial index i as #i: ";

      while(isspace(std::cin.peek())) std::cin.get();

      KLIndex i=UndefKLIndex;
      if (std::cin.peek()=='#' ? std::cin.get(), true : mi.get()==NULL )
	{
	  std::cin >> i;
	  if (i==UndefKLIndex) break;
	  if (i>=n_polynomials)
	    {
	      std::cout << "index too large, limit is " << n_polynomials-1
			<< ".\n";
	      continue;
	    }
	}
      else
	{
	  int c; BlockElt x=UndefBlock,y=UndefBlock,xx;
	  while(ispunct(c=std::cin.peek()) or isspace(c)) std::cin.get();
	  std::cin >> x;
	  if (x==UndefBlock) break;
	  while(ispunct(c=std::cin.peek()) or isspace(c)) std::cin.get();
	  std::cin >> y;
	  if (y==UndefBlock)
	    { std::cout << "failure reading y, try again.\n";
  	      while(isalpha(c=std::cin.peek())) std::cin.get(); // remove text
	      continue;
	    }
	  if (x>=mi->block_size())
	    { std::cout << "first parameter too large, try again.\n";
	      continue;
	    }
	  if (y>=mi->block_size())
	    { std::cout << "second parameter too large, try again.\n";
	      continue;
	    }

	  if (x>y)
	    { std::cout << "Result null by triangularity.\n"; continue; }

	  i=mi->find_pol_nr(x,y,xx);
	  if (xx==UndefBlock)
	    { std::cout << "Result null by raising first argument.\n";
	      continue; }
	  else
	    std::cout << "P_{" << x << ',' << y << "}=P_{" << xx << ',' << y
		      << "}=polynomial #" << i << ':' << std::endl;
	}
      coef_file.seekg(index_begin+5*i,std::ios_base::beg);
      ullong index=read_bytes(5,coef_file);
      ullong next_index=read_bytes(5,coef_file);
      size_t length=(next_index-index)/coefficient_size;
      if (length>=33)
	{
	  std::cout << "Degree found too large: " << length-1
		    << ".\n";
	  continue;
	}
      unsigned long* coefficients = new unsigned long[length];

      coef_file.seekg(coefficients_begin+index,std::ios_base::beg);

      unsigned int sum=0; // sum of coefficients, won't exceed 2^29
      for (size_t i=0; i<length; ++i)
	sum+=coefficients[i]=read_bytes(coefficient_size,coef_file);

      if (not coef_file.good())
	{
	  std::cout << "Input error reading coefficients.\n";
	  continue;
	}

      bool first=true;
      for (size_t i=length; i-->0;)
	if (coefficients[i]!=0)
	  {
	    if (first) first=false; else std::cout << " + ";
	    if (coefficients[i]!=1 or i==0) std::cout << coefficients[i];
	    std::cout << (i==0? "" : "q");
	    if (i>1)
	      if (i<10) std::cout << '^' << i;
	      else std::cout << "^{" << i << '}';
	  }
      if (length==0) std::cout << 0;
      else if (length>1) std::cout << "; value at q=1: " << sum;
      std::cout << '.' << std::endl;

      delete[] coefficients;

    } // while(true)


}
