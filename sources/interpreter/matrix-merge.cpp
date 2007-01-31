#include  <vector>
#include  "../utilities/hashtable.h"
#include  <iostream>
#include  <stdexcept>
#include  "../utilities/bitmap.h"
#include  <string>
#include  <fstream>
#include  <sstream>
#include  "../utilities/arithmetic.h"
#include  <iomanip>


template<unsigned int n>
  class tuple_entry
  { unsigned int comp[n];
  public:
    tuple_entry() // default constructor builds 0
    { for (unsigned int i=0; i<n; ++i) comp[i]=0; }
    tuple_entry(const std::vector<unsigned int>& comps)
    { for (unsigned int i=0; i<n; ++i) comp[i]=comps[i]; }
    tuple_entry(const tuple_entry& other) // copy constructor
    { for (unsigned int i=0; i<n; ++i) comp[i]=other.comp[i]; }
  
    typedef std::vector<tuple_entry> Pooltype;
    size_t hashCode(size_t modulus) const;
    bool operator!= (const tuple_entry& other) const
    { for (unsigned int i=0; i<n; ++i)
        if (comp[i]!=other.comp[i]) return true;
      return false;
    }
  
    unsigned int operator[] (unsigned int i) const { return comp[i]; }
  };

typedef std::vector<std::pair<unsigned int,unsigned int> > coord_vector;

const std::ios_base::openmode binary_out=
			    std::ios_base::out
			  | std::ios_base::trunc
			  | std::ios_base::binary;

const std::ios_base::openmode binary_in=
			    std::ios_base::in
			  | std::ios_base::binary;

template<unsigned int n>
  size_t tuple_entry<n>::hashCode(size_t modulus) const
  { size_t h=0;
    for (unsigned int i=0; i<n; ++i) h=((h<<2)+comp[i])&(modulus-1);
    return h;
  }


unsigned int read_int(std::istream& in)
{ char c; unsigned char uc; unsigned int ui,result;
in.get(c); result=uc=c;
in.get(c); ui=uc=c; result += ui<<8;
in.get(c); ui=uc=c; result += ui<<16;
in.get(c); ui=uc=c; return result + (ui<<24);
}

void write_int(unsigned int n, std::ostream& out)
{
  out.put(char(n&0xFF)); n>>=8;
  out.put(char(n&0xFF)); n>>=8;
  out.put(char(n&0xFF)); n>>=8;
  out.put(char(n));
}
template<unsigned int n>
  void combine_rows
    (unsigned int y,
     atlas::hashtable::HashTable<tuple_entry<n>,unsigned int>& hash,
     std::vector<std::istream*>in, std::ostream& out,
     std::vector<unsigned int>& lim,
     coord_vector* first_use)
  { for (unsigned int i=0; i<n; ++i)
      if (read_int(*in[i])!=y) 
                            { std::cerr << "y=" << y << ", i=" << i << ":\n";
                              throw std::runtime_error("Wrong alignment in source file");
                            }
    unsigned int nr_prim=read_int(*in[0]);
    for (unsigned int i=1; i<n; ++i)
      if (read_int(*in[i])!=nr_prim)
        
        { std::cerr << "y=" << y << ", i=" << i << ":\n";
          throw std::runtime_error("Primitive count mismatch in source files");
        }
    write_int(y,out); write_int(nr_prim,out); // reproduce info in output

    
    typedef atlas::bitmap::BitMap bit_map;
    
    std::vector<bit_map> in_map(n,bit_map(nr_prim));
    bit_map out_map(nr_prim);
    
    for (size_t i=0; i<nr_prim; i+=32)
    { unsigned int b=0;
      for (unsigned int j=0; j<n; ++j)
      { unsigned int bj=read_int(*in[j]); in_map[j].setRange(i,32,bj);
        b |= bj;
      }
      out_map.setRange(i,32,b); write_int(b,out);
    }
    
    for (bit_map::iterator it=out_map.begin(); it(); ++it)
    { std::vector<unsigned int> tuple(n,0);
      for (unsigned int j=0; j<n; ++j)
        if (in_map[j].isMember(*it))
        { tuple[j]=read_int(*in[j]);
          if (tuple[j]>=lim[j]) lim[j]=tuple[j]+1;
            // keep track of limit for modular numbers
        }
        else {} // |tuple[j]| stays 0; and nothing is read from |*in[j]|
      tuple_entry<n> e(tuple);
      unsigned int size=hash.size();
      unsigned int code=hash.match(e);
      if (first_use!=NULL and code==size)
        first_use->push_back(std::make_pair(*it,y));
      write_int(code,out);
    }
  }

template<unsigned int n>
void do_work
  (std::string name_base,
   std::vector<unsigned int>& modulus,
   coord_vector* first_use)
{ 
  std::vector<std::ifstream*>in_file(n,NULL);
    std::vector<std::istream*>in_stream(n,NULL);
    for (unsigned int i=0; i<n; ++i)
    { std::ostringstream name;
      name << name_base << "-mod" << modulus[i];
      in_file[i]=new std::ifstream(name.str().c_str(),binary_in);
      if (in_file[i]->is_open())
        in_stream[i]=in_file[i]; // get stream underlying file stream
      else
      { std::cerr << "Could not open file '" << name.str() << "'.\n";
        exit(1);
      }
    }
  
    unsigned long out_modulus=modulus[0];
    for (unsigned int i=1; i<n; ++i)
      
      out_modulus= atlas::arithmetic::lcm(out_modulus, modulus[i]);
  
    std::ostringstream name;
    name << name_base << "-mod" << out_modulus;
    
    { bool write_protect=false;
      for (ulong i=0; i<n ; ++i)
        if (out_modulus==modulus[i]) write_protect=true;
      if (write_protect) name << '+'; // avoid overwriting file for one modulus
    }
    std::ofstream out_file(name.str().c_str(),binary_out);
    if (out_file.is_open())
      std::cout << "Output to file: " << name.str() << '\n';
    else
    { std::cerr << "Could not open output file '" << name.str() << "'.\n";
        exit(1);
      }

  std::vector<tuple_entry<n> > pool;
  atlas::hashtable::HashTable<tuple_entry<n>,unsigned int> hash(pool);
  hash.match(tuple_entry<n>()); // insert index of Zero, it does not occur!
  std::vector<unsigned int> limits(n,1); // limit of modular sequence numbers

  for (unsigned int y=0; in_stream[0]->peek()!=EOF; ++y)
     // something remains to read
  { std::cerr << y << '\r';
    combine_rows<n>(y,hash,in_stream,out_file,limits,first_use);
  }
  std::cerr << "\ndone!\n";
  for (unsigned int i=0; i<n; ++i) delete in_file[i]; // close files

  
  { std::cout << "Numbers of different polynomials found:\n";
    for (unsigned int i=0; i<n; ++i)
      std::cout << "Mod " << std::setw(10) << modulus[i]
                << ": " << limits[i] << ",\n";
    std::cout << "Mod " << std::setw(10) << out_modulus
              << ": " << hash.size() << ".\n";
  }
  
  for (unsigned int i=0; i<n; ++i)
  { std::ostringstream name;
    name << name_base << "-renumbering-mod" << modulus[i];
    std::ofstream out_file(name.str().c_str(),binary_out);
    if (out_file.is_open())
      std::cout << "Renumbering output to file: " << name.str() << '\n';
    else
      { std::cerr << "Could not open file '" << name.str() << "'.\n";
        exit(1);
      }
  
    for (unsigned int k=0; k<pool.size(); ++k)
      write_int(pool[k][i],out_file);
  }
}

int main(int argc,char** argv)
{ atlas::constants::initConstants(); // don't forget!

  --argc; ++argv; // skip program name

  coord_vector* uses=NULL;
  std::ofstream uses_out;
  if (argc>1 and std::string(*argv)=="-uses-to")
    { uses=new coord_vector;
      uses_out.open(argv[1]);
      if (not uses_out.is_open())
      {
        std::cerr << "Could not open file '" << argv[1] << "' for writing.n";
        exit(1);
      }
      argc-=2; argv+=2;
    }

  std::string base;
  if (argc>0) { base=*argv++; --argc; }
  else
  { std::cout << "File name base (up to '-mod'): " ;
    std::cin >> base;
  }

  std::vector<unsigned int> moduli;

  if (argc==0)
  { std::cout << "Give moduli used, or 0 to terminate.\n";
    while(true)
    { std::cout << "Modulus: ";
      unsigned int m=0; std::cin >> m;
      if (m==0) break;
      moduli.push_back(m);
    }
  }
  else
    while (argc-->0)
    { std::istringstream in(*argv++);
      unsigned int m=0; in >> m;
      if (m!=0) moduli.push_back(m);
      else
      { std::cout << "Illegal modulus argument: " << in.str() << "\n";
	exit(1);
      }
    }


  switch (moduli.size())
  { case 1: do_work<1>(base,moduli,uses); break;
    case 2: do_work<2>(base,moduli,uses); break;
    case 3: do_work<3>(base,moduli,uses); break;
    case 4: do_work<4>(base,moduli,uses); break;
    case 5: do_work<5>(base,moduli,uses); break;
    case 6: do_work<6>(base,moduli,uses); break;
    default: std::cout << "I cannot handle " << moduli.size()
		       << " moduli, sorry.\n";
  }

  std::cerr << "Writing uses of polynomials to file... ";
  if (uses!=NULL)
    for (unsigned int i=0; i<uses->size(); ++i)
      uses_out << i+1 << ": " << (*uses)[i].first
                      << ", " << (*uses)[i].second << ".\n";
  std::cerr << "done.\n";

}


