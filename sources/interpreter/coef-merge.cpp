#include  <iostream>
#include  <vector>
#include  <fstream>
#include  <sstream>


typedef unsigned long int ulong;

class ChineseBox
{ ulong a,b; // the base moduli
  ulong gcd, lcm; // greatest common divisor and least common multiple
  ulong m, mm;
   // multiples of a the are congruent to |gcd| and |-gcd|, respectively

  public:
  ChineseBox(ulong a, ulong b);

  // accessors
  ulong get_gcd() const { return gcd; }
  ulong get_lcm() const { return lcm; }
  ulong lift_remainders(ulong s, ulong t) const;
};

class modulus_info
{ ulong modulus;
  ulong nr_polynomials; // number of polynomials for this modulus
  std::streamoff index_begin;
  std::streamoff coefficients_begin;
  ulong nr_coefficients;
  ulong coefficient_size; // 1 for original files, but may be more
  std::vector<unsigned int> renumber;
  std::ifstream* coefficient_file; // owned file pointer

public:
  modulus_info(ulong mod, std::ifstream* ren_file, std::ifstream* coef_file);
  ~modulus_info();

  ulong length(ulong i) const; // length (degree+1) of polynomial |i|
  std::vector<ulong> coefficients(ulong i) const;
    // coefficients of polynomial |i|
  const std::vector<unsigned int>& renumber_vector() const
    { return renumber; }
};

const std::ios_base::openmode binary_out=
			    std::ios_base::out
			  | std::ios_base::trunc
			  | std::ios_base::binary;

const std::ios_base::openmode binary_in=
			    std::ios_base::in
			  | std::ios_base::binary;

ulong extended_gcd(ulong a,ulong b, ulong&lcm, ulong& m)
{ ulong d0=a, m0=a, d1=b,m1=0;
  while(d0!=0)
  // invariant $d_0\cong+m_0\pmod{b}$, \ and $d_1\cong-m_1\pmod{b}$
  { m1+=(d1/d0)*m0; d1%=d0; // i.e., |d1-=(d1/d0)*d0|
    if (d1==0) break;
    // invariant holds, and |d0>d1>0|
    m0+=(d0/d1)*m1; d0%=d1; // i.e., |d0-=(d0/d1)*d1|
  } // invariant holds, and |0<d0<d1|

  if (d1==0) { m=m0; lcm=m1; return d0; }
  else { m=m0-m1; lcm=m0; return d1; } // |d0==0|
}

ChineseBox::ChineseBox(ulong aa, ulong bb) : a(aa),b(bb)
{ gcd=extended_gcd(a,b,lcm,m); mm=lcm-m; }

ulong ChineseBox::lift_remainders(ulong s, ulong t) const
{ if (gcd==1)
    return (s+ (s<=t ? (t-s)*m : (s-t)*mm))%lcm;
  if ((s<=t ? t-s : s-t)%gcd==0)
     return (s+ (s<=t ? (t-s)/gcd*m : (s-t)/gcd*mm))%lcm;
  std::cerr << "Incompatible remainders " 
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

ulong read_bytes(ulong n, std::istream& in)
{
  if (n==0) return 0;
  char c; in.get(c);  unsigned char low=c; // make unsigned for arithmetic
  return low+(read_bytes(n-1,in)<<8);
}

inline void write_bytes(ulong val, ulong n, std::ostream& out)
{ while (n-->1) { out.put(char(val&0xFF)); val>>=8; } // write |n-1| bytes
  out.put(char(val)); // and one more
}

void read_renumbering_table
  (std::ifstream& in, std::vector<unsigned int>& table)
{ in.seekg(0,std::ios_base::end);
  ulong file_size = in.tellg();
  in.seekg(0,std::ios_base::beg); // return to start (do not collect \$200)

  table.resize(file_size/4);
  for (ulong i=0; i<table.size(); ++i) table[i]=read_bytes(4,in);
}

modulus_info::modulus_info
  (ulong mod, std::ifstream* ren_file, std::ifstream* coef_file)
  : modulus(mod), coefficient_size(1), coefficient_file(coef_file)
{ coefficient_file->seekg(0,std::ios_base::beg); // begin at the beginning
  nr_polynomials=read_bytes(4,*coefficient_file);
  index_begin=coefficient_file->tellg();
  coefficient_file->seekg(5*nr_polynomials,std::ios_base::cur);
  nr_coefficients =read_bytes(5,*coefficient_file);
  coefficients_begin=coefficient_file->tellg();

  --mod; // make largest remainder
  while ((mod>>=8) != 0) ++coefficient_size;
  read_renumbering_table(*ren_file,renumber);
  delete ren_file; // close file when table is read in
}

modulus_info::~modulus_info() {}

ulong modulus_info::length (ulong i) const
{ coefficient_file->seekg(index_begin+5*renumber[i],std::ios_base::beg);
    // locate index in file
  ulong index=read_bytes(5,*coefficient_file);
  ulong next_index=read_bytes(5,*coefficient_file);
  return (next_index-index)/coefficient_size;
}

std::vector<ulong> modulus_info::coefficients (ulong i) const
{ coefficient_file->seekg(index_begin+5*renumber[i],std::ios_base::beg);
  ulong index=read_bytes(5,*coefficient_file);
  ulong next_index=read_bytes(5,*coefficient_file);

  coefficient_file->seekg(coefficients_begin+index,std::ios_base::beg);
  std::vector<ulong> result ((next_index-index)/coefficient_size);
  for (ulong i=0; i<result.size(); ++i)
    result[i]=read_bytes(coefficient_size,*coefficient_file);
  return result;
}

ulong write_indices
 (ulong coefficient_size,
  const std::vector<modulus_info>& mod_info,
  std::ostream& out)
// return value is size of (yet unwritten) coefficient part of output file
{ ulong nr_pol=mod_info[0].renumber_vector().size();
   // number of new polynomials
  write_bytes(nr_pol,4,out);
  for (ulong j=1; j<mod_info.size(); ++j)
    if (mod_info[j].renumber_vector().size()!=nr_pol)
    { std::cerr << "Conflicting numbers of polynomials in renumbering files: "
	      << nr_pol << "!=" << mod_info[j].renumber_vector().size()
		<< " (modulus nrs O, " << j << ").\n";
      exit(1);
    }

  ulong index=0; // index of current polynomial to be written
  for (ulong i=0; i<nr_pol; ++i)
  { ulong len=0; // maximum of degree+1 of polynomials selected
    for (ulong j=0; j<mod_info.size(); ++j)
    { ulong new_len = mod_info[j].length(i);
      if (new_len>len) len=new_len;
    }
    write_bytes(index,5,out); // write index for polynomial
    index+=len*coefficient_size;
    // and advance by its number of coefficient bytes
  }
  write_bytes(index,5,out); // write final index
  return index;
}

ulong write_coefficients
 (ulong coefficient_size,
  const std::vector<modulus_info>& mod_info,
  const std::vector<ChineseBox>& box,
  std::ostream& out)
// return value is maximum of lifted coefficients
{ ulong nr_pol=mod_info[0].renumber_vector().size();
  ulong max=0;
  for (ulong i=0; i<nr_pol; ++i)
  { ulong len=0; // maximum of degree+1 of polynomials selected
    std::vector<std::vector<ulong> > modular_pol;
    for (ulong j=0; j<mod_info.size(); ++j)
    { std::vector<ulong> p=mod_info[j].coefficients(i);
      modular_pol.push_back(p);
      if (modular_pol.back().size()>len) len=modular_pol.back().size();
    }
    std::vector<ulong> lifted_pol(len);
    // the polynomial modulo the lcm of the moduli

    for (ulong d=0; d<len; ++d)
    { ulong c=modular_pol[0][d]; // coefficient for initial modulus
       try
       { for (ulong j=1; j<mod_info.size(); ++j)
         { ulong new_c= d>=modular_pol[j].size() ? 0 : modular_pol[j][d];
           c=box[j-1].lift_remainders(c,new_c); // it happens here!
	 }
	 if (c>max) max=c;
         write_bytes(c, coefficient_size, out);
       }
       catch (bool)
       // incompatibility found during lift; details are already printed
       { std::cerr << "In coefficient " << d
		   << " of polynomial " << i << ".\n";
         exit(1);
       }
    }
  }
  return max;
}

void test()
{ std::vector<unsigned int> moduli;
  std::cout << "Give moduli used, or 0 to terminate.\n";
  while(true)
  { std::cout << "Modulus: ";
    ulong m=0; std::cin >> m;
    if (m==0) break;
    moduli.push_back(m);
  }

  if (moduli.size()<2)
  { std::cerr << "Too few moduli.\n"; exit(1);
  }
  std::vector<ChineseBox> box(1,ChineseBox(moduli[0],moduli[1]));
  for (ulong i=2; i<moduli.size(); ++i)
    box.push_back(ChineseBox(box[i-2].get_lcm(),moduli[i]));
    // defines |box[i-1]|
  ulong lcm=box.back().get_lcm();


  
  std::vector<ulong> remainder(moduli.size());
  while(true)
  { std::cout << "Give remainders mod " << moduli[0];
    for (ulong i=1; i<moduli.size(); ++i) std::cout << ", " << moduli[i];
    std::cout << ": ";
    for (ulong i=0; i<moduli.size(); ++i)
    { remainder[i]=~0ul; std::cin >> remainder[i];
      if (remainder[i]==~0ul)
      { std::cout << "Bye.\n"; exit(0); }
    }
  
    ulong rem=remainder[0];
    try
    { for (ulong i=1; i<moduli.size(); ++i)
        rem=box[i-1].lift_remainders(rem,remainder[i]);
       std::cout << "Solution: " << rem << " (mod " << lcm << ").\n";
    }
    catch (bool)
    { std::cout << "No solution.\n"; }
  }
}

int main(int argc, char** argv)
{ --argc; ++argv; // skip program name
  std::string mat_base,coef_base;
  // base names for renumbering and coefficient files

  if (argc>=2) { mat_base=*argv++; --argc; coef_base=*argv++; --argc; }
  else test();
  // if two names are not given, go into interactive mode

  if (argc<2) { std::cerr<< "Too few moduli"; exit(1); }
  std::vector<ulong> moduli;
  
  while (argc-->0)
  { std::istringstream in(*argv++);
    unsigned int m=0; in >> m;
    if (m!=0) moduli.push_back(m);
    else
    { std::cerr << "Illegal modulus argument: " << in.str() << "\n";
       exit(1);
    }
  }

  std::vector<ChineseBox> box(1,ChineseBox(moduli[0],moduli[1]));
  for (ulong i=2; i<moduli.size(); ++i)
    box.push_back(ChineseBox(box[i-2].get_lcm(),moduli[i]));
    // defines |box[i-1]|
  ulong lcm=box.back().get_lcm();

  ulong coefficient_size=1, rem=lcm-1; // maximal remainder
  while ((rem>>=8)!=0) ++coefficient_size;

  std::vector<modulus_info> mod_info;
   // among other things this will hold the input files
  std::ofstream coefficient_file;

{ for (ulong i=0; i<moduli.size(); ++i)
  { std::ostringstream name0,name1;
    name0 << mat_base << "-renumbering-mod" << moduli[i];
    std::ifstream* renumber_file=
      new std::ifstream(name0.str().c_str(),binary_in);
    if (not renumber_file->is_open())
    { std::cerr << "Could not open file '" << name0.str() << "'.\n";
      exit(1);
    }

    name1 << coef_base << "-mod" << moduli[i];
    std::ifstream* coef_file=new std::ifstream(name1.str().c_str(),binary_in);
    if (not coef_file->is_open())
    { std::cerr << "Could not open file '" << name1.str() << "'.\n";
      exit(1);
    }
    mod_info.push_back(modulus_info(moduli[i],renumber_file,coef_file));
  }

  std::ostringstream name;
  name << coef_base << "-mod" << lcm;
  
  { bool write_protect=false;
    for (ulong i=0; i<moduli.size() ; ++i)
      if (lcm==moduli[i]) write_protect=true;
    if (write_protect) name << '+'; // avoid overwriting file for one modulus
  }
  coefficient_file.open(name.str().c_str(),binary_out);
  if (coefficient_file.is_open())
    std::cout << "Output to file: " << name.str() << '\n';
  else
  { std::cerr << "Could not open output file '" << name.str() << "'.\n";
      exit(1);
    }
}

  ulong nr_c=write_indices(coefficient_size,mod_info,coefficient_file);
  std::cout << "Done writing indices, will now write "
            << nr_c << " coefficient bytes.\n";
  ulong max_coef=
     write_coefficients(coefficient_size,mod_info,box,coefficient_file);
  std::cout << "Done!\nMaximal coefficient found: "
            << max_coef << ".\n";
}


