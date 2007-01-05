#include  <iostream>
#include  <vector>
#include  <fstream>
#include  <iomanip>
#include  <sstream>


typedef unsigned long int ulong;

class ChineseBox
{ protected:
  ulong a,b; // the base moduli
  ulong gcd, lcm; // greatest common divisor and least common multiple
  ulong m, mm;
   // multiples of |b| congruent to |gcd| and |-gcd|, respectively, modulo~|a|

  public:
  ChineseBox(ulong a, ulong b);
  virtual ~ChineseBox()  {}

  // accessors
  ulong get_gcd() const { return gcd; }
  ulong get_lcm() const { return lcm; }
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

class PrimeChineseBox : public ChineseBox
{ // no additional data
public: // constructor from basic |ChineseBox|
  PrimeChineseBox(const ChineseBox& cb) : ChineseBox(cb)
  { if (gcd!=1)
    { std::cerr << "Non relatively prime numbers, gcd=" << gcd << ".\n";
      exit(1); // we should never get here
    }
  }
  virtual ~PrimeChineseBox()  {}

  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

class TabledChineseBox : public ChineseBox
{
protected:
  std::vector<ulong> m_table;
  std::vector<ulong>::const_reverse_iterator mm_table;
public: // constructor from basic |ChineseBox|
  TabledChineseBox(const ChineseBox& cb);
  virtual ~TabledChineseBox()  {}  // destroys table

  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

class PrimeTabledChineseBox : public TabledChineseBox
{ // no additional data
public: // constructor from basic |ChineseBox|
  PrimeTabledChineseBox(const ChineseBox& cb) : TabledChineseBox(cb)
  { if (gcd!=1)
    { std::cerr << "Non relatively prime numbers, gcd=" << gcd << ".\n";
      exit(1); // we should never get here
    }
  }
  virtual ~PrimeTabledChineseBox()  {}

  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

class DoubleTabledChineseBox : public ChineseBox
{
protected:
  std::vector<ulong> m_table,mm_table;
public: // constructor from basic |ChineseBox|
  DoubleTabledChineseBox(const ChineseBox& cb);
  virtual ~DoubleTabledChineseBox()  {}  // destroys tables

  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

class PrimeDoubleTabledChineseBox : public DoubleTabledChineseBox
{ // no additional data
public: // constructor from basic |ChineseBox|
  PrimeDoubleTabledChineseBox(const ChineseBox& cb)
  : DoubleTabledChineseBox(cb)
  { if (gcd!=1)
    { std::cerr << "Non relatively prime numbers, gcd=" << gcd << ".\n";
      exit(1); // we should never get here
    }
  }
  virtual ~PrimeDoubleTabledChineseBox()  {}

  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

class modulus_info
{ ulong modulus;
  ulong nr_polynomials; // number of polynomials for this modulus
  std::streamoff index_begin;
  std::streamoff coefficients_begin;
  ulong nr_coefficients;
  ulong coefficient_size; // 1 for original files, but may be more
  std::vector<unsigned int> renumber;
  bool using_renumber;
  std::ifstream& coefficient_file; // owned file reference

public:
  modulus_info(ulong mod, std::ifstream* ren_file, std::ifstream* coef_file);
  ~modulus_info();

  ulong length(ulong i) const; // length (degree+1) of polynomial |i|
  std::vector<ulong> coefficients(ulong i) const;
    // coefficients of polynomial |i|
  ulong nr_pol() const { return nr_polynomials; }
};

const std::ios_base::openmode binary_out=
			    std::ios_base::out
			  | std::ios_base::trunc
			  | std::ios_base::binary;

const std::ios_base::openmode binary_in=
			    std::ios_base::in
			  | std::ios_base::binary;

ulong extended_gcd(ulong a,ulong b, ulong&lcm, ulong& m)
{ ulong d0=a, m0=0, d1=b,m1=b;
  while(d0!=0)
  // invariant: $d_0\cong-m_0\pmod{a}$, \ and $d_1\cong+m_1\pmod{a}$
  { m1+=(d1/d0)*m0; d1%=d0; // i.e., |d1-=(d1/d0)*d0|
    if (d1==0) break;
    // invariant holds, and |d0>d1>0|
    m0+=(d0/d1)*m1; d0%=d1; // i.e., |d0-=(d0/d1)*d1|
  } // invariant holds, and |0<d0<d1|

  if (d1==0) { lcm=m1; m=m1-m0; return d0; }
  else { lcm=m0; m=m1; return d1; } // |d0==0|
}

ChineseBox::ChineseBox(ulong aa, ulong bb) : a(aa),b(bb)
{ gcd=extended_gcd(a,b,lcm,m); mm=lcm-m; }

ulong ChineseBox::lift_remainders(ulong s, ulong t) const
{ if ((s<=t ? t-s : s-t)%gcd==0)
     return (t+ (s>=t ? (s-t)/gcd*m : (t-s)/gcd*mm))%lcm;
  std::cerr << "Incompatible remainders " 
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

ulong PrimeChineseBox::lift_remainders(ulong s, ulong t) const
{ return (t+ (s>=t ? (s-t)*m : (t-s)*mm))%lcm; }

TabledChineseBox::TabledChineseBox(const ChineseBox& cb)
: ChineseBox(cb), m_table(a/gcd+1), mm_table(m_table.rbegin())
{ ulong last = m_table[0]=0;
  for (ulong i=1; i<m_table.size(); ++i)
    m_table[i]= last<mm ? last+=m : last-=mm;
  assert(last==0);
}

ulong TabledChineseBox::lift_remainders(ulong s, ulong t) const
{ if ((s>=t ? s-t : t-s)%gcd==0)
  { ulong d=s>=t ? m_table[(s-t)/gcd] : mm_table[(t-s)%a/gcd];
    return t<lcm-d ? t+d : t-(lcm-d);
  }
  std::cerr << "Incompatible remainders " 
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

ulong PrimeTabledChineseBox::lift_remainders(ulong s, ulong t) const
{ ulong d=s>=t ? m_table[s-t] : mm_table[(t-s)%a];
return t<lcm-d ? t+d : t-(lcm-d);
}

DoubleTabledChineseBox::DoubleTabledChineseBox(const ChineseBox& cb)
: ChineseBox(cb), m_table(a/gcd), mm_table(b/gcd)
{ ulong last = m_table[0]=0;
  for (ulong i=1; i<m_table.size(); ++i)
    m_table[i]= last<mm ? last+=m : last-=mm;
  assert(last==mm);
  last = mm_table[0]=0;
  for (ulong i=1; i<mm_table.size(); ++i)
    mm_table[i]= last<m ? last+=mm : last-=m;
}

ulong DoubleTabledChineseBox::lift_remainders(ulong s, ulong t) const
{ if ((s>=t ? s-t : t-s)%gcd==0)
  { ulong d=s>=t ? m_table[(s-t)/gcd] : mm_table[(t-s)/gcd];
    return t<lcm-d ? t+d : t-(lcm-d);
  }
  std::cerr << "Incompatible remainders " 
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

ulong PrimeDoubleTabledChineseBox::lift_remainders(ulong s, ulong t) const
{ ulong d=s>=t ? m_table[s-t] : mm_table[t-s];
return t<lcm-d ? t+d : t-(lcm-d);
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
  : modulus(mod), nr_polynomials(), index_begin(), coefficients_begin()
  , coefficient_size(1), renumber(), using_renumber(false)
  , coefficient_file(*coef_file)
{ coefficient_file.seekg(0,std::ios_base::beg); // begin at the beginning
  nr_polynomials=read_bytes(4,coefficient_file);
  index_begin=coefficient_file.tellg();
  coefficient_file.seekg(5*nr_polynomials,std::ios_base::cur);
  nr_coefficients =read_bytes(5,coefficient_file);
  coefficients_begin=coefficient_file.tellg();

  --mod; // make largest remainder
  while ((mod>>=8) != 0)
    ++coefficient_size; // compute number of bytes required
  if (ren_file!=NULL)
  { using_renumber=true; read_renumbering_table(*ren_file,renumber); }
  delete ren_file; // close file when table is read in
}

modulus_info::~modulus_info() { delete &coefficient_file;}

ulong modulus_info::length (ulong i) const
{ if (using_renumber) i=renumber[i];
  coefficient_file.seekg(index_begin+5*i,std::ios_base::beg);
    // locate index in file
  ulong index=read_bytes(5,coefficient_file);
  ulong next_index=read_bytes(5,coefficient_file);
  return (next_index-index)/coefficient_size;
}

std::vector<ulong> modulus_info::coefficients (ulong i) const
{ if (using_renumber) i=renumber[i];
  coefficient_file.seekg(index_begin+5*i,std::ios_base::beg);
  ulong index=read_bytes(5,coefficient_file);
  ulong next_index=read_bytes(5,coefficient_file);

  coefficient_file.seekg(coefficients_begin+index,std::ios_base::beg);
  std::vector<ulong> result ((next_index-index)/coefficient_size);
  for (ulong i=0; i<result.size(); ++i)
    result[i]=read_bytes(coefficient_size,coefficient_file);
  return result;
}

ulong write_indices
 (ulong coefficient_size,
  const std::vector<modulus_info*>& mod_info,
  std::ostream& out,
  bool verbose)
// return value is size of (yet unwritten) coefficient part of output file
{ ulong nr_pols=mod_info[0]->nr_pol(); // number of new polynomials
  write_bytes(nr_pols,4,out);
  for (ulong j=1; j<mod_info.size(); ++j)
    if (mod_info[j]->nr_pol()!=nr_pols)
    { std::cerr << "Conflicting numbers of polynomials in renumbering files: "
	      << nr_pols << "!=" << mod_info[j]->nr_pol()
		<< " (modulus nrs O, " << j << ").\n";
      exit(1);
    }

  ulong index=0; // index of current polynomial to be written
  for (ulong i=0; i<nr_pols; ++i)
  { if (verbose and (i&0xFFF)==0)
      std::cerr << "Polynomial: " << std::setw(10) << i << '\r';
    ulong len=0; // maximum of degree+1 of polynomials selected
    for (ulong j=0; j<mod_info.size(); ++j)
    { ulong new_len = mod_info[j]->length(i);
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
  const std::vector<modulus_info*>& mod_info,
  const std::vector<ChineseBox*>& box,
  std::ostream& out,
  bool verbose)
// return value is maximum of lifted coefficients
{ ulong nr_pols=mod_info[0]->nr_pol();
  ulong n=mod_info.size(); // number of moduli
  ulong max=0;
  std::vector<ulong> rem(2*n-1);
      // remainders for |n| original and |n-1| derived moduli

  for (ulong i=0; i<nr_pols; ++i)
  { if (verbose and (i&0xFFF)==0)
      std::cerr << "Polynomial: " << std::setw(10) << i << '\r';
    ulong len=0; // maximum of degree+1 of polynomials selected
    std::vector<std::vector<ulong> > modular_pol;

    
    for (ulong j=0; j<mod_info.size(); ++j)
    { std::vector<ulong> p=mod_info[j]->coefficients(i);
      modular_pol.push_back(p);
      if (modular_pol.back().size()>len) len=modular_pol.back().size();
    }
    std::vector<ulong> lifted_pol(len);
    // the polynomial modulo the lcm of the moduli

    for (ulong d=0; d<len; ++d)
    { for (ulong j=0; j<n; ++j) // install original remainders
        rem[j]= d>=modular_pol[j].size() ? 0 : modular_pol[j][d];
       try
       { for (ulong j=0; j<n-1; ++j)
           rem[n+j]=box[j]->lift_remainders(rem[2*j],rem[2*j+1]);
           // it happens here!
         ulong c=rem.back();
	 
	 if (c>max)
	 { max=c;
	   std::cerr << (verbose ? "\t\t\tm" : "M") << "aximal coefficient so far: " 
	             << max << ", in polynomial " << i << '\r';
	 }
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

void test(std::vector<ulong>& moduli,std::vector<ChineseBox*>& box)
{ ulong n=moduli.size();
  ulong lcm=box.back()->get_lcm();

  std::vector<ulong> remainder(2*n-1);
  // remainders for |n| original and |n-1| derived moduli

  while(true)
  { std::cout << "Give remainders mod " << moduli[0];
    for (ulong i=1; i<n; ++i) std::cout << ", " << moduli[i];
    std::cout << ": ";
    for (ulong i=0; i<n; ++i)
    { remainder[i]=~0ul; std::cin >> remainder[i];
      if (remainder[i]==~0ul)
      { std::cout << "Bye.\n"; exit(0); }
      remainder[i]%=moduli[i];
    }

    try
    { for (ulong i=0; i<n-1; ++i)
      { remainder[n+i]=
          box[i]->lift_remainders(remainder[2*i],remainder[2*i+1]);
        std::cout << "Lifted to " << remainder[n+i]
		  << " (mod " << box[i]->get_lcm() << ").\n";
      }
      std::cout << "Solution found: " << remainder.back()
		<< " (mod " << lcm << ").\n";
      for (ulong i=0; i<n; ++i)
        if (remainder.back()%moduli[i]!=remainder[i])
	{ std::cout << "Solution does not pass test.\n"; throw true; }
      std::cout << "Solution checks correctly.\n";
    }
    catch (bool b)
    { if (not b) std::cout << "No solution.\n"; }
  }

}

int main(int argc, char** argv)
{ --argc; ++argv; // skip program name
  bool verbose=true, double_tables=false;
  if (argc>0 and std::string(*argv)=="-q") { verbose=false; --argc; ++argv;}
  if (argc>0 and std::string(*argv)=="-d")
    { double_tables=true; --argc; ++argv;}

  std::string mat_base,coef_base;
  // base names for renumbering and coefficient files
  std::vector<ulong> moduli;
  bool interactive= argc<1;

  
  if (interactive) argc=0; // ignore
  else  { mat_base=*argv++; --argc; coef_base=*argv++; --argc; }
  
  if (argc==0) 
             { std::cout << "Give moduli used, or 0 to terminate.\n";
               while(true)
               { std::cout << "Modulus: ";
                 ulong m=0; std::cin >> m;
                 if (m==0) break;
                 moduli.push_back(m);
               }
             }
  else while (argc-->0)
  { std::istringstream in(*argv++);
    unsigned int m=0; in >> m;
    if (m!=0) moduli.push_back(m);
    else
    { std::cerr << "Illegal modulus argument: " << in.str() << "\n";
       exit(1);
    }
  }
  if (moduli.size()<1) { std::cerr << "Too few moduli.\n"; exit(1); }

  ulong n=moduli.size();
  std::vector<ChineseBox*> box(n-1,NULL);
  ulong lcm;
  
  { std::vector<ulong> mod(moduli);
    for (ulong i=0; i<n-1; ++i)
    { ChineseBox b(mod[2*i],mod[2*i+1]);
      mod.push_back(b.get_lcm()); // defines |mod[n+i]|
      if (double_tables)
        box[i]= b.get_gcd()!=1 ? new DoubleTabledChineseBox(b)
  			     : new PrimeDoubleTabledChineseBox(b);
      else
        box[i]= b.get_gcd()!=1 ? new TabledChineseBox(b)
  			     : new PrimeTabledChineseBox(b);
    }
    lcm=mod.back();
  }

  ulong coefficient_size=1, rem=lcm-1; // maximal remainder
  while ((rem>>=8)!=0) ++coefficient_size;

  std::vector<modulus_info*> mod_info;
   // among other things this will hold the input files
  std::ofstream coefficient_file;

  try
  { if (interactive) test(moduli,box); // does not return

    
    { for (ulong i=0; i<moduli.size(); ++i)
      { std::ostringstream name0,name1;
        name0 << mat_base << "-renumbering-mod" << moduli[i];
        std::ifstream* renumber_file=
          new std::ifstream(name0.str().c_str(),binary_in);
        if (not renumber_file->is_open())
        { std::cout << "Assuming canonical order for modulus " << moduli[i]
                    << ".\n";
          renumber_file=NULL; // signal absence of renumbering file
        }
    
        name1 << coef_base << "-mod" << moduli[i];
        std::ifstream* coef_file=new std::ifstream(name1.str().c_str(),binary_in);
        if (not coef_file->is_open())
        { std::cerr << "Could not open file '" << name1.str() << "'.\n";
          exit(1);
        }
        mod_info.push_back(new modulus_info(moduli[i],renumber_file,coef_file));
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
    ulong nr_c=
      write_indices(coefficient_size,mod_info,coefficient_file,verbose);
    std::cout << "\nDone writing indices, will now write "
              << nr_c << " coefficient bytes.\n";
    ulong max_coef=
       write_coefficients
         (coefficient_size,mod_info,box,coefficient_file,verbose);
    std::cout << "\nMaximal coefficient found: "
              << max_coef << ".\n";
  }
  catch (...)
  { for (ulong i=0; i<mod_info.size(); ++i) delete mod_info[i];
    for (ulong i=0; i<box.size(); ++i) delete box[i];
    throw;
  }
}


