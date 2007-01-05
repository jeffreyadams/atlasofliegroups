#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
/*1:*//*2:*/
#line 61 "coef-merge.w"
typedef unsigned long int ulong;/*:2*//*4:*/
#line 125 "coef-merge.w"
class ChineseBox
{protected:
ulong a,b;
ulong gcd,lcm;
ulong m,mm;


public:
ChineseBox(ulong a,ulong b);
virtual~ChineseBox(){}


ulong get_gcd()const{return gcd;}
ulong get_lcm()const{return lcm;}
virtual ulong lift_remainders(ulong s,ulong t)const;
};/*:4*//*7:*/
#line 168 "coef-merge.w"
class PrimeChineseBox:public ChineseBox
{
public:
PrimeChineseBox(const ChineseBox&cb):ChineseBox(cb)
{if(gcd!=1)
{std::cerr<<"Non relatively prime numbers, gcd="<<gcd<<".\n";
exit(1);
}
}
virtual~PrimeChineseBox(){}


virtual ulong lift_remainders(ulong s,ulong t)const;
};/*:7*//*9:*/
#line 212 "coef-merge.w"
class TabledChineseBox:public ChineseBox
{
protected:
std::vector<ulong>m_table;
std::vector<ulong>::const_reverse_iterator mm_table;
public:
TabledChineseBox(const ChineseBox&cb);
virtual~TabledChineseBox(){}


virtual ulong lift_remainders(ulong s,ulong t)const;
};/*:9*//*12:*/
#line 279 "coef-merge.w"
class PrimeTabledChineseBox:public TabledChineseBox
{
public:
PrimeTabledChineseBox(const ChineseBox&cb):TabledChineseBox(cb)
{if(gcd!=1)
{std::cerr<<"Non relatively prime numbers, gcd="<<gcd<<".\n";
exit(1);
}
}
virtual~PrimeTabledChineseBox(){}


virtual ulong lift_remainders(ulong s,ulong t)const;
};/*:12*//*16:*/
#line 353 "coef-merge.w"
class modulus_info
{ulong modulus;
ulong nr_polynomials;
std::streamoff index_begin;
std::streamoff coefficients_begin;
ulong nr_coefficients;
ulong coefficient_size;
std::vector<unsigned int>renumber;
std::ifstream&coefficient_file;

public:
modulus_info(ulong mod,std::ifstream*ren_file,std::ifstream*coef_file);
~modulus_info();

ulong length(ulong i)const;
std::vector<ulong>coefficients(ulong i)const;

const std::vector<unsigned int> &renumber_vector()const
{return renumber;}
};/*:16*//*28:*/
#line 666 "coef-merge.w"
const std::ios_base::openmode binary_out=
std::ios_base::out
|std::ios_base::trunc
|std::ios_base::binary;

const std::ios_base::openmode binary_in=
std::ios_base::in
|std::ios_base::binary;/*:28*//*3:*/
#line 93 "coef-merge.w"
ulong extended_gcd(ulong a,ulong b,ulong&lcm,ulong&m)
{ulong d0=a,m0=0,d1=b,m1=b;
while(d0!=0)

{m1+=(d1/d0)*m0;d1%=d0;
if(d1==0)break;

m0+=(d0/d1)*m1;d0%=d1;
}

if(d1==0){lcm=m1;m=m1-m0;return d0;}
else{lcm=m0;m=m1;return d1;}
}/*:3*//*5:*/
#line 144 "coef-merge.w"
ChineseBox::ChineseBox(ulong aa,ulong bb):a(aa),b(bb)
{gcd=extended_gcd(a,b,lcm,m);mm=lcm-m;}/*:5*//*6:*/
#line 153 "coef-merge.w"
ulong ChineseBox::lift_remainders(ulong s,ulong t)const
{if((s<=t?t-s:s-t)%gcd==0)
return(t+(s>=t?(s-t)/gcd*m:(t-s)/gcd*mm))%lcm;
std::cerr<<"Incompatible remainders "
<<s<<" (mod "<<a<<") and "
<<t<<" (mod "<<b<<").\n";
throw false;
}/*:6*//*8:*/
#line 187 "coef-merge.w"
ulong PrimeChineseBox::lift_remainders(ulong s,ulong t)const
{return(t+(s>=t?(s-t)*m:(t-s)*mm))%lcm;}/*:8*//*10:*/
#line 237 "coef-merge.w"
TabledChineseBox::TabledChineseBox(const ChineseBox&cb)
:ChineseBox(cb),m_table(a/gcd+1),mm_table(m_table.rbegin())
{ulong last=m_table[0]=0;
for(ulong i=1;i<m_table.size();++i)
m_table[i]=last<mm?last+=m:last-=mm;
assert(last==0);
}/*:10*//*11:*/
#line 263 "coef-merge.w"
ulong TabledChineseBox::lift_remainders(ulong s,ulong t)const
{if((s>=t?s-t:t-s)%gcd==0)
{ulong d=s>=t?m_table[(s-t)/gcd]:mm_table[(t-s)%a/gcd];
return t<lcm-d?t+d:t-(lcm-d);
}
std::cerr<<"Incompatible remainders "
<<s<<" (mod "<<a<<") and "
<<t<<" (mod "<<b<<").\n";
throw false;
}/*:11*//*13:*/
#line 297 "coef-merge.w"
ulong PrimeTabledChineseBox::lift_remainders(ulong s,ulong t)const
{ulong d=s>=t?m_table[s-t]:mm_table[(t-s)%a];
return t<lcm-d?t+d:t-(lcm-d);
}/*:13*//*14:*/
#line 315 "coef-merge.w"
ulong read_bytes(ulong n,std::istream&in)
{
if(n==0)return 0;
char c;in.get(c);unsigned char low=c;
return low+(read_bytes(n-1,in)<<8);
}

inline void write_bytes(ulong val,ulong n,std::ostream&out)
{while(n-- >1){out.put(char(val&0xFF));val>>=8;}
out.put(char(val));
}/*:14*//*15:*/
#line 338 "coef-merge.w"
void read_renumbering_table
(std::ifstream&in,std::vector<unsigned int> &table)
{in.seekg(0,std::ios_base::end);
ulong file_size=in.tellg();
in.seekg(0,std::ios_base::beg);

table.resize(file_size/4);
for(ulong i=0;i<table.size();++i)table[i]=read_bytes(4,in);
}/*:15*//*17:*/
#line 381 "coef-merge.w"
modulus_info::modulus_info
(ulong mod,std::ifstream*ren_file,std::ifstream*coef_file)
:modulus(mod),coefficient_size(1),coefficient_file(*coef_file)
{coefficient_file.seekg(0,std::ios_base::beg);
nr_polynomials=read_bytes(4,coefficient_file);
index_begin=coefficient_file.tellg();
coefficient_file.seekg(5*nr_polynomials,std::ios_base::cur);
nr_coefficients=read_bytes(5,coefficient_file);
coefficients_begin=coefficient_file.tellg();

--mod;
while((mod>>=8)!=0)
++coefficient_size;
read_renumbering_table(*ren_file,renumber);
delete ren_file;
}/*:17*//*18:*/
#line 416 "coef-merge.w"
modulus_info::~modulus_info(){delete&coefficient_file;}/*:18*//*19:*/
#line 422 "coef-merge.w"
ulong modulus_info::length(ulong i)const
{coefficient_file.seekg(index_begin+5*renumber[i],std::ios_base::beg);

ulong index=read_bytes(5,coefficient_file);
ulong next_index=read_bytes(5,coefficient_file);
return(next_index-index)/coefficient_size;
}/*:19*//*20:*/
#line 435 "coef-merge.w"
std::vector<ulong>modulus_info::coefficients(ulong i)const
{coefficient_file.seekg(index_begin+5*renumber[i],std::ios_base::beg);
ulong index=read_bytes(5,coefficient_file);
ulong next_index=read_bytes(5,coefficient_file);

coefficient_file.seekg(coefficients_begin+index,std::ios_base::beg);
std::vector<ulong>result((next_index-index)/coefficient_size);
for(ulong i=0;i<result.size();++i)
result[i]=read_bytes(coefficient_size,coefficient_file);
return result;
}/*:20*//*21:*/
#line 458 "coef-merge.w"
ulong write_indices
(ulong coefficient_size,
const std::vector<modulus_info* > &mod_info,
std::ostream&out,
bool verbose)

{ulong nr_pol=mod_info[0]->renumber_vector().size();

write_bytes(nr_pol,4,out);
for(ulong j=1;j<mod_info.size();++j)
if(mod_info[j]->renumber_vector().size()!=nr_pol)
{std::cerr<<"Conflicting numbers of polynomials in renumbering files: "
<<nr_pol<<"!="<<mod_info[j]->renumber_vector().size()
<<" (modulus nrs O, "<<j<<").\n";
exit(1);
}

ulong index=0;
for(ulong i=0;i<nr_pol;++i)
{if(verbose and(i&0xFFF)==0)
std::cerr<<"Polynomial: "<<std::setw(10)<<i<<'\r';
ulong len=0;
for(ulong j=0;j<mod_info.size();++j)
{ulong new_len=mod_info[j]->length(i);
if(new_len>len)len=new_len;
}
write_bytes(index,5,out);
index+=len*coefficient_size;

}
write_bytes(index,5,out);
return index;
}/*:21*//*22:*/
#line 502 "coef-merge.w"
ulong write_coefficients
(ulong coefficient_size,
const std::vector<modulus_info* > &mod_info,
const std::vector<ChineseBox* > &box,
std::ostream&out,
bool verbose)

{ulong nr_pol=mod_info[0]->renumber_vector().size();
ulong n=mod_info.size();
ulong max=0;
std::vector<ulong>rem(2*n-1);


for(ulong i=0;i<nr_pol;++i)
{if(verbose and(i&0xFFF)==0)
std::cerr<<"Polynomial: "<<std::setw(10)<<i<<'\r';
ulong len=0;
std::vector<std::vector<ulong> >modular_pol;/*23:*/
#line 554 "coef-merge.w"
for(ulong j=0;j<mod_info.size();++j)
{std::vector<ulong>p=mod_info[j]->coefficients(i);
modular_pol.push_back(p);
if(modular_pol.back().size()>len)len=modular_pol.back().size();
}/*:23*/
#line 523 "coef-merge.w"
std::vector<ulong>lifted_pol(len);


for(ulong d=0;d<len;++d)
{for(ulong j=0;j<n;++j)
rem[j]=d>=modular_pol[j].size()?0:modular_pol[j][d];
try
{for(ulong j=0;j<n-1;++j)
rem[n+j]=box[j]->lift_remainders(rem[2*j],rem[2*j+1]);

ulong c=rem.back();/*24:*/
#line 565 "coef-merge.w"
if(c>max)
{max=c;
std::cerr<<(verbose?"\t\t\tm":"M")<<"aximal coefficient so far: "
<<max<<", in polynomial "<<i<<'\r';
}/*:24*/
#line 535 "coef-merge.w"
write_bytes(c,coefficient_size,out);
}
catch(bool)

{std::cerr<<"In coefficient "<<d
<<" of polynomial "<<i<<".\n";
exit(1);
}
}
}
return max;
}/*:22*//*32:*/
#line 744 "coef-merge.w"
void test(std::vector<ulong> &moduli,std::vector<ChineseBox* > &box)
{ulong n=moduli.size();
ulong lcm=box.back()->get_lcm();

std::vector<ulong>remainder(2*n-1);


while(true)
{std::cout<<"Give remainders mod "<<moduli[0];
for(ulong i=1;i<n;++i)std::cout<<", "<<moduli[i];
std::cout<<": ";
for(ulong i=0;i<n;++i)
{remainder[i]=~0ul;std::cin>>remainder[i];
if(remainder[i]==~0ul)
{std::cout<<"Bye.\n";exit(0);}
remainder[i]%=moduli[i];
}

try
{for(ulong i=0;i<n-1;++i)
{remainder[n+i]=
box[i]->lift_remainders(remainder[2*i],remainder[2*i+1]);
std::cout<<"Lifted to "<<remainder[n+i]
<<" (mod "<<box[i]->get_lcm()<<").\n";
}
std::cout<<"Solution found: "<<remainder.back()
<<" (mod "<<lcm<<").\n";
for(ulong i=0;i<n;++i)
if(remainder.back()%moduli[i]!=remainder[i])
{std::cout<<"Solution does not pass test.\n";throw true;}
std::cout<<"Solution checks correctly.\n";
}
catch(bool b)
{if(not b)std::cout<<"No solution.\n";}
}

}/*:32*//*25:*/
#line 581 "coef-merge.w"
int main(int argc,char* *argv)
{--argc;++argv;
bool verbose=true;
if(argc>0 and std::string(*argv)=="-q"){verbose=false;--argc;++argv;}

std::string mat_base,coef_base;

std::vector<ulong>moduli;
bool interactive=argc<1;/*26:*/
#line 632 "coef-merge.w"
if(interactive)argc=0;
else{mat_base= *argv++;--argc;coef_base= *argv++;--argc;}

if(argc==0)/*31:*/
#line 728 "coef-merge.w"
{std::cout<<"Give moduli used, or 0 to terminate.\n";
while(true)
{std::cout<<"Modulus: ";
ulong m=0;std::cin>>m;
if(m==0)break;
moduli.push_back(m);
}
}/*:31*/
#line 636 "coef-merge.w"
else while(argc-- >0)
{std::istringstream in(*argv++);
unsigned int m=0;in>>m;
if(m!=0)moduli.push_back(m);
else
{std::cerr<<"Illegal modulus argument: "<<in.str()<<"\n";
exit(1);
}
}
if(moduli.size()<1){std::cerr<<"Too few moduli.\n";exit(1);}/*:26*/
#line 592 "coef-merge.w"

ulong n=moduli.size();
std::vector<ChineseBox* >box(n-1,NULL);
ulong lcm;/*27:*/
#line 653 "coef-merge.w"
{std::vector<ulong>mod(moduli);
for(ulong i=0;i<n-1;++i)
{ChineseBox b(mod[2*i],mod[2*i+1]);
mod.push_back(b.get_lcm());
box[i]=b.get_gcd()!=1?new TabledChineseBox(b)
:new PrimeTabledChineseBox(b);
}
lcm=mod.back();
}/*:27*/
#line 597 "coef-merge.w"

ulong coefficient_size=1,rem=lcm-1;
while((rem>>=8)!=0)++coefficient_size;

std::vector<modulus_info* >mod_info;

std::ofstream coefficient_file;

try
{if(interactive)test(moduli,box);/*29:*/
#line 678 "coef-merge.w"
{for(ulong i=0;i<moduli.size();++i)
{std::ostringstream name0,name1;
name0<<mat_base<<"-renumbering-mod"<<moduli[i];
std::ifstream*renumber_file=
new std::ifstream(name0.str().c_str(),binary_in);
if(not renumber_file->is_open())
{std::cerr<<"Could not open file '"<<name0.str()<<"'.\n";
exit(1);
}

name1<<coef_base<<"-mod"<<moduli[i];
std::ifstream*coef_file=new std::ifstream(name1.str().c_str(),binary_in);
if(not coef_file->is_open())
{std::cerr<<"Could not open file '"<<name1.str()<<"'.\n";
exit(1);
}
mod_info.push_back(new modulus_info(moduli[i],renumber_file,coef_file));
}

std::ostringstream name;
name<<coef_base<<"-mod"<<lcm;/*30:*/
#line 715 "coef-merge.w"
{bool write_protect=false;
for(ulong i=0;i<moduli.size();++i)
if(lcm==moduli[i])write_protect=true;
if(write_protect)name<<'+';
}/*:30*/
#line 700 "coef-merge.w"
coefficient_file.open(name.str().c_str(),binary_out);
if(coefficient_file.is_open())
std::cout<<"Output to file: "<<name.str()<<'\n';
else
{std::cerr<<"Could not open output file '"<<name.str()<<"'.\n";
exit(1);
}
}/*:29*/
#line 610 "coef-merge.w"
ulong nr_c=
write_indices(coefficient_size,mod_info,coefficient_file,verbose);
std::cout<<"\nDone writing indices, will now write "
<<nr_c<<" coefficient bytes.\n";
ulong max_coef=
write_coefficients
(coefficient_size,mod_info,box,coefficient_file,verbose);
std::cout<<"\nMaximal coefficient found: "
<<max_coef<<".\n";
}
catch(...)
{for(ulong i=0;i<mod_info.size();++i)delete mod_info[i];
for(ulong i=0;i<box.size();++i)delete box[i];
throw;
}
}/*:25*//*:1*/
