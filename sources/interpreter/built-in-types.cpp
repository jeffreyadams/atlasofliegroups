#include "built-in-types.h"
#include <stdexcept>
#include <sstream>
#include "lietype.h"
#include "constants.h"
#include <cctype>
#include "basic_io.h"
#include "prerootdata.h"
#include "dynkin.h"
#include "smithnormal.h"
#include "arithmetic.h"
#include "lattice.h"
#include "latticetypes.h"
#include "matrix.h"
#include "tags.h"
#include "tori.h"
#include "setutils.h"
#include "dynkin.h"
#include "poset.h"
#include "prettyprint.h"
#include "blocks.h"
#include "block_io.h"
#include "realredgp_io.h"
#include "kgb.h"
#include "kgb_io.h"
#include "kl.h"
#include "klsupport.h"
#include "kl_io.h"
#include "wgraph.h"
#include "wgraph_io.h"
/*1:*/
#line 14 "built-in-types.w"
namespace atlas{namespace interpreter{
namespace{/*8:*/
#line 88 "built-in-types.w"
std::string num(size_t n)
{std::ostringstream s;s<<n;return s.str();}/*:8*//*10:*/
#line 141 "built-in-types.w"
inline void skip_punctuation(const char* &p)
{while(std::ispunct(*p)||std::isspace(*p))++p;}

void Lie_type_wrapper()throw(std::bad_alloc,std::runtime_error)
{string_value*s=get_string();
std::auto_ptr<Lie_type_value>result(new Lie_type_value);
size_t total_rank=0;
const char*p=s->value.c_str();skip_punctuation(p);
while(std::isalpha(*p))
{char c= *p++;
if(!std::isdigit(*p)){--p;break;}
size_t rank= *p++ -'0';
while(std::isdigit(*p))rank=10*rank+(*p++ -'0');
skip_punctuation(p);
result->add_simple_factor(c,rank);

if((total_rank+=rank)>constants::RANK_MAX)
throw std::runtime_error
("Total rank exceeds implementation limit "+num(constants::RANK_MAX));

}
if(*p!='\0')
throw std::runtime_error
("Error in type string '"+s->value+"' for Lie type");
push_value(result.release());
}/*:10*//*16:*/
#line 226 "built-in-types.w"
void Cartan_matrix_wrapper()
{std::auto_ptr<Lie_type_value>t(get_Lie_type());
std::auto_ptr<latmat_value>
result(new latmat_value(latticetypes::LatticeMatrix()));
prerootdata::cartanMatrix(result->value,t->value);
push_value(result.release());
}/*:16*//*17:*/
#line 242 "built-in-types.w"
void type_of_Cartan_matrix_wrapper()
{std::auto_ptr<latmat_value>m(get_mat());
lietype::LieType t;
dynkin::lieType(t,m->value);
push_value(new Lie_type_value(t));
}/*:17*//*19:*/
#line 275 "built-in-types.w"
void smithBasis(latticetypes::CoeffList&invf,latticetypes::WeightList&b,
const lietype::LieType&lt)

{matrix::initBasis(b,lietype::rank(lt));
latticetypes::WeightList::iterator bp=b.begin();

for(size_t j=0;j<lt.size();++j)

{size_t r=lietype::rank(lt[j]);

if(lietype::type(lt[j])=='T')

invf.insert(invf.end(),r,latticetypes::ZeroCoeff);
else{
latticetypes::LatticeMatrix ms;
prerootdata::cartanMatrix(ms,lt[j]);
ms.transpose();
smithnormal::smithNormal(invf,bp,ms);

if(lietype::type(lt[j])=='D' and r%2==0)

latticetypes::operator+=(bp[r-2],bp[r-1]);
}
bp+=r;
}
}/*:19*//*20:*/
#line 309 "built-in-types.w"
void Smith_Cartan_wrapper()
{std::auto_ptr<Lie_type_value>t(get_Lie_type());
std::auto_ptr<latmat_value>
m(new latmat_value(latticetypes::LatticeMatrix()));
std::auto_ptr<weight_value>inv_factors
(new weight_value(latticetypes::CoeffList(0)));
latticetypes::WeightList b;
smithBasis(inv_factors->value,b,t->value);
latticetypes::LatticeMatrix basis(b);
m->value.swap(basis);
push_value(m.release());push_value(inv_factors.release());wrap_tuple(2);
}/*:20*//*21:*/
#line 343 "built-in-types.w"
void filter_units_wrapper()
{tuple_value*t=get_tuple();push_value(t);
if(t->value.size()!=2)throw std::logic_error("Argument is not a pair");
execution_stack.push_back(t->value[0]);
latmat_value*basis=get_mat();
execution_stack.push_back(t->value[1]);
weight_value*inv_f=get_vec();
if(inv_f->value.size()!=basis->value.numColumns())
throw std::runtime_error("Size mismatch "+
num(inv_f->value.size())+':'+num(basis->value.numColumns())+
" in filter_units");

size_t i=0;
while(i<inv_f->value.size())
if(inv_f->value[i]!=1)++i;
else
{inv_f->value.erase(inv_f->value.begin()+i);
basis->value.eraseColumn(i);
}
}/*:21*//*22:*/
#line 396 "built-in-types.w"
latticetypes::LatticeMatrix
annihilator_modulo
(const latticetypes::LatticeMatrix&M,
latticetypes::LatticeCoeff denominator)

{const size_t m=M.numRows();

latticetypes::WeightList b;matrix::initBasis(b,m);


latticetypes::CoeffList lambda;
smithnormal::smithNormal(lambda,b.begin(),M);


latticetypes::LatticeMatrix A(b);A.invert();A.transpose();

for(size_t j=0;j<lambda.size();++j)
{unsigned long f=(lambda[j]);
unsigned long c=denominator/arithmetic::gcd(f,denominator);
for(size_t i=0;i<m;++i)A(i,j)*=c;
}
return A;
}/*:22*//*23:*/
#line 427 "built-in-types.w"
void ann_mod_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>d(get_int());
std::auto_ptr<latmat_value>m(get_mat());

latticetypes::LatticeMatrix A=
annihilator_modulo(m->value,d->value);
m->value.swap(A);
push_value(m.release());
}/*:23*//*24:*/
#line 451 "built-in-types.w"
void replace_gen_wrapper()
{push_tuple_components();
std::auto_ptr<latmat_value>new_generators(get_mat());
push_tuple_components();
std::auto_ptr<weight_value>inv_f(get_vec());
std::auto_ptr<latmat_value>old_generators(get_mat());

if(new_generators->value.numRows()!=old_generators->value.numRows())
throw std::runtime_error("Column lengths do not match in replace_gen");

if(inv_f->value.size()!=old_generators->value.numColumns())
throw std::runtime_error("Number of columns mismatch in replace_gen");


for(size_t j=0,k=0;j<inv_f->value.size();++j)
if(inv_f->value[j]!=1)

{if(k>=new_generators->value.numColumns())
throw std::runtime_error
("Not enough replacement columns in replace_gen");

for(size_t i=0;i<old_generators->value.numRows();++i)
old_generators->value(i,j)=new_generators->value(i,k);
++k;
}
push_value(old_generators.release());
}/*:24*//*25:*/
#line 499 "built-in-types.w"
lietype::InnerClassType transform_inner_class_type
(const char*s,const lietype::LieType&lt)
throw(std::bad_alloc,std::runtime_error)
{static const std::string types=atlas::lietype::innerClassLetters;

lietype::InnerClassType result(0);
size_t i=0;
for(;skip_punctuation(s),*s!='\0';++s)/*26:*/
#line 520 "built-in-types.w"
{if(types.find(*s)==std::string::npos)throw std::runtime_error
(std::string("Unknown inner class symbol `")+ *s+"'");
if(i>=lt.size())throw std::runtime_error("Too many inner class symbols");
if(*s=='C')
{if(i+1>=lt.size()or lt[i+1]!=lt[i])throw std::runtime_error
("Complex inner class needs two identical consecutive types");
result.push_back('C');
i+=2;
}
else if(*s=='s'){result.push_back('s');++i;}
else if(*s=='c' or*s=='e')
{result.push_back('c');++i;}
else if(*s=='u')
{/*27:*/
#line 545 "built-in-types.w"
{lietype::TypeLetter t=lietype::type(lt[i]);size_t r=lietype::rank(lt[i]);
if(t=='D')result.push_back(r%2==0?'u':'s');
else if(t=='A' and r>=2 or t=='E' and r==6)result.push_back('s');
else throw std::runtime_error
(std::string("unequal rank class is meaningless for type ")+t+num(r));
}/*:27*/
#line 534 "built-in-types.w"
++i;
}
}/*:26*/
#line 510 "built-in-types.w"
if(i<lt.size())throw std::runtime_error("Too few inner class symbols");
return result;
}/*:25*//*28:*/
#line 560 "built-in-types.w"
void basic_involution_wrapper()
{push_tuple_components();
std::auto_ptr<string_value>str(get_string());
std::auto_ptr<Lie_type_value>t(get_Lie_type());
std::auto_ptr<latmat_value>m
(new latmat_value(latticetypes::LatticeMatrix()));
lietype::involution(m->value,t->value
,transform_inner_class_type(str->value.c_str(),t->value));
push_value(m.release());
}/*:28*//*29:*/
#line 590 "built-in-types.w"
void based_involution_wrapper()
{push_tuple_components();
std::auto_ptr<string_value>str(get_string());
std::auto_ptr<latmat_value>basis(get_mat());
std::auto_ptr<Lie_type_value>type(get_Lie_type());

size_t r=type->rank();
if(basis->value.numRows()!=r or basis->value.numRows()!=r)
throw std::runtime_error
("lattice matrix should be "+num(r)+'x'+num(r)+" for this type");
std::auto_ptr<latmat_value>m
(new latmat_value(latticetypes::LatticeMatrix()));
lietype::involution
(m->value,type->value
,transform_inner_class_type(str->value.c_str(),type->value));
m->value*=basis->value;
latticetypes::LatticeCoeff d;basis->value.invert(d);
matrix::leftProd(m->value,basis->value);
if(d==0 or!m->value.divisible(d))throw std::runtime_error
("inner class is not compatible with given lattice");
m->value/=d;push_value(m.release());
}/*:29*//*33:*/
#line 654 "built-in-types.w"
lietype::LieType type_of_datum(const rootdata::RootDatum&rd)
{latticetypes::LatticeMatrix Cartan;rootdata::cartanMatrix(Cartan,rd);
lietype::LieType t;dynkin::lieType(t,Cartan);
if(!rd.isSemisimple())
for(size_t i=rd.semisimpleRank();i<rd.rank();++i)
t.push_back(lietype::SimpleLieType('T',1));
return t;
}/*:33*//*35:*/
#line 677 "built-in-types.w"
void type_of_root_datum_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
push_value(new Lie_type_value(type_of_datum(rd->value)));
}/*:35*//*36:*/
#line 688 "built-in-types.w"
latticetypes::WeightList columns(const latticetypes::LatticeMatrix M)
{latticetypes::WeightList result(M.numColumns());
for(size_t i=0;i<M.numColumns();++i)
M.column(result[i],i);
return result;
}/*:36*//*37:*/
#line 703 "built-in-types.w"
void root_datum_wrapper()
{push_tuple_components();
std::auto_ptr<latmat_value>lattice(get_mat());
std::auto_ptr<Lie_type_value>type(get_Lie_type());
if(lattice->value.numRows()!=lattice->value.numColumns()or
lattice->value.numRows()!=type->rank())
throw std::runtime_error
("lattice matrix should be "+num(type->rank())+'x'+num(type->rank())+
" for this type");/*38:*/
#line 738 "built-in-types.w"
{latticetypes::LatticeMatrix M(lattice->value);
latticetypes::LatticeCoeff d;
M.invert(d);
if(d==0)throw std::runtime_error
("Lattice matrix has dependent columns; in root_datum");
latticetypes::LatticeMatrix tC;
prerootdata::cartanMatrix(tC,type->value);tC.transpose();
M*=tC;

if(!M.divisible(d))
throw std::runtime_error
("Given lattice does not contain the root lattice; in root_datum");
}/*:38*/
#line 714 "built-in-types.w"
prerootdata::PreRootDatum pre(type->value,columns(lattice->value));
push_value(new root_datum_value(rootdata::RootDatum(pre)));
}/*:37*//*39:*/
#line 781 "built-in-types.w"
void quotient_basis_wrapper()
{push_tuple_components();
std::auto_ptr<latmat_value>M(get_mat());

Smith_Cartan_wrapper();
tuple_value*S=get_tuple();push_value(S);
push_value(S->clone());
filter_units_wrapper();
push_tuple_components();
std::auto_ptr<weight_value>v(get_vec());
std::auto_ptr<latmat_value>C(get_mat());
size_t d=1;
for(size_t i=0;i<v->value.size();++i)
if(v->value[i]!=0)d=arithmetic::lcm(d,v->value[i]);/*40:*/
#line 809 "built-in-types.w"
{if(M->value.isEmpty())M->value.resize(v->value.size(),0);
if(M->value.numRows()!=v->value.size())
throw std::runtime_error("Number "+num(M->value.numRows())+
" of rows does not match number "+num(v->value.size())+
" of kernel generators");
for(size_t i=0;i<v->value.size();++i)
if(v->value[i]!=0)
for(size_t j=0;j<M->value.numColumns();++j)
M->value(i,j)*=d/v->value[i];
}/*:40*/
#line 796 "built-in-types.w"
push_value(C.release());
push_value(M.release());
push_value(new int_value(d));
wrap_tuple(2);ann_mod_wrapper();
wrap_tuple(2);mm_prod_wrapper();
wrap_tuple(2);replace_gen_wrapper();
}/*:39*//*41:*/
#line 824 "built-in-types.w"
void quotient_datum_wrapper()
{std::auto_ptr<tuple_value>args(get_tuple());
push_value(args->value[0]->clone());
push_value(args.release());quotient_basis_wrapper();
wrap_tuple(2);root_datum_wrapper();
}/*:41*//*42:*/
#line 838 "built-in-types.w"
void simply_connected_datum_wrapper()
{Lie_type_value*type=get_Lie_type();push_value(type);
push_value(new int_value(type->rank()));
id_mat_wrapper();
wrap_tuple(2);root_datum_wrapper();
}

void adjoint_datum_wrapper()
{Lie_type_value*type=get_Lie_type();push_value(type);
push_value(type->clone());
Cartan_matrix_wrapper();transpose_mat_wrapper();
latmat_value*M=get_mat();
for(size_t i=0;i<type->rank();++i)
if(M->value(i,i)==0)M->value(i,i)=1;
push_value(M);
wrap_tuple(2);root_datum_wrapper();
}/*:42*//*43:*/
#line 871 "built-in-types.w"
void SL_wrapper()
{std::auto_ptr<int_value>n(get_int());
if(n->value<1)throw std::runtime_error("Non positive argument for GL");
const size_t r=n->value-1;
std::auto_ptr<Lie_type_value>type(new Lie_type_value());
if(r>0)type->add_simple_factor('A',r);
push_value(type.release());
std::auto_ptr<latmat_value>lattice
(new latmat_value(latticetypes::LatticeMatrix()));
identityMatrix(lattice->value,r);
for(size_t i=0;i<r-1;++i)lattice->value(i,i+1)= -1;
push_value(lattice.release());
wrap_tuple(2);root_datum_wrapper();
}

void GL_wrapper()
{std::auto_ptr<int_value>n(get_int());
if(n->value<1)throw std::runtime_error("Non positive argument for GL");
const size_t r=n->value-1;
std::auto_ptr<Lie_type_value>type(new Lie_type_value());
if(r>0)type->add_simple_factor('A',r);
type->add_simple_factor('T',1);
push_value(type.release());
std::auto_ptr<latmat_value>lattice
(new latmat_value(latticetypes::LatticeMatrix()));
identityMatrix(lattice->value,r+1);
for(size_t i=0;i<r;++i)
{lattice->value(n->value-1,i)=1;lattice->value(i,i+1)= -1;}
push_value(lattice.release());
wrap_tuple(2);root_datum_wrapper();
}/*:43*//*46:*/
#line 924 "built-in-types.w"
void simple_roots_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
latticetypes::WeightList l
(rd->value.beginSimpleRoot(),rd->value.endSimpleRoot());
latmat_value*result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
push_value(result);
}

void simple_coroots_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
latticetypes::WeightList l
(rd->value.beginSimpleCoroot(),rd->value.endSimpleCoroot());
latmat_value*result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
push_value(result);
}

void datum_Cartan_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
latticetypes::LatticeMatrix M;
rootdata::cartanMatrix(M,rd->value);
latmat_value*result=new latmat_value(M);
push_value(result);
}/*:46*//*47:*/
#line 952 "built-in-types.w"
void roots_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
latticetypes::WeightList l
(rd->value.beginRoot(),rd->value.endRoot());
latmat_value*result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
push_value(result);
}

void coroots_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
latticetypes::WeightList l
(rd->value.beginCoroot(),rd->value.endCoroot());
latmat_value*result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
push_value(result);
}/*:47*//*48:*/
#line 972 "built-in-types.w"
void root_coradical_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
latticetypes::WeightList l
(rd->value.beginSimpleRoot(),rd->value.endSimpleRoot());
l.insert(l.end(),rd->value.beginCoradical(),rd->value.endCoradical());
latmat_value*result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
push_value(result);
}

void coroot_radical_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
latticetypes::WeightList l
(rd->value.beginSimpleCoroot(),rd->value.endSimpleCoroot());
l.insert(l.end(),rd->value.beginRadical(),rd->value.endRadical());
latmat_value*result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
push_value(result);
}/*:48*//*49:*/
#line 994 "built-in-types.w"
void dual_datum_wrapper()
{std::auto_ptr<root_datum_value>rd(get_root_datum());
rootdata::RootDatum dual(rd->value,tags::DualTag());
rd->value.swap(dual);
push_value(rd.release());
}/*:49*//*52:*/
#line 1073 "built-in-types.w"
std::pair<size_t,size_t>classify_involution
(const latticetypes::LatticeMatrix&M,size_t r)
throw(std::bad_alloc,std::runtime_error)
{/*53:*/
#line 1088 "built-in-types.w"
{if(M.numRows()!=r or M.numColumns()!=r)throw std::runtime_error
("involution should be a "+num(r)+"x"+num(r)+" matrix, got a "
+num(M.numRows())+"x"+num(M.numColumns())+" matrix");
latticetypes::LatticeMatrix I,Q(M);
identityMatrix(I,r);Q*=M;
if(!(Q==I))throw std::runtime_error
("given transformation is not an involution");
}/*:53*/
#line 1077 "built-in-types.w"
tori::RealTorus T(M);
return std::make_pair(T.compactRank(),T.splitRank());
}/*:52*//*54:*/
#line 1101 "built-in-types.w"
void classify_wrapper()
{std::auto_ptr<latmat_value>M(get_mat());
size_t r=M->value.numRows();
std::pair<size_t,size_t>p=classify_involution(M->value,r);
push_value(new int_value(p.first));
push_value(new int_value((r-p.first-p.second)/2));
push_value(new int_value(p.second));
wrap_tuple(3);
}/*:54*//*55:*/
#line 1130 "built-in-types.w"
std::pair<lietype::LieType,lietype::InnerClassType>check_involution
(const latticetypes::LatticeMatrix&M,const rootdata::RootDatum&rd)
throw(std::bad_alloc,std::runtime_error)
{size_t r=rd.rank(),s=rd.semisimpleRank();/*53:*/
#line 1088 "built-in-types.w"
{if(M.numRows()!=r or M.numColumns()!=r)throw std::runtime_error
("involution should be a "+num(r)+"x"+num(r)+" matrix, got a "
+num(M.numRows())+"x"+num(M.numColumns())+" matrix");
latticetypes::LatticeMatrix I,Q(M);
identityMatrix(I,r);Q*=M;
if(!(Q==I))throw std::runtime_error
("given transformation is not an involution");
}/*:53*/
#line 1135 "built-in-types.w"
setutils::Permutation p(s);/*56:*/
#line 1152 "built-in-types.w"
{rootdata::WRootIterator first=rd.beginSimpleRoot(),last=rd.endSimpleRoot();
for(size_t i=0;i<s;++i)
{rootdata::Root alpha(r);M.apply(alpha,rd.simpleRoot(i));
size_t p_i=std::find(first,last,alpha)-first;
if(p_i<s)p[i]=p_i;
else throw std::runtime_error
("given transformation does not permute simple roots");
}
for(size_t i=0;i<s;++i)
for(size_t j=0;j<s;++j)
if(rd.cartan(p[i],p[j])!=rd.cartan(i,j))throw std::runtime_error
("given transformation is not a root datum automorphism");
}/*:56*/
#line 1138 "built-in-types.w"
std::pair<lietype::LieType,lietype::InnerClassType>result;
lietype::LieType&type=result.first;
lietype::InnerClassType&inner_class=result.second;/*57:*/
#line 1176 "built-in-types.w"
{latticetypes::LatticeMatrix(C);rootdata::cartanMatrix(C,rd);
dynkin::DynkinDiagram diagr(C);
lietype::LieType type0;dynkin::lieType(type0,C);
bitset::RankFlagsList comps;dynkin::components(comps,diagr);
std::vector<char>letter(comps.size());
size_t nr_Complex_letters=0;
for(size_t i=0;i<comps.size();++i)
{bool equal_rank=true;
for(bitset::RankFlags::iterator j=comps[i].begin();j();++j)
if(p[*j]!= *j){equal_rank=false;break;}
if(equal_rank)letter[i]='c';

else if(!comps[i][p[comps[i].firstBit()]])

{++nr_Complex_letters;letter[i]='C';}
else letter[i]=type0[i].first=='D' and type0[i].second%2==0?'u':'s';

}/*58:*/
#line 1205 "built-in-types.w"
{std::vector<size_t>pos(nr_Complex_letters);
std::vector<size_t>first(nr_Complex_letters);
for(size_t l=0,i=0;l<comps.size();++l)
if(letter[l]=='C')
{pos[i]=l;first[i]=comps[l].firstBit();++i;}
std::vector<size_t>buddy(nr_Complex_letters,nr_Complex_letters);
for(size_t i=0;i<nr_Complex_letters;++i)
if(buddy[i]==nr_Complex_letters)
{size_t b=diagr.component(p[comps[pos[i]].firstBit()]).firstBit();
for(size_t j=i+1;j<nr_Complex_letters;++j)
if(first[j]==b){buddy[i]=j;buddy[j]=i;break;}
}
type.resize(comps.size()+(r-s));
for(size_t l=0,k=0,i=0;l<comps.size();++l)
if(letter[l]!='C')
{type[k++]=type0[l];inner_class.push_back(letter[l]);}
else if(buddy[i]<i)++i;
else
{type[k++]=type0[l];type[k++]=type0[pos[buddy[i]]];
inner_class.push_back('C');++i;
}
if(r>s)/*59:*/
#line 1245 "built-in-types.w"
{using latticetypes::WeightList;using latticetypes::LatticeMatrix;
for(size_t k=comps.size();k<comps.size()+(r-s);++k)
type[k]=lietype::SimpleLieType('T',1);
WeightList b;matrix::initBasis(b,r);
WeightList simple_roots(rd.beginSimpleRoot(),rd.endSimpleRoot());
latticetypes::CoeffList ivf;
smithnormal::smithNormal(ivf,b.begin(),simple_roots);
LatticeMatrix inv(LatticeMatrix(M,b),s,s,r,r);

std::pair<size_t,size_t>cl=classify_involution(inv,r-s);
size_t&compact_rank=cl.first;
size_t&split_rank=cl.second;
size_t Complex_rank=(r-s-compact_rank-split_rank)/2;
while(compact_rank-- >0)inner_class.push_back('c');
while(Complex_rank-- >0)inner_class.push_back('C');
while(split_rank-- >0)inner_class.push_back('s');
}/*:59*/
#line 1228 "built-in-types.w"
}/*:58*/
#line 1196 "built-in-types.w"
}/*:57*/
#line 1142 "built-in-types.w"
return result;
}/*:55*//*66:*/
#line 1396 "built-in-types.w"
void fix_involution_wrapper()
{push_tuple_components();
std::auto_ptr<latmat_value>M(get_mat());
std::auto_ptr<root_datum_value>rd(get_root_datum());
std::pair<lietype::LieType,lietype::InnerClassType>cl
=check_involution(M->value,rd->value);
if(verbosity>0)/*67:*/
#line 1416 "built-in-types.w"
{Lie_type_value t(cl.first);
*output_stream<<"Found "<<t<<", and inner class '";
for(size_t i=0;i<cl.second.size();++i)
*output_stream<<cl.second[i];
*output_stream<<"'.\n";
}/*:67*/
#line 1403 "built-in-types.w"
rootdata::RootDatum*rdp=new rootdata::RootDatum(rd->value);
std::auto_ptr<complexredgp::ComplexReductiveGroup>
G(new complexredgp::ComplexReductiveGroup(rdp,M->value));
inner_class_value*result=new
inner_class_value(G.get(),cl.first,cl.second);
G.release();push_value(result);
}/*:66*//*68:*/
#line 1431 "built-in-types.w"
void set_type_wrapper()
{push_tuple_components();
std::auto_ptr<string_value>ict(get_string());
std::auto_ptr<latmat_value>gen(get_mat());
Lie_type_wrapper();
std::auto_ptr<Lie_type_value>t(get_Lie_type());

push_value(t->clone());push_value(gen.release());
wrap_tuple(2);quotient_basis_wrapper();
std::auto_ptr<latmat_value>basis(get_mat());

push_value(t->clone());push_value(basis->clone());
wrap_tuple(2);root_datum_wrapper();
push_value(t.release());push_value(basis.release());
push_value(ict.release());
wrap_tuple(3);based_involution_wrapper();
wrap_tuple(2);fix_involution_wrapper();
}/*:68*//*69:*/
#line 1466 "built-in-types.w"
void set_inner_class_wrapper()
{push_tuple_components();
std::auto_ptr<string_value>ict(get_string());
std::auto_ptr<root_datum_value>rd(get_root_datum());

push_value(rd->clone());type_of_root_datum_wrapper();
std::auto_ptr<Lie_type_value>t(get_Lie_type());

push_value(rd->clone());
coroot_radical_wrapper();transpose_mat_wrapper();
std::auto_ptr<latmat_value>basis(get_mat());

push_value(rd.release());
push_value(t.release());push_value(basis.release());
push_value(ict.release());
wrap_tuple(3);based_involution_wrapper();
wrap_tuple(2);fix_involution_wrapper();
}/*:69*//*72:*/
#line 1507 "built-in-types.w"
void distinguished_involution_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
push_value(new latmat_value(G->value.distinguished()));
}

void root_datum_of_inner_class_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
push_value(new root_datum_value(G->value.rootDatum()));
}

void dual_inner_class_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
push_value(new inner_class_value(*G,tags::DualTag()));
}/*:72*//*73:*/
#line 1528 "built-in-types.w"
void push_name_list(const realform_io::Interface&interface)
{std::auto_ptr<row_value>result
(new row_value(std::vector<value_ptr>()));
for(size_t i=0;i<interface.numRealForms();++i)
result->value.push_back(new string_value(interface.typeName(i)));
push_value(result.release());
}

void form_names_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
push_name_list(G->interface);
}

void dual_form_names_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
push_name_list(G->dual_interface);
}/*:73*//*74:*/
#line 1555 "built-in-types.w"
void block_sizes_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
G->value.fillCartan();
std::auto_ptr<latmat_value>
M(new latmat_value
(latticetypes::LatticeMatrix(G->value.numRealForms()
,G->value.numDualRealForms())
));
for(size_t i=0;i<M->value.numRows();++i)
for(size_t j=0;j<M->value.numColumns();++j)
M->value(i,j)=
G->value.blockSize(G->interface.in(i),G->dual_interface.in(j));
push_value(M.release());
}/*:74*//*75:*/
#line 1575 "built-in-types.w"
void occurrence_matrix_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
G->value.fillCartan();
size_t nr=G->value.numRealForms();
size_t nc=G->value.numCartanClasses();
std::auto_ptr<latmat_value>
M(new latmat_value(latticetypes::LatticeMatrix(nr,nc)));
for(size_t i=0;i<nr;++i)
{bitmap::BitMap b=G->value.cartanSet(G->interface.in(i));
for(size_t j=0;j<nc;++j)
M->value(i,j)=b.isMember(j)?1:0;
}
push_value(M.release());
}/*:75*//*76:*/
#line 1595 "built-in-types.w"
void dual_occurrence_matrix_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
G->value.fillCartan();
size_t nr=G->value.numDualRealForms();
size_t nc=G->value.numCartanClasses();
std::auto_ptr<latmat_value>
M(new latmat_value(latticetypes::LatticeMatrix(nr,nc)));
for(size_t i=0;i<nr;++i)
{bitmap::BitMap b=G->value.dualCartanSet(G->dual_interface.in(i));
for(size_t j=0;j<nc;++j)
M->value(i,j)=b.isMember(j)?1:0;
}
push_value(M.release());
}/*:76*//*81:*/
#line 1711 "built-in-types.w"
void real_form_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>i(get_int());
std::auto_ptr<inner_class_value>G(get_complexredgp());
if(i->value<0||size_t(i->value)>=G->value.numRealForms())
throw std::runtime_error("illegal real form number: "+num(i->value));
push_value(new real_form_value(*G,G->interface.in(i->value)));
}

void quasisplit_form_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
push_value(new real_form_value(*G,G->value.quasisplit()));
}/*:81*//*84:*/
#line 1747 "built-in-types.w"
void components_rank_wrapper()
{std::auto_ptr<real_form_value>R(get_real_form());
const latticetypes::ComponentList c=R->value.componentReps();
push_value(new int_value(c.size()));
}/*:84*//*85:*/
#line 1758 "built-in-types.w"
void count_Cartans_wrapper()
{std::auto_ptr<real_form_value>rf(get_real_form());
rf->value.fillCartan();
push_value(new int_value(rf->value.numCartan()));
}/*:85*//*86:*/
#line 1768 "built-in-types.w"
void KGB_size_wrapper()
{std::auto_ptr<real_form_value>rf(get_real_form());
rf->value.fillCartan();
push_value(new int_value(rf->value.kgbSize()));
}/*:86*//*87:*/
#line 1780 "built-in-types.w"
void Cartan_order_matrix_wrapper()
{std::auto_ptr<real_form_value>rf(get_real_form());
rf->value.fillCartan();
size_t n=rf->value.numCartan();
latmat_value*M=
new latmat_value(latticetypes::LatticeMatrix(n,n,0));
const poset::Poset&p=rf->value.cartanOrdering();
for(size_t i=0;i<n;++i)
for(size_t j=i;j<n;++j)
if(p.lesseq(i,j))M->value(i,j)=1;

push_value(M);
}/*:87*//*91:*/
#line 1846 "built-in-types.w"
void dual_real_form_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>i(get_int());
std::auto_ptr<inner_class_value>G(get_complexredgp());
if(i->value<0||size_t(i->value)>=G->value.numDualRealForms())
throw std::runtime_error("illegal dual real form number: "+num(i->value));
push_value(new dual_real_form_value(*G,G->dual_interface.in(i->value)));
}

void dual_quasisplit_form_wrapper()
{std::auto_ptr<inner_class_value>G(get_complexredgp());
push_value(new dual_real_form_value(*G,G->dual.quasisplit()));
}/*:91*//*94:*/
#line 1881 "built-in-types.w"
void real_form_from_dual_wrapper()
{std::auto_ptr<dual_real_form_value>d(get_dual_real_form());
push_value(new real_form_value
(inner_class_value(d->parent,tags::DualTag())
,d->value.realForm()));
}/*:94*//*99:*/
#line 1962 "built-in-types.w"
void Cartan_class_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>i(get_int());
std::auto_ptr<real_form_value>rf(get_real_form());
rf->value.fillCartan();
if(i->value<0||size_t(i->value)>=rf->value.numCartan())
throw std::runtime_error
("illegal Cartan class number: "+num(i->value)
+", this real form only has "+num(rf->value.numCartan())+" of them");
bitmap::BitMap cs=rf->value.cartanSet();
push_value(new Cartan_class_value(rf->parent,cs.n_th(i->value)));
}/*:99*//*100:*/
#line 1983 "built-in-types.w"
void most_split_Cartan_wrapper()
{std::auto_ptr<real_form_value>rf(get_real_form());
rf->value.fillCartan();
push_value(new Cartan_class_value(rf->parent,rf->value.mostSplit()));
}/*:100*//*103:*/
#line 2023 "built-in-types.w"
void print_Cartan_info_wrapper()
{using basic_io::operator<<;

std::auto_ptr<Cartan_class_value>cc(get_Cartan_class());

const rootdata::RootDatum&rd=cc->parent.value.rootDatum();

prettyprint::printTorusType(*output_stream,cc->value.fiber().torus())
<<std::endl;

*output_stream<<"twisted involution orbit size: "<<cc->value.orbitSize()
<<std::endl;


lietype::LieType ilt,rlt,clt;
rootdata::lieType(ilt,cc->value.simpleImaginary(),rd);

if(ilt.size()==0)
*output_stream<<"imaginary root system is empty"<<std::endl;
else
*output_stream<<"imaginary root system: "<<ilt<<std::endl;


rootdata::lieType(rlt,cc->value.simpleReal(),rd);

if(rlt.size()==0)
*output_stream<<"real root system is empty"<<std::endl;
else
*output_stream<<"real root system: "<<rlt<<std::endl;


rootdata::lieType(clt,cc->value.simpleComplex(),rd);

if(clt.size()==0)
*output_stream<<"complex factor is empty"<<std::endl;
else
*output_stream<<"complex factor: "<<clt<<std::endl;

wrap_tuple(0);
}/*:103*//*104:*/
#line 2079 "built-in-types.w"
void fiber_part_wrapper()
{push_tuple_components();
std::auto_ptr<real_form_value>rf(get_real_form());
std::auto_ptr<Cartan_class_value>cc(get_Cartan_class());
if(&rf->parent.value!= &cc->parent.value)
throw std::runtime_error
("inner class mismatch between real form and Cartan class");
bitmap::BitMap b(cc->parent.value.cartanSet(rf->value.realForm()));
if(!b.isMember(cc->number))
throw std::runtime_error
("fiber_part: inner class not defined for this real form");

const partition::Partition&pi=cc->value.fiber().weakReal();
const realform::RealFormList rf_nr=
cc->parent.value.realFormLabels(cc->number);

std::auto_ptr<row_value>result
(new row_value(std::vector<value_ptr>()));
for(size_t i=0;i<pi.size();++i)
if(rf_nr[pi(i)]==rf->value.realForm())
result->value.push_back(new int_value(i));
push_value(result.release());
}/*:104*//*105:*/
#line 2114 "built-in-types.w"
void print_gradings_wrapper()
{push_tuple_components();
std::auto_ptr<real_form_value>rf(get_real_form());
std::auto_ptr<Cartan_class_value>cc(get_Cartan_class());
if(&rf->parent.value!= &cc->parent.value)
throw std::runtime_error
("inner class mismatch between real form and Cartan class");
bitmap::BitMap b(cc->parent.value.cartanSet(rf->value.realForm()));
if(!b.isMember(cc->number))
throw std::runtime_error
("fiber_part: inner class not defined for this real form");

const partition::Partition&pi=cc->value.fiber().weakReal();
const realform::RealFormList rf_nr=
cc->parent.value.realFormLabels(cc->number);


const rootdata::RootList&si=cc->value.fiber().simpleImaginary();


latticetypes::LatticeMatrix cm;
setutils::Permutation sigma;/*106:*/
#line 2150 "built-in-types.w"
{rootdata::cartanMatrix(cm,si,cc->parent.value.rootDatum());
dynkin::DynkinDiagram d(cm);dynkin::bourbaki(sigma,d);
}/*:106*//*107:*/
#line 2158 "built-in-types.w"
{*output_stream<<"Imaginary root system is ";
if(si.size()==0)*output_stream<<"empty.\n";
else
{using basic_io::operator<<;
lietype::LieType t;dynkin::lieType(t,cm);

*output_stream<<"of type "<<t<<", with simple root"
<<(si.size()==1?" ":"s ");
for(size_t i=0;i<si.size();++i)
*output_stream<<si[sigma[i]]<<(i<si.size()-1?",":".\n");
}
}/*:107*//*108:*/
#line 2181 "built-in-types.w"
{bool first=true;
for(size_t i=0;i<pi.size();++i)
if(rf_nr[pi(i)]==rf->value.realForm())
{*output_stream<<(first?first=false,'[':',');
gradings::Grading gr;cc->value.fiber().grading(gr,i);
gr.permute(sigma);
prettyprint::prettyPrint(*output_stream,gr,si.size());
}
*output_stream<<"]"<<std::endl;
}/*:108*/
#line 2143 "built-in-types.w"
wrap_tuple(0);
}/*:105*//*109:*/
#line 2207 "built-in-types.w"
void print_block_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_block: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
block_io::printBlock(*output_stream,block);

wrap_tuple(0);
}/*:109*//*110:*/
#line 2232 "built-in-types.w"
void print_blockd_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_blockd: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
block_io::printBlockD(*output_stream,block);

wrap_tuple(0);
}


void print_blocku_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_blocku: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
block_io::printBlockU(*output_stream,block);

wrap_tuple(0);
}/*:110*//*111:*/
#line 2292 "built-in-types.w"
void print_blockstabilizer_wrapper()
{push_tuple_components();
std::auto_ptr<Cartan_class_value>cc(get_Cartan_class());
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

if(&rf->parent.value!= &drf->parent.value or
&rf->parent.value!= &cc->parent.value)
throw std::runtime_error
("blockstabilizer: inner class mismatch between arguments");
bitmap::BitMap b(rf->parent.value.cartanSet(rf->value.realForm()));
b&=bitmap::BitMap(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(not b.isMember(cc->number))
throw std::runtime_error
("blockstabilizer: Cartan class not defined for both real forms");


realredgp_io::printBlockStabilizer
(*output_stream,rf->value,cc->number,drf->value.realForm());

wrap_tuple(0);
}/*:111*//*112:*/
#line 2321 "built-in-types.w"
void print_KGB_wrapper()
{std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();

*output_stream
<<"kgbsize: "<<rf->value.kgbSize()<<std::endl;
kgb::KGB kgb(rf->value);
kgb_io::printKGB(*output_stream,kgb);

wrap_tuple(0);
}/*:112*//*113:*/
#line 2341 "built-in-types.w"
void print_KL_basis_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_KL_basis: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
klsupport::KLSupport kls(block);kls.fill();

kl::KLContext klc(kls);klc.fill();

*output_stream
<<"Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:\n\n";
kl_io::printAllKL(*output_stream,klc);

wrap_tuple(0);
}/*:113*//*114:*/
#line 2371 "built-in-types.w"
void print_prim_KL_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_prim_KL: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
klsupport::KLSupport kls(block);kls.fill();

kl::KLContext klc(kls);klc.fill();

*output_stream
<<"Non-zero Kazhdan-Lusztig-Vogan polynomials for primitive pairs:\n\n";
kl_io::printPrimitiveKL(*output_stream,klc);

wrap_tuple(0);
}/*:114*//*115:*/
#line 2402 "built-in-types.w"
void print_KL_list_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_KL_list: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
klsupport::KLSupport kls(block);kls.fill();

kl::KLContext klc(kls);klc.fill();

kl_io::printKLList(*output_stream,klc);

wrap_tuple(0);
}/*:115*//*116:*/
#line 2435 "built-in-types.w"
void print_W_cells_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_W_cells: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
klsupport::KLSupport kls(block);kls.fill();

kl::KLContext klc(kls);klc.fill();

wgraph::WGraph wg(klc.rank());kl::wGraph(wg,klc);

wgraph_io::printCells(*output_stream,wg);

wrap_tuple(0);
}/*:116*//*117:*/
#line 2466 "built-in-types.w"
void print_W_graph_wrapper()
{push_tuple_components();
std::auto_ptr<dual_real_form_value>drf(get_dual_real_form());
std::auto_ptr<real_form_value>rf(get_real_form());

rf->value.fillCartan();
if(&rf->parent.value!= &drf->parent.value)
throw std::runtime_error
("inner class mismatch between real form and dual real form");
bitmap::BitMap b(rf->parent.value.dualCartanSet(drf->value.realForm()));
if(!b.isMember(rf->value.mostSplit()))
throw std::runtime_error
("print_W_graph: real form and dual real form are incompatible");

blocks::Block block(rf->parent.value
,rf->value.realForm(),drf->value.realForm());
klsupport::KLSupport kls(block);kls.fill();

kl::KLContext klc(kls);klc.fill();

wgraph::WGraph wg(klc.rank());kl::wGraph(wg,klc);

wgraph_io::printWGraph(*output_stream,wg);

wrap_tuple(0);
}/*:117*/
#line 15 "built-in-types.w"
}/*4:*/
#line 46 "built-in-types.w"
void initialise_builtin_types()
{/*11:*/
#line 171 "built-in-types.w"
install_function(Lie_type_wrapper,"Lie_type","(string->LieType)");/*:11*//*18:*/
#line 250 "built-in-types.w"
install_function(Cartan_matrix_wrapper,"Cartan_matrix","(LieType->mat)");
install_function(type_of_Cartan_matrix_wrapper
,"type_of_Cartan_matrix","(mat->LieType)");/*:18*//*30:*/
#line 617 "built-in-types.w"
install_function(Smith_Cartan_wrapper,"Smith_Cartan","(LieType->mat,vec)");
install_function(filter_units_wrapper,"filter_units","(mat,vec->mat,vec)");
install_function(ann_mod_wrapper,"ann_mod","(mat,int->mat)");
install_function(replace_gen_wrapper,"replace_gen",
"((mat,vec),mat->mat)");
install_function(basic_involution_wrapper,"basic_involution",
"(LieType,string->mat)");
install_function(based_involution_wrapper,"based_involution",
"(LieType,mat,string->mat)");/*:30*//*50:*/
#line 1004 "built-in-types.w"
install_function(type_of_root_datum_wrapper,"type_of_root_datum"
,"(RootDatum->LieType)");
install_function(root_datum_wrapper,"root_datum","(LieType,mat->RootDatum)");
install_function(quotient_basis_wrapper
,"quotient_basis","(LieType,mat->mat)");
install_function(quotient_datum_wrapper
,"quotient_datum","(LieType,mat->RootDatum)");
install_function(simply_connected_datum_wrapper
,"simply_connected_datum","(LieType->RootDatum)");
install_function(adjoint_datum_wrapper,"adjoint_datum","(LieType->RootDatum)");
install_function(SL_wrapper,"SL","(int->RootDatum)");
install_function(GL_wrapper,"GL","(int->RootDatum)");
install_function(simple_roots_wrapper,"simple_roots","(RootDatum->mat)");
install_function(simple_coroots_wrapper,"simple_coroots","(RootDatum->mat)");
install_function(datum_Cartan_wrapper,"Cartan_matrix_of_datum"
,"(RootDatum->mat)");
install_function(roots_wrapper,"roots","(RootDatum->mat)");
install_function(coroots_wrapper,"coroots","(RootDatum->mat)");
install_function(root_coradical_wrapper,"root_coradical","(RootDatum->mat)");
install_function(coroot_radical_wrapper,"coroot_radical","(RootDatum->mat)");
install_function(dual_datum_wrapper,"dual_datum","(RootDatum->RootDatum)");/*:50*//*77:*/
#line 1612 "built-in-types.w"
install_function(classify_wrapper,"classify_involution"
,"(mat->int,int,int)");
install_function(fix_involution_wrapper,"fix_involution"
,"(RootDatum,mat->InnerClass)");
install_function(set_type_wrapper,"set_type"
,"(string,mat,string->InnerClass)");
install_function(set_inner_class_wrapper,"set_inner_class"
,"(RootDatum,string->InnerClass)");
install_function(distinguished_involution_wrapper,"distinguished_involution"
,"(InnerClass->mat)");
install_function(root_datum_of_inner_class_wrapper,"root_datum_of_inner_class"
,"(InnerClass->RootDatum)");
install_function(dual_inner_class_wrapper,"dual_inner_class"
,"(InnerClass->InnerClass)");
install_function(form_names_wrapper,"form_names"
,"(InnerClass->[string])");
install_function(dual_form_names_wrapper,"dual_form_names"
,"(InnerClass->[string])");
install_function(block_sizes_wrapper,"block_sizes"
,"(InnerClass->mat)");
install_function(occurrence_matrix_wrapper,"occurrence_matrix"
,"(InnerClass->mat)");
install_function(dual_occurrence_matrix_wrapper,"dual_occurrence_matrix"
,"(InnerClass->mat)");/*:77*//*118:*/
#line 2496 "built-in-types.w"
install_function(real_form_wrapper,"real_form","(InnerClass,int->RealForm)");
install_function(quasisplit_form_wrapper,"quasisplit_form"
,"(InnerClass->RealForm)");
install_function(components_rank_wrapper,"components_rank","(RealForm->int)");
install_function(count_Cartans_wrapper,"count_Cartans","(RealForm->int)");
install_function(KGB_size_wrapper,"KGB_size","(RealForm->int)");
install_function(Cartan_order_matrix_wrapper,"Cartan_order_matrix"
,"(RealForm->mat)");
install_function(dual_real_form_wrapper,"dual_real_form"
,"(InnerClass,int->DualRealForm)");
install_function(dual_quasisplit_form_wrapper,"dual_quasisplit_form"
,"(InnerClass->DualRealForm)");
install_function(real_form_from_dual_wrapper,"real_form_from_dual"
,"(DualRealForm->RealForm)");
install_function(Cartan_class_wrapper,"Cartan_class"
,"(RealForm,int->CartanClass)");
install_function(most_split_Cartan_wrapper,"most_split_Cartan"
,"(RealForm->CartanClass)");
install_function(print_Cartan_info_wrapper,"print_Cartan_info"
,"(CartanClass->)");
install_function(fiber_part_wrapper,"fiber_part"
,"(CartanClass,RealForm->[int])");
install_function(print_gradings_wrapper,"print_gradings"
,"(CartanClass,RealForm->)");
install_function(print_block_wrapper,"print_block"
,"(RealForm,DualRealForm->)");
install_function(print_blocku_wrapper,"print_blocku"
,"(RealForm,DualRealForm->)");
install_function(print_blockd_wrapper,"print_blockd"
,"(RealForm,DualRealForm->)");
install_function(print_blockstabilizer_wrapper,"print_blockstabilizer"
,"(RealForm,DualRealForm,CartanClass->)");
install_function(print_KGB_wrapper,"print_KGB"
,"(RealForm->)");
install_function(print_KL_basis_wrapper,"print_KL_basis"
,"(RealForm,DualRealForm->)");
install_function(print_prim_KL_wrapper,"print_prim_KL"
,"(RealForm,DualRealForm->)");
install_function(print_KL_list_wrapper,"print_KL_list"
,"(RealForm,DualRealForm->)");
install_function(print_W_cells_wrapper,"print_W_cells"
,"(RealForm,DualRealForm->)");
install_function(print_W_graph_wrapper,"print_W_graph"
,"(RealForm,DualRealForm->)");/*:118*/
#line 48 "built-in-types.w"
}/*:4*//*9:*/
#line 108 "built-in-types.w"
void Lie_type_value::add_simple_factor(char c,size_t rank)
throw(std::bad_alloc,std::runtime_error)
{using std::runtime_error;
static const std::string types=atlas::lietype::typeLetters;
size_t t=types.find(c);
if(t==std::string::npos)
throw runtime_error(std::string("Invalid type letter '")+c+'\'');
static const size_t lwb[]={1,2,2,4,6,4,2,0};
static const size_t r=constants::RANK_MAX;
static const size_t upb[]={r,r,r,r,8,4,2,r};
if(rank<lwb[t])
throw runtime_error("Too small rank "+num(rank)+" for Lie type "+c);

if(rank>upb[t])
if(upb[t]!=r)
throw runtime_error("Too large rank "+num(rank)+" for Lie type "+c);

else
throw runtime_error
("Rank "+num(rank)+" exceeds implementation limit "+num(r));

if(c=='T')
while(rank-- >0)value.push_back(lietype::SimpleLieType('T',1));
else
value.push_back(lietype::SimpleLieType(c,rank));
}/*:9*//*12:*/
#line 180 "built-in-types.w"
void Lie_type_value::print(std::ostream&out)const
{if(value.empty())out<<"empty Lie type";
else
{using basic_io::operator<<;
out<<"Lie type '"<<value<<'\'';
}
}/*:12*//*13:*/
#line 193 "built-in-types.w"
size_t Lie_type_value::rank()const
{size_t r=0;
for(size_t i=0;i<value.size();++i)r+=value[i].second;
return r;
}

size_t Lie_type_value::semisimple_rank()const
{size_t r=0;
for(size_t i=0;i<value.size();++i)
if(value[i].first!='T')r+=value[i].second;
return r;
}/*:13*//*15:*/
#line 213 "built-in-types.w"
Lie_type_value*get_Lie_type()throw(std::logic_error)
{value_ptr p=pop_arg();
Lie_type_value*result=dynamic_cast<Lie_type_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a Lie type");}

return result;
}/*:15*//*34:*/
#line 668 "built-in-types.w"
void root_datum_value::print(std::ostream&out)const
{Lie_type_value type(type_of_datum(value));
out<<(value.isSimplyConnected()?"simply connected ":"")
<<(value.isAdjoint()?"adjoint ":"")
<<"root datum of "<<type;
}/*:34*//*45:*/
#line 911 "built-in-types.w"
root_datum_value*get_root_datum()throw(std::logic_error)
{value_ptr p=pop_arg();
root_datum_value*result=dynamic_cast<root_datum_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a root_datum");}

return result;
}/*:45*//*62:*/
#line 1329 "built-in-types.w"
inner_class_value::inner_class_value(const inner_class_value&v)
:value(v.value),dual(v.dual),ref_count(v.ref_count)
,interface(v.interface),dual_interface(v.dual_interface)
{++ref_count;}

inner_class_value::~inner_class_value()
{if(--ref_count==0){delete&value;delete&dual;delete&ref_count;}}/*:62*//*63:*/
#line 1346 "built-in-types.w"
inner_class_value::inner_class_value
(complexredgp::ComplexReductiveGroup*g
,lietype::LieType lt,lietype::InnerClassType ict)
:value(*g)
,dual(*new complexredgp::ComplexReductiveGroup(*g,tags::DualTag()))
,ref_count(*new size_t(1))
,interface(*g,layout::Layout(lt,ict))
,dual_interface(*g,layout::Layout(lt,ict),tags::DualTag())
{}/*:63*//*64:*/
#line 1365 "built-in-types.w"
inner_class_value::inner_class_value(const inner_class_value&v,tags::DualTag)
:value(v.dual),dual(v.value)
,ref_count(v.ref_count)
,interface(v.dual_interface),dual_interface(v.interface)
{++ref_count;}/*:64*//*65:*/
#line 1378 "built-in-types.w"
void inner_class_value::print(std::ostream&out)const
{out<<"Complex reductive group equipped with an involution,\n"
"defining an inner class of "
<<value.numRealForms()<<" real "
<<(value.numRealForms()==1?"form":"forms")<<" and "
<<value.numDualRealForms()<<" dual real "
<<(value.numDualRealForms()==1?"form":"forms");
}/*:65*//*71:*/
#line 1493 "built-in-types.w"
inner_class_value*get_complexredgp()throw(std::logic_error)
{value_ptr p=pop_arg();
inner_class_value*result=dynamic_cast<inner_class_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a complex group");}

return result;
}/*:71*//*80:*/
#line 1695 "built-in-types.w"
void real_form_value::print(std::ostream&out)const
{out<<(value.isQuasisplit()?"quasisplit ":"")
<<"real form '"
<<parent.interface.typeName(parent.interface.out(value.realForm()))
<<"', defining a"
<<(value.isConnected()?" connected":"")
<<(value.isSplit()?" split":"")
<<" real group";
}/*:80*//*83:*/
#line 1733 "built-in-types.w"
real_form_value*get_real_form()throw(std::logic_error)
{value_ptr p=pop_arg();
real_form_value*result=dynamic_cast<real_form_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a real form");}

return result;
}/*:83*//*88:*/
#line 1807 "built-in-types.w"
real_form_value::real_form_value(inner_class_value p,realform::RealForm f
,tags::DualTag)
:parent(p),value(p.dual,f)
{}/*:88*//*90:*/
#line 1833 "built-in-types.w"
void dual_real_form_value::print(std::ostream&out)const
{out<<(value.isQuasisplit()?"quasisplit ":"")
<<"dual real form '"
<<parent.dual_interface.typeName
(parent.dual_interface.out(value.realForm()))
<<"'";
}/*:90*//*93:*/
#line 1867 "built-in-types.w"
dual_real_form_value*get_dual_real_form()throw(std::logic_error)
{value_ptr p=pop_arg();
dual_real_form_value*result=dynamic_cast<dual_real_form_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a dual real form");}

return result;
}/*:93*//*97:*/
#line 1940 "built-in-types.w"
Cartan_class_value::Cartan_class_value(const inner_class_value&p,size_t cn)
:parent(p),number(cn),value(p.value.cartan(cn))
{if(cn>=p.value.numCartanClasses())throw std::runtime_error
(std::string("Cartan class number ")+num(cn)+" does not (currently) exist");
}/*:97*//*98:*/
#line 1950 "built-in-types.w"
void Cartan_class_value::print(std::ostream&out)const
{out<<"Cartan class #"<<number<<", occurring for "
<<value.numRealForms()<<" real "
<<(value.numRealForms()==1?"form":"forms")<<" and for "
<<value.numDualRealForms()<<" dual real "
<<(value.numDualRealForms()==1?"form":"forms");
}/*:98*//*102:*/
#line 1998 "built-in-types.w"
Cartan_class_value*get_Cartan_class()throw(std::logic_error)
{value_ptr p=pop_arg();
Cartan_class_value*result=dynamic_cast<Cartan_class_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a Cartan class");}

return result;
}/*:102*/
#line 16 "built-in-types.w"

}}/*:1*/
