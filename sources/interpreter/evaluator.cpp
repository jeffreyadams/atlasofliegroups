#include "evaluator.h"
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "smithnormal.h"
#include <set>
#include <map>
/*1:*/
#line 10 "evaluator.w"
namespace atlas{namespace interpreter{/*14:*/
#line 259 "evaluator.w"
const char*prim_names[]=
{"int","string","bool",/*41:*/
#line 735 "evaluator.w"
"*",/*:41*//*51:*/
#line 908 "evaluator.w"
"vec","mat",/*:51*//*131:*/
#line 2194 "evaluator.w"
"LieType","RootDatum","InnerClass","RealForm","DualRealForm",
"CartanClass",/*:131*/
#line 260 "evaluator.w"
};/*:14*//*92:*/
#line 1410 "evaluator.w"
Id_table*global_id_table=NULL;/*:92*//*96:*/
#line 1474 "evaluator.w"
std::vector<value_ptr>execution_stack;/*:96*//*130:*/
#line 2184 "evaluator.w"
int verbosity=0;/*:130*//*137:*/
#line 2322 "evaluator.w"
std::ostream*output_stream= &std::cout;/*:137*//*35:*/
#line 617 "evaluator.w"
type_ptr find_type(expr e)throw(std::bad_alloc,program_error);/*:35*//*38:*/
#line 674 "evaluator.w"
void check_type(const type_declarator&t,expr&e)
throw(std::bad_alloc,program_error);/*:38*//*47:*/
#line 822 "evaluator.w"
latticetypes::Weight cast_intlist_to_weight(const value_ptr);
latticetypes::LatticeMatrix cast_intlistlist_to_matrix(const value_ptr);/*:47*//*61:*/
#line 1035 "evaluator.w"
type_ptr make_tuple_type(type_list_ptr l);/*:61*//*65:*/
#line 1076 "evaluator.w"
type_list_ptr find_type_list(expr_list);/*:65*//*5:*/
#line 89 "evaluator.w"
type_ptr copy(const type_declarator*t)
{return type_ptr(new type_declarator(*t));}/*:5*//*6:*/
#line 109 "evaluator.w"
type_node::type_node(const type_node&n)
{type_ptr head=copy(n.t);
next=n.next==NULL?NULL:new type_node(*n.next);
t=head.release();
}/*:6*//*7:*/
#line 121 "evaluator.w"
type_node::~type_node()
{delete t;delete next;}/*:7*//*8:*/
#line 157 "evaluator.w"
type_list_ptr make_type_list(type_ptr t,type_list_ptr l)
{type_node*p=new type_node(t.get(),l.get());
t.release(),l.release();return type_list_ptr(p);
}
type_list_ptr make_type_singleton(type_ptr t)
{type_node*p=new type_node(t.get(),NULL);t.release();
return type_list_ptr(p);
}/*:8*//*10:*/
#line 209 "evaluator.w"
type_declarator::type_declarator(const type_declarator&t):kind(t.kind)
{switch(kind)
{case primitive_type:prim=t.prim;break;
case row_type:
comp_type=new type_declarator(*t.comp_type);
break;/*59:*/
#line 1021 "evaluator.w"
case tuple_type:
tuple=t.tuple==NULL?NULL:new type_node(*t.tuple);
break;/*:59*//*77:*/
#line 1233 "evaluator.w"
case function_type:
{type_ptr a=copy(t.func->arg_type)
,r=copy(t.func->result_type);
func=new func_type(a.get(),r.get());a.release(),r.release();
}
break;/*:77*/
#line 216 "evaluator.w"
}
}

type_declarator::~type_declarator()
{switch(kind)
{case primitive_type:break;
case row_type:delete comp_type;break;/*60:*/
#line 1028 "evaluator.w"
case tuple_type:
delete tuple;
break;/*:60*//*78:*/
#line 1247 "evaluator.w"
case function_type:
delete func->arg_type;delete func->result_type;
delete func;
break;/*:78*/
#line 224 "evaluator.w"
}
}/*:10*//*11:*/
#line 232 "evaluator.w"
type_ptr make_prim_type(primitive p)
{return type_ptr(new type_declarator(p));}

type_ptr make_row_type(type_ptr p)
{type_declarator*t=new type_declarator(p.get());p.release();
return type_ptr(t);
}/*:11*//*16:*/
#line 276 "evaluator.w"
std::ostream&operator<<(std::ostream&out,const type_declarator&t)
{switch(t.kind)
{case primitive_type:out<<prim_names[t.prim];break;
case row_type:out<<"["<< *t.comp_type<<"]";break;/*64:*/
#line 1063 "evaluator.w"
case tuple_type:
out<<'(';
for(type_list l=t.tuple;l!=NULL;l=l->next)
out<< *l->t<<(l->next!=NULL?",":"");
out<<')';
break;/*:64*//*81:*/
#line 1276 "evaluator.w"
case function_type:
out<<'(';
if(t.func->arg_type->kind!=tuple_type)
out<< *t.func->arg_type;
else
for(type_list l=t.func->arg_type->tuple;
l!=NULL;l=l->next)
out<< *l->t<<(l->next!=NULL?",":"");
out<<"->";
if(t.func->result_type->kind!=tuple_type)
out<< *t.func->result_type;
else
for(type_list l=t.func->result_type->tuple;
l!=NULL;l=l->next)
out<< *l->t<<(l->next!=NULL?",":"");
out<<")";
break;/*:81*/
#line 281 "evaluator.w"
}
return out;
}/*:16*//*18:*/
#line 299 "evaluator.w"
bool operator!=(const type_declarator&x,const type_declarator&y)
{if(x.kind!=y.kind)return true;
switch(x.kind)
{case primitive_type:return x.prim!=y.prim;
case row_type:
return*x.comp_type!= *y.comp_type;/*63:*/
#line 1053 "evaluator.w"
case tuple_type:
for(type_list l0=x.tuple,
l1=y.tuple;
l0!=NULL||l1!=NULL;l0=l0->next,l1=l1->next)
if(l0==NULL||l1==NULL|| *l0->t!= *l1->t)return true;
return false;/*:63*//*80:*/
#line 1268 "evaluator.w"
case function_type:
return*x.func->arg_type!= *y.func->arg_type
|| *x.func->result_type!= *y.func->result_type;/*:80*/
#line 307 "evaluator.w"
}
return false;
}/*:18*//*22:*/
#line 358 "evaluator.w"
std::ostream&operator<<(std::ostream&out,const value_base&v)
{v.print(out);return out;}

value_base::value_base(const value_base&x)
{std::cerr<<"Copying "<<x<<std::endl;}
value_base&value_base::operator=(const value_base&x)
{std::cerr<<"Assigning "<< *this<<"<-"<<x<<std::endl;
return*this;
}/*:22*//*26:*/
#line 435 "evaluator.w"
row_value::row_value(const row_value&v):value(v.value)
{for(std::vector<value_ptr>::iterator p=value.begin();p!=value.end();++p)
*p=(*p)->clone();
}

row_value::~row_value()
{for(std::vector<value_ptr>::iterator p=value.begin();p!=value.end();++p)
delete*p;
}/*:26*//*27:*/
#line 451 "evaluator.w"
void row_value::print(std::ostream&out)const
{if(value.empty())out<<"[]";
else
{out<<'[';
std::vector<value_ptr>::const_iterator p=value.begin();
do{(*p)->print(out);++p;out<<(p==value.end()?']':',');}
while(p!=value.end());
}
}/*:27*//*30:*/
#line 489 "evaluator.w"
value_ptr evaluate(expr e)
throw(std::bad_alloc,std::logic_error,std::runtime_error)
{value_ptr result=NULL;
switch(e.kind)
{case integer_denotation:
result=new int_value(e.e.int_denotation_variant);
break;
case string_denotation:
result=new string_value(e.e.str_denotation_variant);
break;
case boolean_denotation:
result=new bool_value(e.e.int_denotation_variant);
break;
case list_display:/*31:*/
#line 519 "evaluator.w"
{size_t length=0;
static const value_ptr nil=NULL;
for(expr_list l=e.e.sublist;l!=NULL;l=l->next)++length;
row_ptr v(new
row_value(std::vector<value_ptr>(length,nil)));
expr_list l=e.e.sublist;
for(size_t i=0;i<length;++i,l=l->next)
v->value[i]=evaluate(l->e);
result=v.release();
}/*:31*/
#line 505 "evaluator.w"
break;/*70:*/
#line 1126 "evaluator.w"
case tuple_display:
{size_t length=0;
static const value_ptr nil=NULL;
for(expr_list l=e.e.sublist;l!=NULL;l=l->next)++length;
row_ptr v(new
tuple_value(std::vector<value_ptr>(length,nil)));
expr_list l=e.e.sublist;
for(size_t i=0;i<length;++i,l=l->next)
v->value[i]=evaluate(l->e);
result=v.release();
}
break;/*:70*//*101:*/
#line 1548 "evaluator.w"
case function_call:
{push_value(evaluate(e.e.call_variant->arg));

std::string name=main_hash_table->name_of(e.e.call_variant->fun);

value_ptr f_val=global_id_table->value_of(e.e.call_variant->fun);
if(f_val==NULL)throw
std::logic_error("Built-in function absent: "+name);

builtin_value*b=dynamic_cast<builtin_value* >(f_val);
if(b==NULL)throw
std::logic_error("Built-in not a function: "+name);

b->value();
result=pop_arg();
}
break;/*:101*//*121:*/
#line 2012 "evaluator.w"
case applied_identifier:
{value_ptr p=global_id_table->value_of(e.e.identifier_variant);
if(p==NULL)throw std::logic_error
("Identifier without value:"
+main_hash_table->name_of(e.e.identifier_variant));

result=p->clone();
}
break;/*:121*/
#line 507 "evaluator.w"
}
return result;
}/*:30*//*34:*/
#line 605 "evaluator.w"
type_error::type_error(const type_error&e)
:program_error(e),offender(e.offender)
{type_ptr p=copy(e.required);
actual=copy(e.actual).release();required=p.release();
}/*:34*//*36:*/
#line 622 "evaluator.w"
type_ptr find_type(expr e)throw(std::bad_alloc,program_error)
{switch(e.kind)
{case integer_denotation:return make_prim_type(integral_type);
case string_denotation:return make_prim_type(string_type);
case boolean_denotation:return make_prim_type(boolean_type);
case list_display:/*37:*/
#line 642 "evaluator.w"
{if(e.e.sublist==NULL)return make_type("[int]");
type_ptr comp_type=find_type(e.e.sublist->e);
for(expr_list l=e.e.sublist->next;l!=NULL;l=l->next)
{type_ptr c2=find_type(l->e);
if(*c2!= *comp_type)
throw type_error(l->e,*c2.release(),*comp_type.release());

}
return make_row_type(comp_type);
}/*:37*//*66:*/
#line 1080 "evaluator.w"
case tuple_display:return make_tuple_type(find_type_list(e.e.sublist));/*:66*//*93:*/
#line 1421 "evaluator.w"
case function_call:
{type_declarator*f_type=global_id_table->type_of(e.e.call_variant->fun);
if(f_type==NULL||f_type->kind!=function_type)
throw program_error("Call of "
+std::string(f_type==NULL?"unknown ":"non-")
+"function '"
+main_hash_table->name_of(e.e.call_variant->fun)
+"'");


check_type(*f_type->func->arg_type,e.e.call_variant->arg);
return copy(f_type->func->result_type);
}/*:93*//*122:*/
#line 2026 "evaluator.w"
case applied_identifier:
{type_declarator*t=global_id_table->type_of(e.e.identifier_variant);
if(t==NULL)throw program_error
("Undefined identifier "
+main_hash_table->name_of(e.e.identifier_variant));

return copy(t);
}/*:122*/
#line 631 "evaluator.w"
}
return type_ptr(NULL);
}/*:36*//*39:*/
#line 691 "evaluator.w"
void check_type(const type_declarator&t,expr&e)
throw(std::bad_alloc,program_error)
{static type_declarator vect(vector_type);
static type_declarator matr(matrix_type);

type_ptr actual(NULL);
switch(e.kind)
{case integer_denotation:
if(t.kind!=primitive_type||t.prim!=integral_type)
actual=make_prim_type(integral_type);
break;
case string_denotation:
if(t.kind!=primitive_type||t.prim!=string_type)
actual=make_prim_type(string_type);
break;
case boolean_denotation:
if(t.kind!=primitive_type||t.prim!=boolean_type)
actual=make_prim_type(boolean_type);
break;/*42:*/
#line 744 "evaluator.w"
case list_display:
if(t.kind==row_type)
for(expr_list l=e.e.sublist;l!=NULL;l=l->next)
check_type(*t.comp_type,l->e);
else
{/*43:*/
#line 765 "evaluator.w"
{expr_list l=e.e.sublist;
type_ptr comp(NULL);
if(t==vect)
{/*44:*/
#line 787 "evaluator.w"
{expr_list arg=make_exprlist_node(e,null_expr_list);
e=make_application_node
(main_hash_table->match_literal(">vec<[int]:"),arg);
}/*:44*/
#line 769 "evaluator.w"
comp=make_type("int");
}
else if(t==matr)
{/*45:*/
#line 795 "evaluator.w"
{expr_list arg=make_exprlist_node(e,null_expr_list);
e=make_application_node
(main_hash_table->match_literal(">mat<[vec]:"),arg);
}/*:45*/
#line 773 "evaluator.w"
comp=make_type("vec");
}
if(comp.get()!=NULL)
{for(;l!=NULL;l=l->next)check_type(*comp,l->e);
break;
}
}/*:43*/
#line 752 "evaluator.w"
actual=make_type("[*]");
}
break;/*:42*//*71:*/
#line 1145 "evaluator.w"
case tuple_display:
if(t.kind==tuple_type)/*72:*/
#line 1166 "evaluator.w"
{type_list l=t.tuple;
for(expr_list a=e.e.sublist;a!=NULL||l!=NULL;a=a->next,l=l->next)
if(a==NULL||l==NULL)
throw program_error("Too "+std::string(a==NULL?"few":"many")
+" components in tuple");


else check_type(*l->t,a->e);
}/*:72*/
#line 1149 "evaluator.w"
else
{type_list_ptr tl(NULL);
static type_declarator unknown(undetermined_type);
for(expr_list l=e.e.sublist;l!=NULL;l=l->next)
tl=make_type_list(copy(&unknown),tl);
actual=make_tuple_type(tl);
}
break;/*:71*//*94:*/
#line 1440 "evaluator.w"
case function_call:
{type_declarator*f_type=global_id_table->type_of(e.e.call_variant->fun);
if(f_type==NULL||f_type->kind!=function_type)
throw program_error("Call of "
+std::string(f_type==NULL?"unknown ":"non-")
+"function '"
+main_hash_table->name_of(e.e.call_variant->fun)
+"'");


if(*f_type->func->result_type!=t)
actual=copy(f_type->func->result_type);
else check_type(*f_type->func->arg_type,e.e.call_variant->arg);
}
break;/*:94*//*123:*/
#line 2043 "evaluator.w"
case applied_identifier:
{type_declarator*it=global_id_table->type_of(e.e.identifier_variant);
if(it==NULL)throw program_error
("Undefined identifier "
+main_hash_table->name_of(e.e.identifier_variant));

if(*it!=t)
{static type_declarator row_int= *make_type("[int]").release();
static type_declarator row_row_int= *make_type("[[int]]").release();
static type_declarator row_vec= *make_type("[vec]").release();
if(t==vect&& *it==row_int)/*44:*/
#line 787 "evaluator.w"
{expr_list arg=make_exprlist_node(e,null_expr_list);
e=make_application_node
(main_hash_table->match_literal(">vec<[int]:"),arg);
}/*:44*/
#line 2055 "evaluator.w"
else if(t==matr&&(*it==row_row_int|| *it==row_vec))/*45:*/
#line 795 "evaluator.w"
{expr_list arg=make_exprlist_node(e,null_expr_list);
e=make_application_node
(main_hash_table->match_literal(">mat<[vec]:"),arg);
}/*:45*/
#line 2057 "evaluator.w"
else actual=copy(it);
}
}
break;/*:123*/
#line 711 "evaluator.w"
}
if(actual.get()!=NULL)
{type_ptr p=copy(&t);
throw type_error(e,*actual.release(),*p.release());

}
}/*:39*//*48:*/
#line 841 "evaluator.w"
latticetypes::Weight row_to_weight(const row_value*r)
{latticetypes::Weight result(r->value.size());
for(size_t i=0;i<r->value.size();++i)
{int_value*n=dynamic_cast<int_value* >(r->value[i]);
if(n==NULL)
throw std::logic_error("Row display entry failed to be integral");

result[i]=n->value;
}
return result;
}

latticetypes::Weight cast_intlist_to_weight(const value_ptr v)
{row_ptr r(dynamic_cast<row_value* >(v));
if(r.get()==NULL)
throw std::logic_error("Row display failed to return a row");

return row_to_weight(r.get());
}/*:48*//*49:*/
#line 869 "evaluator.w"
latticetypes::LatticeMatrix cast_intlistlist_to_latmat(const value_ptr v)
{row_ptr rr(dynamic_cast<row_value* >(v));

if(rr.get()==NULL)
throw std::logic_error("Matrix display failed to return row");

latticetypes::WeightList res_vec(rr->value.size());
size_t depth=0;
for(size_t i=0;i<rr->value.size();++i)
{row_value*r=dynamic_cast<row_value* >(rr->value[i]);
if(r!=NULL)res_vec[i]=row_to_weight(r);
else
{weight_value*lambda=dynamic_cast<weight_value* >(rr->value[i]);
if(lambda==NULL)
throw std::logic_error("Matrix column failed to be a weight");

res_vec[i]=lambda->value;
}
if(res_vec[i].size()>depth)depth=res_vec[i].size();
}
for(size_t i=0;i<res_vec.size();++i)
if(res_vec[i].size()<depth)
{size_t j=res_vec[i].size();
res_vec[i].resize(depth);
for(;j<depth;++j)res_vec[i][j]=0;
}
return latticetypes::LatticeMatrix(res_vec);
}/*:49*//*53:*/
#line 940 "evaluator.w"
void weight_value::print(std::ostream&out)const
{size_t l=value.size(),w=0;std::vector<std::string>tmp(l);
for(size_t i=0;i<l;++i)
{std::ostringstream s;s<<value[i];tmp[i]=s.str();
if(tmp[i].length()>w)w=tmp[i].length();
}
if(l==0)out<<"[ ]";
else
{w+=1;out<<std::right<<'[';
for(size_t i=0;i<l;++i)
out<<std::setw(w)<<tmp[i]<<(i<l-1?",":" ]");
}
}/*:53*//*54:*/
#line 956 "evaluator.w"
void latmat_value::print(std::ostream&out)const
{size_t k=value.numRows(),l=value.numColumns();
std::vector<size_t>w(l,0);
for(size_t i=0;i<k;++i)
for(size_t j=0;j<l;++j)
{std::ostringstream s;s<<value(i,j);size_t len=s.str().length();
if(len>w[j])w[j]=len;
}
out<<std::endl<<std::right;
for(size_t i=0;i<k;++i)
{out<<'|';
for(size_t j=0;j<l;++j)
out<<std::setw(w[j]+1)<<value(i,j)<<(j<l-1?',':' ');
out<<'|'<<std::endl;
}
}/*:54*//*62:*/
#line 1041 "evaluator.w"
type_ptr make_tuple_type(type_list_ptr l)
{type_declarator*p=new type_declarator(l.get());l.release();
return type_ptr(p);
}/*:62*//*67:*/
#line 1084 "evaluator.w"
type_list_ptr find_type_list(expr_list l)
{if(l==NULL)return type_list_ptr(NULL);
type_ptr t=find_type(l->e);
type_list_ptr tail=find_type_list(l->next);
type_node*result=new type_node(t.get(),tail.get());
t.release(),tail.release();
return type_list_ptr(result);
}/*:67*//*69:*/
#line 1113 "evaluator.w"
void tuple_value::print(std::ostream&out)const
{if(value.empty())out<<"()";
else
{out<<'(';
std::vector<value_ptr>::const_iterator p=value.begin();
do{(*p)->print(out);++p;out<<(p==value.end()?')':',');}
while(p!=value.end());
}
}/*:69*//*79:*/
#line 1256 "evaluator.w"
type_ptr make_function_type(type_ptr a,type_ptr r)
{type_declarator*p=new type_declarator(a.get(),r.get());
a.release(),r.release();
return type_ptr(p);
}/*:79*//*86:*/
#line 1346 "evaluator.w"
Id_table::~Id_table()
{for(map_type::const_iterator p=table.begin();p!=table.end();++p)
{delete p->second.value;delete p->second.type;}
}/*:86*//*87:*/
#line 1359 "evaluator.w"
void Id_table::add(Hash_table::id_type id,value_ptr v,type_ptr t)
{id_data data(v,t.get());std::auto_ptr<value_base>safe(v);
std::pair<map_type::iterator,bool>trial
=table.insert(std::make_pair(id,data));
safe.release(),t.release();
if(!trial.second)
{delete trial.first->second.value;delete trial.first->second.type;
trial.first->second=data;
}
}/*:87*//*88:*/
#line 1374 "evaluator.w"
type_declarator*Id_table::type_of(Hash_table::id_type id)const
{map_type::const_iterator p=table.find(id);
return p==table.end()?NULL:p->second.type;
}
value_ptr Id_table::value_of(Hash_table::id_type id)const
{map_type::const_iterator p=table.find(id);
return p==table.end()?NULL:p->second.value;
}/*:88*//*89:*/
#line 1386 "evaluator.w"
void Id_table::print(std::ostream&out)const
{for(map_type::const_iterator p=table.begin();p!=table.end();++p)
out<<main_hash_table->name_of(p->first)<<": "
<< *p->second.type<<": "<< *p->second.value<<std::endl;
}

std::ostream&operator<<(std::ostream&out,const Id_table&p)
{p.print(out);return out;}/*:89*//*98:*/
#line 1494 "evaluator.w"
void clear_execution_stack()
{if(!execution_stack.empty())
{std::cerr<<"Discarding from execution stack:"<<std::endl;
do
{value_ptr v=execution_stack.back();
std::cerr<< *v<<std::endl;
delete v;execution_stack.pop_back();
}
while(!execution_stack.empty());
}
}/*:98*//*105:*/
#line 1660 "evaluator.w"
void push_tuple_components()
{std::auto_ptr<tuple_value>tuple(get_tuple());
for(size_t i=0;i<tuple->length();++i)
{push_value(tuple->value[i]);
tuple->value[i]=NULL;
}
}

void wrap_tuple(size_t n)
{tuple_value*result=new tuple_value(std::vector<value_ptr>(n));
while(n-- >0)
result->value[n]=pop_arg();
push_value(result);
}/*:105*//*106:*/
#line 1680 "evaluator.w"
void intlist_to_weight_wrapper()
{push_value(new weight_value(cast_intlist_to_weight(pop_arg())));
}

void intlistlist_to_latmat_wrapper()
{push_value(new latmat_value(cast_intlistlist_to_latmat(pop_arg())));
}

void id_wrapper(){}/*103:*/
#line 1589 "evaluator.w"
int_value*get_int()throw(std::logic_error)
{value_ptr p=pop_arg();
int_value*result=dynamic_cast<int_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not an integer");}
return result;
}

string_value*get_string()throw(std::logic_error)
{value_ptr p=pop_arg();
string_value*result=dynamic_cast<string_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a string");}
return result;
}

bool_value*get_bool()throw(std::logic_error)
{value_ptr p=pop_arg();
bool_value*result=dynamic_cast<bool_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a boolean");}
return result;
}


weight_value*get_vec()throw(std::logic_error)
{value_ptr p=pop_arg();
weight_value*result=dynamic_cast<weight_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a vector");}
return result;
}

latmat_value*get_mat()throw(std::logic_error)
{value_ptr p=pop_arg();
latmat_value*result=dynamic_cast<latmat_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a matrix");}
return result;
}

row_value*get_row()throw(std::logic_error)
{value_ptr p=pop_arg();
row_value*result=dynamic_cast<row_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a row");}
return result;
}

tuple_value*get_tuple()throw(std::logic_error)
{value_ptr p=pop_arg();
tuple_value*result=dynamic_cast<tuple_value* >(p);
if(result==NULL)
{delete p;throw std::logic_error("Argument is not a tuple");}
return result;
}/*:103*//*111:*/
#line 1771 "evaluator.w"
void plus_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>j(get_int());std::auto_ptr<int_value>i(get_int());
i->value+=j->value;
push_value(i.release());
}
void minus_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>j(get_int());std::auto_ptr<int_value>i(get_int());
i->value-=j->value;
push_value(i.release());
}

void times_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>j(get_int());std::auto_ptr<int_value>i(get_int());
i->value*=j->value;
push_value(i.release());
}

void divide_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>j(get_int());std::auto_ptr<int_value>i(get_int());
if(j->value==0)throw std::runtime_error("Division by zero");
i->value/=j->value;
push_value(i.release());
}

void modulo_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>j(get_int());std::auto_ptr<int_value>i(get_int());
if(j->value==0)throw std::runtime_error("Modulo zero");
i->value%=j->value;
push_value(i.release());
}

void unary_minus_wrapper()
{int_value*i=get_int();i->value= -i->value;push_value(i);}

void divmod_wrapper()
{push_tuple_components();
std::auto_ptr<int_value>j(get_int());std::auto_ptr<int_value>i(get_int());
if(j->value==0)throw std::runtime_error("DivMod by zero");
int mod=i->value%j->value;
i->value/=j->value;j->value=mod;
push_value(i.release());push_value(j.release());wrap_tuple(2);
}/*:111*//*112:*/
#line 1823 "evaluator.w"
void print_wrapper()
{std::auto_ptr<string_value>s(get_string());
*output_stream<<s->value<<std::endl;
wrap_tuple(0);
}/*:112*//*115:*/
#line 1862 "evaluator.w"
void id_mat_wrapper()
{std::auto_ptr<int_value>i(get_int());
std::auto_ptr<latmat_value>m
(new latmat_value(latticetypes::LatticeMatrix()));
identityMatrix(m->value,std::abs(i->value));push_value(m.release());
}

void transpose_mat_wrapper()
{std::auto_ptr<latmat_value>m(get_mat());
m->value.transpose();push_value(m.release());
}

void diagonal_wrapper()
{std::auto_ptr<weight_value>d(get_vec());
size_t n=d->value.size();
std::auto_ptr<latmat_value>m
(new latmat_value(latticetypes::LatticeMatrix(n,n,0)));
for(size_t i=0;i<n;++i)m->value(i,i)=d->value[i];
push_value(m.release());
}/*:115*//*117:*/
#line 1900 "evaluator.w"
void mv_prod_wrapper()
{push_tuple_components();
std::auto_ptr<weight_value>v(get_vec());
std::auto_ptr<latmat_value>m(get_mat());
if(m->value.numColumns()!=v->value.size())
{std::ostringstream s;
s<<"Size mismatch "<<m->value.numColumns()<<":"<<v->value.size()
<<" in mv_prod";
throw std::runtime_error(s.str());
}
std::auto_ptr<weight_value>w
(new weight_value(latticetypes::Weight(m->value.numRows())));
m->value.apply(w->value,v->value);
push_value(w.release());
}

void mm_prod_wrapper()
{push_tuple_components();
std::auto_ptr<latmat_value>r(get_mat());
std::auto_ptr<latmat_value>l(get_mat());
if(l->value.numColumns()!=r->value.numRows())
{std::ostringstream s;
s<<"Size mismatch "<<l->value.numColumns()<<":"<<r->value.numRows()
<<" in mm_prod";
throw std::runtime_error(s.str());
}
l->value*=r->value;
push_value(l.release());
}/*:117*//*118:*/
#line 1938 "evaluator.w"
void invfact_wrapper()
{std::auto_ptr<latmat_value>m(get_mat());
size_t nr=m->value.numRows();
latticetypes::WeightList b;matrix::initBasis(b,nr);
std::auto_ptr<weight_value>inv_factors
(new weight_value(latticetypes::Weight(0)));
smithnormal::smithNormal(inv_factors->value,b.begin(),m->value);
push_value(inv_factors.release());
}
void Smith_basis_wrapper()
{std::auto_ptr<latmat_value>m(get_mat());
size_t nr=m->value.numRows();
latticetypes::WeightList b;matrix::initBasis(b,nr);
latticetypes::Weight inv_factors(0);
smithnormal::smithNormal(inv_factors,b.begin(),m->value);
latticetypes::LatticeMatrix new_basis(b);
m->value.swap(new_basis);
push_value(m.release());
}

void Smith_wrapper()
{std::auto_ptr<latmat_value>m(get_mat());
size_t nr=m->value.numRows();
latticetypes::WeightList b;matrix::initBasis(b,nr);
std::auto_ptr<weight_value>inv_factors
(new weight_value(latticetypes::Weight(0)));
smithnormal::smithNormal(inv_factors->value,b.begin(),m->value);
latticetypes::LatticeMatrix new_basis(b);
m->value.swap(new_basis);
push_value(m.release());push_value(inv_factors.release());wrap_tuple(2);
}/*:118*//*119:*/
#line 1975 "evaluator.w"
void invert_wrapper()
{std::auto_ptr<latmat_value>m(get_mat());
if(m->value.numRows()!=m->value.numColumns())
{std::ostringstream s;
s<<"Cannot invert a "
<<m->value.numRows()<<"x"<<m->value.numColumns()<<" matrix";
throw std::runtime_error(s.str());
}
std::auto_ptr<int_value>denom(new int_value(0));
m->value.invert(denom->value);
push_value(m.release());push_value(denom.release());wrap_tuple(2);
}/*:119*//*:106*//*108:*/
#line 1709 "evaluator.w"
void initialise_evaluator()
{execution_stack.reserve(16);

install_function(intlist_to_weight_wrapper,">vec<[int]:","([int]->vec)");
install_function(intlistlist_to_latmat_wrapper
,">mat<[vec]:","([vec]->mat)");
install_function(id_wrapper,"vec","(vec->vec)");
install_function(id_wrapper,"mat","(mat->mat)");/*113:*/
#line 1833 "evaluator.w"
install_function(plus_wrapper,"+","(int,int->int)");
install_function(minus_wrapper,"-","(int,int->int)");
install_function(times_wrapper,"*","(int,int->int)");
install_function(divide_wrapper,"/","(int,int->int)");
install_function(modulo_wrapper,"%","(int,int->int)");
install_function(unary_minus_wrapper,"-u","(int->int)");
install_function(divmod_wrapper,"/%","(int,int->int,int)");
install_function(print_wrapper,"print","(string->)");/*:113*//*120:*/
#line 1990 "evaluator.w"
install_function(id_mat_wrapper,"id_mat","(int->mat)");
install_function(transpose_mat_wrapper,"transpose_mat","(mat->mat)");
install_function(diagonal_wrapper,"diagonal_mat","(vec->mat)");
install_function(mv_prod_wrapper,"mv_prod","(mat,vec->vec)");
install_function(mm_prod_wrapper,"mm_prod","(mat,mat->mat)");
install_function(invfact_wrapper,"inv_fact","(mat->vec)");
install_function(Smith_basis_wrapper,"Smith_basis","(mat->mat)");
install_function(Smith_wrapper,"Smith","(mat->mat,vec)");
install_function(invert_wrapper,"invert","(mat->mat,int)");/*:120*/
#line 1718 "evaluator.w"
}/*:108*//*110:*/
#line 1742 "evaluator.w"
type_ptr analyse_types(expr&e)throw(std::bad_alloc,std::runtime_error)
{try{return find_type(e);}
catch(type_error&err)
{std::cerr<<err.what()<<std::endl<<
"Subexpression "<<err.offender<<" has wrong type: found "
<< *err.actual<<" while "<< *err.required<<" was needed.\n";

}
catch(program_error&err)
{std::cerr<<err.what()<<
" in expression '"<<e<<"'\n";
}
throw std::runtime_error("Type check failed");

}/*:110*//*124:*/
#line 2071 "evaluator.w"
type_ptr scan_type(const char* &s);
type_list_ptr scan_type_list(const char* &s);

type_ptr make_type(const char*s)
{try{return scan_type(s);}
catch(std::logic_error e)
{std::cerr<<e.what()<<"; text remaining: "<<s<<std::endl;;
return make_prim_type(string_type);
}
}

type_ptr scan_type(const char* &s)
{if(*s=='[')/*125:*/
#line 2096 "evaluator.w"
{type_ptr c=scan_type(++s);
if(*s++ !=']')throw std::logic_error("Missing ']' in type");
return make_row_type(c);
}/*:125*/
#line 2085 "evaluator.w"
if(*s=='(')/*126:*/
#line 2121 "evaluator.w"
{type_list_ptr l0=scan_type_list(++s),l1(NULL);
bool is_tuple= *s==')';
if(*s=='-'&& * ++s=='>')l1=scan_type_list(++s);
if(*s++ !=')')throw std::logic_error("Missing ')' in type");
type_ptr a(NULL);
if(l0.get()!=NULL&&l0->next==NULL){a=type_ptr(l0->t);l0->t=NULL;}
else a=make_tuple_type(l0);
if(is_tuple)return a;
type_ptr r(NULL);
if(l1.get()!=NULL&&l1->next==NULL){r=type_ptr(l1->t);l1->t=NULL;}
else r=make_tuple_type(l1);
return make_function_type(a,r);
}/*:126*//*128:*/
#line 2152 "evaluator.w"
{for(size_t i=0;i<nr_of_primitive_types;++i)
{std::string name=prim_names[i];
if(name.compare(0,name.length(),s,name.length())==0)
{s+=name.length();
return make_prim_type(static_cast<primitive>(i));
}
}
throw std::logic_error("Type unrecognised");
}/*:128*/
#line 2089 "evaluator.w"
}/*:124*//*127:*/
#line 2140 "evaluator.w"

type_list_ptr scan_type_list(const char* &s)
{if(*s==')' or*s=='-')return type_list_ptr(NULL);
type_ptr head=scan_type(s);
if(*s!=',')return make_type_singleton(head);
return make_type_list(head,scan_type_list(++s));
}/*:127*//*133:*/
#line 2230 "evaluator.w"
extern "C"
void global_set_identifier(expr_list ids,expr e)
{using namespace atlas::interpreter;using namespace std;
try
{type_ptr t=analyse_types(e);
if(ids->next!=NULL)/*134:*/
#line 2263 "evaluator.w"
{if(t->kind!=tuple_type)
throw runtime_error("Multiple assignment requires tuple value");

set<Hash_table::id_type>seen;
type_list tl=t->tuple;
for(expr_list l=ids;l!=NULL||tl!=NULL;l=l->next,tl=tl->next)
{if(l==NULL||tl==NULL)
throw runtime_error
("Right hand side has too "+string(l==NULL?"many":"few")
+" components");


if(!seen.insert(l->e.e.identifier_variant).second)

throw runtime_error
("Repeated identifier "
+main_hash_table->name_of(l->e.e.identifier_variant)
+" in multiple assignment");
}
}/*:134*/
#line 2238 "evaluator.w"
value_ptr v=evaluate(e);
if(ids->next==NULL)
{cout<<"Identifier "<<ids->e<<": "<< *t<<std::endl;
global_id_table->add(ids->e.e.identifier_variant,v,t);
}
else/*135:*/
#line 2294 "evaluator.w"
{auto_ptr<tuple_value>tv(dynamic_cast<tuple_value* >(v));
if(tv.get()==NULL)throw logic_error("Non-tuple value assigned");

cout<<"Identifiers ";
size_t i=0;type_list tl=t->tuple;
for(expr_list l=ids;l!=NULL;l=l->next,++i,tl=tl->next)
{cout<<l->e<<": "<< *tl->t<<(l->next!=NULL?", ":".\n");
global_id_table->
add(l->e.e.identifier_variant,tv->value[i],copy(tl->t));
tv->value[i]=NULL;
}
}/*:135*/
#line 2244 "evaluator.w"
}
catch(runtime_error&err)
{cerr<<err.what()<<", identifier"<<(ids->next!=NULL?"s ":" ");
for(expr_list l=ids;l!=NULL;l=l->next)
cerr<<main_hash_table->name_of(l->e.e.identifier_variant)
<<(l->next!=NULL?",":"");
cerr<<" not assigned to.\n";
clear_execution_stack();
}
}/*:133*//*138:*/
#line 2329 "evaluator.w"
extern "C"
void type_of_expr(expr e)
{try{*output_stream<<"type: "<< *analyse_types(e)<<std::endl;}
catch(std::runtime_error&err){std::cerr<<err.what()<<std::endl;}
}/*:138*//*139:*/
#line 2342 "evaluator.w"
extern "C"
void show_ids()
{*output_stream<< *global_id_table;
}/*:139*/
#line 13 "evaluator.w"

}}/*:1*/
