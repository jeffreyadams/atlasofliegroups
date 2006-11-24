#include "parsetree.h"
#include <cstdio>
#include "lexer.h"
/*2:*/
#line 37 "parsetree.w"
extern "C"{/*7:*/
#line 93 "parsetree.w"
void destroy_expr(expr e)
{switch(e.kind)
{/*11:*/
#line 135 "parsetree.w"
case integer_denotation:case boolean_denotation:break;
case string_denotation:delete[]e.e.str_denotation_variant;break;/*:11*//*18:*/
#line 191 "parsetree.w"
case applied_identifier:break;/*:18*//*28:*/
#line 263 "parsetree.w"
case list_display:destroy_exprlist(e.e.sublist);
break;/*:28*//*33:*/
#line 329 "parsetree.w"
case tuple_display:destroy_exprlist(e.e.sublist);
break;/*:33*//*41:*/
#line 387 "parsetree.w"
case function_call:
destroy_expr(e.e.call_variant->arg);
delete e.e.call_variant;
break;/*:41*/
#line 96 "parsetree.w"
}
}/*:7*//*13:*/
#line 151 "parsetree.w"
expr make_int_denotation(int val)
{expr result;result.kind=integer_denotation;
result.e.int_denotation_variant=val;return result;
}
expr make_string_denotation(char*val)
{expr result;result.kind=string_denotation;
result.e.str_denotation_variant=val;return result;
}
expr make_bool_denotation(int val)
{expr result;result.kind=boolean_denotation;
result.e.int_denotation_variant=val;return result;
}/*:13*//*20:*/
#line 202 "parsetree.w"
expr make_applied_identifier(id_type id)
{expr result;result.kind=applied_identifier;
result.e.identifier_variant=id;return result;
}/*:20*//*27:*/
#line 255 "parsetree.w"
void destroy_exprlist(expr_list l)
{while(l!=NULL)
{destroy_expr(l->e);expr_list this_node=l;l=l->next;delete this_node;}
}/*:27*//*30:*/
#line 291 "parsetree.w"
const expr_list null_expr_list=NULL;
expr_list make_exprlist_node(expr e,expr_list l)
{expr_list n=new exprlist_node;n->e=e;n->next=l;return n;}
expr_list reverse_expr_list(expr_list l)
{expr_list r=NULL;
while(l!=NULL){expr_list t=l;l=t->next;t->next=r;r=t;}
return r;
}
expr wrap_list_display(expr_list l)
{expr result;result.kind=list_display;result.e.sublist=l;return result;
}/*:30*//*35:*/
#line 338 "parsetree.w"
expr wrap_tuple_display(expr_list l)
{expr result;result.kind=tuple_display;result.e.sublist=l;return result;
}/*:35*//*43:*/
#line 410 "parsetree.w"
expr make_application_node(id_type f,expr_list args)
{app a=new application_node;a->fun=f;
if(args!=NULL&&args->next==NULL)
{a->arg=args->e;delete args;}
else a->arg=wrap_tuple_display(args);
expr result;result.kind=function_call;result.e.call_variant=a;
return result;
}/*:43*//*45:*/
#line 431 "parsetree.w"
id_type lookup_identifier(const char*name)
{return atlas::interpreter::main_hash_table->match_literal(name);}/*:45*//*46:*/
#line 438 "parsetree.w"
void include_file()
{
atlas::interpreter::main_input_buffer->push_file
(atlas::interpreter::lex->scanned_file_name());
}/*:46*/
#line 38 "parsetree.w"

}
namespace atlas
{namespace interpreter
{/*5:*/
#line 74 "parsetree.w"
std::ostream&operator<<(std::ostream&out,expr e)
{switch(e.kind)
{/*10:*/
#line 125 "parsetree.w"
case integer_denotation:out<<e.e.int_denotation_variant;break;
case string_denotation:
out<<'"'<<e.e.str_denotation_variant<<'"';break;
case boolean_denotation:
out<<(e.e.int_denotation_variant?"true":"false");break;/*:10*//*17:*/
#line 184 "parsetree.w"
case applied_identifier:
out<<main_hash_table->name_of(e.e.identifier_variant);
break;/*:17*//*25:*/
#line 233 "parsetree.w"
case list_display:
{expr_list l=e.e.sublist;
if(l==NULL)out<<"[]";
else
{out<<'[';
do{out<<l->e;l=l->next;out<<(l==NULL?']':',');}
while(l!=NULL);
}
}
break;/*:25*//*32:*/
#line 315 "parsetree.w"
case tuple_display:
{expr_list l=e.e.sublist;
if(l==NULL)out<<"()";
else
{out<<'(';
do{out<<l->e;l=l->next;out<<(l==NULL?')':',');}
while(l!=NULL);
}
}
break;/*:32*//*40:*/
#line 374 "parsetree.w"
case function_call:
{app a=e.e.call_variant;
out<<main_hash_table->name_of(a->fun);
expr arg=a->arg;
if(arg.kind==tuple_display)out<<arg;
else out<<'('<<arg<<')';
}
break;/*:40*/
#line 77 "parsetree.w"
}
return out;
}/*:5*/
#line 43 "parsetree.w"

}
}/*:2*/
