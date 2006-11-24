#include "parsetree.h"
#include "parser.tab.h"
#include "buffer.h"
#include "lexer.h"
#include <iostream>
#include <fstream>
#include <readline/readline.h>
#include <readline/history.h>
#include "built-in-types.h"
#include "constants.h"
#include <stdexcept>
#include "evaluator.h"
#include <cstring>
/*1:*//*2:*/
#line 33 "main.w"
extern "C"
{int yylex(YYSTYPE*,YYLTYPE*);
int yyparse(expr*parsed_expr,int*verbosity);
}/*:2*/
#line 10 "main.w"

namespace{/*3:*/
#line 46 "main.w"
const char*keywords[]=
{"quit","true","false","quiet","verbose","whattype","showall",NULL};/*:3*//*4:*//*:4*/
#line 11 "main.w"
}/*5:*/
#line 68 "main.w"
extern "C"
int yylex(YYSTYPE*valp,YYLTYPE*locp)
{return atlas::interpreter::lex->get_token(valp,locp);}

extern "C"
void yyerror(YYLTYPE*locp,expr*parsed_expr,int*verbosity,char const*s)
{atlas::interpreter::main_input_buffer->show_range(std::cerr,
locp->first_line,locp->first_column,
locp->last_line,locp->last_column);
std::cerr<<s<<std::endl;
atlas::interpreter::main_input_buffer->close_includes();
}/*:5*//*6:*/
#line 98 "main.w"
int main(int argc,char* *argv)
{using namespace std;using namespace atlas::interpreter;/*9:*/
#line 196 "main.w"
bool use_readline=argc<2 or std::strcmp(argv[1],"-nr")!=0;/*:9*/
#line 103 "main.w"
BufferedInput input_buffer("expr> "
,use_readline?readline:NULL
,use_readline?add_history:NULL);
main_input_buffer= &input_buffer;
Hash_table hash;main_hash_table= &hash;
Lexical_analyser ana(input_buffer,hash,keywords);lex= &ana;
Id_table main_table;global_id_table= &main_table;/*7:*/
#line 148 "main.w"
using_history();
rl_completion_entry_function=id_completion_func;

atlas::constants::initConstants();

initialise_evaluator();initialise_builtin_types();/*:7*/
#line 112 "main.w"

cout<<"Enter expressions:\n";
while(ana.reset())
{expr expression;
int old_verbosity=verbosity;
ofstream redirect;
if(yyparse(&expression,&verbosity))
continue;
if(verbosity!=0)
{if(verbosity<0)break;
if(verbosity==2 or verbosity==3)

{/*10:*/
#line 205 "main.w"
{redirect.open(ana.scanned_file_name(),ios_base::out|
(verbosity==2?ios_base::trunc:ios_base::app));
if(redirect.is_open())output_stream= &redirect;
else
{cerr<<"Failed to open "<<ana.scanned_file_name()<<endl;
continue;
}
}/*:10*/
#line 126 "main.w"
verbosity=old_verbosity;
}
if(verbosity==1)
cout<<"Expression before type analysis: "<<expression<<endl;
}/*8:*/
#line 164 "main.w"
{bool type_OK=false;
try
{type_ptr type=analyse_types(expression);
type_OK=true;
if(verbosity==1)
*output_stream<<"Type found: "<< *type<<endl
<<"Expression after type analysis: "<<expression<<endl;
value_ptr v=evaluate(expression);
static type_declarator empty= *make_type("()").release();
if(*type!=empty)*output_stream<<"Value: "<< *v<<endl;
destroy_expr(expression);delete v;
}
catch(runtime_error&err)
{if(type_OK)cerr<<"Runtime error: ";
cerr<<err.what()<<", evaluation aborted.\n";
clear_execution_stack();main_input_buffer->close_includes();
}
catch(logic_error&err)
{cerr<<"Unexpected error: "<<err.what()<<", evaluation aborted.\n";
clear_execution_stack();main_input_buffer->close_includes();
}
catch(exception&err)
{cerr<<err.what()<<", evaluation aborted.\n";
clear_execution_stack();main_input_buffer->close_includes();
}
}/*:8*/
#line 133 "main.w"
output_stream= &cout;
}
clear_history();

cout<<"Bye.\n";
}/*:6*//*:1*/
