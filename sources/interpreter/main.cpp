#define axis_version "1.0"
#include "parse_types.h"
#include "parser.tab.h"
#include <iostream>
#include <fstream>
#include "buffer.h"
#include "lexer.h"
#include "../version.h"
#include <csignal>
#include "atlas-types.h"
#include "global.h"
#include <cstring>
#include <string>
#include <stdexcept>
#include "axis.h"
#include "parsetree.h"
#include <cstdlib>
#include "sl_list.h"
/*2:*//*6:*/
#line 233 "main.w"
#ifdef NREADLINE
#define readline nullptr
#define add_history nullptr
#define clear_history()
#else
#include <readline/readline.h>
#include <readline/history.h>
#endif
#ifdef NOT_UNIX
#define isatty(x) true 
#define chdir(x) (-1) 
#else
#include <unistd.h> 
#endif/*:6*//*3:*/
#line 179 "main.w"
int yylex(YYSTYPE*,YYLTYPE*);
int yyparse(atlas::interpreter::expr_p*parsed_expr,int*verbosity);/*:3*/
#line 147 "main.w"

namespace{/*7:*/
#line 257 "main.w"
const char*keywords[]=
{"quit"
,"set","let","in","begin","end"
,"if","then","else","elif","fi"
,"and","or","not"
,"next","do","dont","from","downto","while","for","od"
,"case","esac","rec_fun"
,"true","false","die","break","return"
,"set_type"
,"whattype","showall","forget"
,nullptr};/*:7*//*12:*/
#line 374 "main.w"
static atlas::interpreter::shared_share
input_path_pointer,prelude_log_pointer;
#ifndef NOT_UNIX
enum{cwd_size=0x1000};
char cwd_buffer[cwd_size];
const char*working_directory_name=nullptr;
const char*first_path=nullptr;
#endif/*:12*//*23:*/
#line 725 "main.w"
char lexical_break_chars[]=" \t\n=<>+-*/\\%,;:.()[]{}#!?$@\"|~";/*:23*/
#line 148 "main.w"
}/*4:*/
#line 191 "main.w"
int yylex(YYSTYPE*valp,YYLTYPE*locp)
{return atlas::interpreter::lex->get_token(valp,locp);}


void yyerror(YYLTYPE*locp,atlas::interpreter::expr_p*,int*,char const*s)
{atlas::interpreter::main_input_buffer->show_range
(std::cerr,
locp->first_line,locp->first_column,
locp->last_line,locp->last_column);
std::cerr<<s<<std::endl;
atlas::interpreter::main_input_buffer->close_includes();
atlas::interpreter::clean=false;
}/*:4*//*5:*/
#line 210 "main.w"
extern "C" void sigint_handler(int)
{atlas::interpreter::interrupt_flag=1;}/*:5*//*22:*/
#line 687 "main.w"
extern "C"
char*id_completion_func(const char*text,int state)
{using namespace atlas;
char*result=nullptr;
static containers::sl_list<const char* >comps;
if(state==0)
comps=interpreter::completions(text);
if(comps.empty())
return result;
{result=static_cast<char* >(std::malloc(std::strlen(comps.front())+1));
if(result!=nullptr)
std::strcpy(result,comps.front());
}
comps.pop_front();
return result;
}/*:22*//*24:*/
#line 750 "main.w"
#ifndef NREADLINE
extern "C" char* *do_completion(const char*text,int start,int end)
{
if(atlas::interpreter::lex->is_initial()and start>0)

{int i;char c;bool need_file=false;bool in=true;
for(i=0;i<start;++i)
if(std::isspace(c=rl_line_buffer[i]))
continue;
else if(c=='<')
need_file=true;
else if(c=='>')
need_file=true,in=false;
else
break;

if(working_directory_name!=nullptr and need_file and i==start)

{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
chdir(in and first_path!=nullptr?first_path:working_directory_name);

#pragma GCC diagnostic pop
return nullptr;
}
}
rl_attempted_completion_over=true;

return rl_completion_matches(text,id_completion_func);

}
#endif/*:24*/
#line 149 "main.w"

namespace atlas
{namespace interpreter
{/*15:*/
#line 443 "main.w"
unsigned int input_path_size()
{const row_value*path=force<row_value>(input_path_pointer->get());
return path->val.size();
}
const std::string&input_path_component(unsigned int i)
{const row_value*path=force<row_value>(input_path_pointer->get());
const string_value*dir=force<string_value>(path->val[i].get());
return dir->val;
}/*:15*/
#line 152 "main.w"
}
}/*10:*/
#line 313 "main.w"
int main(int argc,char* *argv)
{using namespace atlas::interpreter;


Hash_table hash;main_hash_table= &hash;
Id_table main_table;global_id_table= &main_table;
overload_table main_overload_table;
global_overload_table= &main_overload_table;

bool use_readline=true,do_prompting=isatty(STDIN_FILENO);/*11:*/
#line 357 "main.w"
std::vector<const char* >paths,prelude_filenames;
paths.reserve(argc-1);prelude_filenames.reserve(argc-1);/*:11*//*13:*/
#line 393 "main.w"
while(* ++argv!=nullptr)
{static const char*const path_opt="--path=";
static const size_t pol=std::strlen(path_opt);
std::string arg(*argv);
if(arg=="--no-readline")
{use_readline=false;continue;}
if(arg.substr(0,pol)==path_opt)
paths.push_back(&(*argv)[pol]);
else prelude_filenames.push_back(*argv);
}/*:13*/
#line 325 "main.w"
BufferedInput input_buffer(do_prompting?"atlas> ":nullptr
,use_readline?readline:nullptr
,use_readline?add_history:nullptr);
main_input_buffer= &input_buffer;
Lexical_analyser ana(input_buffer,hash,keywords,prim_names);lex= &ana;/*8:*/
#line 282 "main.w"
main_hash_table->match_literal("quiet");
main_hash_table->match_literal("verbose");

ana.set_comment_delims('{','}');/*:8*//*9:*//*25:*/
#line 796 "main.w"
#ifndef NOT_UNIX
working_directory_name=getcwd(cwd_buffer,cwd_size);
#endif/*:25*//*27:*/
#line 821 "main.w"
#ifndef NREADLINE
rl_readline_name="axis";
using_history();
rl_completer_word_break_characters=lexical_break_chars;
rl_attempted_completion_function=do_completion;

#endif/*:27*/
#line 299 "main.w"
signal(SIGINT,sigint_handler);
initialise_evaluator();
initialise_builtin_types();/*:9*//*14:*/
#line 411 "main.w"
{/*16:*/
#line 460 "main.w"
#ifndef NOT_UNIX
if(paths.size()>0)
{if(working_directory_name!=nullptr and chdir(paths[0])==0)

{if(chdir(working_directory_name)!=0)
return EXIT_FAILURE;
first_path=paths[0];
}
else
std::cerr<<"Warning: specified path '"<<paths[0]
<<"' is not a directory."<<std::endl;
}
#endif/*:16*/
#line 412 "main.w"
own_value input_path=std::make_shared<row_value>(paths.size());
id_type ip_id=main_hash_table->match_literal("input_path");
auto oit=force<row_value>(input_path.get())->val.begin();
for(auto it=paths.begin();it!=paths.end();++it)
*oit++=std::make_shared<string_value>(std::string(*it)+'/');
global_id_table->add(ip_id
,input_path,mk_type_expr("[string]"),false);
input_path_pointer=global_id_table->address_of(ip_id);

own_value prelude_log=std::make_shared<row_value>(0);
id_type pl_id=main_hash_table->match_literal("prelude_log");
auto&logs=force<row_value>(input_path.get())->val;
logs.reserve(prelude_filenames.size());
global_id_table->add(pl_id
,prelude_log,mk_type_expr("[string]"),true);
prelude_log_pointer=global_id_table->address_of(pl_id);

own_value back_trace=std::make_shared<row_value>(0);
id_type bt_id=main_hash_table->match_literal("back_trace");
global_id_table->add(bt_id
,back_trace,mk_type_expr("[string]"),false);
back_trace_pointer=global_id_table->address_of(bt_id);
}/*:14*//*19:*/
#line 567 "main.w"
for(auto it=prelude_filenames.begin();it!=prelude_filenames.end();++it)
{std::ostringstream log_stream;output_stream= &log_stream;
main_input_buffer->push_file(*it,true);

while(main_input_buffer->include_depth()>0)
{if(not ana.reset())
{std::cerr<<"Internal error, getline fails reading "<< *it
<<std::endl;
return EXIT_FAILURE;
}
expr_p parse_tree;
if(yyparse(&parse_tree,&verbosity)!=0)
continue;
if(verbosity!=0)
{std::cerr<<"Cannot "
<<(verbosity<0?"quit":
verbosity==1?"set verbose":"redirect output")
<<" during prelude.\n";
verbosity=0;main_input_buffer->close_includes();
}
else
{try
{expression_ptr e;type_expr found_type=analyse_types(*parse_tree,e);
e->evaluate(expression_base::single_value);
if(found_type!=void_type)
log_stream<<"Value: "<< *pop_value()<<'\n';
else
pop_value();
}
catch(error_base&err)
{std::cerr<<err.message<<"\nEvaluation aborted.\n";
clean=false;
reset_evaluator();main_input_buffer->close_includes();
}
catch(std::exception&err)
{
std::cerr<<err.what()<<"\nEvaluation aborted.\n";
clean=false;
reset_evaluator();main_input_buffer->close_includes();
}
}
destroy_expr(parse_tree);
}
row_value*logs=uniquify<row_value>(*prelude_log_pointer);
logs->val.emplace_back(std::make_shared<string_value>(log_stream.str()));
output_stream= &std::cout;
}/*:19*/
#line 335 "main.w"
if(do_prompting)
std::cout<<"This is 'atlas' (version "
<<atlas::version::VERSION<<", axis language version "
axis_version "),\n"<<atlas::version::NAME<<
" interpreter,\ncompiled on "<<atlas::version::COMPILEDATE
<<(atlas::version::debugging?", with debugging":"")
<<", readline "
<<(atlas::version::with_readline?"enabled":"disabled")
<<".\nhttp://www.liegroups.org/\n";/*17:*/
#line 485 "main.w"
last_value=shared_value(new tuple_value(0));
last_type=void_type.copy();

while(ana.reset())
{/*26:*/
#line 806 "main.w"
#ifndef NOT_UNIX
if(working_directory_name!=nullptr)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
chdir(working_directory_name);
#pragma GCC diagnostic pop
#endif/*:26*/
#line 490 "main.w"
expr_p parse_tree;
int old_verbosity=verbosity;
std::ofstream redirect;
try
{if(yyparse(&parse_tree,&verbosity)!=0)

continue;
if(verbosity!=0)
{if(verbosity<0)
break;
if(verbosity==2 or verbosity==3)

{auto write_mode=
verbosity==2?std::ios_base::trunc:std::ios_base::app;
verbosity=old_verbosity;/*20:*/
#line 621 "main.w"
{redirect.open(ana.scanned_file_name(),std::ios_base::out|write_mode);
if(redirect.is_open())
output_stream= &redirect;
else
{std::cerr<<"Failed to open "<<ana.scanned_file_name()<<std::endl;
continue;
}
}/*:20*/
#line 508 "main.w"
}
if(verbosity==1)
std::cout<<"Expression before type analysis: "<< *parse_tree
<<std::endl;
}/*18:*/
#line 528 "main.w"
{expression_ptr e;
type_expr found_type=analyse_types(*parse_tree,e);
if(verbosity>0)
std::cout<<"Type found: "<<found_type<<std::endl
<<"Converted expression: "<< *e<<std::endl;
if(found_type==void_type)
e->void_eval();
else
{e->eval();
last_type=std::move(found_type);
last_value=pop_value();
*output_stream<<"Value: "<< *last_value<<std::endl;
}
destroy_expr(parse_tree);
}/*:18*/
#line 514 "main.w"
}/*21:*/
#line 637 "main.w"
catch(error_base&err)
{if(dynamic_cast<runtime_error* >(&err)!=nullptr)
std::cerr<<"Runtime error:\n  ";
else if(dynamic_cast<logic_error* >(&err)!=nullptr)
std::cerr<<"Internal error: ";
std::cerr<<err.message<<"\n";

set_back_trace(err.back_trace);

std::cerr<<"Evaluation aborted.\n";
clean=false;
reset_evaluator();main_input_buffer->close_includes();
}
catch(std::exception&err)
{std::cerr<<err.what()<<"\nEvaluation aborted.\n";
clean=false;
reset_evaluator();main_input_buffer->close_includes();
}/*:21*/
#line 516 "main.w"
output_stream= &std::cout;
}/*:17*/
#line 346 "main.w"
signal(SIGINT,SIG_DFL);/*28:*/
#line 833 "main.w"
#ifndef NREADLINE
clear_history();

#endif/*:28*/
#line 348 "main.w"
std::cout<<"Bye.\n";
return clean?EXIT_SUCCESS:EXIT_FAILURE;
}/*:10*//*:2*/
