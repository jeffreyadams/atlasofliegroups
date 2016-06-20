% Copyright (C) 2006-2015 Marc van Leeuwen
% This file is part of the Atlas of Lie Groups and Representations (the Atlas)

% This program is made available under the terms stated in the GNU
% General Public License (GPL), see http://www.gnu.org/licences/licence.html

% The Atlas is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% The Atlas is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with the Atlas; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


\def\emph#1{{\it#1\/}}
\def\point{\item{$\bullet$}}

\def\axis.{\.{axis}}

@* Introduction to the {\tt atlas} interpreter.
%
While the Atlas software initially consisted of a single executable
program \.{atlas} that essentially allows running all of the functionality
provided by the software library separately, a second executable program
initially called \.{realex} has been developed since the year~$2006$, which
program also gives access to the functionality provided by the software
library, but in a manner that provides the user with means to capture the
results from previous computations and use them in subsequent computations. It
started out as an interface giving the user just assignable variables and an
expression language for calling library functions, but has grown out to a
complete programming language, which at the time of writing is still being
extended in many directions. The name \.{realex} stood for
Redesigned Expression-based Atlas of Lie groups EXecutable.

At the end of the year~$2015$ it was decided, in order to promote the use of
the programming language version as primary user interface to the Atlas
software, to give the name \.{atlas} to the interpreter program, while
renaming the old \.{atlas} program to \.{Fokko} in honour of the creator of
the Atlas software, Fokko du~Cloux. At the same time it was decided to give
the name \axis. to the programming language (independent of the Atlas
library) that had been developed as a major part of the \.{realex} program
(and which continues to develop). Though no pure \axis. interpreter is
compiled inside the Atlas software, one could fairly easily be extracted.

The sources specific to the interpreter are situated in the
current (\.{sources/interpreter/}) subdirectory, and almost all its names are
placed in the |atlas::interpreter| namespace. This part of the program has
been split into several modules, only some of which are centred around
specific \Cpp\ classes (as is the case for most modules of the Atlas library).
In general these interpreter modules are more interdependent
than those of the main Atlas library, and the subdivision is more arbitrary
and subject to change. We list here the source files for these modules, and
their dependencies notably at the level of their header files.

\point The file \.{buffer.w} defines the classes |BufferedInput| providing an
interface to input streams, and |Hash_table| storing and providing a
translation from \axis. identifiers to small integers |id_type|.

\point The file \.{parsetree.w} defines the |expr| structure representing
parsed \axis. expressions, and many related types, and node constructing
functions for use by the parser. Its header file includes \.{buffer.h}, and
declares some pointer types to types that will be defined in \.{axis-types.h}.

\point The file \.{parser.y} is the source for \.{bison}-generated the parser
file \.{parser.tab.c}. The header file \.{parser.tab.h} includes nothing, but
contains definitions depending on \.{parse\_types.h} having been included.

\point The file \.{lexer.w} defines the lexical analyser class
|Lexical_analyser| and the readline completion function |id_completion_func|.
Its header file includes \.{buffer.h}, \.{parse\_types.h} and \.{parser.tab.h}.

\point The file \.{axis-types.w} defines the main base classes for the \axis.
evaluator, |type_expr| for representing \axis. types, |value_base| for dynamic
values, |shared_context| for dynamic evaluation contexts, |expression_base|
for ``compiled'' expressions, |program_error| for exceptions, and numerous
types related to these. Its header file includes \.{parsetree.h} (which is
needed for the error classes only).

\point The file \.{global.w} defines primitive \axis. types like integers,
rationals, matrices. Also some global aspects of the interpreter like
operating the global identifier tables. Its header file
includes \.{axis-types.h} and some headers from the Atlas library.

\point The file \.{atlas-types.w} defines types primitive to \.{atlas} which
encapsulate types of the Atlas library, and their interface functions. Its
header file includes \.{types.h} and many headers from the Atlas library.

\point The file \.{axis.w} Defines the \axis. type-checker and evaluator. It
defines many classes derived from |expression_base|. Its header file
includes \.{lexer.h} and \.{global.h}.

\point The file \.{main.w}, the current module, brings everything together and
defining the main program. It has no header file.


@* Main program. This file defines a small main program that binds together
the \.{atlas} program. It is written in \Cpp, but is it mainly concerned with
interfacing to the parser that is generated by~\.{bison}, and which is
therefore written in~\Cee. However, we now compile the generated parser file
\.{parser.tab.c} using a \Cpp\ compiler, which means we can write our code
here as an ordinary \Cpp\ program, including use of namespaces.

Since depending on the readline libraries still gives difficulties on some
platforms, we arrange for the possibility of compiling this program in the
absence of that library. This means that we should refrain from any reference
to its header files, and so the corresponding \&{\#include} statements cannot
be given in the usual way, which would cause their inclusion unconditionally.
Like for the \.{Fokko} program, the compile time flag |NREADLINE|, if defined
by setting \.{-DNREADLINE} as a flag to the compiler, will prevent any
dependency on the readline library.

@q axis_version "0.9" @>
@q axis_version "0.9.1"  multiple assignments, multi-overload set command @>
@d axis_version "0.9.2" @q balancing, break, return, ++reach in while loops @>
 // numbering from 0.5 (on 27/11/2010); last change May 6, 2015

@c

@< Conditionally include some header files @>

@< Declaration of interface to the parser @>@;
namespace { @< Local static data @>@; }@;
@< Definitions of global namespace functions @>@;
namespace atlas { namespace interpreter {@< Definitions of other functions @>@;
}@;}@;
@/
@< Main program @>

@ Since the file \.{parser.y} declares \.{\%pure-parser} and \.{\%locations},
the prototype of the lexical analyser (wrapper) function |yylex| is the one
below. Curiously, the program~\.{bison} does not write this prototype to
\.{parser.tab.h}, but it does write the definitions of the types |YYSTYPE| and
|YYLTYPE| there; these require that \.{parse\_types.h} be included first. We
also declare ``{\tt\%parse-param \char`\{} |int* verbosity, expr_p*
parsed_expr@;| {\tt\char`\}}'' in~\.{parser.y}, so that the parser itself,
|yyparse|, takes an integer pointer as parameter, which it uses to signal
special requests from the user (such as verbose output but also termination or
output redirection), and a pointer to an expression, in which it writes the
result of parsing.

The definitions below used to start with |extern "C"|, but no longer do so
since the parser is now compiled as a \Cpp\ program.

@h "parse_types.h"
@h "parser.tab.h"

@< Declaration of interface to the parser @>=

int yylex (YYSTYPE *, YYLTYPE *);
@/int yyparse( atlas::interpreter::expr_p* parsed_expr, int* verbosity );

@ Here are the wrapper function for the lexical analyser and the error
reporting function, which are necessary because the parser cannot directly
call a class method. The prototypes are imposed, in particular the second and
third arguments to |yyerror| are those passed to |yyparse|, even though they
are not used in |yyerror|. In |yyerror| we close any open include files, as
continuing to execute their commands is undesirable.

@< Definitions of global namespace functions @>=

int yylex(YYSTYPE *valp, YYLTYPE *locp)
{@; return atlas::interpreter::lex->get_token(valp,locp); }
@)

void yyerror (YYLTYPE* locp, atlas::interpreter::expr_p* ,int* ,char const *s)
{ atlas::interpreter::main_input_buffer->show_range@|
  (std::cerr,
   locp->first_line, locp->first_column,
   locp->last_line,  locp->last_column);
  std::cerr << s << std::endl;
  atlas::interpreter::main_input_buffer->close_includes();
@/atlas::interpreter::clean=false;
}

@ We have a user interrupt handler that simple raises |interrupt_flag| and
returns. It must be declared |extern "C"|.

@< Definitions of global namespace functions @>=

extern "C" void sigint_handler (int)
@+{@; atlas::interpreter::interrupt_flag=1; }

@ Here are some header files which need to be included for this main program.
As we discussed above, the inclusion of header files for the readline
libraries is made dependent on the flag |NREADLINE|. In case the flag is set,
we define the few symbols used from the readline library as macros, so that
the code using them can be compiled without needing additional \&{\#ifdef}
lines. It turns out that |getline| and |add_history| are not used in calls,
but rather passed as function pointers to the |BufferedInput| constructor, so
the appropriate expansion for these macros is the null pointer.

Similarly we make the inclusion of \.{unistd.h} dependent on the |NOT_UNIX|
not being defined.

@h <iostream>
@h <fstream>

@h "buffer.h"
@h "lexer.h"
@h "../version.h"

@< Conditionally include some header files @>=
#ifdef NREADLINE
#define readline nullptr
#define add_history nullptr
#define clear_history()
#else
#include <readline/readline.h>
#include <readline/history.h>
#endif
#ifdef NOT_UNIX
#define isatty() true // if we cannot find out, guess it is ``yes''
#define chdir()
#else
#include <unistd.h>
#endif


@ Here is an array that declares the keywords that the lexical scanner is to
recognise, terminated by a null pointer. Currently the lexical analyser adds
the offset of the keyword in this list to |QUIT|, so the recognition depends
on the fact that |"quit"| is the first keyword, and that they are listed below
in the same order as in the \.{\%token} declarations in \.{parser.y}.

@< Local static data @>=

const char* keywords[] =
 {"quit"
 ,"set","let","in","begin","end"
 ,"if","then","else","elif","fi"
 ,"and","or","not"
 ,"next","do","dont","from","downto","while","for","od"
 ,"case","esac", "rec_fun"
 ,"true","false", "die", "break", "return"
 ,"whattype","showall","forget"
 ,nullptr};

@~After installing keywords in the lexical analyser, some more preparation is
needed. The identifiers |quiet| and |verbose| that used to be keywords are now
instead recognised only in the special commands \.{set quiet} and \.{set
verbose}; to this end the parser uses their numeric identifier codes. To
ensure that they are respectively at offsets $0,1$ of
|ana.first_identifier()|, we look up these names before any other identifiers
are introduced, notably before |initialise_evaluator| and
|initialise_builtin_types| are called to define built-in operators and
functions.

@< Prepare the lexical analyser... @>=
main_hash_table->match_literal("quiet");
main_hash_table->match_literal("verbose");
// these must be the very first identifiers
ana.set_comment_delims('{','}');

@ Here are several calls necessary to get various parts of this program off to
a good start, starting with the history and readline libraries, and setting a
comment convention. Initialising the constants in the Atlas library is no
longer necessary, as it is done automatically before |main| is called. Our own
compilation units do require explicit calling of their initialisation
functions.

@h <csignal>
@h "atlas-types.h"

@< Initialise various parts of the program @>=
  @< Initialise the \.{readline} library interface @>
  signal(SIGINT,sigint_handler); // install handler for user interrupt
  initialise_evaluator();
  initialise_builtin_types();


@ Our main program constructs unique instances
for various classes of the interpreter, and sets pointers to them so that
various compilation units can access them. Then it processes command line
arguments and does some more declarations that can depend on them.
Finally it executes the main command
loop, from which normally only the \.{quit} command will make it exit.

@< Main program @>=

int main(int argc, char** argv)
{ using namespace atlas::interpreter;
@)

@/Hash_table hash; main_hash_table= &hash;
@/Id_table main_table; @+ global_id_table=&main_table;
@/overload_table main_overload_table;
 @+ global_overload_table=&main_overload_table;
@)
  bool use_readline=true;
  @< Other local variables of |main| @>
  @< Handle command line arguments @>
@/BufferedInput input_buffer(isatty(STDIN_FILENO) ? "atlas> " : nullptr
                            ,use_readline ? readline : nullptr
			    ,use_readline ? add_history : nullptr);
  main_input_buffer= &input_buffer;
@/Lexical_analyser ana(input_buffer,hash,keywords,prim_names); lex=&ana;
  @< Prepare the lexical analyser |ana| after construction and before use @>
@/@< Initialise various parts of the program @>
  @< Enter system variables into |global_id_table| @>
@)
  @< Silently read in the files from |prelude_filenames| @>
  std::cout << "This is 'atlas' (version " @|
       << atlas::version::VERSION @| << ", axis language version " @|
       axis_version "),\n" @| << atlas::version::NAME << @|
       " interpreter,\ncompiled on " @|  << atlas::version::COMPILEDATE
       << ".   http://www.liegroups.org/\n";
@)
  @< Enter the main command loop @>
  signal(SIGINT,SIG_DFL); // reinstall default signal handler
  @< Finalise  various parts of the program @>
  std::cout << "Bye.\n";
  return clean ? EXIT_SUCCESS : EXIT_FAILURE ;
}

@ When reading command line arguments, some options may specify a search path,
any non-option arguments will considered to be file names (forming the
``prelude''). We shall gather both in lists of \Cee-strings.

@< Other local variables of |main| @>=
std::vector<const char*> paths,prelude_filenames;
paths.reserve(argc-1); prelude_filenames.reserve(argc-1);

@ The strings in |paths| will initialise a ``system variable'' created below
(a variable the user can assign to, and which is inspected whenever files are
opened). We shall use a static variable that gives access to its value. It
continues to do so even if the user should manage to forget or hide the user
variable introduced below to hold it. We shall also keep a pointer to the
initial working directory name and to the first path specified, in order to
be able to perform the somewhat messy filename completion manoeuvres below.

@h "global.h" // defines |shared_share|
@< Local static data @>=
static atlas::interpreter::shared_share input_path_pointer,prelude_log_pointer;
#ifndef NOT_UNIX
enum { cwd_size = 0x1000 };
char cwd_buffer[cwd_size]; // since only |getcwd| is standard, fix a buffer
const char* working_directory_name=nullptr;
const char* first_path = nullptr;
#endif

@ Apart from the \.{--no-readline} option to switch off the readline functions
(which might be useful when input comes from a file), the program accepts
options that set the search path for scripts, and a number of scripts that
form the ``prelude''. The readline option must be read early to influence the
constructor of the lexical analyser, but the other options are just stored
away here for later processing.

@h <cstring>

@< Handle command line arguments @>=
while (*++argv!=nullptr)
{ static const char* const path_opt = "--path=";
  static const size_t pol = std::strlen(path_opt);
  std::string arg(*argv);
  if (arg=="--no-readline")
    {@; use_readline = false; continue; }
  if (arg.substr(0,pol)==path_opt)
     paths.push_back(&(*argv)[pol]);
  else prelude_filenames.push_back(*argv);
}

@ Here we create the system variables called |input_path| and |prelude_log|;
both are lists of strings, the latter a constant one. This is also a
convenient time to test and record the first path specified (which has to wait
until other initialisations have been done), but we defer the details.

@< Enter system variables into |global_id_table| @>=
{ @< Record the first specified path, if it was a directory @>
  own_value input_path = std::make_shared<row_value>(paths.size());
  id_type ip_id = main_hash_table->match_literal("input_path");
  auto oit = force<row_value>(input_path.get())->val.begin();
  for (auto it=paths.begin(); it!=paths.end(); ++it)
    *oit++ = std::make_shared<string_value>(std::string(*it)+'/');
@/global_id_table->add@|(ip_id
                       ,input_path, mk_type_expr("[string]"), false);
  input_path_pointer = global_id_table->address_of(ip_id);
@)
  own_value prelude_log = std::make_shared<row_value>(0); // start out empty
  id_type pl_id = main_hash_table->match_literal("prelude_log");
  auto& logs = force<row_value>(input_path.get())->val;
  logs.reserve(prelude_filenames.size());
@/global_id_table->add@|(pl_id
                       ,prelude_log, mk_type_expr("[string]"), true);
  prelude_log_pointer = global_id_table->address_of(pl_id);
}

@ We can now define the functions that are used in \.{buffer.w} to access the
input path.

@< Definitions of other functions @>=

unsigned int input_path_size()
{ const row_value* path = force<row_value>(input_path_pointer->get());
@/return path->val.size();
}
const std::string& input_path_component(unsigned int i)
{ const row_value* path = force<row_value>(input_path_pointer->get());
  const string_value* dir = force<string_value>(path->val[i].get());
@/return dir->val;
}

@ At start up we recorded the working directory, so that we can now easily
test whether a first specified path is actually a directory, by just trying to
change directory to it (which is all we shall do with the value later); if
this succeeds we immediately set it back, but we record the path in this case
for later use.

@< Record the first specified path, if it was a directory @>=
#ifndef NOT_UNIX
if (paths.size()>0)
{ if (working_directory_name!=nullptr and chdir(paths[0])==0)
   // see if we can \.{cd} there
  { if (chdir(working_directory_name)!=0) // but if so switch back
      return EXIT_FAILURE; // in the unlikely case that we cannot, just give up
    first_path=paths[0]; // record that this is a valid directory
  }
  else
    std::cerr << "Warning: specified path '" << paths[0] @|
              << "' is not a directory." << std::endl;
}
#endif

@ The command loop maintains two global variables that were defined
in \.{axis.w}, namely |last_type| and |last_value|; these start off in a
neutral state. In a loop we call the parser until it sets |verbosity<0|, which
is done upon seeing the \.{quit} command. We call the |reset| method of the
lexical scanner before calling the parser, which will discard any input that
is left by a possible previous erroneous input. This also already fetches a
new line of input, or abandons the program in case none can be obtained.

@< Enter the main command loop @>=
last_value = shared_value (new tuple_value(0));
last_type = void_type.copy();
 // |last_type| is a |type_ptr| defined in \.{axis.w}
while (ana.reset()) // get a fresh line for lexical analyser, or quit
{ @< Undo temporary trickery aimed at |readline| filename completion @>
  expr_p parse_tree;
  int old_verbosity=verbosity;
  std::ofstream redirect; // if opened, this will be closed at end of loop
  if (yyparse(&parse_tree,&verbosity)!=0)
     // syntax error (inputs are closed) or non-expression
    continue;
  if (verbosity!=0) // then some special action was requested
  { if (verbosity<0)
      break; // \.{quit} command
    if (verbosity==2 or verbosity==3)
      // indicates output redirection was requested
    { auto write_mode =
        verbosity==2 ? std::ios_base::trunc : std::ios_base::@;app;
      verbosity=old_verbosity; // ensure that verbosity change remains temporary
      @< Open |redirect| to specified file, and if successful make
      |output_stream| point to it; otherwise |continue| @>
    }
    if (verbosity==1) //
      std::cout << "Expression before type analysis: " << *parse_tree
                << std::endl;
  }
  interrupt_flag=0; // clear interrupt before starting evaluation
  @< Analyse types and then evaluate and print, or catch runtime or other
     errors @>
  output_stream= &std::cout; // reset output stream if it was changed
}

@ If a type error is detected by |analyse_types|, then it will have signalled
it and thrown a |program_error|; if not, then evaluation may instead produce a
|runtime_error|. Therefore the manipulation below of |type_OK| to see whether
we passed the analysis phase successfully should be redundant. If the result
of evaluation is an empty tuple, we suppress printing of the uninteresting
value.

@h <stdexcept>
@h "axis.h"

@< Analyse types and then evaluate and print... @>=
{ bool type_OK=false;
  try
  { expression_ptr e;
    type_expr found_type=analyse_types(*parse_tree,e);
    type_OK=true;
    if (verbosity>0)
      std::cout << "Type found: " << found_type << std::endl @|
	        << "Converted expression: " << *e << std::endl;
    e->evaluate(expression_base::single_value);
@)  // now that evaluation did not |throw|, we can record the predicted type
    last_type = std::move(found_type);
    last_value=pop_value();
    if (last_type!=void_type)
      *output_stream << "Value: " << *last_value << std::endl;
    destroy_expr(parse_tree);
  }
  @< Various |catch| phrases for the main loop @>
}

@*1 Reading in the prelude files.
%
Before entering that main loop, we do a simplified version of the command
loop to read in the prelude files. We do not accept anything that changes
|verbosity| (like trying to call \.{quit}) during the prelude, do not maintain
a last value, and in case of type or runtime errors complain more succinctly
than in the main loop. All non-error output goes to a component of the
|prelude_log| user variable. This is mainly due to the assignment
|output_stream = &log_stream;| and the fact that functions in \.{global.w} use
this pointer; the \.{set} commands that form the essence of prelude files
return a nonzero value from |yyparse|, so in practice most of the code below
is rarely executed.

Errors will break from the inner loop simply by popping the open include
file(s) from |main_input_buffer|. We do not attempt to break from the outer
loop upon an error, as this circumstance is rather hard to detect: any error
has already been reported on |std::cerr|, and leaves a situation not very
different from successfully completing reading the file.

@< Silently read in the files from |prelude_filenames| @>=
for (auto it=prelude_filenames.begin(); it!=prelude_filenames.end(); ++it )
{ std::ostringstream log_stream; output_stream = &log_stream;
  main_input_buffer->push_file(*it,true);
    // set up to read |fname|, unless already done
  while (main_input_buffer->include_depth()>0) // go on until file ends
  { if (not ana.reset())
    { std::cerr << "Internal error, getline fails reading " << *it
                  << std::endl;
      return EXIT_FAILURE;
    }
    expr_p parse_tree;
    if (yyparse(&parse_tree,&verbosity)!=0)
      continue; // if a syntax error was signalled input has been closed
    if (verbosity!=0)
    { std::cerr << "Cannot "
                << (verbosity<0 ? "quit" :
                    verbosity==1 ? "set verbose" : "redirect output")
                         << " during prelude.\n";
      verbosity=0; main_input_buffer->close_includes();
    }
    else
    { try
      { expression_ptr e; type_expr found_type=analyse_types(*parse_tree,e);
        e->evaluate(expression_base::single_value);
        if (found_type!=void_type)
          log_stream << "Value: " << *pop_value() << '\n';
        else
          pop_value(); // don't forget to cast away that void value
      }
      catch (std::exception& err)
      { std::cerr << err.what() << std::endl;
        reset_evaluator(); main_input_buffer->close_includes();
      @/clean=false;
      }
    }
    destroy_expr(parse_tree);
  }
  value pl_val = uniquify(*prelude_log_pointer);
  row_value* logs = force<row_value>(pl_val);
  logs->val.emplace_back(std::make_shared<string_value>(log_stream.str()));
  output_stream = &std::cout;
}

@ The |std::ofstream| object was already created earlier in the main loop,
but it will only be opened if we come here. If this fails then we report it
directly and |continue| to the next iteration of the main loop, which is more
practical at this point than throwing and catching an error.

@< Open |redirect| to specified file... @>=
{ redirect.open(ana.scanned_file_name() ,std::ios_base::out | write_mode);
  if (redirect.is_open())
    output_stream = &redirect;
  else
  {@; std::cerr << "Failed to open " << ana.scanned_file_name() << std::endl;
    continue;
  }
}

@ We distinguish runtime errors (which are normal) from internal errors (which
should not happen), and also |catch| and report any other error derived from
|error_base| that could be thrown. After any of these errors we close all
open auxiliary input files; reporting where we were reading is done by the
method |close_includes| defined in \.{buffer.w}.

@< Various |catch| phrases for the main loop @>=
catch (const runtime_error& err)
{ assert(type_OK);
  std::cerr << "Runtime error:\n  " << err.what() << "\nEvaluation aborted."
            << std::endl;
@/clean=false;
  reset_evaluator(); main_input_buffer->close_includes();
}
catch (const program_error& err)
{ assert(not type_OK);
  std::cerr << err.what() << std::endl;
@/clean=false;
  reset_evaluator(); main_input_buffer->close_includes();
}
catch (const logic_error& err)
{ std::cerr << "Internal error: " << err.what() << std::endl;
@/clean=false;
  reset_evaluator(); main_input_buffer->close_includes();
}
catch (const std::exception& err)
{ std::cerr << err.what() << "\nEvaluation aborted.\n";
@/clean=false;
  reset_evaluator(); main_input_buffer->close_includes();
}

@*1 Identifier completion.
%
We define a completion function |id_completion_func| that will be used by the
\.{readline} library. The completion function will be called from the
|readline| function, which is itself probably called by
|BufferedInput::getline|, if the user asks for it by hitting the ``tab'' key.
The function prototype is dictated by the \.{readline} library.

The completion function will find identifiers matching the given prefix that
are currently either known in the overload or global identifier tables, or
have such a low code that they are keywords or predefined types. All such
identifiers are known in |main_hash_table|, so that is what our primary loop
is over, but upon finding a partial match we test the stated conditions. Those
additional tests prevent names of local identifiers of functions loaded and
accidentally typed erroneous identifiers to ``pollute the completion space''.

The completion function has some strange characteristics that are dictated by
the \.{readline} library. It must perform a loop that is actually
started \emph{outside} the function body, so the only possible way to keep
track of the loop state is using static variables. The |state| parameter
signals (by being~|0|) when a new loop starts, so in that case it is time to
(re-)initialise the static variables. A part of the loop can be picked up
inside the function body, namely the search for the next match to the supplied
prefix |text|. A tricky point is that once a partial match is found, we must
increment the static iterator before returning, since that |return| jumps out
of our local loop. For this reason we increment |i| right away while picking
its identifier from the hash table. Otherwise there are no other
complications, except having to produce a string allocated by |malloc|; this
used to be done by calling |strdup|, but since that function appears not to be
part of standard \Cpp\ at all, we do the duplication explicitly. If our local
loop terminates normally there are no more matches and we return |nullptr| to
indicate that circumstance.

@h <cstdlib>
@< Definitions of global namespace functions @>=
extern "C"
char* id_completion_func(const char* text, int state)
{ using namespace atlas::interpreter;
  static size_t l; static id_type i,n;
  if (state==0)
  { i=0; n=main_hash_table->nr_entries();
    l=std::strlen(text); // fix length for during search
  }
  while (i<n)
    // |i| is initialised above when |state==0|, and incremented below
  { id_type id=i++; // take next identifier code and increment
    const char* s=main_hash_table->name_of(id);
      // get stored identifier and increment loop
    if (std::strncmp(text,s,l) == 0 // is |text| a prefix of |s|?
      and (id<lex->first_identifier() @| or
          not global_overload_table->variants(id).empty() @| or
          global_id_table->present(id)))
    { char* result=static_cast<char*>(std::malloc(std::strlen(s)+1));
      if (result==nullptr)
        @[throw std::bad_alloc()@];
      return std::strcpy(result,s);
    }
  }
  return nullptr; /* if loop terminates, report failure */
}

@ The readline library needs to know where to break the input into words that
may be completed. The value |lexical_break_chars| reflects what our lexical
analyser considers separating characters: the list contains those
non-alphanumeric characters that are valid input characters, as implicitly
defined by the |get_token| method of |Lexical_analyser|. This choice will
affect file name completion as well, because |readline| will not include any
of those characters in the partial match it wants to complete; even if we
decide to ask for file name completion after all, it will not extend its
partial match to the left with characters like `\.-' that are not special in
file names, though it will in that case consider any characters it finds in
actual file names when it performs the completion (to the right). Therefore
typing more can be worse: adding characters that are actually part of the
intended file name can in fact thwart file name completion. This is awkward,
but the only way to avoid this would be to build our knowledge about when
completion should be for file names rather than for identifiers \emph{into}
the function bound to the completion key (tab), rather than leave it at its
default binding |rl_complete| as we do. Once |readline| calls one of our
functions (in our case |do_completion| defined below), it is too late to
adjust the |text| it selected to be completed.

@< Local static data @>=
char lexical_break_chars[] = " \t\n=<>+-*/\\%,;:()[]{}#!?$@@\"|~";

@ The above function |id_completion_func| will not be plugged directly into
the readline completion mechanism, but instead we provide an alternative
function for generating matches, which may pass our function to
|rl_completion_matches| when it deems the situation appropriate, or else
returns |nullptr| to indicate that the default function, completing on file
names, should be used instead.

This function is called with |text| pointing to the partial string to be
completed, which is at positions from |start| up to |end| in the buffer
|rl_line_buffer|. We decide by looking at the part before |start| whether we
should interpret that part as a partial file name, and if so in which
directory we should look: if only characters \.<, \.> and spaces precede the
string, then we ask for file name completion, and since input will most likely
come from the directory in |first_path|, if any was specified, we temporarily
change directory there if only and at least one characters \.< was present. In
the more common case where file name completion is not called for, we call the
readline function |rl_completion_matches| with |text| to get us a list of
possible completions, which it does by calling our |id_completion_func|
repeatedly, and we pass the pointer to the list of completions back to our
caller (which is probably some |readline| action function).

@< Definitions of global namespace functions @>=
#ifndef NREADLINE
extern "C" char** do_completion(const char* text, int start, int end)
{
  if (atlas::interpreter::lex->is_initial() and start>0)
   // could be a file name
  { int i; char c; bool need_file=false; bool in=true;
    for (i=0; i<start; ++i)
      if (std::isspace(c=rl_line_buffer[i]))
        continue; // ignore space characters
      else if (c=='<')
        need_file=true;
      else if (c=='>')
        need_file=true,in=false;
      else
        break; // any other preceding characters make it not a file name

    if (working_directory_name!=nullptr and need_file and i==start)
       // the sub-string is preceded by one or more copies of \.<, \.>
    { static_cast<void> // pretend we inspect the result, in fact ignore failure
      (chdir(in and first_path!=nullptr ? first_path : working_directory_name));
         // temporary \.{cd}
      return nullptr; // and signal that file name completion should be used
    }
  }
  rl_attempted_completion_over = true;
    // don't try file name completion if we get here
  return rl_completion_matches(text,id_completion_func);
    // but make |readline| call our function
}
#endif

@ In order to trick |readline| into doing file name completion for a
non-current directory, we apply a truly horrible trick: we save the true
current working directory name in order to temporarily change directory and
afterwards change back.

At start up, we capture the current directory name. Since the only
sufficiently standard function |getcwd| to get the directory name will only
work with a fixed size buffer, our program has to cope with the remote
possibility that the call below sets |working_directory_name =nullptr| for
lack of space in the buffer.

@< Initialise the \.{readline} library interface @>=
#ifndef NOT_UNIX
working_directory_name = getcwd(cwd_buffer,cwd_size);
#endif

@ After the lexical scanner is reset (and therefore |readline| certainly has
returned) and before parsing starts, we ensure any temporary change to the
current working directory are undone. Thus none of the true I/O actions will
notice the temporary change of working directory.

@< Undo temporary trickery aimed at |readline| filename completion @>=
#ifndef NOT_UNIX
if (working_directory_name!=nullptr)
  static_cast<void> // pretend we inspect the result, in fact ignore failure
      (chdir(working_directory_name)); // reset to true working directory
#endif


@ The \.{readline} library requires setting the variable
|rl_completer_word_break_characters|, and by setting the function pointer
|rl_attempted_completion_function|, we hook our functionality into the
|readline| machinery.

@< Initialise the \.{readline} library interface @>=
#ifndef NREADLINE
  rl_readline_name="axis"; // for users who want to customise |readline|
  using_history();
  rl_completer_word_break_characters = lexical_break_chars;
  rl_attempted_completion_function = do_completion; // set up input completion

#endif

@ Just to be proper, we clear out our history and |malloc|ed variable when
terminating the program.

@< Finalise  various parts of the program @>=
#ifndef NREADLINE
  clear_history();
   // clean up (presumably disposes of the lines stored in history)
#endif

@* Index.

% Local IspellDict: british
