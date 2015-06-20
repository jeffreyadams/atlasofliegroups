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

@* Introduction.
This file describes a lexical scanner that provides a layer situated between
the class |BufferedInput| that provides lines of input, and the parser that is
to receive a sequence of tokens. The functionality provided is the recognition
of tokens, in particular keywords and constants.

@ As usual the external interface is written to the header file associated to
this file. In this module all include files are needed in the header file.

@( lexer.h @>=
#ifndef LEXER_H
#define LEXER_H
@p@;
namespace atlas { namespace interpreter {

@< Class declarations @>@;
@< Declarations of static variables @>@;

}@; }@;
#endif

@ Here we include only our own header file \.{lexer.h}; all other include
directives will go to that file, as indicated in the previous section (this
might be a bit wasteful, but only system header files are concerned).

@c
#include "lexer.h"
namespace atlas
{ namespace interpreter
   {
@< Definitions of class members @>@;
@< Definitions of static variables @>@;
   }@;
}@;


@ Since there is one central hash table, and other parts of the program must
have access to it, for instance to look up names of identifiers, we define a
static variable with a pointer to this hash table.

@< Declarations of static variables @>=
extern Hash_table* main_hash_table;

@~We initialise this variable to the null pointer; the main program will make
it point to the main hash table once it is allocated.

@< Definitions of static variables @>=
Hash_table* main_hash_table=nullptr;

@* The lexical analyser class.
%
We now come to the lexical analyser proper. Although only one lexical analyser
is envisaged, we shall define a class for it. The module \.{buffer} is used
for the class |BufferedInput| as well as for the type |id_type| defined within
that class. The file \.{parser.tab.h} contains definitions of |YYSTYPE| and
|YYLTYPE| defined by the parser and used in the code below, but on its turn it
uses (for other purposes) types defined in \.{parse\_types.h}, which therefore
has to be loaded before it (we would like to have put an \&{\#include} of the
file \.{parsetree.h} into \.{parser.tab.h} so that it need no be mentioned
here, but we do not know if or how this could be arranged).

@h "buffer.h"
@h "parse_types.h"
@h "parser.tab.h"

@< Class declarations @>=
class Lexical_analyser
{ enum states @+ { initial, normal, ended };
@)BufferedInput& input;
  Hash_table& id_table;
  id_type keyword_limit; // first non-keyword identifier
  id_type type_limit; // first non-type identifier
  int nesting; // number of pending opening symbols
  char prevent_termination; // either |'\0'| or character requiring more input
  int comment_start, comment_end; // characters that start/end a comment
  states state; // to trigger special behaviour
  std::string file_name; // stores a file name for I/O redirection
public:
  Lexical_analyser
    (BufferedInput&, Hash_table&, const char**, const char** type_names);
  int get_token(YYSTYPE *valp, YYLTYPE* locp);
  bool reset(); // get clean slate, return |false| if |getline()| fails
  void set_comment_delims (char c, char d)
          {@; assert(c!='\0' and d!='\0'); comment_start=c; comment_end=d; }
  const char* scanned_file_name() const @+{@; return file_name.c_str(); }
  id_type first_identifier() const @+{@; return type_limit; }
private:
  void skip_space() const;
  bool becomes_follows();
  std::string scan_quoted_string() const;
};

@ Since there is one lexical analyser object, and other parts of the program
must have access to it, we define a static variable with a pointer to it.

@< Declarations of static variables @>=
extern Lexical_analyser* lex;

@~We initialise this variable to the null pointer; the main program will make
it point to the main hash table once it is allocated.

@< Definitions of static variables @>=
Lexical_analyser* lex=nullptr;


@ Here is the constructor for the lexical analyser, which assumes that a
buffered input object and an empty hash table object have been previously
constructed, and are passed by reference. Currently it is called with a list
of keyword strings, and a list of predefined type names; both will be
installed into the hash table and determine the values of |keyword_limit| and
|type_limit|. There will probably be a need to further parametrise the lexical
analyser, if we do not want to hard-code all lexical details into it (there is
nothing wrong with that as long as there is only one object of this class, but
the class concept invites us to envision some more flexible use). One such
parametrisation is via the |comment_start| and |command_end| characters, that
if set using |set_comment_delims| will automatically skip text enclosed
between them (they may or may not be equal). In the unset state they are set
to an integer that cannot match any |char| value; we would have like to use
|EOF| defined in \.{ctype.h} here, but it is only guaranteed to be
non-|(unsigned char)|, and since using the type |(unsigned char*)| is
unwieldy, we use another value.

@< Definitions of class members @>=
Lexical_analyser::Lexical_analyser
  (BufferedInput& source, Hash_table& hash,
   const char** keywords, const char** type_names)
: input(source),id_table(hash),nesting(0)
 ,prevent_termination('\0'),state(initial)
{ @< Install |keywords| and |type_names| into |id_table| @>
  comment_start=comment_end=256; // a non-|char| value
}

@ Keywords are identified by sequence number in the order by which they are
entered into |id_table|; therefore the order in the list |keyword| passed to
the constructor should match the numeric \.{\%token} values defined in
\.{parser.y}. The actual code transmitted for keywords will be obtained by
adding the constant |QUIT| to the value returned from the hash table look-up.
Type names are next in |id_table|, but they all will return the token |TYPE|,
while recording which names was entered in the semantic value.

@< Install |keywords| and |type_names| into |id_table| @>=
{ for (size_t i=0; keywords[i]!=0; ++i)
    id_table.match_literal(keywords[i]);
  keyword_limit=id_table.nr_entries();
  for (size_t i=0; type_names[i]!=0; ++i)
    id_table.match_literal(type_names[i]);
  type_limit=id_table.nr_entries();
}

@ The member function |reset| can be called to reset the lexical analyser,
discarding any remaining input on the current line and clearing the |nesting|
level. It is safe to call |reset| after a successfully executed command, which
will fetch a new line without getting any tokens, but |reset| should not be
called twice in succession, as this will discard the newly fetched line. The
result returned tells whether a fresh line was successfully obtained.

@< Definitions of class members @>=
bool Lexical_analyser::reset()
{@; nesting=0; state=initial; input.reset(); return input.getline(); }

@ Skipping spaces is a rather common activity during scanning; it is performed
by |skip_space|. When it is called, the first potential space character has
been shifted in, and when it returns a non-space character is left in
shifted-in position. Skipping comments is also performed by |skip_space|. In
spite of its name, this function does not unconditionally skip all spaces, so
as to be able to catch some runaway constructions. Also the subtle decision
mechanism to determine which newlines terminate a command kicks in here, in
the |std::isspace(c)| branch below, although the actual control of this
decision is distributed in the various pieces of code that maintain of the
fields |prevent_termination| and |nesting|.

In case end of input occurs one obtains |shift()=='\0'|, and for the end of an
included file one obtains |shift()=='\f'|, a form-feed. If the former happens
in the main loop of |skip_space|, then we break out of it like we do for any
non-space character (the predicate |is_space| does not hold for the null
character), and although the following |input.unshift()| does nothing, the
value |'\0'| should reappear at the next call of |shift|.

@h<cctype>
@h<iostream>

@s level x

@< Definitions of class members @>=
void Lexical_analyser::skip_space() const
{ if (prevent_termination!='\0') input.push_prompt(prevent_termination);
  do
  { char c=input.shift();
    if (std::isspace(c))
     // ignore unless file ends, or a newline where a command could end
  @/{@; if (c=='\f' or
            c=='\n' and prevent_termination=='\0' and nesting==0)
        break;
    }
    else if (c==comment_start) @< Skip comment, possibly nested @>
    else break;
      // non-space and non-comment character (possibly |'\0'|): terminate
  } while(true);
  if (prevent_termination!='\0') input.pop_prompt();
  input.unshift(); // prepare to re-read character that ended space
}

@ If |comment_start| is set and |comment_end!='\n'|, then newlines will be
skipped inside comments, but if |comment_end=='\n'| then a closing newline
will be considered for the possibility of ending input. In the former case we
change the input prompt by pushing the comment character to warn the user that
no action has been performed. In case |prevent_termination| is set, that
character is also pushed into the prompt for the duration of skipping spaces.

If end of input or a form-feed occurs during a comment, then the comment is
terminated (with an error message) due to the explicit test, and the character
that provoked this will be reconsidered so the end of file will not be masked
by the comment it occurred inside.

@< Skip comment, possibly nested @>=
{ int line, column;
  input.locate(input.point(),line,column); // maybe needed for reporting
  if (comment_end=='\n' or comment_end==comment_start)
     // non nestable comment
  { do c=input.shift();
 @/ while (c!=comment_end and c!='\0' and c!='\f');
    if (c!=comment_end)
    { std::cerr << "Comment that started on line " << line
                << ", column " << column @| << " is never closed.\n";
    @/input.unshift(); // reconsider end of file
    }
    else if (comment_end=='\n')
      input.unshift(); // reconsider newline character
  }
  else // skip comment, allowing for nested comment groups
  { int level=1;
    input.push_prompt(comment_start);
    do
    {
      do c=input.shift(); @/ @q forced break fixes broken cweavex @>
      while (c!= comment_start and c!=comment_end and c!='\0' and c!='\f');
      if (c=='\0' or c=='\f')
      { std::cerr << "Comment that started on line " << line
                  << ", column " << column @| << " is never closed.\n";
      @/ input.unshift(); break;
      // force out of comment loop, reconsider character
      }
      if (c==comment_start)
      {@; ++level;
        input.push_prompt(comment_start);
      }
      else
      {@; --level;
        input.pop_prompt();
      } // |c==comment_end|
    }
    while (level>0);
  }
} // after skipping the comment, skipping spaces (and comments) continues

@ Another auxiliary method is used for scanning a quoted string; it should be
called when an initial double-quote character has been recognised, and after
scanning the string returns it (with quotes and escapes removed) as a
|std::string| value. We currently use a simple model for strings. They should
be contained in a single line, and the only escapes used are the doubling of
double-quote characters. We copy the string while reducing doubled
double-quote characters to single ones.

@h <string>

@< Definitions of class members @>=
std::string Lexical_analyser::scan_quoted_string() const
{ const char* start=input.point(); bool broken=false;
  std::string result; char c;
  while (true)
  { for (c=input.shift(); c!='"' and c!='\n' and c!='\0'; c=input.shift())
      result.push_back(c);
    if (c!='"')
      {@; broken=true; break; }
    else if ((c=input.shift())!='"')
      break; // normal ending of string
    result.push_back(c); // doubled quote; insert one copy and continue
  }
  input.unshift();
  if (broken)
  { const char* end=input.point();
    int l0,c0,l1,c1;
    input.locate(start,l0,c0); input.locate(end,l1,c1);
    input.show_range(std::cerr,l0,c0,l1,c1);
    std::cerr << "Closing string denotation.\n";
  }
  return result;
}

@*1 The main scanning routine.
%
Now we come to the function |get_token| that will actually be called by the
parser and return a token to it. The token is an integral value that is either
a character code in the case of single character tokens, or a token value
defined in \.{parser.tab.h}. The null token value has the special meaning of
signalling the end of the input to be recognised by the parser, which it wants
to see before returning successfully. Since our parser is written to recognise
individual commands rather than the full input source, we take measures to
send a null token after each command, even though the input buffer does
not signal end of input.

If it should happen that the input buffer does signal end of input by making
|input.shift()=='\0'|, then this will be transmitted to the parser as the end
of a command by the code below, by the general mechanism that characters that
are not specifically recognised are transmitted with a token value equal to
the character value. But this is really a marginal case, since the input
buffer guarantees that end of input cannot occur in the middle of a line;
normally the previous line will therefore be completely processed, and the
end-of-input condition is signalled by a failure of the call of the |reset|
method of the lexical analyser rather than by |input.shift()=='\0'| occurring
in its |getline| method. The only way the latter can happen is if the
preceding newline character was ignored by |skip_space|, due to
|prevent_termination| or |nesting|.

The way in which we arrange to signal the end of a command from |get_token| is
by sending \emph{two} successive tokens, a |'\n'| followed by a null token.
This strange setup is a work-around of a nasty property of the parsers
generated by \.{bison}, namely that on one hand the required terminal null
token is not representable in syntax rules (it should follow input matching
the target symbol, but is not part of its production), while on the other hand
the parser will perform any unambiguous reduction of the input it has seen,
even if the lookahead token is wrong for that reduction (so a syntax error
will be reported \emph{after} the reduction is performed). This means that if
we would define the syntax of commands normally, and if such a command would
occur correctly in the input but followed by extra tokens, then the reduction
recognising and executing the command would be performed \emph{before}
flagging a syntax error for the spurious tokens, which is undesirable.
Therefore we include into the production for any command the |'\n'| token that
is guaranteed by our lexical analyser to precede the terminating null token,
thus ensuring that the command will only be executed if it constitutes the
complete user input.

(One may imagine that premature reduction could cause problems in other
positions as well, but this is by far the most irritating case, and one where
no grammar rewriting can help.)

The implementation of this two-token termination reporting is by setting
|state=ended| when the token |'\n'| is returned, and testing this condition as
the first action in |get_token|, sending the following null token if it holds.
The |state| variable is also used to allow the scanner to behave differently
while scanning the very first token of a command, when |state==initial| will
hold; once a token is scanned, |state| is set to |normal| unless scanning had
set it to |ended|.

Besides returning a token code, each token defines a value for
|prevent_termination|. Since the previous value is only used in |skip_space|,
we set this variable to its most common value |'\0'| after the call to that
function, so that only for those tokens that cannot end a command we have to
set that value explicitly.

The code below also deals with setting the fields of the |locp| to delimit the
token scanned.

@< Definitions of class members @>=
int Lexical_analyser::get_token(YYSTYPE *valp, YYLTYPE* locp)
{ if (state==ended)
    {@; state=initial; return 0; } // send end of input
  skip_space(); prevent_termination='\0';
  input.locate(input.point(),locp->first_line,locp->first_column);
  int code; char c=input.shift();
  if (std::isalpha(c) or c=='_')
    @< Scan an identifier or a keyword @>
  else if (std::isdigit(c))
    @< Scan a number @>
  else
    @< Scan a token starting with a non alpha-numeric character @>
  input.locate(input.point(),locp->last_line,locp->last_column);
@/if (state==initial)
   state=normal;
  return code;
}

@ Everything that looks like an identifier is either that, or a keyword, or a
type name. In any case we start with looking it up in the |id_table|, and then
the numeric value of the code returned, which is determined by the order in
which names were first entered into |id_table|, will allow us to discriminate
the possibilities. All keywords get distinct code numbers, determined by their
offset from the first keyword; this implies that we cannot handle multiple
distinct keywords that scan as the same category because of a similar
syntactic role. However type names for primitive types all get the same
category |TYPE|, with the actual name stored in the token value
|valp->type_code|. While user defined type names (abbreviations) are
equivalent to primitive ones at the syntactic level, we give them a different
category |TYPE_ID|, because the nature of the associated token value is
different: the identifier code is stored, and finding the designated type
will require looking it up in |global_id_table|. Distinguishing type
identifiers from other identifiers at lexical analysis is also uses
|global_id_table|: type identifiers are known in that table, but have a null
pointer as value (unlike global identifiers that have been declared without
initial value: there the associated values is a shared pointer to an empty
slot (holding a null pointer) instead.

@h "global.h" // need to inspect |global_id_table|

@< Scan an identifier or a keyword @>=
{ const char* p=input.point@[()@]-1; // start of token
  do
    c=input.shift();
  while(std::isalnum(c) || c=='_');
  input.unshift();
  id_type id_code=id_table.match(p,input.point()-p);
  if (id_code>=type_limit)
  { valp->id_code=id_code;
    if (global_id_table->is_defined_type(id_code))
      code=TYPE_ID;
    else
      code=IDENT;
  }
  else if (id_code>=keyword_limit)
  {@; valp->type_code=id_code-keyword_limit; code=TYPE; }
  else // we have |id_code<keyword_limit|, so it is a keyword
  { code=QUIT+id_code;
    switch(code)
    {
      case LET: ++nesting; input.push_prompt('L'); break;
      case BEGIN:
      case IF:
      case WHILE:
      case FOR: ++nesting; input.push_prompt('G'); break;
      case IN: if (input.top_prompt()=='L')
      @/{@; --nesting; input.pop_prompt(); prevent_termination='I'; }
      break;
      case END:
      case FI:
      case OD: --nesting; input.pop_prompt(); break;
      case AND: case OR: case NOT: prevent_termination='~'; break;
      case WHATTYPE: prevent_termination='W'; break;
    }
  }
}

@ Numbers do not include a sign, and are always scanned as |unsigned long|
values.

@< Scan a number @>=
{ const char* p=input.point@[()@]-1; // start of token
  do
    c=input.shift();
  while(std::isdigit(c));
  input.unshift();
  const char* end=input.point();
  unsigned long val=*p++-'0'; @+  while (p<end) val=10*val+(*p++-'0');
  valp->val=val; code=INT;
}

@ For reasons of limited look-ahead in the parser, certain operator symbols
will be contracted with a possible following |":="| into a single token during
lexical analysis. The function |becomes_follows| will be called to perform the
necessary look-ahead, possibly include the token |":="| and report whether it
did find it. While unlikely, we do allow spaces or comments to intervene
between the operator and the becomes, so as to cover our tracks for this
lexical hack.

@< Definitions of class members @>=
bool Lexical_analyser::becomes_follows ()
{ skip_space(); // allow intervening space, comments
  const char* p = input.point();
   // take a peek into the input buffer without advancing
  if (*p==':' and *++p=='=') // then indeed |":="| follows our operator
  { input.shift(); input.shift(); // gobble up the symbol |":="|
    prevent_termination = ':'; return true;
  }
  return false;
}

@ After splitting off the alphanumeric characters, the scanner plunges into a
large switch on the value of the look-ahead character~|c|. We increase and
decrease nesting on obvious grouping characters, and for certain characters we
store them into |prevent_termination| to prevent a newline to be interpreted
as command termination. If in the |initial| state we find a |'<'| or |'>'|
character, then we prepare for input file inclusion or output redirection;
both tokens can be doubled to designate forced inclusion (if file was already
included before) respectively appending output redirection.

@< Scan a token... @>=
{ switch(c)
  {      case '"': @< Scan a string denotation @> @+
  break; case '(':
         case '{':
         case '[': ++nesting; input.push_prompt(c); code=c;
  break; case ')':
         case '}':
         case ']': --nesting; input.pop_prompt(); code=c;
  break; case ';': code=prevent_termination=c;
  break; case '<':
         case '>':
         if (state==initial)
         { code= c=='<'
           ? input.shift()=='<' ? FORCEFROMFILE:  (input.unshift(),FROMFILE)
           : input.shift()=='>' ? ADDTOFILE : (input.unshift(),TOFILE) ;
	   @< Read in |file_name| @>
           @+ break;
         }
         prevent_termination=c;
         valp->oper.priority=2;
         if (input.shift()=='=')
           valp->oper.id=id_table.match_literal(c=='<' ? "<=" : ">=");
         else
           {@; input.unshift();
               valp->oper.id=id_table.match_literal(c=='<' ? "<" : ">");
           }
         code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
  break; case ':': prevent_termination=c;
    code = input.shift()=='=' ? BECOMES  : (input.unshift(),c);
  break; case '=':
         valp->oper.id = id_table.match_literal("=");
         valp->oper.priority = 2; // in case
	 prevent_termination=c;
	 code = becomes_follows() ? OPERATOR_BECOMES : '=';
  break; case '!':
         if (input.shift()=='=')
         { valp->oper.id = id_table.match_literal("!=");
           valp->oper.priority = 2;
           prevent_termination='=';
           code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
         }
         else
           code= (input.unshift(),c);
         // currently unused; might some day be factorial operator
  break; case '~': @< Handle the |'~'| case, involving some look-ahead @>
  break; @/@< Cases of arithmetic operators, ending with |break| @>
         case '\n': state=ended; // and {\bf fall through}.
         default: code=c;
  }
}

@ The tilde character was long unused, but was introduced in subscriptions,
slices and loops, where it can be placed before certain punctuation tokens or
the keywords \&{do} or \&{od} to indicate various kinds of reversals to be
incorporated in the semantics of these expressions. While innocent looking,
these syntactic extensions cause difficulties due to the limited one-token
look-ahead the parser allows: in general it is clear what use of the tilde is
intended by looking at the next token, but the parser needs to know this to
know what is can allows for the \emph{preceding} expression, and cannot wait
until seeing the next token. The basic problem involves the $a\,{\sim}[i]$
syntax which could occur at the position of $j$ inside a $b[i:j{\sim}]$ slice;
the presence of |'['| after the tilde influences the parsing rules for~$a$. We
solve the dilemma by requiring |'~['| to appear as single unit (no space) and
have the lexer transmit the combination as a single token |TLSUB| to the
parser.

Given that in these uses |'~'| is never followed by a token that could start a
sub-formula, it seems reasonable to allow the character to, independently
from these uses, also be used as operator symbol (where it always precedes a
sub-formula; the only caveat to the user will be to separate by a space in
case the operator should be applied to a list display). We can allow that, but
unless the distinction between the two uses is made by the lexer, we again run
into parsing problems, similar to the previous ones. So to find out whether
this is a use of tilde as special symbol rather than as operator, we need to
look at the next token from within the lexer, and see it if is one of |':'|,
|']'|, |','|, or the keywords \&{do} or \&{od}. But we need to do that within
the lexer, which can be done even though it is quite ugly. We use the fact that
the input buffer allows us get a pointer~|p| to the next characters to be
read, so we just test the various possibilities manually. The keyword cases
are the hardest (we don't want to call |get_token| recursively to scan them,
especially since we cannot contribute the token yet), but fortunately the
keywords in question are not very long.

@< Handle the |'~'| case, involving some look-ahead @>=
if (input.shift()=='[') // recognise combination for parse reason
{@; code = TLSUB; ++nesting; input.push_prompt('['); }
else
 // now see if, skipping spaces and comments, next can follow |'~'|
{ input.unshift(); prevent_termination='~'; skip_space();
   // for these cases we do allow intervening space, comments
  const char* p = input.point();
   // take a peek into the input buffer without advancing
  if (*p==':' or *p==']' or *p==',' or
 @|  (*p=='d' and p[1]=='o' or *p=='o' and p[1]=='d')
   and @| not (p[2]=='_' or std::isalnum((unsigned char)p[2])))
  @/ code=c; // return |'~'|
  else
  { valp->oper.id = id_table.match_literal("~");
    valp->oper.priority = 8; // same precedence as \.{\#}
    code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
  }
}

@ Here are some cases split off to avoid the module with the |switch| on |c|
getting too long.

@< Cases of arithmetic operators... @>=
    case '+': prevent_termination=code=c;
       valp->oper.id = id_table.match_literal("+");
       valp->oper.priority = 4;
       code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
break; case '-': prevent_termination=c;
       if (input.shift()=='>')
          code= ARROW;
       else
       { input.unshift();
         valp->oper.id = id_table.match_literal("-");
         valp->oper.priority = 4;
         code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
       }
break; case '*': prevent_termination=c;
       valp->oper.id =
          id_table.match_literal(c=='*' ? "*" : c=='%' ? "%" : "/");
       valp->oper.priority = 6;
       code = becomes_follows() ? OPERATOR_BECOMES : c;
break; case '%': case '/': prevent_termination=c;
       valp->oper.id =
          id_table.match_literal(c=='%' ? "%" : "/");
       valp->oper.priority = 6;
       code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
break; case '\\':
       valp->oper.priority = 6;
       if (input.shift()=='%')
       @/{@; prevent_termination='%';
         valp->oper.id = id_table.match_literal("\\%");
       }
       else
       {@; input.unshift();
         valp->oper.id = id_table.match_literal("\\");
       }
       code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
break; case '^': prevent_termination=c;
       valp->oper.id = id_table.match_literal("^");
       valp->oper.priority = 7; // exponentiation is right-associative
       code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
break; case '#': prevent_termination=c;
       valp->oper.id = id_table.match_literal("#");
       valp->oper.priority = 8; // basic meaning is non-associative
       code = becomes_follows() ? OPERATOR_BECOMES : OPERATOR;
break;

@ We hand to the parser a pointer to a dynamic variable move-constructed from
the |std::string| value returned by |scan_quoted_string()|. For technical
reasons having a |std::string| itself as token value is not practical (it
would define a union member with non-trivial destructor, and though \Cpp\ will
allow that if certain provisions are made, the parser generator which
is \Cpp-agnostic does not make those provisions). There is only one parser
action for |STRING| tokens, and it will turn the token value into a string
denotation expression; the code below used to circumvent to problem by already
performing this conversion inside the scanner, and letting the mentioned
parser action be empty. The current solution, though a doing a bit of extra
work, is cleaner. It was chosen so that the parser can pass a complete
location description when the string denotation expression is built; when the
current code is executed the start of the token is recorded in |locp|, but not
yet its end. Also this is allows for possible future addition of parser rules
involving |STRING|.

@< Scan a string... @>=
{@; valp->str = new std::string(scan_quoted_string()); code=STRING; }

@ Since file names have a different lexical structure than identifiers, they
are treated separately in the scanner; moreover since at most one file name
can occur per command, we store it in a field |file_name| of the lexical
scanner reserved for that purpose. We allow a file name to be either a
sequence of non-space characters not starting with a quote, or a quoted
string.

@< Read in |file_name| @>=
if ((skip_space(),c=input.shift())=='"')
  file_name=scan_quoted_string();
else
{ file_name="";
  while (!std::isspace(c))
  {@; file_name.push_back(c); c=input.shift(); }
  input.unshift();
}

@* Index.

% Local IspellDict: british
