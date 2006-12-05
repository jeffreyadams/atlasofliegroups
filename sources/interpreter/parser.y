%{
  /*
   Copyright (C) 2006 Marc van Leeuwen
   This file is part of the Atlas of Reductive Lie Groups software (the Atlas)

   This program is made available under the terms stated in the GNU
   General Public License (GPL), see http://www.gnu.org/licences/licence.html

   The Atlas is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   The Atlas is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the Atlas; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  */


  #include <stdio.h>
  #include "parsetree.h"  /* types and functions used to construct nodes */
%}
%union {
  int   val;        /* For integral constants.  */
  short id_code;    /* For identifier codes  */
  expr	expression; /* For generic expressions */
  expr_list expression_list; /* For any list of expressions */
}

%locations
%parse-param { expr* parsed_expr}
%parse-param { int* verbosity }
%pure-parser
%error-verbose

%token QUIT TRUE FALSE QUIET VERBOSE WHATTYPE SHOWALL
%token DIVMOD
%token <val> INT
%token <expression> STRING
%token <id_code> IDENT
%token TOFILE ADDTOFILE FROMFILE
%type  <expression> exp primary
%destructor { destroy_expr ($$); } exp primary
%type  <expression_list> commalist commabarlist idlist
%destructor { destroy_exprlist($$); } commalist commabarlist idlist
%nonassoc DIVMOD
%left '-' '+'
%left '*' '/' '%'
%left NEG     /* negation--unary minus */

%{
  int yylex (YYSTYPE *, YYLTYPE *);
  void yyerror (YYLTYPE *, expr*, int*,const char *);
%}
%%

input:  '\n'			{ YYABORT } /* null input, skip evaluator */
        | exp    '\n'		{ *parsed_expr=$1; }
        | idlist '=' exp '\n'
		{ global_set_identifier(reverse_expr_list($1),$3); YYABORT }
        | QUIT	'\n'		{ *verbosity =-1; } /* causes immediate exit */
        | QUIET	'\n'		{ *verbosity =0; YYABORT } /* quiet mode */
        | VERBOSE '\n'		{ *verbosity =1; YYABORT } /* verbose mode */
	| TOFILE exp '\n'	{ *parsed_expr=$2; *verbosity=2; }
	| ADDTOFILE exp '\n'	{ *parsed_expr=$2; *verbosity=3; }
        | FROMFILE '\n'		{ include_file(); YYABORT } /* include file */
	| WHATTYPE exp '\n'     { type_of_expr($2); YYABORT } /* print type */
        | SHOWALL '\n'          { show_ids(); YYABORT } /* print id table */
;

exp     : exp '+' exp
        { $$ = make_application_node
	       (lookup_identifier("+")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | exp '-' exp
        { $$ = make_application_node
	       (lookup_identifier("-")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | exp '*' exp
        { $$ = make_application_node
	       (lookup_identifier("*")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | exp '/' exp
        { $$ = make_application_node
	       (lookup_identifier("/")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | exp '%' exp
        { $$ = make_application_node
	       (lookup_identifier("%")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | '-' exp  %prec NEG
        { $$ = make_application_node
	       (lookup_identifier("-u")
	       ,make_exprlist_node($2,null_expr_list)
	       );
        }
	| exp DIVMOD exp
        { $$ = make_application_node
	       (lookup_identifier("/%")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | primary
;

primary:  primary '[' exp ']' { $$ = make_subscription_node($1,$3); }
        | primary'[' commalist ',' exp ']'
          { $$=make_subscription_node
               ($1,wrap_tuple_display
                    (reverse_expr_list(make_exprlist_node($5,$3)))
               ) ;
	  }
	| IDENT '(' commalist ')'
		{ $$=make_application_node($1,reverse_expr_list($3)); }
        | INT { $$ = make_int_denotation($1); }
        | TRUE { $$ = make_bool_denotation(1); }
        | FALSE { $$ = make_bool_denotation(0); }
	| STRING
        | IDENT { $$=make_applied_identifier($1); }
        | '(' exp ')'          { $$=$2; }
        | '[' commalist ']'    { $$=wrap_list_display(reverse_expr_list($2)); }
	| '[' commabarlist ']'
          { $$=make_application_node
	        (lookup_identifier("transpose_mat")
                ,make_exprlist_node(wrap_list_display(reverse_expr_list($2))
				   ,null_expr_list)
                );
	  }
	| '(' ')' 		{ $$=wrap_tuple_display(NULL); }
        | '(' commalist ',' exp ')'
	{ $$=wrap_tuple_display(reverse_expr_list(make_exprlist_node($4,$2)));
        }
;

commalist:  /* empty */  { $$=null_expr_list; }
	| exp            { $$=make_exprlist_node($1,null_expr_list); }
	| commalist ',' exp { $$=make_exprlist_node($3,$1); }
;

commabarlist: commalist '|' commalist
        { $$ = make_exprlist_node(wrap_list_display(reverse_expr_list($3))
		 ,make_exprlist_node(wrap_list_display(reverse_expr_list($1))
		   ,null_expr_list));
        }
	| commabarlist '|' commalist
        { $$=make_exprlist_node(wrap_list_display(reverse_expr_list($3)),$1); }
;

idlist	: IDENT
        { $$=make_exprlist_node(make_applied_identifier($1),null_expr_list); }
	| idlist ',' IDENT
	{ $$=make_exprlist_node(make_applied_identifier($3),$1); }
	;


%%
