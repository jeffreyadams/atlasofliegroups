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
  let_list decls; /* declarations in a LET expression */
}

%locations
%parse-param { expr* parsed_expr}
%parse-param { int* verbosity }
%pure-parser
%error-verbose

%token QUIT LET IN TRUE FALSE QUIET VERBOSE WHATTYPE SHOWALL
%token DIVMOD
%token <val> INT
%token <expression> STRING
%token <id_code> IDENT
%token TOFILE ADDTOFILE FROMFILE
%type  <expression> exp quaternary lettail formula primary
%destructor { destroy_expr ($$); } exp quaternary lettail formula primary
%type  <expression_list> commalist commalistopt commabarlist idlist
%destructor { destroy_exprlist($$); } commalist commalistopt commabarlist idlist
%type <decls> declarations
%destructor { destroy_letlist($$); } declarations
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

exp     : quaternary ;


quaternary: LET lettail { $$=$2; }
        | formula
;

lettail : declarations IN quaternary { $$ = make_let_expr_node($1,$3); }
        | declarations ';' lettail  { $$ = make_let_expr_node($1,$3); }
;

formula : formula '+' formula
        { $$ = make_application_node
	       (lookup_identifier("+")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '-' formula
        { $$ = make_application_node
	       (lookup_identifier("-")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '*' formula
        { $$ = make_application_node
	       (lookup_identifier("*")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '/' formula
        { $$ = make_application_node
	       (lookup_identifier("/")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '%' formula
        { $$ = make_application_node
	       (lookup_identifier("%")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | '-' formula  %prec NEG
        { $$ = make_application_node
	       (lookup_identifier("-u")
	       ,make_exprlist_node($2,null_expr_list)
	       );
        }
	| formula DIVMOD formula
        { $$ = make_application_node
	       (lookup_identifier("/%")
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | primary
;

declarations: declarations ',' IDENT '=' exp { $$ = add_let_node($1,$3,$5); }
        | IDENT '=' exp { $$ = add_let_node(NULL,$1,$3); }
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
        | '[' commalistopt ']' { $$=wrap_list_display(reverse_expr_list($2)); }
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

commalistopt: commalist
        | /* empty */       { $$=null_expr_list; }
;

commalist: exp              { $$=make_exprlist_node($1,null_expr_list); }
	| commalist ',' exp { $$=make_exprlist_node($3,$1); }
;

commabarlist: commalistopt '|' commalistopt
        { $$ = make_exprlist_node(wrap_list_display(reverse_expr_list($3))
		 ,make_exprlist_node(wrap_list_display(reverse_expr_list($1))
		   ,null_expr_list));
        }
	| commabarlist '|' commalistopt
        { $$=make_exprlist_node(wrap_list_display(reverse_expr_list($3)),$1); }
;

idlist	: IDENT
        { $$=make_exprlist_node(make_applied_identifier($1),null_expr_list); }
	| idlist ',' IDENT
	{ $$=make_exprlist_node(make_applied_identifier($3),$1); }
	;


%%
