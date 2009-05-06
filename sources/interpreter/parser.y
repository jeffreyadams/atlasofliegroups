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
  unsigned short type_code; /* For type names */
  expr	expression; /* For generic expressions */
  expr_list expression_list; /* For any list of expressions */
  let_list decls; /* declarations in a LET expression */
  struct id_pat ip;
  struct { ptr typel; patlist patl; } id_sp;
  patlist pl;
  ptr gen_ptr;
}

%locations
%parse-param { expr* parsed_expr}
%parse-param { int* verbosity }
%pure-parser
%error-verbose

%token QUIT LET IN TRUE FALSE QUIET VERBOSE WHATTYPE SHOWALL
%token DIVMOD "/%"
%token <val> INT
%token <expression> STRING
%token <id_code> IDENT
%token TOFILE ADDTOFILE FROMFILE
%type  <expression> exp quaternary tertiary lettail formula primary comprim subscription
%destructor { destroy_expr ($$); } exp quaternary tertiary lettail formula primary comprim subscription
%type  <expression_list> commalist commalist_opt commabarlist idlist
%destructor { destroy_exprlist($$); } commalist commalist_opt commabarlist idlist
%type <decls> declarations
%destructor { destroy_letlist($$); } declarations
%nonassoc DIVMOD
%left '-' '+'
%left '*' '/' '%'
%left NEG     /* negation--unary minus */

%type <ip> pattern pattern_opt
%destructor { destroy_id_pat(&$$); } pattern pattern_opt
%type <pl> pat_list
%destructor { destroy_pattern($$); } pat_list
%token <type_code> TYPE
%type <gen_ptr> type types types_opt
%destructor { destroy_type($$); } type
%destructor { destroy_type_list($$); } types types_opt
%type <id_sp> id_specs
%destructor { destroy_type_list($$.typel);destroy_pattern($$.patl); } id_specs

%token ARROW BECOMES

%{
  int yylex (YYSTYPE *, YYLTYPE *);
  void yyerror (YYLTYPE *, expr*, int*,const char *);
%}
%%

input:  '\n'			{ YYABORT } /* null input, skip evaluator */
        | exp    '\n'	        { *parsed_expr=$1; }
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

exp     : type ':' exp      { $$ = make_cast($1,$3); }
        | quaternary
;


quaternary: tertiary ';' quaternary { $$=make_sequence($1,$3); }
        | tertiary
;

tertiary: LET lettail { $$=$2; }
        | '(' ')' ':' tertiary { $$=make_lambda_node(NULL,NULL,$4); }
        | '(' id_specs ')' ':' tertiary
          { $$=make_lambda_node($2.patl,$2.typel,$5); }
        | '(' ')' type ':' tertiary
          { $$=make_lambda_node(NULL,NULL,make_cast($3,$5)); }
        | '(' id_specs ')' type ':' tertiary
          { $$=make_lambda_node($2.patl,$2.typel,make_cast($4,$6)); }
        | IDENT BECOMES tertiary { $$ = make_assignment($1,$3); }
        | subscription BECOMES tertiary { $$ = make_comp_ass($1,$3); }
        | formula
;

lettail : declarations IN tertiary { $$ = make_let_expr_node($1,$3); }
        | declarations ';' lettail  { $$ = make_let_expr_node($1,$3); }
;

declarations: declarations ',' pattern '=' tertiary
          { $$ = add_let_node($1,$3,$5); }
        | pattern '=' tertiary { $$ = add_let_node(NULL,$1,$3); }
;

formula : formula '+' formula
        { $$ = make_application_node
	       (make_applied_identifier(lookup_identifier("+"))
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '-' formula
        { $$ = make_application_node
	       (make_applied_identifier(lookup_identifier("-"))
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '*' formula
        { $$ = make_application_node
	       (make_applied_identifier(lookup_identifier("*"))
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '/' formula
        { $$ = make_application_node
	       (make_applied_identifier(lookup_identifier("/"))
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | formula '%' formula
        { $$ = make_application_node
	       (make_applied_identifier(lookup_identifier("%"))
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
        | '-' formula  %prec NEG
        { $$ = make_application_node
	       (make_applied_identifier(lookup_identifier("-u"))
	       ,make_exprlist_node($2,null_expr_list)
	       );
        }
	| formula DIVMOD formula
        { $$ = make_application_node
	       (make_applied_identifier(lookup_identifier("/%"))
	       ,make_exprlist_node($1,make_exprlist_node($3,null_expr_list))
	       );
        }
	| '(' ')' /* don't allow this as first part in subscription or call */
          { $$=wrap_tuple_display(NULL); }
        | primary
;

primary: comprim
        | IDENT { $$=make_applied_identifier($1); }
;
comprim: subscription
	| comprim '[' exp ']' { $$ = make_subscription_node($1,$3); }
        | comprim '[' commalist ',' exp ']'
          { $$=make_subscription_node
               ($1,wrap_tuple_display
                    (reverse_expr_list(make_exprlist_node($5,$3)))
               ) ;
	  }
	| primary '(' commalist_opt ')'
		{ $$=make_application_node($1,reverse_expr_list($3)); }
        | INT { $$ = make_int_denotation($1); }
        | TRUE { $$ = make_bool_denotation(1); }
        | FALSE { $$ = make_bool_denotation(0); }
	| STRING
        | '(' exp ')'          { $$=$2; }
        | '[' commalist_opt ']'
                { $$=wrap_list_display(reverse_expr_list($2)); }
	| '[' commabarlist ']'
          { $$=make_application_node
	        (make_applied_identifier(lookup_identifier("transpose_mat"))
                ,make_exprlist_node(wrap_list_display(reverse_expr_list($2))
				   ,null_expr_list)
                );
	  }
        | '(' commalist ',' exp ')'
	{ $$=wrap_tuple_display(reverse_expr_list(make_exprlist_node($4,$2)));
        }
;

subscription: IDENT '[' exp ']'
          { $$ = make_subscription_node(make_applied_identifier($1),$3); }
        | IDENT '[' commalist ',' exp ']'
          { $$=make_subscription_node
               (make_applied_identifier($1),
                wrap_tuple_display
                    (reverse_expr_list(make_exprlist_node($5,$3)))
               ) ;
	  }
;

pattern : IDENT             { $$.kind=0x1; $$.name=$1; }
        | '(' pat_list ')'  { $$.kind=0x2; $$.sublist=reverse_patlist($2); }
        | IDENT ':' '(' pat_list ')'
            { $$.kind=0x3; $$.name=$1; $$.sublist=reverse_patlist($4);}
;

pattern_opt :/* empty */ { $$.kind=0x0; }
        | pattern
;

pat_list: pattern_opt ',' pattern_opt
          { $$=make_pattern_node(make_pattern_node(NULL,&$1),&$3); }
        | pat_list ',' pattern_opt { $$=make_pattern_node($1,&$3); }
;

id_specs: type pattern
        { $$.typel=mk_type_singleton($1);
          $$.patl=make_pattern_node(NULL,&$2);
        }
        | type pattern ',' id_specs
        { $$.typel=mk_type_list($1,$4.typel);
          $$.patl=make_pattern_node($4.patl,&$2);
        }
;

type    : TYPE                    { $$=mk_prim_type($1); }
        | '(' type ARROW type ')' { $$=mk_function_type($2,$4); }
        | '(' types_opt ARROW type ')'
                          { $$=mk_function_type(mk_tuple_type($2),$4); }
        | '(' type ARROW types_opt ')'
                          { $$=mk_function_type($2,mk_tuple_type($4)); }
        | '(' types_opt ARROW types_opt ')'
           { $$=mk_function_type(mk_tuple_type($2),mk_tuple_type($4)); }
        | '[' type ']'            { $$=mk_row_type($2); }
        | '(' types ')'           { $$=mk_tuple_type($2); }
;

types   : type ',' type  { $$=mk_type_list($1,mk_type_singleton($3)); }
        | type ',' types { $$=mk_type_list($1,$3); }
;

types_opt : /* empty */ { $$=NULL; }
        | types
;

commalist_opt: /* empty */       { $$=null_expr_list; }
        | commalist
;

commalist: exp              { $$=make_exprlist_node($1,null_expr_list); }
	| commalist ',' exp { $$=make_exprlist_node($3,$1); }
;

commabarlist: commalist_opt '|' commalist_opt
        { $$ = make_exprlist_node(wrap_list_display(reverse_expr_list($3))
		 ,make_exprlist_node(wrap_list_display(reverse_expr_list($1))
		   ,null_expr_list));
        }
	| commabarlist '|' commalist_opt
        { $$=make_exprlist_node(wrap_list_display(reverse_expr_list($3)),$1); }
;

idlist	: IDENT
        { $$=make_exprlist_node(make_applied_identifier($1),null_expr_list); }
	| idlist ',' IDENT
	{ $$=make_exprlist_node(make_applied_identifier($3),$1); }
	;


%%
