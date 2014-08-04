%{
  /*
   Copyright (C) 2006-2012 Marc van Leeuwen
   This file is part of the Atlas of Lie Groups and Representations (the Atlas)

   This program is made available under the terms stated in the GNU
   General Public License (GPL), see http://www.gnu.org/licences/licence.html

   The Atlas is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   The Atlas is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the Atlas; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  */


  #include <stdio.h>
#include "parsetree.h"  // types and functions used to construct nodes
#include "evaluator.h" // action functions invoked from within the parser
  using namespace atlas::interpreter; // to allow simplifying action code
%}
%union {
  int	val;	    /* For integral constants.	*/
  short id_code;    /* For identifier codes  */
  struct { short id, priority; } oper; /* for operator symbols */
  atlas::interpreter::raw_form_stack ini_form;
  unsigned short type_code; /* For type names */
  atlas::interpreter::expr_p    expression; /* For generic expressions */
  atlas::interpreter::raw_expr_list expression_list; /* list of expressions */
  atlas::interpreter::raw_let_list decls; /* declarations in a LET expression */
  atlas::interpreter::raw_id_pat ip;
  struct {
    atlas::interpreter::raw_type_list typel;
    atlas::interpreter::raw_patlist patl;
  } id_sp;
  atlas::interpreter::raw_patlist pl;
  atlas::interpreter::type_p type_pt;
  atlas::interpreter::raw_type_list type_l;
}

%locations
%parse-param { atlas::interpreter::expr_p* parsed_expr}
%parse-param { int* verbosity }
%pure-parser
%error-verbose

%token QUIT SET LET IN BEGIN END IF THEN ELSE ELIF FI AND OR NOT
%token WHILE DO OD NEXT FOR FROM DOWNTO
%token TRUE FALSE QUIET VERBOSE WHATTYPE SHOWALL FORGET
%token <oper> OPERATOR '='
%token <val> INT
%token <expression> STRING '$'
%token <id_code> IDENT
%token TOFILE ADDTOFILE FROMFILE FORCEFROMFILE

%token <type_code> TYPE
%token ARROW "->"
%token BECOMES ":="

%type <expression> exp quaternary tertiary lettail or_expr and_expr not_expr
%type <expression>  formula operand secondary primary iftail
%type <expression> comprim subscription
%type <ini_form> formula_start
%type <oper> operator
%destructor { destroy_expr ($$); } exp quaternary tertiary lettail or_expr
%destructor { destroy_expr ($$); } and_expr not_expr formula operand iftail
%destructor { destroy_expr ($$); } secondary primary comprim subscription
%destructor { destroy_formula($$); } formula_start
%type  <expression_list> commalist commalist_opt commabarlist
%destructor { destroy_exprlist($$); } commalist commalist_opt commabarlist
%type <decls> declarations declaration
%destructor { destroy_letlist($$); } declarations declaration

%left OPERATOR

%type <ip> pattern pattern_opt
%destructor { destroy_id_pat($$); } pattern pattern_opt
%type <pl> pat_list
%destructor { destroy_pattern($$); } pat_list
%type <type_pt> type
%destructor { destroy_type($$); } type
%type <type_l> types types_opt
%destructor { destroy_type_list($$); } types types_opt
%type <id_sp> id_specs id_specs_opt
%destructor { destroy_type_list($$.typel);destroy_pattern($$.patl); } id_specs id_specs_opt

%{
  int yylex (YYSTYPE *, YYLTYPE *);
  void yyerror (YYLTYPE *, atlas::interpreter::expr_p*, int*,const char *);
%}
%%

input:	'\n'			{ YYABORT; } /* null input, skip evaluator */
	| '\f'	   { YYABORT; } /* allow form feed as well at command level */
	| exp '\n'		{ *parsed_expr=$1; }
	| tertiary ';' '\n'
	  { *parsed_expr=make_sequence($1,wrap_tuple_display(NULL),1); }
	| SET pattern '=' exp '\n' { global_set_identifier($2,$4,1); YYABORT; } 
	| SET IDENT '(' id_specs_opt ')' '=' exp '\n'
	  { struct raw_id_pat id; id.kind=0x1; id.name=$2;
	    global_set_identifier(id,make_lambda_node($4.patl,$4.typel,$7),1);
	    YYABORT;
	  }
	| FORGET IDENT '\n'	{ global_forget_identifier($2); YYABORT; }
	| SET operator '(' id_specs ')' '=' exp '\n'
	  { struct raw_id_pat id; id.kind=0x1; id.name=$2.id;
	    global_set_identifier(id,make_lambda_node($4.patl,$4.typel,$7),2);
	    YYABORT;
	  }
	| SET operator '=' exp '\n'
	  { struct raw_id_pat id; id.kind=0x1; id.name=$2.id;
	    global_set_identifier(id,$4,2); YYABORT;
	  }
	| FORGET IDENT '@' type '\n'
	  { global_forget_overload($2,$4); YYABORT;  }
	| FORGET operator '@' type '\n'
	  { global_forget_overload($2.id,$4); YYABORT; }
	| IDENT ':' exp '\n'
		{ struct raw_id_pat id; id.kind=0x1; id.name=$1;
		  global_set_identifier(id,$3,0); YYABORT; }
	| IDENT ':' type '\n'	{ global_declare_identifier($1,$3); YYABORT; }
	| QUIT	'\n'		{ *verbosity =-1; } /* causes immediate exit */
	| QUIET '\n'		{ *verbosity =0; YYABORT; } /* quiet mode */
	| VERBOSE '\n'		{ *verbosity =1; YYABORT; } /* verbose mode */
	| TOFILE exp '\n'	{ *parsed_expr=$2; *verbosity=2; }
	| ADDTOFILE exp '\n'	{ *parsed_expr=$2; *verbosity=3; }
	| FROMFILE '\n'		{ include_file(1); YYABORT; } /* include file */
	| FORCEFROMFILE '\n'	{ include_file(0); YYABORT; } /* force include */
	| WHATTYPE exp '\n'	{ type_of_expr($2); YYABORT; } /* print type */
	| WHATTYPE operator '?' '\n'
			     { show_overloads($2.id); YYABORT; } /* show types */
	| WHATTYPE IDENT '?' '\n'
				{ show_overloads($2); YYABORT; } /* show types */
	| SHOWALL '\n'		{ show_ids(); YYABORT; } /* print id table */
;

exp: LET lettail { $$=$2; }
	| '(' ')' ':' exp { $$=make_lambda_node(NULL,raw_type_list(NULL),$4); }
	| '(' id_specs ')' ':' exp
	  { $$=make_lambda_node($2.patl,$2.typel,$5); }
	| '(' ')' type ':' exp
	  { $$=make_lambda_node(NULL,raw_type_list(NULL),make_cast($3,$5)); }
	| '(' id_specs ')' type ':' exp
	  { $$=make_lambda_node($2.patl,$2.typel,make_cast($4,$6)); }
        | type ':' exp	 { $$ = make_cast($1,$3); }
	| quaternary NEXT exp { $$=make_sequence($1,$3,0); }
        | quaternary
;

lettail : declarations IN exp { $$ = make_let_expr_node($1,$3); }
	| declarations THEN lettail  { $$ = make_let_expr_node($1,$3); }
;

declarations: declarations ',' declaration { $$ = append_let_node($1,$3); }
        | declaration
;

declaration: pattern '=' exp { $$ = make_let_node($1,$3); }
        | IDENT '(' id_specs_opt ')' '=' exp
	  { struct raw_id_pat p; p.kind=0x1; p.name=$1;
	    $$ = make_let_node(p,make_lambda_node($3.patl,$3.typel,$6));
	  }
;

quaternary: tertiary ';' quaternary { $$=make_sequence($1,$3,1); }
	| tertiary
;

tertiary: IDENT BECOMES tertiary { $$ = make_assignment($1,$3); }
	| subscription BECOMES tertiary { $$ = make_comp_ass($1,$3); }
	| IDENT operator BECOMES tertiary
	{ $$ = make_assignment($1,
		make_binary_call($2.id,make_applied_identifier($1),$4)); }
	| or_expr
;

or_expr : and_expr OR or_expr
	  { $$ = make_conditional_node($1,make_bool_denotation(true),$3); }
	| and_expr
;

and_expr: not_expr AND and_expr
	  { $$ = make_conditional_node($1,$3,make_bool_denotation(false)); }
	| not_expr
;

not_expr: NOT secondary
	  { $$ = make_conditional_node($2,make_bool_denotation(false),
					  make_bool_denotation(true)); }
	| secondary
;

secondary : formula
	| '(' ')' /* don't allow this as first part in subscription or call */
	  { $$=wrap_tuple_display(NULL); }
	| primary
;

formula : formula_start operand { $$=end_formula($1,$2); }
;
formula_start : operator       { $$=start_unary_formula($1.id,$1.priority); }
	| comprim operator     { $$=start_formula($1,$2.id,$2.priority); }
	| IDENT operator
	  { $$=start_formula(make_applied_identifier($1),$2.id,$2.priority); }
	| formula_start operand operator
	  { $$=extend_formula($1,$2,$3.id,$3.priority); }
;


operator : OPERATOR | '=';

operand : operator operand { $$=make_unary_call($1.id,$2); }
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
	| TRUE { $$ = make_bool_denotation(true); }
	| FALSE { $$ = make_bool_denotation(false); }
	| STRING
        | '$' { $$=make_dollar(); }
	| IF iftail { $$=$2; }
	| WHILE exp DO exp OD { $$=make_while_node($2,$4); }
	| FOR pattern IN exp DO exp OD
	  { struct raw_id_pat p,x; p.kind=0x2; x.kind=0x0;
	    p.sublist=make_pattern_node(make_pattern_node(NULL,$2),x);
	    $$=make_for_node(p,$4,$6);
	  }
	| FOR pattern '@' IDENT IN exp DO exp OD
	  { struct raw_id_pat p,i; p.kind=0x2; i.kind=0x1; i.name=$4;
	    p.sublist=make_pattern_node(make_pattern_node(NULL,$2),i);
	    $$=make_for_node(p,$6,$8);
	  }
	| FOR IDENT ':' exp FROM exp DO exp OD
	  { $$=make_cfor_node($2,$4,$6,true,$8); }
	| FOR IDENT ':' exp DOWNTO exp DO exp OD
	  { $$=make_cfor_node($2,$4,$6,false,$8); }
	| FOR IDENT ':' exp DO exp OD
	  { $$=make_cfor_node($2,$4,wrap_tuple_display(NULL),true,$6); }
	| '(' exp ')'	       { $$=$2; }
	| BEGIN exp END	       { $$=$2; }
	| '[' commalist_opt ']'
		{ $$=wrap_list_display(reverse_expr_list($2)); }
	| '[' commabarlist ']'
	  { $$=make_unary_call
		(lookup_identifier("^"),
                 make_cast
                 (make_prim_type(5) /* |matrix_type| */
		 ,wrap_list_display(reverse_expr_list($2))));
	  }
	| '(' commalist ',' exp ')'
	{ $$=wrap_tuple_display(reverse_expr_list(make_exprlist_node($4,$2)));
	}
	| operator '@' type { $$=make_op_cast($1.id,$3); }
	| IDENT '@' type    { $$=make_op_cast($1,$3); }
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

iftail	: exp THEN exp ELSE exp FI { $$=make_conditional_node($1,$3,$5); }
	| exp THEN exp ELIF iftail { $$=make_conditional_node($1,$3,$5); }
	| exp THEN exp FI
	  { $$=make_conditional_node($1,$3,wrap_tuple_display(NULL)); }
;

pattern : IDENT		    { $$.kind=0x1; $$.name=$1; }
	| '(' pat_list ')'  { $$.kind=0x2; $$.sublist=reverse_patlist($2); }
	| '(' pat_list ')' ':' IDENT
	    { $$.kind=0x3; $$.name=$5; $$.sublist=reverse_patlist($2);}
	| '(' ')' { $$.kind=0x2; $$.sublist=0; } /* allow throw-away value */
;

pattern_opt :/* empty */ { $$.kind=0x0; }
	| pattern
;

pat_list: pattern_opt ',' pattern_opt
	  { $$=make_pattern_node(make_pattern_node(NULL,$1),$3); }
	| pat_list ',' pattern_opt { $$=make_pattern_node($1,$3); }
;

id_specs: type pattern
	{ $$.typel=make_type_singleton($1);
	  $$.patl=make_pattern_node(NULL,$2);
	}
	| type pattern ',' id_specs
	{ $$.typel=make_type_list($1,$4.typel);
	  $$.patl=make_pattern_node($4.patl,$2);
	}
;

id_specs_opt: id_specs
	| /* empty */ { $$.typel=NULL; $$.patl=NULL; }
;

type	: TYPE			  { $$=make_prim_type($1); }
        | '(' type ')'            { $$=$2; }
	| '(' type ARROW type ')' { $$=make_function_type($2,$4); }
	| '(' types_opt ARROW type ')'
			  { $$=make_function_type(make_tuple_type($2),$4); }
	| '(' type ARROW types_opt ')'
			  { $$=make_function_type($2,make_tuple_type($4)); }
	| '(' types_opt ARROW types_opt ')'
	   { $$=make_function_type(make_tuple_type($2),make_tuple_type($4)); }
	| '[' type ']'		  { $$=make_row_type($2); }
	| '(' types ')'		  { $$=make_tuple_type($2); }
;

types	: type ',' type	 { $$=make_type_list($1,make_type_singleton($3)); }
	| type ',' types { $$=make_type_list($1,$3); }
;

types_opt : /* empty */ { $$=NULL; }
	| types
;

commalist_opt: /* empty */	 { $$=raw_expr_list(nullptr); }
	| commalist
;

commalist: exp		    { $$=make_exprlist_node($1,raw_expr_list(nullptr)); }
	| commalist ',' exp { $$=make_exprlist_node($3,$1); }
;

commabarlist: commalist_opt '|' commalist_opt
	{ $$ = make_exprlist_node(wrap_list_display(reverse_expr_list($3))
		 ,make_exprlist_node(wrap_list_display(reverse_expr_list($1))
		   ,raw_expr_list(nullptr)));
	}
	| commabarlist '|' commalist_opt
	{ $$=make_exprlist_node(wrap_list_display(reverse_expr_list($3)),$1); }
;


%%
