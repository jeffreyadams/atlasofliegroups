%{
  /*
   Copyright (C) 2006-2015 Marc van Leeuwen
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


#include <iostream>
#include "parsetree.h"  // types and functions used to construct nodes
#include "lexer.h"  // pointer |lex| to lexical analyser
#include "evaluator.h" // action functions invoked from within the parser
  using namespace atlas::interpreter; // to allow simplifying action code
%}
%union {
  int	val;	    /* For integral constants.	*/
  short id_code;    /* For identifier codes, or defined types  */
  std::string* str;  // string denotation
  struct { short id, priority; } oper; /* for operator symbols */
  atlas::interpreter::raw_form_stack ini_form;
  unsigned short type_code; /* For type names */
  atlas::interpreter::expr_p    expression; /* For generic expressions */
  atlas::interpreter::raw_expr_list expression_list; /* list of expressions */
  atlas::interpreter::raw_let_list decls; /* declarations in a LET expression */
  atlas::interpreter::raw_id_pat ip;
  struct {
    atlas::interpreter::type_p type_pt;
    atlas::interpreter::raw_id_pat ip;
  } id_sp1;
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
%token WHILE DO OD NEXT FOR FROM DOWNTO CASE ESAC REC_FUN
%token TRUE FALSE DIE WHATTYPE SHOWALL FORGET
%token <oper> OPERATOR OPERATOR_BECOMES '=' '*'
%token <val> INT
%token <str> STRING
%token <id_code> IDENT TYPE_ID
%token TOFILE ADDTOFILE FROMFILE FORCEFROMFILE

%token <type_code> TYPE
%token ARROW "->"
%token BECOMES ":="
%token TLSUB "~["

%type <expression> expr expr_opt tertiary lettail or_expr and_expr not_expr
%type <expression>  formula operand secondary primary iftail
%type <expression> subscription slice comprim assignable_subsn ident_expr
%type <ini_form> formula_start
%type <oper> operator
%type <val> tilde_opt
%destructor { destroy_expr ($$); } expr expr_opt tertiary lettail or_expr
%destructor { destroy_expr ($$); } and_expr not_expr formula operand iftail
%destructor { destroy_expr ($$); } secondary primary comprim subscription slice
%destructor { destroy_expr ($$); } assignable_subsn ident_expr
%destructor { destroy_formula($$); } formula_start
%destructor { delete $$; } STRING
%type  <expression_list> commalist commalist_opt commabarlist
%destructor { destroy_exprlist($$); } commalist commalist_opt commabarlist
%type <decls> declarations declaration
%destructor { destroy_letlist($$); } declarations declaration

%left OPERATOR

%type <ip> pattern pattern_opt
%destructor { destroy_id_pat($$); } pattern pattern_opt
%type <pl> pat_list
%destructor { destroy_pattern($$); } pat_list
%type <type_pt> nostar_type type
%destructor { destroy_type($$); } nostar_type type
%type <type_l> types types_opt
%destructor { destroy_type_list($$); } types types_opt
%type <id_sp1> id_spec
%destructor { destroy_type($$.type_pt);destroy_id_pat($$.ip); } id_spec
%type <id_sp> id_specs id_specs_opt
%destructor { destroy_type_list($$.typel);destroy_pattern($$.patl); } id_specs id_specs_opt

%{
  int yylex (YYSTYPE *, YYLTYPE *);
  void yyerror (YYLTYPE *, atlas::interpreter::expr_p*, int*,const char *);
%}
%%

input:	'\n'			{ YYABORT; } /* null input, skip evaluator */
	| '\f'	   { YYABORT; } /* allow form feed as well at command level */
	| expr '\n'		{ *parsed_expr=$1; }
	| SET declarations '\n' { global_set_identifiers($2); YYABORT; }
	| FORGET IDENT '\n'	{ global_forget_identifier($2); YYABORT; }
	| FORGET TYPE_ID '\n'	{ global_forget_identifier($2); YYABORT; }
	| SET operator '(' id_specs ')' '=' expr '\n'
	  { struct raw_id_pat id; id.kind=0x1; id.name=$2.id;
	    global_set_identifier(id,
				  make_lambda_node($4.patl,$4.typel,$7,@$),
				  2);
	    YYABORT;
	  }
	| SET operator '=' expr '\n'
	  { struct raw_id_pat id; id.kind=0x1; id.name=$2.id;
	    global_set_identifier(id,$4,2); YYABORT;
	  }
	| FORGET IDENT '@' type '\n'
	  { global_forget_overload($2,$4); YYABORT;  }
	| FORGET operator '@' type '\n'
	  { global_forget_overload($2.id,$4); YYABORT; }
	| IDENT ':' expr '\n'
		{ struct raw_id_pat id; id.kind=0x1; id.name=$1;
		  global_set_identifier(id,$3,0); YYABORT; }
	| IDENT ':' type '\n'	{ global_declare_identifier($1,$3); YYABORT; }
	| ':' IDENT '=' type '\n' { type_define_identifier($2,$4); YYABORT; }
	| ':' TYPE_ID '=' type '\n' { type_define_identifier($2,$4); YYABORT; }
	| QUIT	'\n'		{ *verbosity =-1; } /* causes immediate exit */
	| SET IDENT // set an option; option identifiers have lowest codes
          { unsigned n=$2-lex->first_identifier();
            if (n<2)
              *verbosity=n; // |quiet| gives 0, and |verbose| gives 1
            else
	      std::cerr << '\'' << main_hash_table->name_of($2)
                        << "' is not something one can set" << std::endl;
            YYABORT;
	  }
	| TOFILE expr '\n'	{ *parsed_expr=$2; *verbosity=2; }
	| ADDTOFILE expr '\n'	{ *parsed_expr=$2; *verbosity=3; }
	| FROMFILE '\n'		{ include_file(1); YYABORT; } /* include file */
	| FORCEFROMFILE '\n'	{ include_file(0); YYABORT; } // force include
	| WHATTYPE expr '\n'	{ type_of_expr($2); YYABORT; } // print type
	| WHATTYPE operator '?' '\n'
			      { show_overloads($2.id); YYABORT; } // show types
	| WHATTYPE IDENT '?' '\n'
				{ show_overloads($2); YYABORT; } // show types
	| SHOWALL '\n'		{ show_ids(); YYABORT; } /* print id table */
;

expr    : LET lettail { $$=$2; }
	| '@' ':' expr { $$=make_lambda_node(nullptr,nullptr,$3,@$); }
	| '@' type ':' expr
	  { $$=make_lambda_node(nullptr,nullptr,make_cast($2,$4,@$),@$); }
	| '(' id_specs ')' ':' expr
	  { $$=make_lambda_node($2.patl,$2.typel,$5,@$); }
	| '(' id_specs ')' type ':' expr
	  { $$=make_lambda_node($2.patl,$2.typel,make_cast($4,$6,@$),@$); }
        | REC_FUN IDENT '(' id_specs ')' type ':' expr
	  { auto l=make_lambda_node($4.patl,$4.typel,make_cast($6,$8,@8),@$);
            $$ = make_recfun($2,l,@$,@2);
          }
        | type ':' expr { $$ = make_cast($1,$3,@$); }
	| tertiary ';' expr { $$=make_sequence($1,$3,true,@$); }
        | tertiary NEXT expr { $$=make_sequence($1,$3,false,@$); }
	| tertiary
;

lettail : declarations IN expr { $$ = make_let_expr_node($1,$3,@$); }
	| declarations THEN lettail  { $$ = make_let_expr_node($1,$3,@$); }
;

declarations: declarations ',' declaration { $$ = append_let_node($1,$3); }
        | declaration
;

declaration: pattern '=' expr { $$ = make_let_node($1,$3); }
        | IDENT '(' id_specs_opt ')' '=' expr
	  { struct raw_id_pat p; p.kind=0x1; p.name=$1;
	    $$ = make_let_node(p,make_lambda_node($3.patl,$3.typel,$6,@$));
	  }
;

tertiary: IDENT BECOMES tertiary { $$ = make_assignment($1,$3,@$); }
	| assignable_subsn BECOMES tertiary { $$ = make_comp_ass($1,$3,@$); }
	| IDENT OPERATOR_BECOMES tertiary
	  { $$ = make_assignment($1,
		  make_binary_call($2.id,
		    make_applied_identifier($1,@1),$3,@$,@2),@$); }
	| assignable_subsn OPERATOR_BECOMES tertiary
	  { $$ = make_comp_upd_ass($1,$2.id,$3,@$,@2); }
	| or_expr
;

or_expr : or_expr OR and_expr
	  { $$ =
	      make_conditional_node($1,make_bool_denotation(true,@$),$3,@$); }
	| and_expr
;

and_expr: and_expr AND not_expr
	  { $$ =
	      make_conditional_node($1,$3,make_bool_denotation(false,@$),@$); }
	| not_expr
;

not_expr: NOT not_expr { $$ = make_negation($2,@$); }
	| secondary
;

secondary : formula
	| '(' ')' /* don't allow this as first part in subscription or call */
	  { $$=wrap_tuple_display(nullptr,@$); }
	| primary
;

formula : formula_start operand { $$=end_formula($1,$2,@$); }
;
formula_start : operator       { $$=start_unary_formula($1.id,$1.priority,@1); }
	| comprim operator     { $$=start_formula($1,$2.id,$2.priority,@2); }
	| IDENT operator       { $$=start_formula
	      (make_applied_identifier($1,@1),$2.id,$2.priority,@2); }
	| formula_start operand operator
	  { $$=extend_formula($1,$2,$3.id,$3.priority,@3); }
;


operator : OPERATOR | '=' | '*' ;

operand : operator operand { $$=make_unary_call($1.id,$2,@$,@1); }
	| primary
;


tilde_opt : '~' { $$ = 1; }
	| { $$ = 0; }
;

primary: comprim | ident_expr ;
ident_expr : IDENT { $$=make_applied_identifier($1,@1); } ;

comprim: subscription | slice
	| primary '(' commalist_opt ')'
	{ $$=make_application_node($1,reverse_expr_list($3),@$,@2,@4); }
	| INT { $$ = make_int_denotation($1,@$); }
	| TRUE { $$ = make_bool_denotation(true,@$); }
	| FALSE { $$ = make_bool_denotation(false,@$); }
	| STRING { $$ = make_string_denotation($1,@$); }
        | '$' { $$=make_dollar(@$); }
	| IF iftail { $$=$2; }
	| CASE expr IN commalist ESAC
	  { $$=make_int_case_node($2,reverse_expr_list($4),@$); }
	| WHILE expr DO expr OD { $$=make_while_node($2,$4,@$); }
	| FOR pattern_opt IN expr tilde_opt DO expr tilde_opt OD
	  { struct raw_id_pat p,x; p.kind=0x2; x.kind=0x0;
	    p.sublist=make_pattern_node(make_pattern_node(nullptr,$2),x);
	    $$=make_for_node(p,$4,$7,$5+2*$8,@$);
	  }
	| FOR pattern_opt '@' IDENT IN expr tilde_opt DO expr tilde_opt OD
	  { struct raw_id_pat p,i; p.kind=0x2; i.kind=0x1; i.name=$4;
	    p.sublist=make_pattern_node(make_pattern_node(nullptr,$2),i);
	    $$=make_for_node(p,$6,$9,$7+2*$10,@$);
	  }
	| FOR IDENT ':' expr tilde_opt DO expr tilde_opt OD
	  { $$=make_cfor_node($2,$4,wrap_tuple_display(nullptr,@$)
                             ,$7,$5+2*$8,@$); }
	| FOR IDENT ':' expr FROM expr tilde_opt DO expr tilde_opt OD
	  { $$=make_cfor_node($2,$4,$6,$9,$7+2*$10,@$); }
	| FOR ':' expr DO expr tilde_opt OD
	  { $$=make_cfor_node(-1,$3,wrap_tuple_display(nullptr,@$)
                             ,$5,2*$6+4,@$); }
	| FOR ':' expr FROM expr DO expr tilde_opt OD
	  { $$=make_cfor_node(-1,$3,$5,$7,2*$8+4,@$); }
	| FOR IDENT ':' expr DOWNTO expr DO expr OD
	  { $$=make_cfor_node($2,$4,$6,$8,1,@$); }
	| '(' expr ')'	       { $$=$2; }
	| BEGIN expr END	       { $$=$2; }
	| '[' commalist_opt ']'
		{ $$=wrap_list_display(reverse_expr_list($2),@$); }
	| '[' commabarlist ']'
	  { $$=make_unary_call
		(lookup_identifier("transpose "),
                 make_cast
                 (make_prim_type(5) /* |matrix_type| */
		  ,wrap_list_display(reverse_expr_list($2),@$),@$),@$,@1);
	  }
	| '(' commalist ',' expr ')'
	{ $$=wrap_tuple_display
	    (reverse_expr_list(make_exprlist_node($4,$2)),@$);
	}
	| operator '@' type { $$=make_op_cast($1.id,$3,@$); }
	| IDENT '@' type    { $$=make_op_cast($1,$3,@$); }
	| DIE { $$= make_die(@$); }
;

assignable_subsn:
          IDENT '[' expr ']'
	  { $$ = make_subscription_node
		  (make_applied_identifier($1,@1),$3,false,@$); }
	| IDENT TLSUB expr ']'
	  { $$ = make_subscription_node
		  (make_applied_identifier($1,@1),$3,true,@$); }
	| IDENT '[' expr ',' expr ']'
	  { $$=make_subscription_node(make_applied_identifier($1,@1),
		wrap_tuple_display
		(make_exprlist_node($3,
                   make_exprlist_node($5,raw_expr_list(nullptr))),@$)
		,false,@$);
          }
	| IDENT TLSUB expr ',' expr ']'
	  { $$=make_subscription_node(make_applied_identifier($1,@1),
		wrap_tuple_display
		(make_exprlist_node($3,
                   make_exprlist_node($5,raw_expr_list(nullptr))),@$)
		,true,@$);
          }
;

subscription: assignable_subsn
	| comprim '[' expr ']'
	  { $$ = make_subscription_node($1,$3,false,@$); }
	| comprim TLSUB expr ']'
	  { $$ = make_subscription_node($1,$3,true,@$); }
	| comprim '[' expr ',' expr ']'
	  { $$=make_subscription_node($1,
		wrap_tuple_display
		(make_exprlist_node($3,
                   make_exprlist_node($5,raw_expr_list(nullptr))),@$)
		,false,@$);
          }
	| comprim TLSUB expr ',' expr ']'
	  { $$=make_subscription_node($1,
		wrap_tuple_display
		(make_exprlist_node($3,
                   make_exprlist_node($5,raw_expr_list(nullptr))),@$)
		,true,@$);
          }
;

expr_opt : expr | { $$=nullptr; } ;

slice   : IDENT '[' expr_opt tilde_opt ':' expr_opt tilde_opt ']'
	  { unsigned l_rev = $3==nullptr ? 0x0 : $4*0x2;
            unsigned u_rev = $6==nullptr ? 0x4 : $7*0x4;
            $$=make_slice_node(
	         make_applied_identifier($1,@1),
                 $3==nullptr ? make_int_denotation(0,@3) : $3,
		 $6==nullptr ? make_int_denotation(0,@6) : $6,
		 l_rev^u_rev,@$);
          }
        | comprim '[' expr_opt tilde_opt ':' expr_opt tilde_opt ']'
	  { unsigned l_rev = $3==nullptr ? 0x0 : $4*0x2;
            unsigned u_rev = $6==nullptr ? 0x4 : $7*0x4;
            $$=make_slice_node( $1,
                 $3==nullptr ? make_int_denotation(0,@3) : $3,
		 $6==nullptr ? make_int_denotation(0,@6) : $6,
		 l_rev^u_rev,@$);
          }
	| IDENT TLSUB expr_opt tilde_opt ':' expr_opt tilde_opt ']'
	  { unsigned l_rev = $3==nullptr ? 0x0 : $4*0x2;
            unsigned u_rev = $6==nullptr ? 0x4 : $7*0x4;
            $$=make_slice_node(
	         make_applied_identifier($1,@1),
                 $3==nullptr ? make_int_denotation(0,@3) : $3,
		 $6==nullptr ? make_int_denotation(0,@6) : $6,
		 0x1^l_rev^u_rev,@$);
          }
	| comprim TLSUB expr_opt tilde_opt ':' expr_opt tilde_opt ']'
	  { unsigned l_rev = $3==nullptr ? 0x0 : $4*0x2;
            unsigned u_rev = $6==nullptr ? 0x4 : $7*0x4;
            $$=make_slice_node( $1,
                 $3==nullptr ? make_int_denotation(0,@3) : $3,
		 $6==nullptr ? make_int_denotation(0,@6) : $6,
		 0x1^l_rev^u_rev,@$);
          }
	| IDENT '[' expr_opt tilde_opt ':' expr_opt tilde_opt
	        ',' expr_opt tilde_opt ':' expr_opt tilde_opt ']'
	  { unsigned r_l_rev =  $3==nullptr ? 0x0  :  $4*0x2;
	    unsigned r_u_rev =  $6==nullptr ? 0x4  :  $7*0x4;
	    unsigned c_l_rev =  $9==nullptr ? 0x0  : $10*0x10;
	    unsigned c_u_rev = $12==nullptr ? 0x20 : $13*0x20;
	    auto arg = raw_expr_list(nullptr);
	    arg = make_exprlist_node
	      ($12==nullptr ? make_int_denotation(0,@12) : $12,arg);
	    arg = make_exprlist_node
	      ($9==nullptr ? make_int_denotation(0,@9) : $9,arg);
	    arg = make_exprlist_node
	      ($6==nullptr ? make_int_denotation(0,@6) : $6,arg);
	    arg = make_exprlist_node
	      ($3==nullptr ? make_int_denotation(0,@3) : $3,arg);
	    arg = make_exprlist_node
	      (make_cast(make_prim_type(5),make_applied_identifier($1,@1),@1),
	       arg);
	    arg = make_exprlist_node
	      (make_int_denotation(r_l_rev^r_u_rev^c_l_rev^c_u_rev,@$),arg);
	    auto arg_tup=wrap_tuple_display(arg,@$);
	    $$ =
	      make_unary_call(lookup_identifier("matrix slicer"),arg_tup,@$,@2);
	  }
	| comprim '[' expr_opt tilde_opt ':' expr_opt tilde_opt
	          ',' expr_opt tilde_opt ':' expr_opt tilde_opt ']'
	  { unsigned r_l_rev =  $3==nullptr ? 0x0  :  $4*0x2;
	    unsigned r_u_rev =  $6==nullptr ? 0x4  :  $7*0x4;
	    unsigned c_l_rev =  $9==nullptr ? 0x0  : $10*0x10;
	    unsigned c_u_rev = $12==nullptr ? 0x20 : $13*0x20;
	    auto arg = raw_expr_list(nullptr);
	    arg = make_exprlist_node
	      ($12==nullptr ? make_int_denotation(0,@12) : $12,arg);
	    arg = make_exprlist_node
	      ($9==nullptr ? make_int_denotation(0,@9) : $9,arg);
	    arg = make_exprlist_node
	      ($6==nullptr ? make_int_denotation(0,@6) : $6,arg);
	    arg = make_exprlist_node
	      ($3==nullptr ? make_int_denotation(0,@3) : $3,arg);
	    arg = make_exprlist_node(make_cast(make_prim_type(5),$1,@1),arg);
	    arg = make_exprlist_node
	      (make_int_denotation(r_l_rev^r_u_rev^c_l_rev^c_u_rev,@$),arg);
	    auto arg_tup=wrap_tuple_display(arg,@$);
	    $$ =
	      make_unary_call(lookup_identifier("matrix slicer"),arg_tup,@$,@2);
	  }
;

iftail	: expr THEN expr ELSE expr FI { $$=make_conditional_node($1,$3,$5,@$); }
	| expr THEN expr ELIF iftail { $$=make_conditional_node($1,$3,$5,@$); }
	| expr THEN expr FI
	  { $$=make_conditional_node($1,$3,wrap_tuple_display(nullptr,@$),@$); }
;

pattern : IDENT		    { $$.kind=0x1; $$.name=$1; }
	| '!' IDENT         { $$.kind=0x5; $$.name=$2; } // IDENT declared const
	| '(' pat_list ')'  { $$.kind=0x2; $$.sublist=reverse_patlist($2); }
	| '(' pat_list ')' ':' IDENT
	    { $$.kind=0x3; $$.name=$5; $$.sublist=reverse_patlist($2);}
	| '(' pat_list ')' ':' '!' IDENT
	    { $$.kind=0x7; $$.name=$6; $$.sublist=reverse_patlist($2);}
	| '(' ')' { $$.kind=0x2; $$.sublist=0; } /* allow throw-away value */
;

pattern_opt :/* empty */ { $$.kind=0x0; }
	| pattern
;

pat_list: pattern_opt ',' pattern_opt
	  { $$=make_pattern_node(make_pattern_node(nullptr,$1),$3); }
	| pat_list ',' pattern_opt { $$=make_pattern_node($1,$3); }
;

id_spec: nostar_type pattern { $$.type_pt=$1; $$.ip=$2; }
        | '(' id_specs ')'
	{ $$.type_pt=make_tuple_type($2.typel);
          $$.ip.kind=0x2; $$.ip.sublist=reverse_patlist($2.patl);
	}
;

id_specs: id_spec
        { $$.typel=make_type_singleton($1.type_pt);
	  $$.patl=make_pattern_node(nullptr,$1.ip);
	}
	| id_specs ',' id_spec
	{ $$.typel=make_type_list($1.typel,$3.type_pt);
	  $$.patl=make_pattern_node($1.patl,$3.ip);
	}
;

id_specs_opt: id_specs
	| /* empty */ { $$.typel=nullptr; $$.patl=nullptr; }
;

nostar_type : TYPE	{ $$=make_prim_type($1); }
	| TYPE_ID
	  { bool c; $$=acquire(global_id_table->type_of($1,c)).release(); }
        | '(' type ')'	{ $$=$2; }
	| '(' type ARROW type ')' { $$=make_function_type($2,$4); }
	| '(' types_opt ARROW type ')'
			{ $$=make_function_type(make_tuple_type($2),$4); }
	| '(' type ARROW types_opt ')'
			{ $$=make_function_type($2,make_tuple_type($4)); }
	| '(' types_opt ARROW types_opt ')'
	   { $$=make_function_type(make_tuple_type($2),make_tuple_type($4)); }
	| '[' type ']'	{ $$=make_row_type($2); }
	| '(' types ')'	{ $$=make_tuple_type($2); }
;

type : nostar_type
	| '*' { $$=new type_expr;}
;

types	: type ',' type
	  { $$=make_type_list(make_type_singleton($1),$3); }
	| types ',' type { $$=make_type_list($1,$3); }
;

types_opt : /* empty */ { $$=nullptr; }
	| types
;

commalist_opt: /* empty */	 { $$=raw_expr_list(nullptr); }
	| commalist
;

commalist: expr  { $$=make_exprlist_node($1,raw_expr_list(nullptr)); }
	| commalist ',' expr { $$=make_exprlist_node($3,$1); }
;

commabarlist: commalist_opt '|' commalist_opt
	{ $$ = make_exprlist_node(wrap_list_display(reverse_expr_list($3),@$)
		,make_exprlist_node(wrap_list_display(reverse_expr_list($1),@$)
				    ,raw_expr_list(nullptr)));
	}
	| commabarlist '|' commalist_opt
	{ $$=make_exprlist_node
	    (wrap_list_display(reverse_expr_list($3),@$),$1); }
;


%%
