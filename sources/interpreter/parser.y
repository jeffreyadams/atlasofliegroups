%{
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
%parse-param { int* running }
%pure-parser
%error-verbose

%token QUIT TRUE FALSE QUIET VERBOSE
%token DIVMOD
%token <val> INT
%token <expression> STRING
%token <id_code> IDENT
%type  <expression> exp
%destructor { destroy_expr ($$); } exp
%type  <expression_list> commalist commabarlist idlist
%destructor { destroy_exprlist($$); } commalist commabarlist
%nonassoc DIVMOD
%left '-' '+'
%left '*' '/' '%'
%left NEG     /* negation--unary minus */

%{
  int yylex (YYSTYPE *, YYLTYPE *);
  void yyerror (YYLTYPE *, expr*, int*,const char *);
%}
%%

input:  '\n'			{ YYABORT } /* return 1 to skip evaluator */
        | exp    '\n'		{ *parsed_expr=$1; }
        | idlist '=' exp '\n'
		{ global_set_identifier(reverse_expr_list($1),$3); YYABORT }
        | QUIT	'\n'		{ *running =-1; } /* causes immediate exit */
        | QUIET	'\n'		{ *running =0; YYABORT } /* quiet mode */
        | VERBOSE '\n'		{ *running =1; YYABORT } /* verbose mode */
;

exp     : INT { $$ = make_int_denotation($1); }
        | TRUE { $$ = make_bool_denotation(1); }
        | FALSE { $$ = make_bool_denotation(0); }
	| STRING
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
        | exp '+' exp
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
        | IDENT { $$=make_applied_identifier($1); }
	| IDENT '(' commalist ')'
		{ $$=make_application_node($1,reverse_expr_list($3)); }
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

