#define version_string  "x3.6"
#define banner  "This is CWEAVE (Version "version_string")"
#define max_refs  10000 
#define max_toks  10000 
#define max_texts  2500 
#define max_scraps  4000 
#define max_no_of_nodes  300 
#define line_length  80 
#define stack_size  400 
#define sort_stack_size  500 
#define variant  xref_info
#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  "common.h"
#define max_bytes  50000L \
   
#define max_modules  1000 
#define max_idents  5000 
#define max_sections  4000 \
   
#define hash_size  353 
#define buf_size  500 
#define longest_name  1000 \
   
#define long_buf_size  (buf_size+longest_name) 
#define local  static
#define array_size(a) ((int)(sizeof(a)/sizeof(a[0])))
#define false  (boolean) 0
#define true  (boolean) 1
#define ctangle  0
#define cweave  1
#define and_and  04 
#define lt_lt  020 
#define gt_gt  021 
#define plus_plus  013 
#define minus_minus  01 
#define minus_gt  031 
#define not_eq  032 
#define lt_eq  034 
#define gt_eq  035 
#define eq_eq  036 
#define or_or  037 
#define find_char() (loc<=limit || get_line())
#define id_index(p) ((sixteen_bits)((p)-id_table))
#define id_at(i)    (&id_table[i])
#define mod_index(p) ((sixteen_bits)((p)-mod_table))
#define mod_at(i)    (&mod_table[i])
#define name_begin(p) ((p)->byte_start)
#define length(p) ((int)(strlen(name_begin(p))))
#define name_end(p) (name_begin(p)+length(p))
#define complete_name(p) ((p)->byte_start[-1]=='\0')
#define print_mod(p) \
   printf(": <%s%s>",name_begin(p), complete_name(p) ? "" : "..." )
#define spotless  0 
#define harmless_message  1 
#define error_message  2 
#define fatal_message  3 
#define mark_harmless() \
    if (history==spotless) history=harmless_message;  else 
#define mark_error() (history=error_message)
#define overflow(t) fatal("\n! Sorry, %s capacity exceeded",t) \
			
#define confusion(s) fatal("\n! This can't happen: %s",s) 
#define show_banner  flags['b'] 
#define show_happiness  flags['h'] 
#define show_progress  flags['p'] 
#define show_stats  flags['s'] 
#define C_plus_plus  flags['+'] 
#define compatibility_mode  flags['c'] 
#define update_terminal() fflush(stdout) 
#define new_line() putchar('\n')
#define term_write(string,leng) printf("%.*s",(int)(leng),string) 
#define reserved(a) (a->ilk>=type_defined) \
  
#define unindexed(a) (a->ilk>=NULL_like) \
  
#define xref  equiv_or_xref
#define next_xref(x) (&xmem[(x)->next])
#define xnum(i) (xmem[i].num)
#define xlink(i) (xmem[i].next)
#define xref_index(p) ((sixteen_bits)((p)-xmem))
#define cite_flag  0x4000 \
       
#define def_flag  0x8000 
#define num_mask  (cite_flag-1) 
#define make_xref(n,i) \
    /* create \xr. node with |num==n| and successor |xmem[i]| */ \
  if (++xref_ptr >= &xmem[max_refs]) \
    overflow ("cross-reference");  \
  else xref_ptr->num=n, xref_ptr->next=i 
#define no_xref  (!flags['x'])
#define make_xrefs  flags['x'] 
#define file_flag  (cite_flag-1) \
		
#define tok_mem_end   (&tok_mem[max_toks]) 
#define text_mem_end   (&text_mem[max_texts]) 
#define text_index(p) ((sixteen_bits)((p)-text_mem)) \
                 
#define text_at(i) (&text_mem[i]) 
#define text_begin(p) (*(p)) \
              
#define text_end(p) (*(p+1)) \
              
#define code_of(c) ccode[(unsigned char)(c)]
#define compress2(char2,code) if (*loc==char2) return ++loc, code 
#define compress3(char2,char3,code) \
  if (*loc==char2 && loc[1]==char3) return loc+=2, code 
#define comp_ass_op2(code) \
  if (*loc=='=' && !compatibility_mode) return ++loc, code 
#define comp_ass_op3(char2,code) \
  if (*loc==char2 && loc[1]=='=' && !compatibility_mode) return loc+=2,code 
#define shift_and_store(ch) (*id_loc++=ch,c=*++loc)
#define shift() (next_control=get_next())
#define typedef_tracking(b) (typedef_master += b ? 5 : -5)
#define emit_space() out_str ("\\Y"); 
#define scrap_info  scrap_union.scrap_field
#define scrap_info_end   (&scrap_info[max_scraps]) 
#define app(a) (*tok_ptr++ = a)
#define app_tok(a)  \
  if (tok_ptr>tok_mem_end-2) overflow("token");  \
   else app(a) 
#define app_char_tok(c) app_tok((unsigned char)(c))
#define freeze_text() (*++text_ptr = tok_ptr)
#define pack_scrap(c,m) \
 ( scrap_ptr->cat = c, scrap_ptr->trans = text_ptr, freeze_text(), \
  (scrap_ptr++)->mathness = 5*(m) )
#define id_flag  10240U 
#define res_flag  (2*id_flag) 
#define mod_flag  (3*id_flag) 
#define text_flag  (4*id_flag) 
#define inner_text_flag  (5*id_flag) 
#define enter_block(i) save[i].txt=text_ptr, save[i].tok=tok_ptr;
#define leave_block(i) text_ptr=save[i].txt, tok_ptr=save[i].tok;
#define check_toks(n)  if (tok_ptr>tok_mem_end-n) \
   overflow("token");   else 
#define check_text()   if (text_ptr>=text_mem_end-1) \
   overflow("text");   else 
#define check_scrap()  if (scrap_ptr>=scrap_info_end) \
   overflow("scrap");   else check_text() 
#define dangling_tokens() (compatibility_mode && tok_ptr>*text_ptr)
#define max_category  end_expr 
#define min_nonrule_cat  lproc 
#define valid_cat(c) ((c)>0 && (c)<min_nonrule_cat)
#define yes_math  1 
#define no_math  2 
#define maybe_math  0 
#define start_scrap(s,c,m) p=&token_trans[s],p->cat=c, p->mathness=5*(m)
#define end_scrap  p->trans=text_ptr, freeze_text();
#define max_lhs_length  4
#define trie_root  (&trie_node_mem[0])
#define successor(q,c) \
   (&trie_node_mem[(q)->next[(c)-1]]) 
#define no_successor(q,c) ((q)->next[(c)-1]==0)
#define set_successor(q,c,x) ((q)->next[(c)-1]=(sixteen_bits)((x)-trie_node_mem))
#define rule_error  install_failed=true,print \
   
#define left_math(a) (a->mathness>>2)
#define right_math(a) (a->mathness&0x3)
#define set_mode(x) \
if (cur_mathness==maybe_math) cur_mathness=init_mathness=x; \
else if (cur_mathness!=x) { app('$'); cur_mathness=x; } \
 \
else  
#define app_trans(a) app_tok(text_flag+text_index((a)->trans))
#define add_trans(a) \
{ scrap_pointer scr=a; /* so that evaluating |a| may have side effects */ \
  if (left_math(scr)!=maybe_math) \
  { if (cur_mathness==maybe_math) init_mathness=left_math(scr); \
    else if (cur_mathness!=left_math(scr)) app('$'); \
 \
    cur_mathness=right_math(scr); \
  } \
  app_trans(scr); \
}
#define math_char(x) ((x)==yes_math ? '+' : (x)==no_math ? '-' : '?')
#define cwebx 		0x0001 
#define compatibility 	0x0002 
#define only_plus_plus 	0x0004 
#define no_plus_plus 		0x0008 
#define unaligned_braces 	0x0050 
#define aligned_braces 	0x0020 
#define wide_braces 		0x0030 
#define standard_braces 	0x0060 
#define merged_decls 		0x0080 
#define forced_statements 	0x0100 
#define no_forced_statements  0x0600 
#define all_stats_forced 	0x0300 
#define not_all_stats_forced 	0x0400 
#define cur_tok  cur_state.tok_field 
#define cur_end  cur_state.end_field 
#define cur_mode  cur_state.mode_field 
#define cur_state  stack[0] 
#define stack_end  (&stack[stack_size]) 
#define pop_level()  cur_state = *--stack_ptr
#define out_id_full(p) out_id_part(name_begin(p),length(p))
#define out_line  (&out_buf[1]) 
#define out_buf_end   (&out_line[line_length]) 
#define tex_putc(c) putc (c, tex_file)
#define tex_new_line() (putc('\n', tex_file),++out_line_nr)
#define tex_printf(format) fprintf(tex_file, format)
#define output_line_empty() (out_ptr==out_line)
#define out(c) \
  *(out_ptr>=out_buf_end ? (break_out(),out_ptr++) : out_ptr++)=c
#define triple_file_output  flags['t']
#define even_out_pages  flags['e']
#define sort_info  scrap_union.id_list_field
#define sort_info_end  (&sort_info[sort_stack_size])
#define ilink(p) index_link[id_index(p)]
#define infinity  255 
#define report(k,c,m) \
  printf("\t%lu %ss (out of %lu)\n",(unsigned long)(c),k,(unsigned long)(m))


boolean names_match (id_pointer,char*,int,int);
void init_id_name (id_pointer,int);
void init_module_name (mod_pointer);


enum 
{ normal, /* |ilk| of ordinary identifiers */
  roman, /* |ilk| of roman type index entries */
  wildcard, /* |ilk| of user-formatted index entries */
  typewriter, /* |ilk| of typewriter type entries */
  reference, /* |ilk| of identifiers used for explicit \xr.s */
  header_file_name, /* |ilk| of file names seen as included header file */
  type_defined, /* |ilk| of identifiers that are defined by |typedef| */
  TeX_like, NULL_like,
    /* |ilk| of identifiers with user-given control sequences */
  const_like, typedef_like, and_like, not_like, namespace_like, typename_like
     /* special reserved words */
};

typedef struct xref_info
{ sixteen_bits num; /* section number plus a multiple of |cite_flag| */
  sixteen_bits next; /* index of the next \xr. in the list */
} xref_info, *xref_pointer;

typedef sixteen_bits token, * token_pointer, ** text_pointer;

enum 
{ at_sign_image = UCHAR_MAX+1, /* quoted `\.{@}' */
  or, /* \:v */
  mul_assign, div_assign, mod_assign, plus_assign, minus_assign,
  left_assign, right_assign, and_assign, xor_assign, or_assign,
  sh_sh, ellipsis, colon_colon, 
  start_preproc, end_preproc, /* begin and end of a preprocessor directive */
  join, /* \:\& */
  thin_space, /* \:, */
  math_break, /* \:\v */
  line_break, /* \:/ */
  big_line_break, /* \:) */
  no_line_break, /* \:+ */
  backup_line, /* \:\\ */
  pseudo_semi, /* \:; */
  force_expr_open, force_expr_close, /* \:[, \ \:] */
  include_preproc, /* \:p */

  ignore, /* control code of no interest to \.{\me.} */
  constant, string, /* the next five codes should remain in this order */
  identifier, /* any (possibly reserved) word found in \Cee\ text */
  xref_roman, xref_wildcard, xref_typewriter, xref_mark,
    /* \:\^, \ \:?, \ \:., \ \:: */
  refer, /* \:\# */
  TeX_string, /* \:t */
  verbatim, /* \:= */
  ignored_text, /* \:q */
  char_trans, /* \:l */
  ASCII_code, /* \:' */
  begin_comment, end_comment,
  underline, /* \:! */
#ifdef DEBUG
  trace0, trace1, trace2, trace3, /* \:0, \dots, \:3 */
#endif
  format, /* \:f */
  definition, /* \:d */
  header, /* \:h */
  begin_C, /* \:c */
  module_name, /* \:< and \:( */
  new_section /* \:\ , \:\~ and \:* */
};

typedef struct
{ eight_bits cat; /* category code */
  eight_bits mathness; /* whether in math mode at left and right boundary */
  text_pointer trans; /* translation text */
} scrap, *scrap_pointer;

enum 
{ cancel=UCHAR_MAX+1,/* the following 9 items should remain in this order */
  indent, outdent, opt, flush_left, break_space, force, big_force,
  backup, big_backup, 
  relax, 
  space=opt, tilde=flush_left 
};

enum  
{ unop = 1, /* a unary operator like `|!|' */
  binop, /* a binary operator like `|<|' */
  unorbinop, /* an operator that can be either, like `|-|' */
  select, /* structure selection: `|.|' or `|->|' */
  question, /* a question mark operator */
  lbrace, rbrace, lpar, rpar, lbrack, rbrack,
    /* `|{|', \ `|}|', `(', \ `)', `[', \ `]' */
  comma, semi, colon, colcol, magic, /* `,', `;', `:', `$\CC$', \:; */
  subscript, /* an array subscript, like `|[]|' or `|[i++]|' */
  struct_head, /* the beginning of a struct specifier, like `|struct s{|' */
  short_lbrace, short_struct_head, /* from `\.{\{@;}', for one-liners */
  compound_statement, /* a complete compound statement */
  statement, /* a complete statement, possibly compound */
  function, /* a complete function definition */
  function_head, /* a function identifier followed by formal parameters */
  parameters,
  /* parameters in function declaration, or casting operator like `|(int)|' */
  label, /* a statement label */
  if_head, /* `|if|' followed by a (parenthesised) expression */
  if_else_head, /* \lq|if @t\dots@>@; else|', \lq|while(@t\dots@>)|'
                or \lq|switch(@t\dots@>)|' */
  do_head, /* `|do @t\dots@>@; while|' */
  mod_scrap, /* module name */
  declarator, /* abstract declarator, like `|(*)(int,char*[])|' */
  declaration, /* a complete declaration */
  expression, /* an expression, possibly a single identifier */
  while_like, /* `|for|', `|while|', `|switch|' */
  do_like,
  if_like,
  else_like,
  extern_like,
  throw_like, try_like, catch_like, /* single tokens */
  int_like, /* `|int|', `|char|', \dots  */
  case_like, /* `|case|', `|default|' */
  sizeof_like, /* `|sizeof|', `|const|', `\&{new}', `\&{delete} */
  struct_like, /* `|struct|', `|union|', `|enum|',
                  `\&{class}', `\&{typename}' */
  return_like, /* `|return|', `|break|', `|continue|', `|goto|' */
  template_like, langle, rangle,
              /* `\&{template}', `$\langle$', `$\rangle$', for \Cpp */
  templ_params, /* template parameters */
  lproc, /* `\&\#' and following identifier starting preprocessor directive */
  rproc, /* end of a preprocessor directive */
  insert, /* comment or other syntactically inert item */
  begin_expr, end_expr /* \:[ and \:] */
};

typedef struct
{ short id; /* for debugging */
  struct
    { eight_bits category[max_lhs_length]; signed char context,length; } lhs;
  struct { eight_bits category; char *translation; } rhs;
  sixteen_bits mask;
  short displacement;
} reduction;

typedef struct { reduction* rule; sixteen_bits next[min_nonrule_cat-1]; }
  trie_node;

typedef enum { inner, outer } mode;
typedef struct
{ token_pointer tok_field; /* present location within token list */
  token_pointer end_field; /* ending location of token list */
  mode mode_field; /* interpretation of control tokens */
} output_stack_element, *stack_pointer;

typedef struct { id_pointer head; int depth; } sort_node, * sort_pointer;


void phase_one (void);
  /* read all the user's text and store the \xr.s */
void phase_two (void);
  /* read all the text again and translate it to \TeX\ form */
void phase_three (void); 

int scan_comment(int* bal, boolean one_liner);
   

 int get_next (void);

boolean push_header_file(boolean suspend);

void C_xref (boolean); /* make \xr.s within in straight \Cee~text */
void outer_xref (void); /* make \xr.s in \Cee~text with comments */
void mod_check (mod_pointer); 

void do_C (void); /* handle \Cee~text enclosed in `\pb' */
void outer_read (void); /* transform input into scraps */
void finish_C (void); /* finishes a definition or a \Cee~part */
void finish_line(void); /* send out a line of output */
void out_str (char*); /* write multiple characters */
void out_sec_nr (int); /* output a section number */
xref_pointer list_refs (xref_pointer,sixteen_bits);
  /* output module \xr.s */
void footnote(xref_pointer*,sixteen_bits); /* same with heading text */
void app_str(char*); 

text_pointer translate(void);
  /* build formatted text from collected scraps */
void make_output(text_pointer,mode);
  

void install_rule (reduction* rule);

void out_identifier (id_pointer);
void out_keyword (id_pointer);
xref_pointer out_module_name(mod_pointer);

void break_out (void);

void unbucket (eight_bits);

void list_modules (mod_pointer);


xref_info xmem[max_refs]; /* contains \xr. information */
xref_pointer xref_ptr = &xmem[0]; /* the last used position in |xmem| */
sixteen_bits xref_switch = 0, mod_xref_switch = 0;
	

token tok_mem[max_toks]; /* tokens */
token_pointer text_mem[max_texts]; /* directory into |tok_mem| */
token_pointer tok_ptr = tok_mem; /* first unused position in |tok_mem| */
text_pointer text_ptr = text_mem; /* first unused position in |text_mem| */
#ifdef STAT
token_pointer max_tok_ptr = tok_mem; /* largest value of |tok_ptr| */
text_pointer max_text_ptr = text_mem; /* largest value of |text_ptr| */
#endif

int ccode[UCHAR_MAX + 1]; 

id_pointer cur_id; /* identifier or index entry just scanned */
mod_pointer cur_mod; /* module name just scanned */
int preprocessing=0;

boolean change_exists=false; /* has any section changed? */
int next_control; 

local int typedef_master=-5; /* tracking disabled outside \Cee~parts */
local int brace_level, par_level;

union
{ scrap scrap_field[max_scraps]; /* memory array for scraps */
  
  sort_node id_list_field[sort_stack_size];

} scrap_union;
scrap_pointer scrap_base=scrap_info;
	/* beginning of the current scrap sequence */
scrap_pointer scrap_ptr = scrap_info;
	/* points to end of the current scrap sequence */
#ifdef STAT
scrap_pointer max_scr_ptr = scrap_info;
	/* largest value assumed by |scrap_ptr| */
#endif

struct { text_pointer txt; token_pointer tok; } save[2];

scrap token_trans[ignore];

trie_node trie_node_mem[max_no_of_nodes];
int node_no = 1;  /* number of trie nodes allocated */
#ifdef DEBUG
boolean install_failed=false;
#endif

 sixteen_bits rule_mask;

scrap_pointer pp; /* current position for reducing scraps */
scrap_pointer lo_ptr; /* end of sequence of scraps that have been examined */
scrap_pointer hi_ptr; 

#ifdef DEBUG
int tracing; /* how much parsing details to show */
#endif

output_stack_element stack[stack_size]; /* info for non-current levels */
stack_pointer stack_ptr=&stack[0];
  /* first unused location in the output state stack */
#ifdef STAT
stack_pointer max_stack_ptr = stack; /* largest value assumed by |stack_ptr| */
#endif

char out_buf[line_length+1]; /* assembled characters */
char* out_ptr; /* first unused position in |out_buf| */
int out_line_nr=1; 

sort_pointer sort_ptr=sort_info;
#ifdef STAT
sort_pointer max_sort_ptr=sort_info; /* largest value of |sort_ptr| */
#endif

eight_bits collate[UCHAR_MAX-25]; /* collation order */
int end_collate; 

id_pointer index_link[max_idents]; /* links identifiers during sorting */
id_pointer bucket[UCHAR_MAX+1];


int main (int argc,char** argv)
{ program=cweave;
  make_xrefs=true;
  common_init(argc,argv,banner);
  
  xnum(0)=0;  
  
  *text_ptr=tok_ptr;
  
  { unsigned char c=0;
    do ccode[c] = isspace(c) ? new_section : ignore; while(c++!=UCHAR_MAX);
    ccode['@'] = at_sign_image;
  ccode['v'] = ccode['V'] = or;
  ccode['!'] = underline; /* set definition flag */
  ccode['^'] = xref_roman; /* index entry to be typeset normally */
  ccode['?'] = xref_wildcard; /* index entry to be in user format */
  ccode['.'] = xref_typewriter; /* index entry to be in typewriter type */
  ccode[':'] = xref_mark; ccode['#']=refer; /* explicit \xr.s */
  ccode['t'] = ccode['T'] = TeX_string; /* \TeX\ box within \Cee\ text */
  ccode['='] = verbatim;
  ccode['q'] = ccode['Q'] = ignored_text;
  ccode['l'] = ccode['L'] = char_trans;
  ccode['\''] = ASCII_code;
  ccode['&'] = join; /* concatenate two tokens */
  ccode[','] = thin_space;
  ccode['|'] = math_break;
  ccode['/'] = line_break;
  ccode[')'] = big_line_break;
  ccode['\\']= backup_line;
  ccode['+'] = no_line_break;
  ccode[';'] = pseudo_semi;
  ccode['['] = force_expr_open; ccode[']'] = force_expr_close;
  ccode['p'] = ccode['P'] = include_preproc;
  #ifdef DEBUG
    ccode['0'] = trace0; ccode['1'] = trace1;
    ccode['2'] = trace2; ccode['3'] = trace3;
  #endif
  ccode['f'] = ccode['F'] = ccode['s'] = ccode['S'] = format;
  ccode['d'] = ccode['D'] = definition;
  ccode['h'] = ccode['H'] = header;
  ccode['c'] = ccode['C'] = begin_C;
        /* \Cee\ text in unnamed module */
  ccode['<'] = ccode['('] = module_name; /* beginning of a module name */
  ccode['~'] = ccode['*'] = new_section; /* beginning of a new section */
    if (compatibility_mode)
      
      { ccode['h']=ccode['H']=include_preproc; /* \:h means \:p */
        ccode['p']=ccode['P']=begin_C; /* \:p means \:c */
        ccode['#']=big_line_break; /* \:\# means \:) */
        ccode[':']=xref_wildcard; /* \:: means \:? */
      }
  }
  
  { static struct { short tok; eight_bits cat, mathness; char* tr; }
      trans_ini [] = { 
                      { '!', unop,	yes_math, "\\R" },  
                      { '~', unop,	yes_math, "\\CM" },  
                      { '/', binop,	yes_math, "/" }, 
                      { '<', binop,	yes_math, "<" }, 
                      { '>', binop,	yes_math, ">" }, 
                      { '.', select,	yes_math, "." }, 
                      { '=', binop,	yes_math, "\\K" },  
                      { '|', binop,	yes_math, "\\OR" }, 
                         { or, binop,	yes_math, "\\OR" }, 
                      { '^', binop,	yes_math, "\\XOR" },  
                      { '%', binop,	yes_math, "\\MOD" },  
                      { '+', unorbinop,	yes_math, "+" }, 
                      { '-', unorbinop,	yes_math, "-" }, 
                      { '*', unorbinop,	yes_math, "*" }, 
                      { '&', unorbinop,	yes_math, "\\AND" },  
                      { '?', question,	yes_math, "\\?" }, 
                      { '(', lpar,	yes_math, "(" }, 
                      { ')', rpar,	yes_math, ")" }, 
                      { '[', lbrack,	maybe_math, "[" }, 
                      { ']', rbrack,	maybe_math, "]" }, 
                      { '{', lbrace,	yes_math, "\\{" }, 
                      { '}', rbrace,	yes_math, "\\}" }, 
                      { ',', comma,	yes_math, "," }, 
                      { ';', semi,	yes_math, ";" }, 
                      { ':', colon,	maybe_math, ":" }, 
                      { '#', insert,	maybe_math, "\\#" },
                          /* this should occur only in macro definitions */  
                      { at_sign_image, insert, maybe_math, "@" },
                          /* this should not occur in legal \Cee~text */
                      
                      { not_eq,	binop,	yes_math, "\\I" },  
                      { lt_eq,	binop,	yes_math, "\\Z" },  
                      { gt_eq,	binop,	yes_math, "\\G" },  
                      { eq_eq,	binop,	yes_math, "\\E" },  
                      { and_and,	binop,	yes_math, "\\W" },  
                      { or_or,	binop,	yes_math, "\\V" },  
                      { plus_plus,	unop,	yes_math, "\\PP" },  
                      { minus_minus,	unop, 	yes_math, "\\MM" },  
                      { minus_gt,	select,	yes_math, "\\MG" },  
                      { gt_gt,	binop,	yes_math, "\\GG" },  
                      { lt_lt,	binop,	yes_math, "\\LL" },  
                      { mul_assign,	binop,	yes_math, "\\KK*" },  
                      { div_assign,	binop,	yes_math, "\\KK/" }, 
                      { mod_assign,	binop,	yes_math, "\\KK\\MOD" },  
                      { plus_assign,	binop,	yes_math, "\\KK+" }, 
                      { minus_assign,	binop,	yes_math, "\\KK-" }, 
                      { left_assign,	binop,	yes_math, "\\KK\\LL" }, 
                      { right_assign,	binop,	yes_math, "\\KK\\GG" }, 
                      { and_assign,	binop,	yes_math, "\\KK\\AND" },  
                      { xor_assign,	binop,	yes_math, "\\KK\\XOR" },  
                      { or_assign,	binop,	yes_math, "\\KK\\OR" },  
                      { thin_space,	insert, yes_math, "\\," },  
                      { pseudo_semi,	magic,	maybe_math, "" }, 
                      { force_expr_open,  begin_expr,	  maybe_math, "" }, 
                      { force_expr_close, end_expr,	  maybe_math, "" }, 
                      { join,		insert,	no_math, "\\J" },  
                      { ellipsis,   int_like,	yes_math, "\\ldots" }, 
                      { sh_sh,	binop,	yes_math, "\\SS" },  
                      { colon_colon,	colcol, yes_math, "\\CC" }
   };
    int i,n=array_size(trans_ini);
    scrap o={ insert, 5*maybe_math, NULL }; /* completely inert scrap */
  
    for (i=0; i<n; ++i)
    { scrap* p=&token_trans[trans_ini[i].tok];
      p->cat=trans_ini[i].cat; p->mathness=5*trans_ini[i].mathness;
      app_str(trans_ini[i].tr); p->trans=text_ptr; freeze_text();
    }
    if (C_plus_plus) 
                     { token_trans['<'].cat=langle; token_trans['>'].cat=rangle; }
    
    { scrap* p;
      start_scrap(math_break,insert,maybe_math); /* \:\v */
      app(opt), app('0'); end_scrap;
    start_scrap(line_break,insert,no_math); /* \:/ */
      app(force); end_scrap;
    start_scrap(end_preproc,rproc,no_math); /* end of preprocessor directive */
      app(force); end_scrap;
    start_scrap(' ',insert,no_math); /* space within preprocessor directive */
      app(break_space); end_scrap;
    start_scrap(big_line_break,insert,no_math); /* \:) */
      app(big_force); end_scrap;
    start_scrap(backup_line,insert,no_math); /* \:\\ */
      app(backup); end_scrap;
    start_scrap(no_line_break,insert,no_math); /* \:+ */
      app(cancel),app(relax),app(break_space),app(relax),app(cancel); end_scrap;
    start_scrap(include_preproc,insert,yes_math); /* \:p */
      app(force),app_str("\\ATP"),app(force); end_scrap; 
    }
    o.trans=text_ptr; freeze_text(); /* empty translation */
    for (i=0; i<=UCHAR_MAX; ++i) /* clear all remaining tokens */
      if (token_trans[i].cat==0) token_trans[i]=o;
    enter_block(0); /* fix tokens; will be restored after each section */
  }
  
  #ifdef DEBUG
  tracing = flags['d'] ? trace1 : trace0;
  #endif
  
  rule_mask= (compatibility_mode ? 0x0001 : 0x0002)
  	 | (C_plus_plus ? 0x0008 : 0x0004)
  	 | (flags['w'] ? 0x0040 : flags['u'] ? 0x0020 : 0x0010)
  	 | (flags['m'] ? 0x0000 : 0x0080)
  	 | (flags['a'] ? 0x0400 : flags['f'] ? 0x0200 : 0x0100)
  	 ;
  { static reduction rule[] = { 
                               { 1, {{expression, unop}},			{expression, NULL}},	
                               { 2, {{expression, binop, expression}},		{expression, NULL}},	
                               { 2, {{expression, langle, expression}},
                               				{expression, NULL},only_plus_plus},	
                               { 2, {{expression, rangle, expression}},
                               				{expression, NULL},only_plus_plus},	
                               { 3, {{expression, unorbinop, expression}},	{expression, NULL}},	
                               { 4, {{expression, select, expression}},	{expression, NULL}},	
                               { 5, {{expression, select, int_like}},		{expression, "__$_"}},	
                               
                               { 6, {{expression, comma, expression}},		{expression, "__p3_"}},	
                               { 7, {{expression, expression}},		{expression, NULL}},	
                               { 8, {{expression, lpar, rpar}},		{expression, "__,_"}},	
                               { 9, {{expression, subscript}},			{expression, NULL}},	
                               {10, {{lpar, expression, rpar}},		{expression, NULL}},	
                               {11, {{lbrack, expression, rbrack}},		{subscript,  NULL}},	
                               {12, {{lbrack, rbrack}},			{subscript,  "_,_"}},	
                               {13, {{unop, expression}},			{expression, NULL}},	
                               {14, {{unorbinop, expression}},			{expression, "o__"}},	
                               {15, {{unop, int_like}},	      {int_like, NULL},only_plus_plus}, 
                               
                               {20, {{question, expression, colon}},		{binop, "__m_"}},	
                               {21, {{parameters, expression}},		{expression, "_,_"}},	
                               {22, {{sizeof_like, parameters}},		{expression, NULL}},	
                               {23, {{sizeof_like, expression}},		{expression,"_~_"}},	
                               {24, {{sizeof_like, int_like}},	                {expression,"_~_"}
                               						      ,only_plus_plus},	
                               {25, {{int_like, lpar, expression, rpar}},	{expression, NULL}
                               						      ,only_plus_plus},	
                               {26, {{sizeof_like, templ_params}}, {int_like,NULL}   ,only_plus_plus}, 
                               
                               {30, {{int_like, int_like}},		   {int_like, "_~_"}},	
                               {31, {{int_like, expression, semi}},	   {declaration, "_~!__"}},	
                               {32, {{int_like, semi}},		   {declaration, NULL}},	
                               {33, {{int_like, expression, comma}},	   {int_like, "_~!__p8"}},	
                               {34, {{unorbinop, int_like}},		   {int_like, "o__"}},		
                               {35, {{sizeof_like,int_like, int_like},1},
                               				 {int_like, "_~_"},only_plus_plus},	
                               
                               {40, {{magic, int_like}},		   {expression, "_$_"}},	
                               
                               {41, {{expression, parameters}},	   {function_head, "_B_"}},	
                               {42, {{lpar, int_like, expression, rpar}}, {parameters, "_+++_~!_---_"}},
                               {43, {{colon, expression, semi}},          {semi, "m___"}},		
                               
                               {50, {{unorbinop, rpar}, -1},			{declarator, "_,"}},	
                               {51, {{unorbinop, comma},-1},			{declarator, "_,"}},	
                               {52, {{unorbinop, rangle},-1},	{declarator, "_,"},only_plus_plus},	
                               {53, {{int_like, subscript},1},			{declarator, ",_"}},	
                               {54, {{unorbinop, subscript},1},		{declarator, ",_"}},	
                               {55, {{lpar, subscript},1},			{declarator, ",_"}},	
                               
                               {60, {{unorbinop, declarator}},  {declarator, "o__"}},		
                               {61, {{declarator, subscript}},  {declarator, NULL}},		
                               {62, {{declarator, parameters}}, {declarator, NULL}},		
                               {63, {{lpar, declarator, rpar}}, {declarator, NULL}},		
                               
                               {70, {{lpar, int_like, declarator, comma}}, {lpar, "____p5"}},		
                               {71, {{lpar, int_like, comma}},		    {lpar, "___p5"}},		
                               {72, {{lpar, int_like, declarator, rpar}}, {parameters, "_+++__---_"}},	
                               {73, {{lpar, int_like, rpar}},		   {parameters, "_+++_---_"}},	
                               {74, {{lpar, rpar}},			    {parameters, "_,_"}},	
                               
                               {75, {{langle, int_like, declarator, comma}},
                               				{langle, "____p5"},only_plus_plus},	
                               {76, {{langle, int_like, comma}},{langle, "___p5"},only_plus_plus},	
                               {77, {{langle, int_like, declarator, rangle}},
                               			{templ_params, "a___a_"},only_plus_plus},	
                               {77, {{langle, int_like, rangle}},
                               			{templ_params, "a__a_"},only_plus_plus},	
                               {78, {{langle, rangle}}, {templ_params, "a_a_"},only_plus_plus},	
                               {79, {{langle, expression, rangle}},
                               			{templ_params, "a__a_"},only_plus_plus},	
                               {80, {{expression, langle, expression, rangle},1},
                               			{templ_params, "a__a_"},only_plus_plus},	
                               {81, {{int_like,templ_params}},{int_like, NULL},only_plus_plus},	
                               {82, {{expression,templ_params}},{expression, NULL},only_plus_plus},	
                               
                               {90, {{struct_like, lbrace}}, {struct_head, "_ft_"},standard_braces},	
                               {90, {{struct_like, lbrace}}, {struct_head, "_~_"},unaligned_braces},	
                               {90, {{struct_like, lbrace}}, {struct_head, "_f_"},wide_braces},	
                               {91, {{struct_like, expression, lbrace}},
                               		{struct_head, "_~!_ft_"},standard_braces},		
                               {91, {{struct_like, expression, lbrace}},
                               		{struct_head, "_~!_~_"},unaligned_braces},		
                               {91, {{struct_like, expression, lbrace}},
                               		{struct_head, "_~!_f_"},wide_braces},			
                               {92, {{struct_like, int_like, lbrace}},
                               	{struct_head, "_~!$_ft_"},standard_braces|no_plus_plus},	
                               
                               {92, {{struct_like, int_like, lbrace}},
                               	{struct_head, "_~!_ft_"},standard_braces|only_plus_plus},	
                               {92, {{struct_like, int_like, lbrace}},
                               	{struct_head, "_~!$_~_"},unaligned_braces|no_plus_plus},	
                               
                               {92, {{struct_like, int_like, lbrace}},
                               	{struct_head, "_~!_~_"},unaligned_braces|only_plus_plus},	
                               {92, {{struct_like, int_like, lbrace}},
                               	{struct_head, "_~!$_f_"},wide_braces|no_plus_plus},		
                               
                               {92, {{struct_like, int_like, lbrace}},
                               	{struct_head, "_~!_f_"},wide_braces|only_plus_plus},		
                               {93, {{struct_like, expression}}, {int_like, "_~_"}},			
                               {94, {{struct_like, int_like}}, {int_like, "_~$_"},no_plus_plus},	
                               
                               {94, {{struct_like, int_like, semi}},
                               				{declaration, "_~!__"},only_plus_plus},	
                               {95, {{struct_head, declaration, rbrace}},
                               		{int_like, "_+_-f_"},standard_braces},			
                               {95, {{struct_head, declaration, rbrace}},
                               		{int_like, "_+f_-f_"},unaligned_braces & wide_braces},	
                               {96, {{struct_head, function, rbrace}},
                               		{int_like, "_+_-f_"},standard_braces|only_plus_plus},	
                               {96, {{struct_head, function, rbrace}},
                                  {int_like, "_+f_-f_"},(unaligned_braces&wide_braces)|only_plus_plus},
                               {97, {{struct_head, rbrace}}, {int_like,"_B_"}},			
                               {98, {{label, declaration}}, {declaration, "b_f_"},only_plus_plus},     
                               
                               {100, {{extern_like, expression, lbrace},-2},
                               				{struct_like, NULL},only_plus_plus},	
                               {101, {{extern_like, expression}},
                               				{int_like, "_~_"},only_plus_plus},	
                               {102, {{extern_like}},		{int_like, NULL},only_plus_plus},	
                               
                               {105, {{struct_like, lbrace, expression},-1}, {struct_head, "_B_"}},	
                               {106, {{struct_like, expression, lbrace, expression},-1},
                               		{struct_head, "_~_B_"}},				
                               {107, {{struct_head, expression, comma, expression},1},
                               		{expression, "__B!_"}},					
                               {108, {{struct_head, expression, rbrace}}, {int_like, "_~+!_-B_"}},	
                               {109, {{struct_head, expression, comma, rbrace}},
                               					{int_like, "_~+!__-B_"}},	
                               
                               {110, {{struct_like, lbrace, magic}}, {short_struct_head, "_B__+"}},	
                               {111, {{struct_like, expression, lbrace, magic}},
                               		{short_struct_head, "_~!_B__+"}},			
                               {112, {{struct_like, int_like, lbrace, magic}},
                               		{short_struct_head, "_~!$_B__+"}, no_plus_plus},	
                               
                               {112, {{struct_like, int_like, lbrace, magic}},
                               		{short_struct_head, "_~!_B__+"}, only_plus_plus},	
                               {113, {{short_struct_head, declaration}}, {short_struct_head, "_B_"}},	
                               {114, {{short_struct_head, rbrace}}, {int_like, "_-B_"}},		
                               
                               {120, {{expression, semi}},		{statement, NULL}},	
                               {121, {{semi}},				{statement, NULL}},	
                               {122, {{expression, colon}},		{label, "!_h_"}},	
                               {123, {{case_like, expression, colon}},	{label, "_ _h_"}},	
                               {124, {{case_like, colon}},		{label, "_h_"}},	
                               {125, {{label, label}},			{label, "_B_"}},	
                               {126, {{label, statement}}, {statement, "b_B_"},not_all_stats_forced},	
                               {126, {{label, statement}}, {statement, "b_f_"},all_stats_forced},	
                               {127, {{magic}},				{semi, NULL}},		
                               
                               {130, {{declaration, declaration}}, {declaration, "_f_"}},		
                               {131, {{lbrace, declaration, declaration},1},
                               				  {declaration, "_B_"},merged_decls},	
                               {132, {{declaration, statement}},  {statement, "_F_"},no_plus_plus},	
                               {132, {{declaration, statement}},
                                        {statement, "_f_"},only_plus_plus|forced_statements},  	
                               {132, {{declaration, statement}},
                                        {statement, "_B_"},only_plus_plus|no_forced_statements},	
                               {133, {{statement, statement}},{statement, "_f_"},forced_statements},	
                               {133, {{statement, statement}},{statement, "_B_"},no_forced_statements},
                               {134, {{statement, declaration}}, {declaration, "_f_"},only_plus_plus},	
                               {135, {{lbrace, rbrace}},	  {statement, "_,_"}},			
                               {136, {{lbrace, statement, rbrace}},
                               	{compound_statement, "ft_+_-f_"},standard_braces},		
                               {136, {{lbrace, statement, rbrace}},
                               	{compound_statement, "_+f_-f_"},unaligned_braces},		
                               {136, {{lbrace, statement, rbrace}},
                               	{compound_statement, "f_+f_-f_"},wide_braces},			
                               {137, {{lbrace, declaration, rbrace}},
                               	{compound_statement, "ft_+_-f_"},standard_braces},		
                               {137, {{lbrace, declaration, rbrace}},
                               	{compound_statement, "_+f_-f_"},unaligned_braces},		
                               {137, {{lbrace, declaration, rbrace}},
                               	{compound_statement, "f_+f_-f_"},wide_braces},			
                               {138, {{compound_statement}},			{statement, "f_f"}},	
                               {139, {{lbrace, expression, comma, rbrace}},	{expression, "_,__,_"}},
                               {140, {{lbrace, expression, rbrace}},		{expression, "_,_,_"}},	
                               
                               {150, {{lbrace, magic}},		{short_lbrace, "__+"}},	
                               {151, {{short_lbrace, declaration}},	{short_lbrace, "_B_"}},	
                               {152, {{short_lbrace, statement}},	{short_lbrace, "_B_"}},	
                               {153, {{short_lbrace, rbrace}},		{statement, "_-B_"}},	
                               
                               {160, {{if_like, expression}},	      {if_head, "f_~_"}},		
                               {161, {{lbrace,if_like,expression},1}, {if_head, "_~_"},standard_braces},
                               {162, {{if_head, compound_statement, else_like, if_like}},
                               			{if_like, "__f_~_"},aligned_braces},		
                               {162, {{if_head, compound_statement, else_like, if_like}},
                               			{if_like, "_~_~_~_"},unaligned_braces},		
                               {163, {{if_head, statement, else_like, if_like}},
                               			{if_like, "_+B_-f_~_"},not_all_stats_forced},	
                               {163, {{if_head, statement, else_like, if_like}},
                               			{if_like, "_+f_-f_~_"},all_stats_forced},	
                               {164, {{if_head, compound_statement, else_like}},
                               			{if_else_head, "__f_"},aligned_braces},		
                               {164, {{if_head, compound_statement, else_like}},
                               			{if_else_head, "_~_~_"},unaligned_braces},	
                               {165, {{if_head, statement, else_like}},
                               			{if_else_head, "_+B_-f_"},not_all_stats_forced},
                               {165, {{if_head, statement, else_like}},
                               			{if_else_head, "_+f_-f_"},all_stats_forced},	
                               {166, {{if_head, compound_statement}},
                               			{statement, "__f"},aligned_braces},		
                               {166, {{if_head, compound_statement}},
                               			{statement, "_~_f"},unaligned_braces},		
                               {167, {{if_head, statement}},
                               			{statement, "_+B_-f"},not_all_stats_forced},	
                               {167, {{if_head, statement}},
                               			{statement, "_+f_-f"},all_stats_forced},	
                               {168, {{if_else_head, compound_statement}},
                               			{statement, "__f"},aligned_braces},		
                               {168, {{if_else_head, compound_statement}},
                               			{statement, "_~_f"},unaligned_braces},		
                               {169, {{if_else_head, statement}},
                               			{statement, "_+B_-f"},not_all_stats_forced},	
                               {169, {{if_else_head, statement}},
                               			{statement, "_+f_-f"},all_stats_forced},	
                               
                               {170, {{short_lbrace, if_like, expression},1},
                               		{if_head, "_~_"}},				
                               {171, {{short_lbrace, if_head, statement, else_like}},
                               		{short_lbrace, "_B_B_B_"}},			
                               {172, {{short_lbrace, if_head, statement}},
                               		{short_lbrace, "_B_B_"}},			
                               
                               {180, {{while_like, expression}}, {if_else_head, "f_~_"}},		
                               {181, {{lbrace, while_like, expression},1},
                               			{if_else_head, "_~_"},standard_braces},		
                               {182, {{lpar, statement, statement}, 1},
                               			{statement, "_B_"}, forced_statements},		
                               {183, {{lpar, statement, expression, rpar}},	{expression, "__B__"}},	
                               {184, {{lpar, statement, rpar}},	{expression, NULL}},		
                               {185, {{lpar, declaration, statement}, 1},
                               			{statement, "_B_"}, only_plus_plus},		
                               {186, {{do_like, compound_statement, while_like}},
                               			{do_head, "__~_"},standard_braces},		
                               {186, {{do_like, compound_statement, while_like}},
                               			{do_head, "_~_~_"},unaligned_braces},		
                               {186, {{do_like, compound_statement, while_like}},
                               			{do_head, "__f_"},wide_braces},			
                               {187, {{do_like, statement, while_like}},
                               			{do_head, "_+B_-B_"},not_all_stats_forced},	
                               {187, {{do_like, statement, while_like}},
                               			{do_head, "_+f_-f_"},all_stats_forced},		
                               {188, {{do_head, expression, semi}}, {statement, "f_~__f"}},		
                               {189, {{lbrace, do_head, expression, semi},1}, {statement, "_~__f"}},	
                               
                               {200, {{short_lbrace, while_like, expression}},
                               					{short_lbrace, "_B_~_"}},	
                               {201, {{short_lbrace, do_like, statement, while_like},1},
                               					{do_head, "_B_B_"}},		
                               {202, {{short_lbrace, do_head, expression, semi}},
                               					{short_lbrace, "_B_~__"}},	
                               
                               {210, {{return_like, semi}},	  {statement, NULL}},		
                               {211, {{return_like, expression}}, {expression, "_~_"}},	
                               
                               {220, {{function_head, statement}}, {function, "!_f_"}},		
                               {221, {{expression, statement}}, {function, "!_f_"},no_plus_plus},	
                               {221, {{int_like,expression, statement}},
                                                                {function, "!_ _f_"},only_plus_plus},	
                               {222, {{expression, declaration, statement}},
                               			       {function, "!_++f_--f_"},no_plus_plus},	
                               {223, {{int_like, function}},	 {function, "_ _"},no_plus_plus},	
                               {224, {{declaration, function}}, {function, "_F_"},no_plus_plus},	
                               {225, {{function, declaration}}, {declaration, "_F_"},no_plus_plus},	
                               {226, {{function, function}},	 {function, "_F_"},no_plus_plus},	
                               {224, {{function}},		 {declaration, NULL},only_plus_plus},	
                               {227, {{function_head, semi},-1},  {expression, NULL},no_plus_plus},	
                               {228, {{function_head, comma},-1}, {expression, NULL}},			
                               {229, {{function_head, binop},-1}, {expression, NULL}},			
                               {229, {{function_head, unorbinop},-1}, {expression, NULL}},		
                               {229, {{function_head, langle},-1}, {expression, NULL}},		
                               {229, {{function_head, rangle},-1}, {expression, NULL}},		
                               {230, {{function_head, rpar},-1},  {expression, NULL}},			
                               
                               {241, {{mod_scrap}},			{statement, "_f"},cwebx},	
                               {242, {{short_lbrace, mod_scrap},1},	{statement, NULL},cwebx},	
                               {243, {{mod_scrap, magic}},		{declaration, "f__f"},cwebx},	
                               {244, {{lbrace, mod_scrap, magic},1},
                               			{declaration, "__f"},cwebx|standard_braces},	
                               {245, {{short_lbrace, mod_scrap, magic},1}, {declaration, NULL},cwebx},	
                               {246, {{short_struct_head, mod_scrap, magic},1},
                               					    {declaration,NULL},cwebx},	
                               {247, {{mod_scrap, magic, magic}},	{expression, NULL},cwebx},	
                               {248, {{lbrace, mod_scrap, magic, magic},1},
                               			{expression, NULL},cwebx|standard_braces},	
                               {249, {{short_lbrace, mod_scrap, magic, magic},1},
                               			{expression, NULL},cwebx},			
                               
                               {241, {{mod_scrap, semi}}, 	{statement, "__f"},compatibility},	
                               {242, {{mod_scrap, magic}}, 	{statement, "__f"},compatibility},	
                               {243, {{short_lbrace, mod_scrap, semi},1},
                               				{statement, NULL},compatibility},	
                               {244, {{short_lbrace, mod_scrap, magic},1},
                               				{statement, NULL},compatibility},	
                               {245, {{mod_scrap}},        	{expression, NULL},compatibility},	
                               {246, {{statement, function}},	 {function, "_F_"},compatibility},	
                               
                               {250, {{binop, binop}},			{binop,"r__"},compatibility},	
                               {251, {{unorbinop, binop}},		{binop,"r__"},compatibility},	
                               {252, {{lpar, expression, comma}}, {lpar, "___p3"}, compatibility},	
                               
                               {260, {{case_like, binop}},	{expression, "_o{_}"},only_plus_plus},	
                               {260, {{case_like, unorbinop}},	{expression, "_o_"},only_plus_plus},	
                               {260, {{case_like, unop}},	{expression, NULL},only_plus_plus},	
                               {260, {{case_like, langle}},	{expression, "_o_"},only_plus_plus},	
                               {260, {{case_like, rangle}},	{expression, "_o_"},only_plus_plus},	
                               {261, {{case_like, subscript}},	{expression, NULL},only_plus_plus},	
                               {262, {{case_like, lpar, rpar}}, {expression, NULL},only_plus_plus},	
                               {263, {{colcol, expression}},	{expression, NULL},only_plus_plus},	
                               {263, {{colcol, int_like}},	{int_like, NULL},only_plus_plus},	
                               {264, {{expression, colcol, expression}},
                               				{expression, NULL},only_plus_plus},	
                               {264, {{expression, colcol, int_like}},
                               				{int_like, NULL},only_plus_plus},	
                               {265, {{int_like, colcol, expression}},
                               				{expression, NULL},only_plus_plus},	
                               {265, {{int_like, colcol, int_like}},
                               				{int_like, NULL},only_plus_plus},	
                               {266, {{int_like, int_like, colcol, expression},1},
                               				{expression, NULL},only_plus_plus},	
                               {266, {{int_like, int_like, colcol, int_like},1},
                               				{int_like, NULL},only_plus_plus},	
                               
                               {270, {{int_like, binop, expression}},
                               				{int_like, NULL},only_plus_plus},	
                               {271, {{int_like, colon, case_like, int_like}},
                                                               {int_like,"_m__~_"},only_plus_plus },	
                               {272, {{int_like,colon, int_like}}, {int_like,"_m__"},only_plus_plus },	
                               {273, {{struct_like, int_like, comma},-1},
                               				{int_like,"_~_"},only_plus_plus},	
                               {274, {{struct_like, int_like, rangle},-1},
                               				{int_like,"_~_"},only_plus_plus},	
                               {275, {{langle, int_like, expression, comma}},
                               				{langle,"__~!__p5"},only_plus_plus},	
                               {276, {{langle, int_like, expression, rangle}},
                               			{templ_params,"a__~!_a_"},only_plus_plus},	
                               {277, {{struct_like,int_like,expression},-1},
                               				{int_like,"_~_"},only_plus_plus},	
                               {278, {{struct_like,expression,expression},-1},
                               				{int_like,"_~_"},only_plus_plus},	
                               {279, {{template_like,templ_params,declaration}},
                               				{declaration,"__f+_-"},only_plus_plus},	
                               
                               {280, {{int_like, parameters}},	{function_head, "!__"},only_plus_plus},	
                               {281, {{int_like, function_head, semi},-1},
                               				{function_head, "_~_"},only_plus_plus},	
                               {282, {{int_like, function_head, statement},-1},
                               				{function_head, "_~_"},only_plus_plus},	
                               {283, {{function_head, semi}},	{declaration, NULL},only_plus_plus},	
                               {284, {{function_head, int_like}},
                               				{function_head, "_ _"},only_plus_plus},	
                               {285, {{int_like,expression,int_like},1},
                                                               {function_head,"_ _"},only_plus_plus},	
                               {286, {{expression, binop, function_head}},
                               				{expression,NULL},only_plus_plus},	
                               {287, {{return_like, function_head}},
                               				{expression,"_~_"},only_plus_plus},	
                               
                               {288, {{function_head,colon,expression,lbrace},-1},
                               			{function_head,"_+p1m__-"},only_plus_plus},	
                               {289, {{sizeof_like,subscript}},{sizeof_like,NULL},only_plus_plus},	
                               
                               {290, {{throw_like,expression}},  {expression,"_~_"},only_plus_plus},	
                               {291, {{throw_like,function_head}}, {expression,"_~_"},only_plus_plus},	
                               {292, {{throw_like,parameters}},  {int_like,NULL},only_plus_plus},	
                               {293, {{throw_like,semi}},   	  {statement,NULL},only_plus_plus},	
                               {295, {{try_like}},		  {if_else_head,"f_"},only_plus_plus},	
                               {296, {{catch_like,parameters}},  {if_else_head,"f__"},only_plus_plus},	
   };
  int i=array_size(rule);  do install_rule(&rule[--i]); while (i>0);
  #ifdef DEBUG
    if (install_failed) fatal("inconsistent grammar",0);
  #endif
  }
  
  { char* line1=" \\input cwebxmac";
    out_ptr=&out_buf[0];  do *out_ptr++=*line1++; while (*line1!='\0');
    if (compatibility_mode) out_ptr[-4]='c'; /* change to \.{cwebcmac} */
  }
  
  { char *p="_abcdefghijklmnopqrstuvwxyz0123456789";
    int c='\1', k=2;
    collate[0]='\0'; collate[1]=' ';
    do  if (!isalnum(c) && c!='_' && c!=' ') collate[k++]=c; 
    while (++c<=UCHAR_MAX);
    while ((c=*p++)!='\0') collate[k++]=c;
    end_collate=k; /* record the length */
  }
  if (show_banner)
    print("%s, in %s mode.\n",banner,C_plus_plus ? "C++" : "C");
    /* print a ``banner line'' */
  
  { int i; static char* int_likes[]=
      { "auto","char","double","float","int","long","register"
      , "short","signed","static","unsigned","void" };
    static char* defined_types[] =
      { "FILE", "size_t", "ptrdiff_t", "wchar_t"
      , "jmp_buf", "sig_atomic_t", "fpos_t", "div_t", "ldiv_t"
      , "clock_t","time_t"
      , "va_list"
      };
    static char* return_likes[]=
      {"break","continue","goto","return"};
    int int_like_nr=array_size(int_likes),
        defined_type_nr=array_size(defined_types),
        return_like_nr=array_size(return_likes);
  
    for (i=0; i<int_like_nr; ++i) id_lookup(int_likes[i],NULL,int_like);
    for (i=0; i<defined_type_nr; ++i)
      id_lookup(defined_types[i],NULL,type_defined);
    for (i=0; i<return_like_nr; ++i) id_lookup(return_likes[i],NULL,return_like);
  
    id_lookup("case", NULL, case_like);
    id_lookup("const", NULL, const_like);
    id_lookup("constexpr", NULL, const_like);
    id_lookup("default", NULL, case_like);
    id_lookup("decltype", NULL, sizeof_like);
    id_lookup("do", NULL, do_like);
    id_lookup("else", NULL, else_like);
    id_lookup("enum", NULL, struct_like);
    id_lookup("extern", NULL, C_plus_plus ? extern_like : int_like);
    id_lookup("for", NULL, while_like);
    id_lookup("if", NULL, if_like);
    id_lookup("mutable", NULL, const_like);
    id_lookup("noexcept", NULL, const_like);
    id_lookup("sizeof", NULL, sizeof_like);
    id_lookup("struct", NULL, struct_like);
    id_lookup("switch", NULL, while_like);
    id_lookup("typedef", NULL, typedef_like);
    id_lookup("union", NULL, struct_like);
    id_lookup("va_dcl",NULL, declaration);
    id_lookup("volatile", NULL, const_like);
    id_lookup("while", NULL, while_like);
    id_lookup("NULL", NULL, NULL_like);
    id_lookup("nullptr", NULL, NULL_like);
    id_lookup("TeX", NULL, TeX_like);
    if (C_plus_plus) 
                     { int i;
                       static char* cpp_types[] =
                       { "exception", "bad_exception", "bad_cast", "bad_typeid", "logic_error",
                         "domain_error", "invalid_argument", "length_error", "out_of_range",
                         "bad_alloc", "runtime_error", "range_error", "overflow_error",
                         "underflow_error",
                         "string",
                         "iterator", "const_iterator", "reverse_iterator", "const_reverse_iterator",
                         "size_type", "value_type",
                         "ios_base","ios", "istream", "ostream", "iostream",
                         "istringstream", "ifstream", "ostringstream", "ofstream",
                         "stringstream", "fstream",
                         "streambuf", "stringbuf", "filebuf",
                         "streamoff", "streampos",
                         "input_iterator_tag", "output_iterator_tag", "forward_iterator_tag",
                         "bidirectional_iterator_tag", "random_access_iterator_tag",
                         "pair", "auto_ptr", "allocator", "raw_storage_iterator",
                         "vector", "list", "deque",
                         "set", "multiset", "map", "multimap",
                         "stack", "queue", "priority_queue", "bitset",
                         "shared_ptr", "weak_ptr", "unique_ptr"
                       };
                       for (i=0; i<array_size(cpp_types); ++i)
                         id_lookup(cpp_types[i],NULL,type_defined);
                       id_lookup("asm", NULL, int_like);
                       id_lookup("and", NULL, and_like);
                       id_lookup("bool", NULL, int_like);
                       id_lookup("catch", NULL, catch_like);
                       id_lookup("class", NULL, struct_like);
                       id_lookup("delete", NULL, sizeof_like);
                       id_lookup("explicit", NULL, int_like);
                       id_lookup("false", NULL, expression);
                       id_lookup("friend", NULL, int_like);
                       id_lookup("inline", NULL, int_like);
                       id_lookup("namespace", NULL, namespace_like);
                       id_lookup("new", NULL, sizeof_like);
                       id_lookup("not", NULL, not_like);
                       id_lookup("operator", NULL, case_like);
                       id_lookup("or", NULL, and_like);
                       id_lookup("private", NULL, case_like);
                       id_lookup("protected", NULL, case_like);
                       id_lookup("public", NULL, case_like);
                       id_lookup("template", NULL, template_like);
                       id_lookup("this", NULL, expression);
                       id_lookup("throw", NULL, throw_like);
                       id_lookup("true", NULL, expression);
                       id_lookup("try", NULL, try_like);
                       id_lookup("typeid", NULL, sizeof_like);
                       id_lookup("typename", NULL, typename_like);
                       id_lookup("using", NULL, int_like);
                       id_lookup("virtual", NULL, int_like);
                       id_lookup("xor", NULL, and_like);
                       id_lookup("const_cast", NULL, sizeof_like);
                       id_lookup("static_cast", NULL, sizeof_like);
                       id_lookup("dynamic_cast", NULL, sizeof_like);
                       id_lookup("reinterpret_cast", NULL, sizeof_like);
                     }
  }
  phase_one (); /* read all the user's text and store the \xr.s */
  if (history>harmless_message) wrap_up(); /* stop in case of trouble */
  open_output_file();
  phase_two (); /* read all the text again and translate it to \TeX\ form */
  if (history>harmless_message) wrap_up(); /* stop in case of trouble */
  phase_three (); /* output the \xr. index */
  wrap_up (); /* and exit gracefully */
  return 0; /* for completeness---not reached */
}

void new_id_xref (id_pointer p)
{ sixteen_bits f=xref_switch; xref_switch=0;
  if (p->ilk==reference) f=0;
  else if (f==0 && (unindexed(p) || length(p)==1)
        || no_xref || including_header_file) return;
  if ((p->xref->num&num_mask)==section_count) p->xref->num|=f;
  else
  { make_xref(section_count|f,xref_index(p->xref)); p->xref=xref_ptr; }
}

void new_mod_xref (mod_pointer p)
{ sixteen_bits head, *q, m=section_count+mod_xref_switch;
  if (p->xref->num==file_flag) q=&p->xref->next; /* skip |file_flag| */
  else head=xref_index(p->xref),q=&head;
  if (mod_xref_switch!=def_flag)
  { while (xnum(*q)>m)
      q=&xlink(*q); /* skip the $d_i$'s and possibly $c_i$'s */
    if (xnum(*q)==m) return; /* don't duplicate */
  }
  make_xref(m,*q); mod_xref_switch=0;
  if (q==&head) p->xref=xref_ptr;  else *q=xref_index(xref_ptr);
}

void set_file_flag (mod_pointer p)
{ if (p->xref->num!=file_flag)
    { make_xref(file_flag,xref_index(p->xref)); p->xref=xref_ptr; }
}

boolean names_match (id_pointer x, char* q, int l, int ilk)
{ char* p=name_begin(x);
  if ((x->ilk==ilk || ilk==normal && reserved(x)))
  { while (--l>=0) if (*p++!=*q++) return false; return *p=='\0'; }
  else return false;
}

void init_id_name (id_pointer p, int t)
{ p->ilk = t; p->xref = &xmem[0]; }

void init_module_name (mod_pointer p) { p->xref=&xmem[0]; }

void copy_limbo (void) /* copy \TeX\ code until the next section begins */
{ while (loc<=limit || (finish_line(),get_line()))
  { eight_bits c;
    limit[1]='@'; /* place a sentinel */
    while ((c=*loc++)!='@')
      if (!(output_line_empty() && isspace(c))) out(c);
    if (loc<=limit) /* then we have hit a control code */
      switch(code_of(*loc++))
      {
      case new_section: return; /* the only exit, unless no sections exist */
      case ignored_text: get_control_text(); break;
      case format: get_next(); get_next(); break; /* skip two identifiers */
      case char_trans: out_str("\\ATL "); break; 
      default: err_print("! Double @ required in limbo part");
			  
        /* fall through */
      case at_sign_image: out('@');
      }
  }
}

int skip_TeX (void) /* skip past pure \TeX\ code */
{ char c;
  while (find_char())
  { limit[1]='@';
    while ((c=*loc++)!='@' && c!='%')
      if (c=='|') return c;
      else if (c=='\\' && *loc!='@') ++loc;
	/* ignore `\.{\\\%}'  and `\.{\\\v}' */
    if (loc<=limit)
      if (c=='@') return code_of(*loc++);
      else /* ignore remainder of line unless a major control code occurs */
	do 
	  if ((c=*loc++)=='@' && code_of(*loc++)>=format)
	    return code_of(loc[-1]);
	while (loc<limit);
  }
  return new_section;
}

int copy_TeX (void) /* copy pure \TeX\ material */
{ eight_bits c; /* current character being copied */
  while (loc<=limit || (finish_line(),get_line()))
  { limit[1]='@';
    while((c=*loc++)!='@')
    { if (c=='|') return '|';
      if (!(output_line_empty() && isspace(c))) out(c);
      if (c=='%') break;
      if (c=='\\' && *loc!='@') out(*loc++);
	/* copy `\.{\\\%}' and `\.{\\\v}' */
    }
    if (loc<=limit)
      if (c=='@') return code_of(*loc++);
      else /* ignore remainder of line unless a major control code occurs */
	do
	  if ((c=*loc++)=='@' && code_of(*loc++)>=format)
	    return finish_line(),code_of(loc[-1]);
	while(loc<limit);
  }
  return new_section;
}

int scan_comment (int* bal, boolean one_liner)
{ char c; boolean forced_out=false; /* prematurely terminated? */
  while (one_liner ? loc<limit
		   : find_char() && (*loc!='*' || loc[1]!='/' ))
  
  if (including_header_file) ++loc; /* don't process characters here */
  else
  { switch(c=*loc++)
    {
    case '|': return '|'; /* beginning of `\pb' inside comment */
    case '@':
      if (*loc++!='@')
        if (code_of(loc[-1])!=new_section)
  	err_print("! Double @ required in comment");
  		   
        else
        { err_print("! Section ended in mid-comment");
  		     
  	forced_out=true; goto done;
        }
      break;
    case '\\':
      if (*loc!='@') { if (phase==2) app_char_tok(c); c=*loc++; }
      break;
    case '{': ++*bal; break;
    case '}':
      if (*bal>0) --*bal;
      else err_print("! Extra } in comment");
  		     
      break;
    case '/':  if (*loc=='*') err_print("! Nested comment");
  					 
  
    }
    if (phase==2) app_char_tok(c);
  }
  if (input_has_ended)
    forced_out=true,err_print("! Input ended in mid-comment");
			       
  else if (!one_liner) loc+=2; /* move past `\.{*{}/}' */
done:
  if (*bal>0) err_print("! Too few closing braces in comment");
			 
  return forced_out ? new_section : end_comment;
}

int get_next (void) /* produces the next input token */
{ eight_bits c; /* the current character */
restart:
  if (!find_char()) { preprocessing=0; return new_section; }
  
  if (preprocessing>0 && loc==limit)
  { preprocessing=0; return end_preproc; }
  if ((c=*loc++)=='@')
    
    if (including_header_file) goto restart; /* ignore `\.@' in header files */
    else
    { int cc=code_of(*loc++);
      switch (cc)
      { case ignore: goto restart;
        case underline: xref_switch=def_flag; goto restart;
    #ifdef DEBUG
        case trace0: case trace1: case trace2: case trace3: 
          if (phase==2) tracing=cc;  goto restart;
    #endif
        case char_trans:
          err_print("! `@l' only allowed in limbo"); goto restart;
    		 
        case ASCII_code: 
                         { id_first=&mod_text[1]; strncpy(id_first,"@'",2); id_loc=&id_first[2];
                           while ((*id_loc++=c=*loc++)!='\'')
                           { if (c=='\\')
                               *id_loc++=*loc++; /* copy any character following backslash */
                             else if (c=='@' && *loc++!='@')
                             { err_print("! Double @ required in strings"); --loc; }
                         		   
                             if (loc>=limit) { err_print("! ASCII constant didn't end"); break; }
                         				   
                           }
                         }
      return string;
        case module_name:
          
          { boolean file_module=loc[-1]=='(';
            cur_mod=get_module_name();
            if (file_module && phase==1 && cur_mod!=NULL) set_file_flag(cur_mod);
          }
      break;
        case ignored_text: get_control_text(); goto restart;
        case verbatim: case TeX_string: get_control_text(); break;
        case xref_roman: case xref_wildcard: case xref_typewriter:
        case xref_mark: case refer:
          if (get_control_text()) goto restart; /* don't index empty strings */
          if (cc==refer) cur_id=id_lookup(id_first,id_loc,reference);
          else if (phase==1)
    	cur_id=id_lookup(id_first,id_loc,cc-xref_roman+roman);
      }
      return cc;
    }
  if (isspace(c))
    if (preprocessing>0) return ' '; /* keep spaces in preprocessor lines */
    else goto restart; /* ignore other white space */
  if (c=='L' && (*loc=='\'' || *loc=='"'))
  { get_string(); return string; }
  if (isalpha(c) || c=='_' || c>=0x80)
    { 
      {  id_first=--loc; /* mark beginning of identifier */
         do c=*++loc; while (isalnum(c) || c=='_' || c>=0x80);
         cur_id= id_lookup(id_first,loc,normal);
      }
 return identifier; }
  if (isdigit(c) || c=='.' && isdigit((eight_bits)*loc))
    { 
      { id_first=id_loc=&mod_text[1];
      
        if (c=='0' && (isdigit(c=*loc) || tolower(c)=='x')) /* octal or hex */
        { if (isdigit(c)) /* octal constant with at least two digits */
          { *id_loc++ = '~'; /* store `\.\~' in place of leading `\.0' */
            do shift_and_store(c); while (isdigit(c));
      	/* copy second and following digits */
          }
          else /* hex constant */
          { shift_and_store('^'); /* replace `\.{0x}' by `\.\^' */
            while (isxdigit(c)) shift_and_store(c);
          }
        }
        else /* decimal constant */
        { c=*--loc; /* recover first digit or decimal point */
          while (isdigit(c)) shift_and_store(c);
          if (c=='.')  do shift_and_store(c); while (isdigit(c));
          if (tolower(c)== 'e') /* floating point constant with exponent */
          { shift_and_store('_'); /* replace `\.e' by `\.\_' */
            if (c=='+' || c=='-') { *id_loc++ = c; c=*++loc; }
            while (isdigit(c)) shift_and_store(c); /* exponent */
          }
        }
        if (isalpha(c)) /* `\.{U}', `\.{L}', and/or `\.{F}' suffix */
        { *id_loc++ = '$'; 
          do shift_and_store(c); while (isalpha(c)); }
      }
 return constant; }
  if (c=='\'' || c=='"' || (c=='<' && preprocessing==2))
    { get_string(); return string; }
  if (c=='#' && loc==&buffer[1])
  { 
    { while (loc<limit && isspace((eight_bits)*loc)) ++loc;
        /* allow spaces after `\.\#' */
      if (limit-loc>=7 && strncmp(loc,"include",7)==0) /* `\.{\#include}' line */
        if (including_header_file) /* start nested header file */
        { loc+=7; push_header_file(false); goto restart; }
        else preprocessing=2;
      else preprocessing=1;
    }
     return start_preproc; }
  if (c=='\\' && preprocessing>0 && loc==limit)
    { ++loc; /* move past |limit|, so |get_line| will be called */
      goto restart;
    }
  
  switch (c) {
  case '/': compress2('*',begin_comment); 
    if (C_plus_plus) compress2('/',begin_comment);
    comp_ass_op2(div_assign); break;
  case '*': compress2('/',end_comment);	comp_ass_op2(mul_assign);  break;
  case '%': comp_ass_op2(mod_assign);	break;
  case '+': compress2('+',plus_plus);	comp_ass_op2(plus_assign); break;
  case '-': compress2('-',minus_minus);	compress2 ('>', minus_gt);
  	  comp_ass_op2(minus_assign);	break;
  case '=': compress2('=',eq_eq);		break;
  case '>': compress2('=',gt_eq);		comp_ass_op3('>',right_assign);
  	  compress2 ('>',gt_gt);	break;
  case '<': compress2('=', lt_eq);	comp_ass_op3('<',left_assign);
  	  compress2 ('<', lt_lt);	break;
  case '&': compress2('&',and_and);	comp_ass_op2(and_assign);   break;
  case '^': comp_ass_op2(xor_assign);	break;
  case '|': compress2('|',or_or);		comp_ass_op2(or_assign);    break;
  case '!': compress2('=',not_eq);	break;
  case '.': compress3('.','.', ellipsis); break;
  case '#': compress2 ('#', sh_sh);	break;
  case ':':  if (C_plus_plus)		compress2 (':',colon_colon);
  }
  return c;
}

boolean push_header_file(boolean suspend)
{ id_pointer p;
  if (locate_file_name() &&
      (p=id_lookup(id_first,id_loc,header_file_name))->xref==&xmem[0])
  { p->xref=&xmem[1]; /* mark file as seen, to avoid multiple inclusion */
    return push_input_file(true,suspend);
  }
  return false;
}

void phase_one (void)
  /* read all the user's text and store the \xr.s */
{ phase=1; reset_input(); section_count=0;
  
  while (find_char())
  { limit[1]='@'; /* place a sentinel */
    while (*loc++!='@') {}
    if (loc<=limit) /* note that |loc!=limit+1| since |*limit==' '| */
    { int c=code_of(*loc++);
      if (c==new_section) break;
      if (c==format) 
                   if (tolower((eight_bits)loc[-1])=='f')
                     err_print("! Double @ required in limbo part");
                   	     
                   else
                   { id_pointer lhs;
                     if (shift()==identifier && (lhs=cur_id,shift()==identifier))
                       lhs->ilk=cur_id->ilk;
                     else err_print("! Improper format definition");
                   		  
                   }
    }
  }
  while (!input_has_ended)
    
    { if (++section_count==max_sections)
        overflow("section number"); 
      if (loc[-1]=='*') print_section_progress ();
      
      do
        switch (next_control=skip_TeX())
        { case underline: xref_switch=def_flag; break;
          case '|': C_xref(true); break;
          case module_name: case refer: loc-=2; get_next(); break;
          case ignored_text: get_control_text(); break;
          case char_trans: err_print("! `@l' only allowed in limbo"); break;
      				
          case xref_roman: case xref_wildcard: case xref_typewriter:
          case xref_mark: loc-=2; get_next(); new_id_xref(cur_id);
        }
      while (next_control<format);
      
      while (next_control<begin_C) /* |format|, |definition| or |header| */
        if(next_control!=header)
        { xref_switch=def_flag; /* implied \:! for first identifier */
          if (next_control==format) 
                                  { boolean f= tolower((eight_bits)loc[-1])=='f';
                                    id_pointer lhs;
                                    if (shift()==identifier && (lhs=cur_id,shift()==identifier))
                                    { if (f) new_id_xref(lhs);  else xref_switch=0;
                                      lhs->ilk=cur_id->ilk;
                                    }
                                    else err_print("! Improper format definition");
                                  		  
                                  }
          outer_xref(); /* macro definition or comment after format definition */
        }
        else 
             { if (push_header_file(true)) /* prepare for reading header file */
                 including_header_file=true; /* will be reset on closing the file */
               typedef_tracking(true); /* this is what we are doing it for */
               outer_xref();
                 /* |shift()| and  collect typedefs until |next_control>=format| */
               typedef_tracking(false);
             }
      
      
      { if (next_control<new_section) /* |begin_C| or |module_name| */
        { typedef_tracking(true);
          mod_xref_switch= next_control==module_name ? def_flag : 0;
          do
          { if (next_control==module_name && cur_mod!=NULL)
      	new_mod_xref(cur_mod);
            outer_xref();
          } while (next_control<new_section);
          typedef_tracking(false);
        }
      }
      if (section_changed(section_count)) change_exists=true;
     }
  if (change_exists) mark_section_as_changed(section_count);
    /* the index changes if anything does */
  
  mod_check(root);
  
  { id_pointer name; id_pointer *h; /* pointer into |hash| */
    for (h=hash; h<hash_end; h++)
      for (name=*h; name!=NULL; name=name->hash_link)
        if (name->ilk!=header_file_name)
        /* traverse hash lists, except for header file names */
      
      {
        sixteen_bits x=xref_index(name->xref),t,y=0;
          /* index of the sentinel node */
        while (xnum(x)!=0)
         { t=x; x=xlink(t); xlink(t)=y; y=t; }
        name->xref=&xmem[y]; /* don't forget to link in the reversed list */
      }
  }
}

void C_xref (boolean inner)
{ while (next_control<format || next_control==module_name && inner)
  { if (preprocessing==0)
      
      { if (typedef_master==0 &&
            next_control==identifier && cur_id->ilk==typedef_like)
        { typedef_master=2; brace_level=par_level=0; }
        else if (typedef_master>0) switch(next_control)
        { case identifier:
            if (brace_level==0)
      	if (typedef_master==2)
      	{ if (cur_id->ilk==int_like || cur_id->ilk==type_defined)
      	    typedef_master=4;
      	  else if (cur_id->ilk==struct_like) typedef_master=3;
      	}
      	else if (typedef_master==4)
      	{ if(cur_id->ilk==normal||cur_id->ilk==type_defined) /* this is it */
      	    cur_id->ilk=type_defined, typedef_master=1;
      	}
      	else if (typedef_master==3) typedef_master=4;
            break;
          case '{': 
            if (brace_level++==0 && typedef_master==3) typedef_master=4;  break;
          case '}': --brace_level; break;
          case '<':  if (C_plus_plus) brace_level++;  break;
          case '>':  if (C_plus_plus) --brace_level;  break;
          case ',': 
            if (typedef_master==1 && par_level==0) typedef_master=4;  break;
          case '(': ++par_level; break;
          case ')': --par_level; break;
          case ';': 
            if (brace_level==0)
            { if (typedef_master>=2)
                
                { if (including_header_file)
                    print("In file %s:\n",cur_file_name);
                  print("\nUnrecognised typedef at line %d in section %d:\n"
                	   ,cur_line, section_count);
                  mark_harmless();
                }
              typedef_master=0;
            } 
            break;
          case colon_colon:
              if (C_plus_plus && brace_level==0 && typedef_master==4)
             typedef_master=2;  break;
        }
        if (C_plus_plus)
        
        { static int class_seen=0; static id_pointer this_id;
          switch (class_seen)
          { case 0:
            if (next_control==identifier)
              if (cur_id->ilk==struct_like) class_seen=1;
              else if(cur_id->ilk==typename_like)
                class_seen=2;
          break;
            case 1:
            if (next_control==identifier && cur_id->ilk==normal)
              cur_id->ilk=type_defined;
            class_seen=0;
          break;
            case 2:
            if (next_control==identifier && cur_id->ilk==normal)
              { this_id=cur_id; class_seen=3; }
            else class_seen=0;
          break;
            case 3:
            if (next_control!=colon_colon) this_id->ilk=type_defined;
            class_seen=0;
          }
        }
      }
    if (next_control>=identifier && next_control<=xref_mark)
      new_id_xref(cur_id);
    else if (next_control==module_name && cur_mod!=NULL)
      mod_xref_switch=cite_flag,new_mod_xref(cur_mod);
    if (next_control==start_preproc && shift()!=end_preproc
      &&next_control!=identifier)
      err_print("! Identifier should follow `#'");
		 
    else shift();
    if (next_control=='|' && inner
     || next_control==begin_comment || next_control==end_comment)
      return;
  }
}

void outer_xref (void) /* extension of |C_xref| */
{ shift(); /* move past previously processed token */
  while (next_control<format)
    if (next_control!=begin_comment) C_xref(false);
    else
    { boolean one_liner=loc[-1]=='/'; int bal=0; /* brace level in comment */
      typedef_tracking(false);
      while ((next_control=scan_comment(&bal,one_liner))=='|')
	{ C_xref(true); if (next_control!='|') break; }
      typedef_tracking(true);
    }
}

void mod_check (mod_pointer p) /* print anomalies in subtree |p| */
{ if (p != NULL)
  { mod_check (p->llink); /* traverse left subtree */
    { boolean file_module = p->xref->num==file_flag;
      sixteen_bits head, *q, threshold;
        /* lower limit of |num| values of current interest */
      if (file_module) q=&p->xref->next; 
      else head=xref_index(p->xref),q=&head;
      if (!complete_name(p))
      { print("\n! Never completed"); print_mod(p); mark_harmless(); }
		     
      if (xnum(*q)<=(threshold=def_flag))
      { print("\n! Never defined"); print_mod(p); mark_harmless(); }
		     
      else
      
      { sixteen_bits x=*q,y=0,t;
        do { t=xlink(x); xlink(x)=y; y=x; } while (xnum(x=t)>threshold);
        xlink(t=*q)=x; *q=y; q=&xlink(t);
      }
      if (xnum(*q)>(threshold=cite_flag)) 
                                       { sixteen_bits x=*q,y=0,t;
                                         do { t=xlink(x); xlink(x)=y; y=x; } while (xnum(x=t)>threshold);
                                         xlink(t=*q)=x; *q=y; q=&xlink(t);
                                       }
      if (xnum(*q)==(threshold=0))
      { if(!file_module)
	  { print("\n! Never used"); print_mod(p); mark_harmless(); }
		       
      }
      else 
           { sixteen_bits x=*q,y=0,t;
             do { t=xlink(x); xlink(x)=y; y=x; } while (xnum(x=t)>threshold);
             xlink(t=*q)=x; *q=y; q=&xlink(t);
           }
      if (!file_module) p->xref=&xmem[head];
	/* set pointer to possibly modified value */
    }
    mod_check (p->rlink); /* traverse right subtree */
  }
}

void phase_two (void)
   /* read all the text again and translate it to \TeX\ form */
{ phase=2; reset_input ();
  print_progress("\nWriting the output file...");
		  
  section_count=0; copy_limbo(); finish_line();
  tex_new_line(); /* insert a blank line, it looks nice */
  while (!input_has_ended) 
                           { section_count++;
                             
                             { out('\\'); out(loc[-1]=='*' ? 'N' : loc[-1]=='~' ? 'n' : 'M' );
                                 
                               if (loc[-1]=='*')
                               { print_section_progress(); 
                                                           { if (*loc=='*') ++loc,out_str("-1");
                                                             else if (!isdigit((eight_bits)*loc)) out('0');
                                                             else do out(*loc++); while (isdigit((eight_bits)*loc));
                                                             out(' '); /* terminate level by a space */
                                                           }
                              }
                               out_sec_nr(section_count); out_str(". ");
                             }
                             
                             do
                               switch (next_control=copy_TeX())
                               { case '|': typedef_master=0; do_C(); break;
                                 case at_sign_image: out('@'); break;
                                 case thin_space: case math_break: case ASCII_code: case line_break:
                                 case big_line_break: case no_line_break: case join: case pseudo_semi:
                                 case force_expr_open: case force_expr_close:
                                   err_print("! You can't do that in TeX text");
                             		 
                                   break;
                             #ifdef DEBUG
                                 case trace0: case trace1: case trace2: case trace3: tracing=next_control;
                                   break;
                             #endif
                                 case module_name: loc-=2; get_next(); break; /* get module name */
                                 case refer: loc-=2; get_next(); /* get name referred to */
                                   if (cur_id->xref->num==0) err_print("! Undefined reference");
                                   else list_refs(cur_id->xref,0);
                                   break;
                                 case TeX_string: err_print("! TeX string should be in C text only");
                             				
                                 /* fall through */
                                 case xref_roman: case xref_wildcard: case xref_typewriter:
                                 case xref_mark: case ignored_text:
                                   get_control_text(); /* skip to \:> */
                               }
                             while (next_control<format);
                             if (next_control<begin_C)
                             { emit_space();
                               
                               { typedef_tracking(false);
                                 do
                                 { boolean suppressed=false; /* whether output suppressed by \:s */
                                   if (next_control==format) 
                                                           if (tolower((eight_bits)loc[-1])=='s')
                                                           { suppressed=true; shift(); shift(); shift(); }
                                                             /* skip format definition */
                                                           else
                                                           { int saved_code=0,saved_mathness;
                                                             app_str("\\F"); shift(); /* this will produce `\&{format}' */ 
                                                             if (cur_id->ilk!=TeX_like && cur_id->ilk!=NULL_like)
                                                               app(id_flag+id_index(cur_id));
                                                             else 
                                                                  { char* p=name_begin(cur_id);
                                                                    saved_mathness=cur_id->ilk==TeX_like ? no_math : yes_math;
                                                                    saved_code=id_flag+id_index(cur_id);/* save to print afterwards */
                                                                    app_str("\\\\{"); 
                                                                    do { if (*p=='_') app('\\'); app_tok(*p); } while (*++p!='\0');
                                                                    app('}'); check_toks(10);
                                                                  }
                                                             app('~'); pack_scrap(insert,yes_math); shift();
                                                             app((cur_id->ilk==normal || cur_id->ilk==TeX_like || cur_id->ilk==NULL_like
                                                                 ? id_flag : res_flag
                                                                 )+id_index(cur_id));
                                                           check_scrap();
                                                             pack_scrap(insert,cur_id->ilk==TeX_like ? no_math : yes_math);
                                                             shift();
                                                             if (saved_code!=0)
                                                             { app_str("\\quad("); app(saved_code); app(')');
                                                             check_scrap(); pack_scrap(insert,saved_mathness);
                                                             }
                                                           }
                                   else if (next_control==definition) 
                                                                    { if (shift()!=identifier)
                                                                        err_print("! Improper macro definition");
                                                                    	       
                                                                      else
                                                                      { app_str("\\D$"); 
                                                                              /* this will produce \&{\#define} */ 
                                                                        app(id_flag+id_index(cur_id));
                                                                        if (*loc=='(')
                                                                        { shift();
                                                                          do
                                                                          { app_char_tok(next_control);
                                                                    	if (shift()!=identifier) break;
                                                                    	app(id_flag+id_index(cur_id));
                                                                          } while(shift()==',');
                                                                          check_toks(2);
                                                                          if (next_control==')') { app(')'); shift(); }
                                                                          else err_print("! Improper macro definition");
                                                                        }
                                                                        else shift();
                                                                        app('$');
                                                                        app(break_space); pack_scrap(insert,no_math);
                                                                      }
                                                                    }
                                   else 
                                        { app_str("\\h"); /* this will produce \&{\#include} */ 
                                          pack_scrap(insert,no_math);
                                          { int save=preprocessing; preprocessing=2; /* emulate `\.{\#include}' */
                                            while (shift()==' ') {} /* skip spaces and read file name as string */
                                            preprocessing=save;
                                          }
                                        }
                                   if (!suppressed) outer_read(), finish_C();
                                   else if (next_control<format)
                                   { err_print("! Improper stuff after `@s' format definition");
                               		 
                                     if (next_control==begin_comment) loc-=2; /* try to get back in phase */
                                     outer_xref(); /* skip illegal stuff */
                                   }
                                 } while (next_control<begin_C); /* |format|, |definition|, or |header| */
                               }
                             }
                             if (next_control<new_section)
                             { mod_pointer this_module=NULL; /* the current module name */
                               emit_space(); 
                                             { typedef_master=0;
                                               if (next_control==begin_C) shift();
                                               else
                                               { this_module=cur_mod; /* register the name for this module */
                                                 
                                                 { if (shift()=='=' || next_control==eq_eq || next_control==plus_assign)
                                                   { if (next_control!=plus_assign || shift()=='=') shift(); }
                                                       /* accept `\.=', `\.{==}', `\.{+=}' or `\.{+==}' */
                                                   else err_print("! You need an = sign after the module name");
                                                 		  
                                                   if (this_module!=NULL) /* i.e., unless module name was bad */
                                                   { xref_pointer x=this_module->xref;
                                                     if (x->num==file_flag) x=next_xref(x);
                                                     app_str("\\4$"); 
                                                                      /* module name will be flush left */ 
                                                     app(mod_flag+mod_index(this_module));
                                                     if (x->num != section_count+def_flag)
                                                     { app_str("\\PE"); /* module has also been defined before */ 
                                                       this_module = NULL; /* so we won't give \xr. info here */
                                                     }
                                                     else app_str("\\EQ"); /* output a module definition sign */ 
                                                     app_str("{}$"); 
                                                     app(force); pack_scrap(insert,no_math);
                                                       /* this forces a line break unless \:+ follows */
                                                   }
                                                 }
                                               }
                                               do
                                               { outer_read();
                                                 if (next_control==new_section) break;
                                                 if (next_control==module_name) 
                                                                              { if (cur_mod!=NULL)
                                                                                  app(mod_flag+mod_index(cur_mod)), pack_scrap(mod_scrap,yes_math);
                                                                              }
                                                 else err_print("! You can't do that in C text");
                                             		    
                                                   /* |format|, |definition| or |begin_C| */
                                                 shift();
                                               } while (true);
                                               finish_C();
                                             }
                               
                               { if (this_module != NULL)
                                 { xref_pointer foot_ref=this_module->xref;
                                   if (foot_ref->num==file_flag) foot_ref=next_xref(foot_ref);
                                   foot_ref=next_xref(foot_ref); /* don't \xr. to yourself */
                                   footnote(&foot_ref,def_flag);
                                     /* display further defining sections; advance |foot_ref| */
                                   footnote(&foot_ref,cite_flag); /* display any citations */
                                   footnote(&foot_ref,0); /* display uses */
                                 }
                               }
                             }
                             
                             { out_str ("\\fi"); finish_line ();  tex_new_line(); }
                              
                           }
}

void C_read (boolean inner) /* creates scraps from \Cee\ tokens */
{ while (next_control<format || next_control==module_name && inner)
  { 
    { check_scrap(); check_toks(6); /* `\.{\\hbox\{}' */
      switch (next_control)
      { case string: case constant: case verbatim:
          
          { int count = -1; /* characters remaining before string break */
          
            if (next_control==constant) app_str("\\T{"); 
            else if (next_control==string) { count=20; app_str("\\.{"); } 
            else app_str("\\vb{"); 
          
            while (id_first<id_loc)
            { if (count--==0) /* insert a discretionary break in a long string */
                if (id_first[-1]=='\\') count=0; /* no break after backslash */
                else { check_toks(2); app_str("\\)"); count = 20; } 
              if (strchr(" \\#%$^{}~&_",*id_first)!=NULL) app('\\');
               
               
              app_char_tok(*id_first++);
            }
            app('}');
            if (next_control==verbatim) pack_scrap(insert,maybe_math);
            else pack_scrap(expression,yes_math);
          }
      goto done;
        case TeX_string: 
                         { app_str("\\hbox{");
                           while (id_first<id_loc) app_char_tok(*id_first++);
                           app('}');
                           if (!compatibility_mode) pack_scrap(expression,maybe_math);
                         }
      goto done;
        case identifier: 
                         { id_pointer p=cur_id; int cat=p->ilk;
                           
                           { if (typedef_master==0 && cat==typedef_like)
                               typedef_master=2, brace_level=par_level=0;
                             else if (typedef_master>0 && brace_level==0)
                               if (typedef_master==2)
                               { if (cat==int_like || cat==type_defined) typedef_master=4;
                                 else if (cat==struct_like) typedef_master=3;
                               }
                               else if (typedef_master==4 && cat==type_defined) /* this is it */
                                 cat=expression, typedef_master=1;
                               else if (typedef_master==3) typedef_master=4;
                           }
                           if (cat==normal || cat==TeX_like || cat==NULL_like)
                           { app(id_flag+id_index(p));
                             pack_scrap(expression
                                       , cat==TeX_like && !compatibility_mode ? no_math : yes_math);
                           }
                           else
                           { if (cat==and_like || cat==not_like) /* provide text operators with space */
                             { if (cat==and_like) app('~');
                               app(res_flag+id_index(p)); app('~');
                             }
                             else app(res_flag+id_index(p)); /* append reserved word */
                         
                             if (cat==type_defined || cat==const_like || cat==typedef_like)
                               cat=int_like;
                             else if (cat==and_like) cat=binop;
                             else if (cat==not_like) cat=unop;
                             else if (cat==namespace_like || cat==typename_like) cat=struct_like;
                             pack_scrap(cat,maybe_math);
                           }
                         }
      goto done;
        case module_name: 
                          { if (cur_mod!=NULL)
                              app(mod_flag+mod_index(cur_mod)), pack_scrap(mod_scrap,yes_math);
                          }
      goto done;
        case start_preproc:
          
          { app(force); app(flush_left); app_str("\\&\\#");
            if (shift()==identifier)
            { app(res_flag+id_index(cur_id)); pack_scrap(lproc,no_math); }
            else if (next_control==end_preproc)
            { pack_scrap(lproc,no_math);
              check_scrap(); *scrap_ptr++=token_trans[end_preproc];
            }
            else confusion("no identifier after `#'");
          		
          
          }
      goto done;
        case refer:
          err_print("! You can't use `@#' in C text"); /*fall through */
        case ignore: case begin_comment: case end_comment:
        case xref_roman: case xref_wildcard: case xref_typewriter:
        case xref_mark: goto done;
      
      case '{': 
        if (typedef_master>0 && brace_level++==0 && typedef_master==3)
          typedef_master=4;
        break;
      case '}': if (typedef_master>0) --brace_level; break;
      case '<': if (C_plus_plus && typedef_master>0) ++brace_level; break;
      case '>': if (C_plus_plus && typedef_master>0) --brace_level; break;
      case ',': if (typedef_master==1 && par_level==0) typedef_master=4; break;
      case '(': if (typedef_master>0) ++par_level; break;
      case ')': if (typedef_master>0) --par_level; break;
      case ';': if (typedef_master>0 && brace_level==0) typedef_master=0; break;
      case colon_colon:
              if (C_plus_plus && typedef_master>0 && brace_level==0)
                typedef_master=2 ;  break;
        case '|':  if (inner) goto done; /* skip initial `\.\v' of `\pb' */
      }
      *scrap_ptr=token_trans[next_control]; /* fixed scrap for this input token */
      
      { if (dangling_tokens())
        { app_trans(scrap_ptr); scrap_ptr->trans=text_ptr; freeze_text(); }
      }
      ++scrap_ptr; /* incorporate the scrap */
      done: {}
    }
    if (shift()=='|' && inner
      || next_control==begin_comment || next_control==end_comment) return;
  }
}

void app_str(char* s) { while(*s!='\0') app(*s++); }

text_pointer C_translate(void)
{ text_pointer p; scrap_pointer save_base=scrap_base;
  scrap_base=scrap_ptr;
  C_read(true); /* get the scraps together */
  if (next_control != '|') err_print("! Missing `|' after C text");
				      
  p=translate(); /* make the translation */
#ifdef STAT
  if (scrap_ptr>max_scr_ptr) max_scr_ptr=scrap_ptr;
#endif
  scrap_ptr=scrap_base; scrap_base=save_base; return p;
}

void outer_read (void) /* makes scraps from \Cee\ tokens and comments */
{ while (next_control<format)
    if (next_control!=begin_comment) C_read(false);
    else 
         { boolean one_liner=loc[-1]=='/'; int bal=0; /* brace level in comment */
           typedef_tracking(false);
           check_scrap(); check_toks(4);
         #if 0
           if (scrap_ptr==scrap_base || scrap_ptr[-1].cat!=insert)
             app(cancel);
         #endif
           app_str(one_liner ? "\\SHC{" : "\\C{");  
           while ((next_control=scan_comment(&bal,one_liner))=='|')
           { text_pointer p=text_ptr, q=(freeze_text(), C_translate());
             check_toks(7);
             app_tok(text_flag+text_index(p)); /* initial text */
             if (compatibility_mode) app_str("\\PB{"); 
             app(inner_text_flag+text_index(q)); /* text from `\pb' */
             if (compatibility_mode) app('}');
             check_text();
           }
           app_char_tok('}'); app(force); pack_scrap(insert, no_math);
         	    /* the full comment becomes a scrap */
           typedef_tracking(true);
         }
  check_scrap(); check_toks(11); /* `\.{\$\\4$m$\\PE\{\}\$$f$}' */
}

void do_C (void) /* read, translate, and output \Cee~text in `\pb' */
{ enter_block(1);
  if (compatibility_mode) out_str("\\PB{"); 
  make_output(C_translate(),inner); /* output the list */
  if (compatibility_mode) out('}');
#ifdef STAT
  if (text_ptr>max_text_ptr) max_text_ptr = text_ptr;
  if (tok_ptr>max_tok_ptr) max_tok_ptr = tok_ptr;
#endif
  leave_block(1); /* forget the tokens */
}

void finish_C (void)
{ out_str ("\\B"); 
  make_output(translate(),outer);
  out_str("\\par"); finish_line();
#ifdef STAT
  if (text_ptr>max_text_ptr) max_text_ptr=text_ptr;
  if (tok_ptr>max_tok_ptr) max_tok_ptr=tok_ptr;
  if (scrap_ptr>max_scr_ptr) max_scr_ptr=scrap_ptr;
#endif
  leave_block(0); scrap_ptr=scrap_info;
	/* forget the tokens and the scraps */
}

void footnote (xref_pointer* p,sixteen_bits flag)
{ if ((*p)->num<=flag) return;
  finish_line(); out('\\');
  out(flag==0 ? 'U' : flag==cite_flag ? 'Q' : 'A');   
  *p=list_refs(*p,flag);
  out('.');
}

xref_pointer list_refs (xref_pointer x,sixteen_bits flag)
{ xref_pointer q=next_xref(x); /* second element in \xr. list */
  if (q->num>flag) out('s');
     /* use `\.{\\As}', `\.{\\Qs}' or `\.{\\Us}' */
       
  out(' ');
  while (out_sec_nr(x->num&num_mask),x=next_xref(x),x->num>flag)
    if (next_xref(x)->num>flag) out_str(", "); /* |x| is not the last */
    else
    { out_str("\\ET"); /* next number printed will be the last */
      if (x!=q) out('s'); /* `\.{\\ETs}' for the last of more than two */
    }  
  return x;
}


#ifdef DEBUG
void print_cat (int c) /* symbolic printout of a category */
{ static char* cat_name[]=
  { "unop", "binop", "op", "select"
  , "?", "{", "}", "(", ")", "[", "]", ",", ";", ":", "::", "@;"
  , "subscr", "struct_head", "short_{", "short_struct_head"
  , "cmp_stmt", "stmt"
  , "function", "function_head", "params", "label"
  , "if_head", "if_else_head", "do_head"
  , "mod_name", "declarator", "decl", "exp"
  , "for", "do", "if", "else", "extern"
  , "throw", "try", "catch,"
  , "int", "case", "sizeof", "struct", "return"
  , "template", "<", ">", "templ_params"
  , "#{", "#}", "insert", "@[", "@]"
  };
  if (c<=max_category && c>0) printf("%s",cat_name[c-1]);
  else printf ("IMPOSSIBLE");
}
#endif 

 trie_node *get_new_trie_node(void)
{ if (node_no>=max_no_of_nodes)
    overflow("trie node"); 
  trie_node_mem[node_no].rule=NULL;
  return &trie_node_mem[node_no++];
}

void install_rule(reduction *rule)
{ if ((rule_mask & rule->mask)==0)
  { eight_bits* p=rule->lhs.category, i=0;
    while (i<max_lhs_length && p[i]!=0) ++i;
    rule->lhs.length=i;
#ifdef DEBUG
    
    { if (rule->lhs.length<=abs(rule->lhs.context))
        rule_error("\nNo scraps to replace in rule %d.\n", rule->id);
    		
      for(i=0; i<rule->lhs.length; ++i)
        if (!valid_cat(p[i]))
          rule_error("\nUnknown category %d in LHS of rule %d.\n", p[i], rule->id);
    		  
    }
    
    { int c=rule->rhs.category; char* s=rule->rhs.translation;
      if (!valid_cat(c))
          rule_error("\nUnknown category %d in RHS of rule %d.\n", c, rule->id);
    		  
      if (s!=NULL)
      { if (*s=='\0') s=rule->rhs.translation=NULL; /* replace empty string */
        else
        { i=0;
          do
    	if (*s!='p') i+= *s++=='_'; /* count underscores */
    	else if (++s,isdigit((eight_bits)*s)) ++s; /* skip digit and advance */
    	else rule_error("\nDigit should follow 'p' in format of rule %d.\n"
    			 	, rule->id);
          while (*s!='\0');
          if (i!=rule->lhs.length-abs(rule->lhs.context))
    	rule_error("\nCount of '_' not equal to length LHS in rule %d.\n"
    		     , rule->id);
        }
      }
    }
#endif
    
    { trie_node* q=trie_root;
      for (i=0; i<rule->lhs.length; ++i)
      { if (no_successor(q,p[i])) set_successor(q,p[i],get_new_trie_node());
        q=successor(q,p[i]);
      }
    #ifdef DEBUG
      if (q->rule!=NULL)
        rule_error("\nIdentical left-hand sides in rules %d and %d.\n"
    	           , q->rule->id, rule->id);
    #endif
      q->rule=rule;
    }
    
    { int k=rule->lhs.context,d;
      if (k<0) { rule->lhs.length+=k; k=0; } else rule->lhs.length-=k;
      d=k-(max_lhs_length-1); /* this cannot be positive */
      if (rule->lhs.category[k]==rule->rhs.category) /* no category change */
      { ++d;
    #ifdef DEBUG
        if (rule->lhs.length==1)
          rule_error("\nNo categories change in rule %d.\n", rule->id);
    		  
    #endif
      }
      rule->lhs.context=k;
      rule->displacement=d; /* if positive, an error was reported */
    }
  }
}

reduction *match (scrap_pointer p)
{ trie_node* q=trie_root; reduction* rule=NULL; int c;
  while (c=p++->cat,valid_cat(c) && !no_successor(q,c))
    if ((q=successor(q,c))->rule!=NULL) rule=q->rule;
  return rule;
}

void fuse (scrap_pointer s, int n)
{ int cur_mathness=maybe_math, init_mathness=maybe_math; scrap_pointer p=s;
  check_toks(n); check_text();
  do add_trans(p++) while (--n>0); /* gather all the translations */
  s->trans=text_ptr; freeze_text();
  s->mathness=(init_mathness<<2)+cur_mathness;
}


void make_nonreserved (scrap_pointer p)
{ text_pointer q=p->trans; token t;
  while (text_flag<=(t=text_begin(q)[0]) && t<inner_text_flag)
    q=text_at(t-text_flag);
  if (text_end(q)==text_begin(q)+1 && res_flag<=t && t<mod_flag)
    text_begin(q)[0] = t-res_flag+id_flag;
}

id_pointer first_ident(text_pointer p)
{ token_pointer q; token t;
  if (p>=text_ptr) confusion("first_ident"); 
  for (q=text_begin(p); q<text_end(p); ++q)
    if (id_flag<=(t=*q) && t<mod_flag) return id_at(t%id_flag);
    else if (text_flag<=t) /* text or inner text */
    { id_pointer r=first_ident(text_at(t%id_flag));
      if (r!=NULL) return r;
    }
  return NULL;
}

void make_underlined (scrap_pointer p)
  /* underline entry for first identifier in |p->trans| */
{ id_pointer name=first_ident(p->trans); /* name of first identifier */
  if (name==NULL) return;
    /* this happens for unsyntactical things like `|int 3;|' */
  { sixteen_bits head=xref_index(name->xref),* r=&head; int n;
    while ((n=xnum(*r)&num_mask)!=0 && n<section_count) r=&xlink(*r);
    if (n==section_count) xnum(*r)|=def_flag;
    else /* this may happen for one-letter identifiers */
    { make_xref(section_count+def_flag,*r);
      if (r==&head) name->xref=xref_ptr;  else *r=xref_index(xref_ptr);
    }
  }
}

void reduce (reduction* rule)
{ int k=rule->lhs.context, l=rule->lhs.length;
  scrap_pointer s = pp+k, p=s;/* position of the new scrap */
  char* f=rule->rhs.translation; /* format string for translation */

  s->cat=rule->rhs.category;

  if (l>1 || f!=NULL) /* otherwise ready already */
  { if (f==NULL) fuse(s,l), p+=l; /* default translation */
    else 
         { int cur_mathness=maybe_math, init_mathness=maybe_math;
           check_toks(23); check_text();
           do
             switch (*f++) 
             { case '+': app(indent); break;
               case '-': app(outdent); break;
               case 'p': app(opt); app(*f++); break; /* penalty with numeric argument */
               case 'f': set_mode(no_math); app(force); break;
               case 'F': set_mode(no_math); app(big_force); break;
               case 'b': set_mode(no_math); app(backup); break;
               case 'B': set_mode(no_math); app(break_space); break;
               case 't': set_mode(yes_math); app_str("\\a"); break; 
         					/* next item in tab space */
               case ',': set_mode(yes_math); app_str("\\,"); break; 
         					/* thin space */
               case 'h': set_mode(no_math); break; /* force horizontal mode */
               case 'm': set_mode(yes_math); app(' '); break;
         				/* force math mode, avoid `\.{\$\$}' */
               case 'o': set_mode(yes_math); app_str("\\m"); break; 
         					/* make ``mathord'' */
               case 'r': set_mode(yes_math); app_str("\\MRL"); break; 
         					/* make ``mathrel'' */
               case 'a': set_mode(yes_math); app_str("\\ang"); break; 
         		         /* change `\.<' or `\.>' to `$\ang<$' or~$\ang>$'  */
               case '!': make_underlined(p); break;
               case '$': make_nonreserved(p); break; 
               case ' ': set_mode(no_math); app(' '); break;
               case '~': case '{': case '}': app(f[-1]); break;
                        /* insert character literally */
               default: printf("%c: ",f[-1]);
                 confusion("illegal character in format string");
               case '_':	add_trans(p++);
             }
           while (*f!='\0');
           s->trans=text_ptr; freeze_text();
           s->mathness=(init_mathness<<2)+cur_mathness;
         }
    if (l>1) 
            { scrap_pointer q=s+1; /* position after the newly formed scrap */
              while (p<lo_ptr) *q++=*p++;
              lo_ptr=q;
            }
  }

  
  #ifdef DEBUG
  { scrap_pointer k; /* pointer into |scrap_info| */
    if (tracing>=trace2)
    { print("\n%3d:", rule->id);
      for (k=scrap_base; k<lo_ptr; k++)
      { putchar (' ');
        if (tracing==trace3) putchar(math_char(left_math(k)));
        if (k==s) putchar('>'), print_cat(k->cat), putchar('<');
        else print_cat(k->cat);
        if (tracing==trace3) putchar(math_char(right_math(k)));
      }
      print("%s\n", hi_ptr<scrap_ptr ? " ..." : ".");
    }
  }
  #endif
  
  if (pp<scrap_base-rule->displacement) pp=scrap_base; 
  else pp+=rule->displacement;
}

text_pointer translate (void) /* converts a sequence of scraps */
{ pp=lo_ptr=hi_ptr=scrap_base;
  if (scrap_ptr==pp || dangling_tokens()) /* then append dummy scrap */
  { check_scrap(); pack_scrap(insert,no_math); }
  
  #ifdef DEBUG
  { if (tracing>=trace2)
    { print("\nTracing after l.%d:\n", cur_line); 
      if (loc>buffer+50) /* shorten long lines to keep within margins */
      { printf("..."); term_write (loc-50,50); }
      else term_write(buffer,loc-buffer);
      new_line(); /* |term_line_empty| is still valid */
    }
  }
  #endif
  
  do
  { reduction *rule;
    
    { scrap_pointer lo_min = pp+max_lhs_length;
      while (lo_ptr<lo_min && lo_ptr->cat!=0)
        if (hi_ptr>=scrap_ptr) lo_ptr->cat=0;
        else
        { *lo_ptr++ = *hi_ptr++;
          while (hi_ptr<scrap_ptr && hi_ptr->cat==insert)
    	{ *lo_ptr = *hi_ptr++; fuse(lo_ptr-1,2); }
        }
    }
    if ((rule=match(pp))!=NULL) reduce(rule);
    else
    { ++pp;
      
      if (pp->cat==end_expr || pp->cat==rproc)
      { int start=pp->cat-1; /* the opening category matching |pp->cat| */
        scrap_pointer s=pp, p=pp+1;
        while ((--s)->cat!=start && s>scrap_base) {}
        if (s->cat==start) /* if opening symbol is missing, take no action */
        { if (start==begin_expr) s->cat=expression;
          else if (s==scrap_base) s->cat=insert; 
          else --s; /* position of new scrap */
          fuse(s,(int)(p-s));
          
          { scrap_pointer q=s+1; /* position after the newly formed scrap */
            while (p<lo_ptr) *q++=*p++;
            lo_ptr=q;
          }
       /* using values of |p| and |s| */
          pp= s-scrap_base<max_lhs_length ? scrap_base : s+1-max_lhs_length;
        }
      }
    }
  }
  while (pp<lo_ptr);
  
  { scrap_pointer j;
    if (scrap_base->cat==insert && lo_ptr>scrap_base+1)
    { fuse(scrap_base,2); /* merge initial |insert| into the next scrap */
      j=scrap_base; j->cat=j[1].cat; --lo_ptr; while (++j<lo_ptr) *j=j[1];
    }
    
    #ifdef DEBUG
    { if (tracing==trace1 && lo_ptr>=scrap_base+2)
      { print("\nIrreducible scrap sequence at line %d in section %d:\n"
    	     ,cur_line, section_count);
        mark_harmless();
        for (j=scrap_base; j<lo_ptr-1; j++)
          print_cat(j->cat), putchar(' ');
        print_cat(j->cat); new_line(); /* |term_line_empty| is still valid */
      }
    }
    #endif
    check_toks(1);
    for (j=scrap_base; j<lo_ptr; j++)
    { if (j!=scrap_base) app_char_tok(' ');
      if (left_math(j)==yes_math) app('$');
      app_trans(j);
      if (right_math(j)==yes_math) app('$');
    }
    freeze_text();
    return text_ptr-1;
  }
}

void push_level (text_pointer p) /* suspends the current level */
{ if (stack_ptr==stack_end) overflow("stack"); 
  *stack_ptr++=cur_state;
#ifdef STAT
  if (stack_ptr>max_stack_ptr) max_stack_ptr=stack_ptr;
#endif
  cur_tok=text_begin(p); cur_end=text_end(p);
}

void make_output(text_pointer t,mode m) /* output a complete text */
{ int save_next_control=next_control; id_pointer save_cur_id=cur_id;
  stack_pointer stack_bot=stack_ptr; token state=cancel;
  push_level(t); cur_mode=m;
  do
    if (cur_tok==cur_end) pop_level();
    else
    { token a= *cur_tok % id_flag;
      switch (*cur_tok++/id_flag)
      {
      
      case 1: 
              { if (state>=space)
                  if (state<break_space) out(state==space ? ' ' : '~');
                  else
                  { out('\\');
                    if (state<backup) out(state-break_space+'5');
              	/* `\.{\\5}', `\.{\\6}', or `\.{\\7}' */   
                    else out(state-backup+'6'),out_str("\\4");
                      /* `\.{\\6\\4}' or `\.{\\7\\4}' */   
                    finish_line();
                  }
                state=0;
              }
       out_identifier(id_at(a)); break;
      case 2: 
              { if (state>=space)
                  if (state<break_space) out(state==space ? ' ' : '~');
                  else
                  { out('\\');
                    if (state<backup) out(state-break_space+'5');
              	/* `\.{\\5}', `\.{\\6}', or `\.{\\7}' */   
                    else out(state-backup+'6'),out_str("\\4");
                      /* `\.{\\6\\4}' or `\.{\\7\\4}' */   
                    finish_line();
                  }
                state=0;
              }
       out_keyword(id_at(a)); break;
      case 3: 
              { if (state>=space)
                  if (state<break_space) out(state==space ? ' ' : '~');
                  else
                  { out('\\');
                    if (state<backup) out(state-break_space+'5');
              	/* `\.{\\5}', `\.{\\6}', or `\.{\\7}' */   
                    else out(state-backup+'6'),out_str("\\4");
                      /* `\.{\\6\\4}' or `\.{\\7\\4}' */   
                    finish_line();
                  }
                state=0;
              }
       out_module_name(mod_at(a)); break;
        case 4: push_level(text_at(a)); break;
	case 5: push_level(text_at(a)); cur_mode=inner; break;
	case 0: 
	        { switch (a)
	          { case relax: 
	                        { if (state>=space)
	                            if (state<break_space) out(state==space ? ' ' : '~');
	                            else
	                            { out('\\');
	                              if (state<backup) out(state-break_space+'5');
	                        	/* `\.{\\5}', `\.{\\6}', or `\.{\\7}' */   
	                              else out(state-backup+'6'),out_str("\\4");
	                                /* `\.{\\6\\4}' or `\.{\\7\\4}' */   
	                              finish_line();
	                            }
	                          state=0;
	                        }
	         break;
	            case cancel: state=cancel; break;
	            case indent: case outdent:
	              if (cur_mode==outer) { out('\\'); out(a-indent+'1'); }
	              break;    
	            case opt:
	            { int digit=*cur_tok++;
	              if (state==0)
	              { out('\\'); out(cur_mode==outer ? '3' : '0'); out(digit); }
	        					 
	              break;
	            }
	            case flush_left:
	              if (cur_mode==outer)
	              { 
	                { if (state>=space)
	                    if (state<break_space) out(state==space ? ' ' : '~');
	                    else
	                    { out('\\');
	                      if (state<backup) out(state-break_space+'5');
	                	/* `\.{\\5}', `\.{\\6}', or `\.{\\7}' */   
	                      else out(state-backup+'6'),out_str("\\4");
	                        /* `\.{\\6\\4}' or `\.{\\7\\4}' */   
	                      finish_line();
	                    }
	                  state=0;
	                }
	         out_str("\\8"); }  
	              break;
	            case big_force: case backup:
	              if (a+state==big_force+backup) a=big_backup; /* fall through */
	            case break_space: case force: case big_backup:
	              if (cur_mode==inner) a=space;
	            up_state:
	              if (state!=cancel && state<a) state=a;
	              break;
	            case ' ': case '~':
	              if (cur_mode==inner) { a= a==' ' ? space : tilde; goto up_state; }
	              if (state==cancel || state>=break_space) break;
	                /* else fall through */
	            default: 
	                     { if (state>=space)
	                         if (state<break_space) out(state==space ? ' ' : '~');
	                         else
	                         { out('\\');
	                           if (state<backup) out(state-break_space+'5');
	                     	/* `\.{\\5}', `\.{\\6}', or `\.{\\7}' */   
	                           else out(state-backup+'6'),out_str("\\4");
	                             /* `\.{\\6\\4}' or `\.{\\7\\4}' */   
	                           finish_line();
	                         }
	                       state=0;
	                     }
	         out(a);
	          }
	        }
      }
    }
  while(stack_ptr>stack_bot);
  
  { if (cur_mode==outer && (state==big_force || state==big_backup))
      out_str("\\Y");
  }
  next_control=save_next_control; cur_id=save_cur_id;
}

void out_sec_nr (int n) /* output a section number */
{ char s[6];
  sprintf(s,"%d",n); out_str(s);
  if (section_changed(n)) out_str ("\\*"); 
}

void out_id_part (char* s, int l)
{ boolean b=l!=1;
  if (b) out ('{');
  while (--l>=0) { if (*s=='_') out ('\\'); out (*s++); }
  if (b) out('}');
}

void out_index (id_pointer p)
{ char* s=name_begin(p); boolean b= s[1]!='\0';
  if (compatibility_mode) { out_id_full(p); return; }
  if (b) out ('{');
  out_str(s);
  if (b) out('}');
}

void out_keyword(id_pointer p)
{ out_str("\\&"); out_id_full(p); }
 

void out_identifier (id_pointer p)
{ int k=length(p); eight_bits* ch=(eight_bits*)name_begin(p);
enum { ord, indexed, caps, single, indexed_single } kind;
  if (p->ilk==TeX_like || p->ilk==NULL_like)  
  {
    
    { if (p->ilk==TeX_like) out_str("\\\\{");
      out('\\');  do out(*ch=='_' ? 'x' : *ch); while (*++ch!='\0');
      if (p->ilk==TeX_like) out('}');
    }
    return;
  }
  if (!compatibility_mode) /* then search for possibly trailing digits */
  { do --k; while (isdigit((eight_bits)ch[k]));
    /* terminates because |!isdigit(ch[0])| */
    ++k; /* point to end of identifier without its index (if any) */
  }
  
  if (k==1) { out(' '); kind= length(p)==1 ? single : indexed_single; }
  else
  { int i=k; 
    while (--i>=0)  if (!isupper(ch[i])&&!isdigit(ch[i])&&ch[i]!='_') break;
    kind= i<0 ? caps : k<length(p) ? indexed : ord;
  out('\\'); out(kind==caps ? '.' : '\\');	  
  }
  if (kind==indexed || kind==indexed_single)
  { out_id_part(name_begin(p),k); /* main part */
    out ('_'); out_id_part(name_begin(p)+k,length(p)-k); /* subscript */
  }
  else out_id_full (p);
}

xref_pointer out_module_name(mod_pointer name)
{ xref_pointer x=name->xref; boolean file_module= x->num==file_flag;
  if (file_module) x=next_xref(x);
  out_str ("\\X"); 
  if (x->num>=def_flag)
  { out_sec_nr(x->num-def_flag); /* output the defining section number */
    if (phase==3) /* all of them in Phase III */
      while (x=next_xref(x), x->num>=def_flag)
	out_str (", "), out_sec_nr(x->num-def_flag);
  }
  else out ('0'); /* section number `0' means `nowhere defined' */
  out (':'); if (file_module) out_str("\\.{"); 
  
  { char* k=name_begin(name),c;
    while ((c=*k++)!='\0')
    { if (file_module) { if (strchr(" \\#%$^{}~&_",c)!=NULL) out ('\\'); }
      if (c=='@' && *k++!='@') 
                { print("\n! Illegal control code in module name"); print_mod(name);
                	     
                  mark_error();
                }
      if (file_module || c!='|') out(c);
      else
      { char* save_loc=loc, *save_limit=limit;
        
        { char delimiter='\0'; /* |'"'| or |'\''|, or |'\0'| outside of strings */
          next_control=*limit++='|'; loc=limit;
          do
            if ((c=*k++)=='\0') 
                      { print("\n! C text in module name didn't end"); print_mod(name);
                      	   
                        mark_error();
                        if (delimiter!='\0') *limit++=delimiter;
                        *limit++='|';
                        if (limit>&buffer[long_buf_size-2])
                          fatal("fix that first, you sneaky devil");
                               
                        break;
                      }
            else
            { *limit++=c;
            
            { if (c=='@') /* control code; now |*k!='\0'| */
              { if ((*limit++=*k++)=='\'' && delimiter=='\0')	/* copy code, test for \:' */
                  delimiter='\''; /* which behaves like `\.'' */
              }
              else if (c=='\\' && delimiter!='\0')
              { char d=*limit++=*k++; /* escaped character, possibly |delimiter| */
                if (d=='\0') --k,limit-=2; /* remove backslash, error is issued anyway */
              }
              else if (c=='\'' || c=='"')
                if (delimiter=='\0') delimiter=c;
                else if (delimiter==c) delimiter='\0';
            }
            }
          while (c!='|' || delimiter!='\0');
          *limit=' ';
        }
        do_C(); loc=save_loc; *(limit=save_limit)=' ';
      }
    }
  }
  if (file_module) out_str ("}");
  out_str ("\\X");
  return x;
}

void flush_buffer(char* b, boolean percent)
   /* output from |out_line| to |b|, where |b<=out_ptr| */
{ int j=(int)(b-out_line); /* number of characters to be output */
  if (!percent)  while (j>0 && out_line[j-1]==' ') --j;
    /* remove trailing blanks */
  fprintf(tex_file, "%.*s",j,out_line);
  if (percent) tex_putc('%');
  tex_new_line();
  { char* p=out_line;
    while (b<out_ptr) *p++=*b++; /* shift back remainder of line */
    out_ptr=p; /* adjust to end of (possibly empty) shifted part */
  }
}

void finish_line(void) /* do this at the end of a line */
{ if (!output_line_empty()) flush_buffer(out_ptr, false);
  else if (limit==buffer) tex_new_line(); /* copy blank input line */
}

void out_str (char* s) { while (*s!='\0') out (*s++); }

void break_out (void) /* finds a way to break the output line */
{ char* k=out_ptr,c; int count=0; /* number of backslashes seen */
  do  if ((c=*--k)==' ') goto found;  while (c!='\\');
  do ++count; while ((c=*--k)=='\\');
found:
  if (++k>out_line) flush_buffer(k,c!=' ');  else 
                                          if (count==0) flush_buffer(out_ptr-1,true);
                                          else if (count>=2) flush_buffer
                                            (&out_line[count&=~1]==out_buf_end ? out_buf_end-2 : &out_line[count],true);
                                          else
                                          { print("\n! Line had to be broken (output l.%d):\n",out_line_nr);
                                          	   
                                            term_write(out_line,out_ptr-out_line); new_line(); mark_harmless();
                                            flush_buffer(out_ptr,false);
                                          }
}

void phase_three (void) /* output the \xr. index */
{ finish_line();
  if (no_xref) out_str("\\endcodemode\\vfill\\end");
  else
  { phase=3; print_progress("\nWriting the index...");
			     
    typedef_tracking(false); /* during parse of `\pb' in module names */
    if (change_exists)
    { 
      { int k=0; boolean first=true;
        out_str("\\ch ");  /* changes */
        while (k<section_count)
        { do ++k; while (!section_changed(k));
          if (first) first=false;  else out_str(", ");
          out_sec_nr(k);
        }
        out('.');
      }
 finish_line(); }
    if (triple_file_output)
      
      { out_str("\\inx \\input \\jobname.idx\n" 
      	  "\\fin \\input \\jobname.scn\n" 
      	  "\\con");
        if (even_out_pages) out_str("even");
        finish_line(); fclose(tex_file);
        if ((tex_file=fopen(idx_file_name,"w"))==NULL)
          fatal("! Cannot open \"%s\" as output file",idx_file_name);
      }
    else { out_str("\\inx"); finish_line(); }   /* index */
    
    { 
      { id_pointer name;
        eight_bits c=UCHAR_MAX;
        id_pointer *h; /* pointer into |hash| */
      
        do bucket[c]=NULL; while (c--!=0);
        for (h=hash; h<hash_end; h++)
          for (name=*h; name!=NULL; name=name->hash_link)
      	/* traverse all hash lists */
            if (name->ilk!=reference && name->ilk!=header_file_name 
               && name->xref->num!=0)
            	/* leave out non-identifiers and unreferenced names */
            { c=name_begin(name)[0]; c=tolower(c);
      	ilink(name)=bucket[c]; bucket[c]=name;
            }
      }
       /* the first time, entries do not come from |sort_info| */
      unbucket(1); /* pick up first-order bucketed lists */
      while (sort_ptr>sort_info) /* i.e., the stack is not empty */
      { eight_bits depth=(--sort_ptr)->depth;
        id_pointer name=sort_ptr->head;
        if (ilink(name)==NULL || depth==infinity)
            /* singleton or set of look-alikes */
          
          do
          { out_str("\\@");   
            
            switch (name->ilk)
            { case normal: case NULL_like: out('m'); out_identifier(name); break;
              case TeX_like: out('h'); out_identifier(name); break;
              case roman: out('h'); out_index(name); break;
              case wildcard: out_str("h\\9"); out_index(name); break;   
              case typewriter: out_str("h\\."); out_index(name); break;   
              default: out('h'); out_keyword(name);
            }
            
            { xref_pointer x=name->xref;
              do
              { sixteen_bits n=x->num;
                out_str(", ");
                if (n<def_flag) out_sec_nr(n);
                else { out_str("\\["); out_sec_nr(n-def_flag); out(']'); }   
              } while ((x=next_xref(x))->num!=0);
              out('.'); finish_line();
            }
          } while ((name=ilink(name))!= NULL);
        else
        { 
          do
          { eight_bits c= tolower((eight_bits)name_begin(name)[depth]);
            id_pointer next_name=ilink(name); /* save link */
            ilink(name)=bucket[c]; bucket[c]=name; name=next_name;
              /* put into bucket */
          }
          while (name!=NULL);
          unbucket(depth+1);
        }
      }
    }
    if (triple_file_output)
      
      { finish_line(); fclose(tex_file);
        if ((tex_file=fopen(scn_file_name,"w"))==NULL)
          fatal("! Cannot open \"%s\" as output file",scn_file_name);
      }
    else { out_str("\\fin"), finish_line(); }  /* end of index */
    
     list_modules(root);
    if (!triple_file_output)
    { out_str("\\con");  /* table of contents */
      if (even_out_pages) out_str("even"); 
    }
  }
  finish_line(); fclose(tex_file);
  print_progress("\nDone.\n");
  check_complete(); /* was all of the change file used? */
}

void unbucket (eight_bits d) /* empties buckets having depth |d| */
{ int i=end_collate; /* index into |collate| */
  while(--i>=0)	 if (bucket[collate[i]]!=NULL)
  { if (sort_ptr>=sort_info_end)
      overflow("sorting"); 
    sort_ptr->depth= i==0 ? infinity : d;
      /* |infinity| means there is nothing left to compare */
    sort_ptr++->head=bucket[collate[i]];
    bucket[collate[i]]=NULL; /* push and empty bucket */
#ifdef STAT
    if (sort_ptr>max_sort_ptr) max_sort_ptr=sort_ptr;
#endif
  }
}

void list_modules (mod_pointer p) /* print all module names in subtree |p| */
{ if (p != NULL)
  { list_modules(p->llink);
  out_str("\\@$");   
    leave_block(0); scrap_ptr=scrap_info; /* get ready for parsing */
    { xref_pointer x=out_module_name(p); out('$');
      footnote(&x,cite_flag); footnote(&x,0);
    }
    finish_line();
  list_modules(p->rlink);
  }
}

#ifdef STAT
void print_stats()
{ print("\nMemory usage statistics:\n"); 
report("identifier",		id_index(id_ptr),	max_idents);
report("module name",		mod_index(mod_ptr),	max_modules);
report("byte",		byte_ptr-byte_mem,	max_bytes);
report("cross-reference",	xref_ptr-xmem,		max_refs-1);
printf("Parsing:\n");
report("scrap",		max_scr_ptr-scrap_info,	max_scraps);
report("text",		max_text_ptr-text_mem,	max_texts);
report("token",		max_tok_ptr-tok_mem,	max_toks);
report("trie node",		node_no,		max_no_of_nodes);
report("level",		max_stack_ptr-stack,	stack_size);
printf("Sorting:\n");
report("level",		max_sort_ptr-sort_info,	sort_stack_size);
}
#endif

