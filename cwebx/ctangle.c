#define version_string  "x3.6"
#define banner  "This is CTANGLE (version "version_string")"
#define max_toks  150000L
#define max_texts  2500
#define max_files  50
#define stack_size_max  50
#define max_indent  1000
#define variant  text
#define line_output  flags['l']
#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  "common.h"
#define max_bytes  250000L
#define max_modules  1000
#define max_idents  10000
#define max_sections  4000
#define hash_size  353
#define buf_size  500
#define longest_name  1000
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
#define overflow(t) fatal("\n! Sorry, %s capacity exceeded",t)
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
#define tok_begin(p) (p)->tok_start
#define tok_end(p) ((p)+1)->tok_start
#define text_table_end  (&text_table[max_texts])
#define tok_mem_end  (&tok_mem[max_toks])
#define store_byte(c)  \
if (tok_ptr==tok_mem_end) overflow("token");  else *tok_ptr++=c
#define macro_flag  (text_table_end-1)
#define header_flag  text_table_end
#define next_sec(m) ((m)->text_link)
#define equiv  equiv_or_xref
#define verb_quote  0x2
#define join  0x3
#define cur_repl  cur_state.repl_field
#define cur_byte  cur_state.byte_field
#define cur_end  cur_state.end_field
#define cur_sec  cur_state.sec_nr_field
#define cur_ind  cur_state.indent_field
#define cur_state  stack[0]
#define stack_end  (&stack[stack_size_max])
#define stack_empty() (stack_ptr==&stack[0])
#define C_printf(format,x) fprintf(C_file,format,x)
#define C_putc(c) putc(c,C_file)
#define put_indent() \
   (indent_buffer[ind_i=cur_ind]='\0',C_printf("%s",indent_buffer))
#define append_white(c) \
  if (ind_i>=max_indent) overflow("indent buffer"); \
  else indent_buffer[ind_i++]= isspace(c) ? c : ' '
#define trans_limit  9
#define trans_of(c) c_trans[(unsigned char)(c)-0x80]
#define translation_exists(c) (trans_of(c)[0]!='\0')
#define comp_op(op) \
  (C_printf(out_state==operator && line_output ? " %s" : "%s",op) \
  ,out_state=operator)
#define code_of(c) ccode[(unsigned char)(c)]
#define compress(char2,code)  \
  if (*loc==char2) return ++loc,code
#define preproc_directive   0
#define section_body  1
#define report(k,c,m) \
  printf("%lu %ss (out of %lu)\n",(unsigned long)(c),k,(unsigned long)(m))


boolean names_match (id_pointer,char*,int,int);
void init_id_name (id_pointer,int);
void init_module_name (mod_pointer);


typedef struct text
{ eight_bits* tok_start; /* pointer into |tok_mem| */
  struct text* text_link; /* relates replacement texts */
} text,* text_pointer;

typedef struct
{ text_pointer repl_field; /* replacement text of active section */
  eight_bits *byte_field; /* present location within replacement text */
  eight_bits *end_field; /* ending location of replacement text */
  sixteen_bits sec_nr_field; /* section number */
  sixteen_bits indent_field; /* amount of indentation */
} output_state, * stack_pointer;

enum {identifier=0x80, section_start, section_end, line_mark };

enum { no_space, num_or_id, operator, literal};

enum
{ ignore=0x80, /* control code of no interest to \.{\me.}, or \:> */
  id_code, constant, /* token codes but not control codes */
  verbatim, /* control code for \:= */
  at_sign_image, /* control code for \:@ */
  join_code, /* control code for \:\& */
  ord, /* control code for \:' */
  control_text, /* control code for \:t, \:\^, etc. */
  include_preproc, /* control code for \:p */
  char_trans, /* control code for \:l */
  format, /* control code for \:f */
  definition, /* control code for \:d */
  header, /* control code for \:h */
  begin_C, /* control code for \:c */
  module_name, /* control code for \:< or \:( */
  new_section /* control code for \:\ , \:\~, or \:* */
};


void phase_two (void); /* output the contents of the compressed tables */
void output (text_pointer); /* recursively write out modules */
void output_preproc_directives(void); /* write out all \:d and \:h stuff */
void C_newline (void); /* send a newline to the \Cee~file */
void out_char (eight_bits);

 boolean mark_line(void);

void phase_one (void);



text text_table[max_texts];
text_pointer text_ptr=&text_table[0]; /* first unused position in |text_table| */
eight_bits tok_mem[max_toks];
eight_bits *tok_ptr=&tok_mem[0];

text_pointer text_root=NULL;
text_pointer* last_unnamed=&text_root;

boolean atp_seen=false;

mod_pointer output_file[max_files];
int output_file_count=0;

output_state stack[stack_size_max]; /* info for non-current levels */
stack_pointer stack_ptr;

char indent_buffer[max_indent];
   /* white space characters to produce indentation */
sixteen_bits ind_i;


int cur_val;

eight_bits out_state; /* current status of partial output */
boolean protect;

char c_trans[UCHAR_MAX+1-0x80][trans_limit+1];

eight_bits ccode[UCHAR_MAX+1];

id_pointer cur_id; /* identifier just scanned */
mod_pointer cur_mod;

text_pointer cur_text; /* replacement text formed by |scan_repl| */
eight_bits next_control;


int main (int argc, char** argv)
{ program=ctangle;
  line_output=true;
  common_init(argc,argv,banner);

  tok_begin(text_ptr)=tok_ptr;

  if (compatibility_mode)
  { unsigned char c=UCHAR_MAX;
    do sprintf(trans_of(c),"X%X",c); while (--c>=0x80);
  }

  { unsigned char c=0;
    do ccode[c] = isspace (c) ? new_section : ignore; while(c++!=UCHAR_MAX);
  ccode['v']=ccode['V']='|';
  ccode['=']=verbatim;
  ccode['@']=at_sign_image;
  ccode['&']=join_code;
  ccode['\'']=ord;
  ccode['^']=ccode['?']=ccode['.']=ccode[':']=ccode['#']=
      ccode['t']=ccode['T']=ccode['q']=ccode['Q']=control_text;
  ccode['p']=ccode['P']=include_preproc;
  ccode['l']=ccode['L']=char_trans;
  ccode['f']=ccode['F']=ccode['s']=ccode['S']=format;
  ccode['d']=ccode['D']=definition;
  ccode['h']=ccode['H']=header;
  ccode['c']=ccode['C']=begin_C;
  ccode['<']=ccode['(']=module_name;
  ccode['~']=ccode['*']=new_section;
  if (compatibility_mode)

    { ccode['h']=ccode['H']=include_preproc; /* \:h means \:p */
      ccode['p']=ccode['P']=begin_C; /* \:p means \:c */
      ccode['#']=ignore; /* \:\# means \:) */
    }
  }
  if (show_banner)
    print("%s, in %s mode.\n",banner,C_plus_plus ? "C++" : "C");
    /* print a ``banner line'' */
  phase_one(); /* read all the user's text and compress it into |tok_mem| */
  phase_two(); /* output the contents of the compressed tables */
  wrap_up(); /* and exit gracefully */
  return 0;
}

void store_two_bytes (sixteen_bits x)
{ if (tok_ptr+2>tok_mem_end) overflow("token");
  *tok_ptr++ = x >> 8; /* store high byte */
  *tok_ptr++ = x & 0xFF; /* store low byte */
}

boolean names_match (id_pointer x,char* q,int l,int dummy)
{ char* p=name_begin(x); while (--l>=0) if (*p++!=*q++) return false;
  return *p=='\0';
}

void init_module_name(mod_pointer node)
{ node->equiv=NULL; }

void init_id_name (id_pointer dummy,int ilk)  {}

 void phase_two (void)
{ phase=2;
  if (text_root==NULL && output_file_count==0)
  { print("\n! No program text was specified."); mark_harmless(); }

  else
  { if (show_progress)
    { print("\nWriting the output file%s"
	   ,(text_root!=NULL)+output_file_count>1 ? "s" : "");
      if (text_root!=NULL) printf(" (%s):",C_file_name);
      update_terminal();
    }
    if (text_root==NULL) C_file=NULL;
    else
    { open_output_file(); cur_line=1;
      if (!atp_seen) output_preproc_directives();
      output(text_root);
    }

    { int i;
      char output_file_name[longest_name+1]; /* name of the file */
      for (i=0; i<output_file_count; i++)
      { mod_pointer output_module=output_file[i];

        { char* p=output_file_name, *q=name_begin(output_module);
          while (*q!='\0')
            if ((*p++=*q++)=='@')
              if (*q++!='@')
              { print("\n! Illegal control code in file name");

        	print_mod(output_module); err_print("");
              }
          *p='\0';
        }
        if (C_file!=NULL) fclose(C_file);
        if (output_module->equiv==NULL)
        { print("\n! Module not present");

          print_mod(output_module); err_print("");
        }
        else if ((C_file=fopen(output_file_name,"w"))==NULL)
        { print("\n! Cannot open \"%s\" as output file",output_file_name);

          err_print("");
        }
        else
        { if (show_progress) print("\n(%s):",output_file_name);
          cur_line=1;
          output(output_module->equiv);
        }
      }
    }
    print_progress("\nDone.\n");
  }
}

void push_level (mod_pointer p) /* suspends the current level */
{ if (stack_ptr==stack_end) overflow("output stack");
  *stack_ptr++=cur_state;
cur_repl=p->equiv;
  cur_byte=tok_begin(cur_repl); cur_end=tok_end(cur_repl);
cur_ind=ind_i;
}

void continue_or_pop_level (void)
{ if (cur_repl->text_link != NULL)
    /* then link to a continuation, staying on the same level */
  { cur_repl=next_sec(cur_repl);
    cur_byte=tok_begin(cur_repl); cur_end=tok_end(cur_repl);
  }
  else if (--stack_ptr > &stack[0]) cur_state=*stack_ptr;
     /* otherwise |stack_empty()| holds */
  ind_i=cur_ind; /* reset indentation to base of current level */
}

void output (text_pointer repl) /* sends tokens to |out_char| */
{ stack_ptr=&stack[1];
cur_repl=repl; cur_byte=tok_begin(cur_repl); cur_end=tok_end(cur_repl);
  cur_ind=ind_i=0;
  do
    if (cur_byte==cur_end)
    { cur_val=cur_sec;
      continue_or_pop_level();
      out_char(section_end); /* output ending section number comment */
    }
    else

    { int a = *cur_byte++;
      if (a<0x80) out_char(a); /* single byte token */
      else if (a>=0xF8)
        if (a<0xFA) out_char(a==0xF8 ? *cur_byte++ : line_mark);
        else { C_newline(); output_preproc_directives(); } /* |a==0xFA| */
      else
      { cur_val=(((a-=0x80)%0x28)<<8)+*cur_byte++;
        switch (a/0x28)
        { case 0: out_char(identifier); break;
          case 1:
                  { mod_pointer mod_name=mod_at(cur_val);
                    if (mod_name->equiv!=NULL) push_level(mod_name);
                    else
                    { print("\n! Module not present"); print_mod(mod_name); err_print("");

                    }
                  }
      break;
          case 2: cur_sec=cur_val; out_char(section_start);
    	  /* set the correct section number and output comment */
        }
      }
    }
  while(!stack_empty());
  C_newline();
}

void output_preproc_directives (void)
{ text_pointer repl,l;
  eight_bits* p, *end;
  protect=true; /* newlines should be preceded by |'\\'| */
  for (repl=&text_table[0]; repl<text_ptr; repl++)
    if ((l=repl->text_link)==macro_flag || l==header_flag)
    { p=tok_begin(repl); end=tok_end(repl);
      C_printf ("#%se ",l==macro_flag ? "defin" : "includ");
      out_state=no_space;
      while (p<end)

      { int a=*p++;
        if (a<0x80) out_char(a); /* single byte token */
        else if (a>=0xF8)
          if (a==0xF8) out_char(*p++);
          else confusion("`@p' within macro");
        else
        { cur_val=(((a-=0x80)%0x28)<<8)+*p++;
          if (a<0x28) out_char(identifier);
          else confusion("module within macro");
        }
      }
      C_newline(); /* this newline is not escaped */
    }
  protect=false;
}

void C_newline (void) /* writes one line to output file */
{ C_putc('\n');
  if (!line_output) put_indent();
  if (cur_line%100==0 && show_progress)
  { if (cur_line%500!=0) print("."); /* progress report */
    else print(cur_line%2500==0 ? "%u\n" : "%u",cur_line);
    update_terminal();
  }
  ++cur_line;
}

void out_char (eight_bits c)
{ if (out_state==literal)
    if (c==verb_quote) out_state=num_or_id; /* close literal mode */
    else C_putc(c); /* write literal character, possibly 8-bit */
  else if (isalnum(c)&&c<0x80 || c=='_')
    /* single character identifier or number */
  { if (out_state==num_or_id && line_output) C_putc(' ');
    C_putc(c); out_state=num_or_id;
    if (!line_output) append_white(c);
  }
  else switch (c)
  { case verb_quote:
      if (out_state==num_or_id && line_output) C_putc(' ');
      out_state=literal; break;
    case join: out_state=no_space; break;
    case '\n':
           { if (protect) { C_putc(' '); C_putc('\\'); }
             C_newline();
             if (out_state!=literal) out_state=no_space;
           }
  break;
    case identifier:
                     { char* p=name_begin(id_at(cur_val)); int l=0;
                       if (out_state==num_or_id && line_output) C_putc(' ');
                       do
                         if ((unsigned char)(*p)<0x80) { C_putc(*p);++l; }
                         else
                         { char* q=trans_of(*p); do { C_putc(*q); ++l; } while (*++q!='\0'); }
                       while (*++p!='\0');
                       out_state=num_or_id;
                       if (!line_output) do append_white(' '); while (--l>0);
                     }
  break;
    case section_start:
      if (line_output) C_printf("/*%d:*/",cur_val);  else C_newline();
      out_state=no_space; break;
    case section_end:
      if (line_output) C_printf("/*:%d*/",cur_val);  else C_newline();
      out_state=no_space; break;
    case line_mark:
                    { sixteen_bits a;
                      a=(*cur_byte++)<<8; a+=*cur_byte++; /* get the line number */
                      C_newline(); C_printf("#line %u \"",a);
                      a=(*cur_byte++)<<8; a+=*cur_byte++; /* get the file name index */
                      C_printf("%s\"",name_begin(id_at(a)));
                      C_newline(); out_state=no_space;
                    }
  break;

    case '+': case '-': case '*': case '/': case '%': case'?':
    case '<': case '>': case '&': case '|':
      if (out_state==operator && line_output) C_putc(' '); /* fall through */
    case '=': C_putc(c); out_state=operator; break;

    case plus_plus: comp_op("++"); break;
    case minus_minus: comp_op("--"); break;
    case minus_gt: comp_op("->"); break;
    case gt_gt: comp_op(">>"); break;
    case eq_eq: comp_op("=="); break;
    case lt_lt: comp_op("<<"); break;
    case gt_eq: comp_op(">="); break;
    case lt_eq: comp_op("<="); break;
    case not_eq: comp_op("!="); break;
    case and_and: comp_op("&&"); break;
    case or_or: comp_op("||"); break;
    default: C_putc(c); out_state=no_space;
      if (!line_output) append_white(c);
  }
}

eight_bits skip_ahead (void) /* skip to next control code */
{ eight_bits c; /* control code found */
  while (find_char())
  { limit[1]='@'; /* place a sentinel */
    while (*loc++!='@') {}
    if (loc<=limit && (c=code_of(*loc++))!=ignore) return c;
  }
  return new_section;
}

boolean skip_comment (boolean one_liner) /* skips over comments */
{ char c; /* current character */
  do
  { if (loc>=limit)
      if (one_liner) return false;
      else if (get_line ()) return true;
      else
      { err_print("! Input ended in mid-comment"); return false; }

    if ((c=*loc++)=='/' && *loc=='*')
      err_print("! `/*' inside comment, did you forget `*/' before? ");

    if (c=='@')
    { eight_bits cc=code_of(*loc++);
      if (cc==new_section) /* watch out for \:\ , \:\~, and \:* */
      { err_print("! Section ended in mid-comment"); loc-=2; return false; }

      if (cc==module_name) {
                           { boolean file_module= loc[-1]=='(';
                             /* does this module define an output file? */
                             cur_mod=get_module_name();
                             if (file_module && cur_mod!=NULL)

                                { int i=0;
                                  while(i<output_file_count)
                                    if (output_file[i]==cur_mod) break;  else ++i;
                                  if (i==output_file_count) /* if not present, add |cur_mod| */
                                    if (output_file_count==max_files) overflow ("output files");
                                    else output_file[output_file_count++]=cur_mod;
                                }
                           }
 continue; }
    }
    if (!line_output)
    { store_byte(c);
      if (c=='@' && (c=loc[-1])!='@') store_byte(c);
    }
  } while (c!='*' || *loc!='/' || one_liner);
  ++loc;  if (!line_output) store_byte('/');
  return false;
}

eight_bits get_next (void) /* produces the next input token */
{ static boolean preprocessing=false; /* did this line start with `\.\#'? */
  static boolean comment_continues=false; /* were we scanning a comment? */
  eight_bits c; /* the current character */
restart:

  { if (loc>=limit)
    { if (preprocessing && limit>buffer && limit[-1]!='\\')
        preprocessing=false;
      return get_line() ? '\n' : new_section;
    }
    if (comment_continues
     ||(c=*loc++)=='/' && (*loc=='*' || C_plus_plus && *loc=='/'))

    { boolean one_liner=false;
      if (!comment_continues)
      { if (!line_output) { store_byte('/'); store_byte(*loc); }
        one_liner=*loc++=='/';
          /*/ record kind of comment, and advance to comment proper {/}*/
      }
      else if (preprocessing)
      { print("\nWarning: Multi-line comment in preprocessor line");

        preprocessing=false; mark_harmless();
      }
      if (comment_continues=skip_comment(one_liner))
        /* scan to end of comment or newline */
        return '\n'; /* comment contains a newline; |get_line| has been called */
      else goto restart; /* comment complete, get next token */
    }
    if (isspace(c))
      if (line_output)
        if (preprocessing) return ' ';  else goto restart;
        /* ignore other spaces */
      else return c;
        /* when |line_output==false| preserve white space faithfully */
    if (c=='#' && loc==buffer+1) preprocessing=true;
  }
  if (c=='L' && (*loc=='\'' || *loc=='\"'))
  { get_string(); return constant; }
  if (c<0x80 ? isalpha(c) || c=='_' : translation_exists(c))
  {
    { id_first=--loc; /* mark beginning of identifier */
      do c=*++loc;
      while (c<0x80 ? isalnum(c) || c=='_' : translation_exists(c));
      cur_id=
        loc==id_first+1 && (eight_bits)(*id_first)<0x80
          ? NULL : id_lookup(id_first,loc,0);
    }
 return id_code; }
  if (c>=0x80) { err_print("! Illegal 8-bit character"); goto restart; }

  if (isdigit(c) || c=='.' && isdigit((eight_bits)*loc))
  {
    { if (*(id_first=loc-1)=='0' && tolower((eight_bits)*loc)=='x')
        /* hex constant */
        do c=*++loc; while (isxdigit(c));
      else /* octal, decimal or float constant */
      { while (isdigit(c)) c=*loc++;
        if (c=='.')  do c=*loc++; while (isdigit(c));
        if (tolower(c)=='e') /* floating point constant with exponent */
        { if ((c=*loc)=='+' || c=='-') c=*++loc;
          while (isdigit(c)) c=*++loc;
        }
        else --loc; /* back up to first character after constant */
      }
      while (isalpha(c)) c=*++loc;
        /* incorporate any `\.{U}', `\.{L}', or `\.{F}' suffixes */
      id_loc=loc;
    }
 return constant; }
  switch(c)
  { case '\'': case '"': get_string(); return constant;
    case '@':
           { eight_bits cc=code_of(*loc++);
             switch(cc)
             { case ignore: goto restart;
               case control_text: get_control_text(); goto restart;
               case verbatim: if (get_control_text()) goto restart;  else break;
               case ord:
                         id_first=loc; /* first character after opening quote */
                         while (*loc!='\'')
                         { if (*loc++=='\\') loc++; /* accept any character following backslash */
                           if (loc>=limit) { err_print("! ASCII constant didn't end"); break; }

                         }
                         id_loc=loc++;
             break;
               case module_name:

               { boolean file_module= loc[-1]=='(';
                 /* does this module define an output file? */
                 cur_mod=get_module_name();
                 if (file_module && cur_mod!=NULL)

                    { int i=0;
                      while(i<output_file_count)
                        if (output_file[i]==cur_mod) break;  else ++i;
                      if (i==output_file_count) /* if not present, add |cur_mod| */
                        if (output_file_count==max_files) overflow ("output files");
                        else output_file[output_file_count++]=cur_mod;
                    }
               }
             }
             return cc;
           }

  case '+': compress('+',plus_plus); break;
  case '-': compress('-',minus_minus); compress('>',minus_gt); break;
  case '=': compress('=',eq_eq); break;
  case '>': compress('=',gt_eq); compress('>',gt_gt); break;
  case '<': compress('=',lt_eq); compress('<',lt_lt); break;
  case '&': compress('&',and_and); break;
  case '|': compress('|',or_or); break;
  case '!': compress('=',not_eq); break;
  }
  return c;
}

void scan_repl (eight_bits context) /* creates a replacement text */
{ eight_bits a; /* the current token */
  eight_bits* keep=tok_ptr; /* non-discardable stuff up to this point */
  int brace_level=0, par_level=0;
  if (context==section_body)

    { if (mark_line()) { a=new_section; goto done; } }
  do
     { switch (a=get_next())
       {

       case id_code:

         if (cur_id==NULL) store_byte(*id_first); /* single-character identifier */
         else store_two_bytes(0x8000+id_index(cur_id));
         keep=tok_ptr; continue;
       case module_name: if (context==preproc_directive) goto done;
         if (cur_mod!=NULL) /* don't record bad module name */
         { sixteen_bits n = mod_index(cur_mod); /* index of module name */

           { char *p=loc;
             while (*p==' ' && p<limit) ++p;
             if (*p=='+') ++p; /* an optional `\.+' is allowed */
             if (*p=='=')
               err_print
                 ("! Illegal defining occurrence of module name; did you forget `@ '?");

           }
           if (line_output) tok_ptr=keep; /* back up */
           store_two_bytes(0xA800+n); keep=tok_ptr;
             /* store reference to module name */

           { if (mark_line()) { a=new_section; goto done; } }
             /* to get in phase after module insertion */
         }
         continue;
       case constant: case verbatim:

         if (id_loc==id_first+1) store_byte(*id_first); /* single digit number */
         else
         { store_byte(verb_quote); /* opening */
           do
             if (*id_first==verb_quote)
             { ++id_first; store_byte('\\'); store_byte('0'+(verb_quote>>6));
               store_byte('0'+((verb_quote>>3)&7)); store_byte('0'+(verb_quote&7));
             }
             else
             { if ((eight_bits)(*id_first)>=0x80) store_byte(0xF8);
         	/* quote 8-bit characters */
               store_byte(*id_first++);
             }
           while (id_first<id_loc);
           store_byte(verb_quote); /* closing */
           if (context==section_body && print_where)

             { if (mark_line()) { a=new_section; goto done; } }
         }
         keep=tok_ptr;
         continue;
       case ord:

         { int c=*id_first++;
           if (c=='@' && *id_first++ !='@')
           { err_print("! Double `@' should be used in ASCII constant");

               --id_first;
           }
           if (c=='\\')
           { c=*id_first++;
             switch (c)
             { case 't':c='\011';break;
               case 'n':c='\012';break;
               case 'b':c='\010';break;
               case 'f':c='\014';break;
               case 'v':c='\013';break;
               case 'r':c='\015';break;
               case '0':c='\0';break;
               case '\\':c='\134';break;
               case '\'':c='\047'; break;
               case '\"':c='\042'; break;
               default: err_print("! Unrecognised escape sequence");

             }
           }
           else
           { /* at this point |c| should be converted to its \ASCII. code number */
           }
           if (id_first!=id_loc)
             if (id_loc>id_first)
               err_print("! ASCII constant should be single character");

             else { err_print("! Empty ASCII constant"); c=0; }

           store_byte(verb_quote);
           if (c>=100) store_byte('0'+c/100);  if (c>=10) store_byte('0'+(c/10)%10);
           store_byte('0'+c%10);
           store_byte(verb_quote);
         }
         keep=tok_ptr; continue;
       case include_preproc:
         if (context==preproc_directive)
           err_print("! `@p' is forbidden in preprocessor directive");

         else
         { if (line_output) tok_ptr=keep;
           store_byte(0xFA); atp_seen=true; keep=tok_ptr;

           { if (mark_line()) { a=new_section; goto done; } }
         }
         continue;
       case char_trans:
         err_print("! `@l' is only allowed in limbo");

         continue;
       case format: case definition: case header: case begin_C:
         if (context==preproc_directive) goto done;
         err_print
           ("! `@f', `@d', `@h', and `@c' are ignored in section body");

         continue;
       case new_section: goto done;

       case '(': ++par_level; break;
       case ')':
         if (par_level<=0)
         { err_print("! Unmatched closing parenthesis"); continue; }

         --par_level; break;
       case '{': ++brace_level; break;
       case '}':
         if (brace_level<=0)
         { err_print("! Unmatched closing brace"); continue; }
         --brace_level; break;
       case '\n': store_byte('\n');
         if (context==section_body && print_where) /* input file was switched */
         { tok_ptr=keep; /* back up over discardable items */

           { if (mark_line()) { a=new_section; goto done; } }
         }
         continue;
       case join_code: a=join;
       }
       store_byte(a); keep=tok_ptr; /* mark as non-discardable */
     }
  while (true);
done: tok_ptr=keep; /* backup over trailing newlines and \&{\#line} marks */
  next_control=a; /* this control code must be reconsidered */
  if (par_level>0 || brace_level>0)
                                { char *p, *s; int l;
                                  if (par_level>0) l=par_level, s="parenthes", p= l>1 ? "es" : "is";
                                  else l=brace_level, s="brace", p=l>1 ? "s" : "";
                                  print("\n! There %s %d unclosed %s%s"
                                       ,par_level+brace_level>1 ? "are" : "is", l, s, p );

                                  if (par_level>0 && brace_level>0)
                                    print(" and %d unclosed brace%s"
                                         , brace_level, brace_level>1 ? "s" : "");
                                  print(" in the previous ");
                                  err_print(context==preproc_directive ? "macro" : "section");
                                  while (--par_level>=0) store_byte(')');
                                  while (--brace_level>=0) store_byte('}');
                                }

  { cur_text=text_ptr++; /* consolidate the replacement text */
    if (text_ptr>=text_table_end) overflow("text");
    tok_begin(text_ptr)=tok_ptr; /* mark end of replacement text */
  }
}

boolean mark_line(void)
{ while (loc>=limit)
     if (!get_line()) return true; /* skip over newlines */
  print_where=false; /* request is being serviced */
  if (line_output)
  { store_byte(0xF9);
    id_first=changing ? change.name : cur_file_name;
    store_two_bytes(changing ?  change_line : cur_line);
    store_two_bytes(id_index(id_lookup(id_first,NULL,0)));
  }
  return false;
}

void scan_section (void)
{ ++section_count;
  if (loc[-1]=='*') print_section_progress();

  { next_control=ignore; /* clear |new_section| */
    do
    { if (next_control<definition)
      { if ((next_control=skip_ahead())==module_name)
        { loc -= 2; get_next(); }
  	/* back up and scan the module name itself */
      }
      else if (next_control==definition)
      {
        { do next_control=get_next(); /* allow white space before identifier */
          while (line_output ? next_control=='\n'
                : next_control<0x80 && isspace(next_control));
          if (next_control!=id_code)
          { err_print("! Macro definition ignored, must start with identifier");

            continue;
          }

          if (cur_id==NULL) store_byte(*id_first); /* single-character identifier */
          else store_two_bytes(0x8000+id_index(cur_id));
          if (isspace((eight_bits)*loc)) store_byte(' ');
              /* separate the replacement text from a parameterless macro */
        }
        scan_repl(preproc_directive);
        cur_text->text_link=macro_flag; /* characterise as macro */
      }
      else /* |next_control==header| */
      { scan_repl(preproc_directive); cur_text->text_link=header_flag; }
      if (next_control==module_name)

      { eight_bits t=get_next();
        if (t=='+') t=get_next(); /* skip an optional `\.+' */
        if (t!='=' && t!=eq_eq) /* then apparently a cited occurrence, ignore it */
        { next_control=ignore;
          if (t!='|' && !compatibility_mode)
            err_print("! `=' sign missing, module name ignored");

        }
      }
    } while (next_control<begin_C);
  }

  { mod_pointer p; /* module name for the current section */
    switch (next_control)
    { default: return; /* no \Cee\ part present */
      case begin_C: p=NULL; break; /* section contributing to unnamed module */
      case module_name: p=cur_mod; /* section contributing to named module */
    }
    store_two_bytes(0xD000+section_count); /* record section number */
    scan_repl(section_body);
        /* now |cur_text| points to the replacement text */

    { if (p == NULL) /* unnamed section, or bad module name */
      { *last_unnamed = cur_text; last_unnamed = &cur_text->text_link; }
      else
      { text_pointer* q=&p->equiv; /* text for the current module */
        while(*q!=NULL) q=&(*q)->text_link; /* find end of list */
        *q=cur_text; /* add section to module */
      }
      cur_text->text_link=NULL;
      /* end list, also marking replacement text as a non-macro */
    }
  }
}

void phase_one (void)
{ phase=1; section_count=0; reset_input();
  while ((next_control=skip_ahead())!=new_section)
    if (next_control==char_trans)
                                { int c;
                                  while (loc<limit && isspace((eight_bits)*loc)) ++loc; /* skip space */
                                  if (!(isxdigit((eight_bits)loc[0])&&
                                	isxdigit((eight_bits)loc[1])&&
                                	isspace ((eight_bits)loc[2])))
                                    err_print("! Two-digit hex number and space should follow `@l'");

                                  else if (sscanf(loc,"%x",&c),c<0x80)
                                    err_print("! You cannot translate characters < 0x80");

                                  else
                                  { char* p=trans_of(c); int i=0;
                                    loc+=3;
                                    while (find_char() && isspace((eight_bits)*loc)) ++loc; /* skip space */
                                    if (!input_has_ended)
                                      while (isalnum((eight_bits)*loc) || *loc=='_')
                                	if (++i<=trans_limit) *p++=*loc++;  else break;
                                    if (i>0) *p='\0'; /* terminate translation unless none was given */

                                    if (i==0) err_print("! Translation string absent after `@l'");

                                    else if (i>trans_limit) err_print("! Translation string too long");
                                    else if (!isspace((eight_bits)*loc))
                                      err_print("! Translation string not terminated by space");
                                  }
                                }
  while (!input_has_ended) scan_section(); /* load all sections */
  check_complete(); /* verify that change file hasn't got out of sync */
}

#ifdef STAT
void print_stats (void)
{ print("\nMemory usage statistics:\n");
report("identifier", id_index(id_ptr), max_idents);
report("module name", mod_index(mod_ptr), max_modules);
report("byte", byte_ptr-byte_mem, max_bytes);
report("replacement text", text_ptr-text_table, max_texts);
report("token", tok_ptr-tok_mem, max_toks);
}
#endif
