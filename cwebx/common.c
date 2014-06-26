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
#define buf_size  200 
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
#define buffer_end   (&buffer[buf_size-2]) 
#define max_include_depth  32 \
  
#define max_include_paths  16 \
  
#define max_path_length  200 
#define lines_match() \
  (change_limit-change_buffer==limit-buffer \
  && strncmp(buffer, change_buffer, limit-buffer)==0)
#define byte_mem_end   (&byte_mem[max_bytes]) 
#define id_table_end   (&id_table[max_idents]) 
#define mod_table_end  (&mod_table[max_modules]) 
#define copy_char(c)  if (id_loc<mod_text_end) *id_loc++=c;  else (void)(c) 
#include  <stdarg.h> 


boolean names_match (id_pointer,char*,int,int);
void init_id_name (id_pointer,int);
void init_module_name (mod_pointer);


int program, phase;

char buffer[long_buf_size]; /* where each line of input goes */
char *loc=buffer;
    /* points to the next character to be read from the buffer */
char *limit=buffer; 

struct f file[max_include_depth]; /* stack of non-change files */
struct f change; /* change file */
local char web_file_name[max_file_name_length]
         , change_file_name[max_file_name_length]
         , alt_web_file_name[max_file_name_length];
int include_depth; /* current level of nesting */
boolean input_has_ended; /* whether there is no more input */
boolean changing; /* whether the current line is from |change_file| */
boolean web_file_open=false; /* whether the web file is being read */
boolean print_where=false; /* should |CTANGLE| print line and file info? */
local struct { char* name; int length; }
  at_h_path[max_include_paths],at_i_path;
  /* alternative search paths for \:h and \:i */
boolean including_header_file=false; 

local boolean saved_changing; /* were we changing before it was suspended? */
local char* saved_change_limit; /* end of suspended change line */
local int saved_include_depth=0; 

local char change_buffer[buf_size]; /* next line of |change_file| */
local char *change_limit; 

sixteen_bits section_count; /* the current section number */
eight_bits changed_section[(max_sections+7)/8]; 

char byte_mem[max_bytes]; /* characters of names */
char *byte_ptr=&byte_mem[0]; /* first unused position in |byte_mem| */
id_info id_table[max_idents]; /* information about identifiers */
id_pointer id_ptr=&id_table[0]; /* first unused position in |id_table| */
mod_info mod_table[max_modules]; /* information about module names */
mod_pointer mod_ptr=&mod_table[0]; 

id_pointer hash[hash_size]; 

mod_pointer root=NULL;
 

char mod_text[longest_name+1]; /* name being sought for */
char *id_first; /* where the current identifier begins in the buffer */
char *id_loc; 

int history=spotless; 

boolean flags[UCHAR_MAX+1]; /* an option for each character code */
char C_file_name[max_file_name_length]; /* name of |C_file| */
local char tex_file_name[max_file_name_length]; /* name of |tex_file| */
char idx_file_name[max_file_name_length]; /* name of index file */
char scn_file_name[max_file_name_length]; /* name of module names file */
local boolean change_file_explicit=false;
  

FILE *C_file; /* where output of \.{CTANGLE} goes */
FILE *tex_file; 

local boolean term_line_empty=true;
	


enum mod_comparison 
{ less, /* the first name is lexicographically less than,
	   but no prefix of the second */
  equal, /* the first name is equal to the second */
  greater, /* the first name is lexicographically greater than,
	      but no extension of the second */
  prefix, /* the first name is a proper prefix of the second */
  extension /* the first name is a proper extension of the second */
};

local void scan_args (int argc,char** argv);


void common_init (int argc,char** argv,char* version)
{ 
  if (argc==2 && strcmp(argv[1],"--version")==0)
    { print("%s\n",version); exit(0); }
  
  { char* cwebinputs=getenv("CWEBINPUTS");
    at_h_path[0].name=at_i_path.name=NULL; /* defaults */
  #ifdef CWEBHEADERS
    at_h_path[0].name=CWEBHEADERS;
    at_h_path[0].length=(int)strlen(CWEBHEADERS);
  #endif
    if (cwebinputs!=NULL)
    { at_i_path.length=(int)strlen(cwebinputs);
      at_i_path.name=strcpy(byte_ptr,cwebinputs);
      byte_ptr+=at_i_path.length+1;
    }
    else
    {
  #ifdef CWEBINPUTS
      at_i_path.name=CWEBINPUTS; at_i_path.length=(int)strlen(CWEBINPUTS);
  #endif
    }
  }
  
  { int i=hash_size; do hash[--i]=NULL; while(i>0); }
  
  *byte_ptr++='\0'; 
  
  mod_text[0]=' ';
  
  show_banner=show_happiness=show_progress=true;
  scan_args(argc,argv);
}

local boolean input_ln (FILE *f)
  /* copies a line into |buffer| or returns |false| */
{ register int c; /* the character read */
  register char* k=limit=buffer;  /* where next character goes */
  while ((c=getc(f))!='\n' && c!=EOF)
    if (k<=buffer_end) { *k++=c;  if (!isspace(c)) limit=k; }
  if (k>buffer_end)
  { loc=&buffer[0]; /* now |err_print| will display unbroken input line */
    err_print ("! Input line too long");   
    if (limit>buffer_end) limit=buffer_end; /* truncate line */
  }
  if (buffer[0]=='@' && limit>&buffer[1] && strchr("IXYZ",buffer[1])!=NULL)
    buffer[1]=tolower(buffer[1]);
  return c!=EOF || limit>buffer; /* whether anything new has been found */
}

boolean locate_file_name()
{ char delim=' '; /* the character being used to delimit the file name */
  while (loc<limit && (isspace((eight_bits)*loc))) ++loc;
  if (*loc=='"') delim=*loc++; /* file name in quotes */
  else if (*loc=='<') delim='>',++loc; /* file name in angle brackets */
  if (loc>=limit)
  { err_print("! Include file name not given");
	       
    return false;
  }
  id_first = loc;
  while (loc<limit &&(delim==' ' ? !isspace((eight_bits)*loc) : *loc!=delim))
    ++loc;
  id_loc=loc;
  loc=&limit[1]; /* force |input_ln| before next character is read */
  return true;
}

boolean push_input_file(boolean header,boolean suspend)
{ boolean success=false; /* whether a file has been opened */
  boolean non_system = *id_loc!='>';
    /* true unless a system header is asked for */
 if (++include_depth>=max_include_depth)
  { --include_depth;
    err_print("! Too many nested includes");

    print(" (%d)",max_include_depth);
  }
  else
  { 
    { char* k=cur_file_name;
      while (id_first<id_loc)
        if (k==&cur_file_name[max_file_name_length-1])
        { err_print("! Include file name truncated"); break; }
    		     
        else *k++=*id_first++;
      *k='\0';
    }
    if (non_system && (cur_file=fopen(cur_file_name,"r"))!=NULL)
      success=true;
    else 
         { char name_buf[max_path_length+max_file_name_length]; int i;
           if (header) /* \:h include file, or subsidiary \.{\#include} */
             for (i=0; i<max_include_paths; ++i)
               if (at_h_path[i].name==NULL) break;
               else
               { strcpy(name_buf,at_h_path[i].name);
         	strcpy(&name_buf[at_h_path[i].length],cur_file_name);
                 if ((cur_file=fopen(name_buf,"r"))!=NULL) { success=true; break; }
               }
           else if (at_i_path.name!=NULL) /* \:i include file */
           { strcpy(name_buf,at_i_path.name);
             strcpy(&name_buf[at_i_path.length],cur_file_name);
             success= (cur_file=fopen(name_buf,"r"))!=NULL;
           }
         }
    if (success)
    { cur_line=0; print_where=true;
      
      if (suspend)
      { saved_changing=changing; changing=false;
      saved_change_limit=change_limit; change_limit=change_buffer;
      saved_include_depth=include_depth;
      }
    }
    else
    { --include_depth;
      if (non_system) /* don't complain about system header files */
	err_print("! Cannot open include file");
		   
    }
  }
  return success;
}

local boolean get_web_line(void)
{ do
    if (++cur_line,input_ln(cur_file)) /* then a line has been found */
      if (!compatibility_mode
       && limit>&buffer[1] && buffer[0]=='@' && buffer[1]=='i')
      { loc=&buffer[2]; print_where=true;
        if (locate_file_name()) /* expand \:i */
         push_input_file(false,false);
      }
      else return true; /* return the line without further action */
    else if (include_depth==0) /* then end of input has been reached */
    { input_has_ended=true; web_file_open=false; return false; }
    else
    { fclose(cur_file); print_where=true;
      if (include_depth--==saved_include_depth) /* then restore |changing| */
      { changing=saved_changing; change_limit=saved_change_limit;
	saved_include_depth=0; including_header_file=false;
	if (changing) return false; /* fall back into change file */
      }
    }
  while (true);
}

local void prime_the_change_buffer (void)
{ change_limit=change_buffer;
    /* this value is used if the change file ends */
  
  do
  { if (++change_line,!input_ln(change_file)) return;
    if (limit>&buffer[1] && buffer[0]=='@')
      if (buffer[1]=='x') break;
      else if (buffer[1]=='y' || buffer[1]=='z')
      { loc=&buffer[2]; /* point out error after \:y or \:z */
        err_print ("! Where is the matching @x?"); 
      }
      else 
           { if (buffer[1]=='i' && !compatibility_mode)
             { loc=&buffer[2]; err_print ("! No includes allowed in change file"); }
           }				    
  } while (true);
  
  do
    if (++change_line,!input_ln(change_file))
    { loc=&buffer[0]; err_print("! Change file ended after @x"); return; }
  while (limit==buffer); 
  
  { int n=(int)(limit-buffer); change_limit=change_buffer+n;
    strncpy(change_buffer,buffer,n);
  }
}

local void check_change (void)
  /* switches to |change_file| if the buffers match */
{ int n=0; /* the number of discrepancies found */
  if (!lines_match()) return;
  print_where=true; /* indicate interrupted line sequencing */
  do
  { changing=true;
    
    { if (++change_line,!input_ln(change_file))
      { loc=&buffer[0]; err_print("! Change file ended before @y");
    			       
        change_limit=change_buffer; changing=false; return;
      }
      if (limit>&buffer[1] && buffer[0]=='@')
        if (buffer[1]=='y') break;
        else if (buffer[1]=='x' || buffer[1]=='z')
        { loc=&buffer[2]; err_print("! Where is the matching @y?"); }
    				     
        else 
             { if (buffer[1]=='i' && !compatibility_mode)
               { loc=&buffer[2]; err_print ("! No includes allowed in change file"); }
             }				    
      
      { int n=(int)(limit-buffer); change_limit=change_buffer+n;
        strncpy(change_buffer,buffer,n);
      }
    }
    changing=false;
    if (!get_web_line())
    { loc=&buffer[0];
      err_print("! CWEB file ended during a change"); return;
		 
    }
    if (!lines_match()) ++n;
  } while (true);
  if (n>0)
  { loc=&buffer[2];
    print("\n! Hmm... %d of the preceding lines failed to match",n);
	    
    err_print("");
  }
}

void reset_input (void) /* initialise to read the web file and change file */
{ boolean use_change_file= change_file_name[0]!='\0';
  
  { if ((web_file=fopen(web_file_name,"r"))!=NULL)
      strcpy(file[0].name,web_file_name);
    else if ((web_file=fopen(alt_web_file_name,"r"))!=NULL)
      strcpy(file[0].name,alt_web_file_name);
    else fatal("! Cannot open \"%s\" as input file", web_file_name);
  	      
    web_file_open=true;
    if (use_change_file)
      if ((change_file=fopen(change_file_name,"r"))!=NULL)
        strcpy(change.name,change_file_name);
      else if (!change_file_explicit)
        use_change_file=false; /* forget about the change file */
      else fatal("! Cannot open \"%s\" as change file", change_file_name);
  		
  }
  cur_line=0; change_line=0; include_depth=0;
  if (use_change_file) { changing=true; prime_the_change_buffer(); }
    /* prepare change file */
  else change_limit=change_buffer;
    /* emulate that change file that has ended */
  limit=buffer; loc=&buffer[1]; /* now |find_char()| will read a line */
  changing=false; input_has_ended=false;
}

boolean get_line (void) /* inputs the next line */
{
restart:
  if (changing) mark_section_as_changed(section_count);
  else 
       { if (get_web_line()
          && change_limit>change_buffer
          && limit-buffer==change_limit-change_buffer
          && buffer[0]==change_buffer[0]
            ) check_change();
       }
  if (changing)
  { 
    { if (++change_line,!input_ln (change_file))
      { err_print("! Change file ended without @z"); 
      buffer[0]='@'; buffer[1]='z'; limit=&buffer[2];
      }
      if (limit>&buffer[1] && buffer[0]=='@') /* check if the change has ended */
        if (buffer[1]=='z')
        { prime_the_change_buffer(); changing=false; print_where=true; }
        else if (buffer[1]=='x' || buffer[1]=='y')
        { loc=&buffer[2]; err_print("! Where is the matching @z?"); }
    				     
        else 
             { if (buffer[1]=='i' && !compatibility_mode)
               { loc=&buffer[2]; err_print ("! No includes allowed in change file"); }
             }				    
    }
    if (!changing)
    { mark_section_as_changed(section_count); goto restart; }
  }
  loc=&buffer[0]; *limit= ' '; /* place sentinel space */
  if (compatibility_mode && buffer[0]=='@' && buffer[1]=='i')
  { loc+=2; print_where=true;
    if (locate_file_name())
      push_input_file(false,changing);
    goto restart;
  }
  if (limit-buffer>5
   && strncmp(buffer,"#line",5)==0 && isspace((eight_bits)buffer[5]))
    
    { sixteen_bits line=0;
      print_where=true; /* output a \&{\#line} directive soon */
      loc=&buffer[6];  while (loc<limit && isspace((eight_bits)*loc)) ++loc;
      if (isdigit((eight_bits)*loc))
      { do line=10*line + *loc++ -'0'; while (isdigit((eight_bits)*loc));
        while (loc<limit && isspace((eight_bits)*loc)) ++loc;
        if (*loc++=='"')
        { int i=0;  while (&loc[i]<limit && loc[i]!='"') ++i;
          if (loc[i]=='"' && i<max_file_name_length)
          { struct f* cur_f= changing ? &change : &file[include_depth];
            cur_f->line=line-1; /* directive applies to next line, not this one */
            strncpy(cur_f->name,loc,i); cur_f->name[i]='\0';
    	goto restart;
          }
        }
      }
      err_print("! Improper #line directive"); goto restart;
    	    
    }
  return !input_has_ended;
}

void check_complete (void) /* checks that all changes were picked up */
{ if (change_limit!=change_buffer)
  { int l=(int)(change_limit-change_buffer);
    strncpy(buffer,change_buffer,l); limit=&buffer[l];
    changing=true; loc=buffer; web_file_open=true;
      /* prepare unmatched line for display */
    err_print("! Change file entry did not match");
	       
  }
}

char* store_string(char* s, int l)
{ char* dest=byte_ptr;
  if (byte_mem_end-byte_ptr<=l) overflow ("byte memory");
  byte_ptr+=l; *byte_ptr++='\0'; return strncpy(dest,s,l);
}

id_pointer id_lookup (char* first,char* last,int ilk)
  /* look up an identifier */
{ int l,h; /* length and hash code of the given identifier */
  if (last==NULL) last=first+(l=(int)strlen(first)); /* null-terminated string */
  else l=(int)(last-first); /* compute the length */
  
  { char* p=first;
    h=*p; while (++p<last) h=((h<<1)+*p)%hash_size;
  }
  
  { id_pointer p=hash[h]; /* the head of the hash list */
    while (p!=NULL && !names_match(p,first,l,ilk)) p=p->hash_link;
    if (p==NULL) /* we haven't seen this identifier before */
    
    { p=id_ptr; /* this is where the new name entry will be created */
      if (id_ptr++>=id_table_end) overflow ("identifier");
      name_begin(p)=store_string(first,l);
      if (program==cweave) init_id_name(p,ilk);
      p->hash_link=hash[h]; hash[h]=p; /* insert |p| at beginning of hash list */
    }
    return p;
  }
}

local enum mod_comparison mod_name_cmp
    (char* p, int l1, char* q, int l2)
{ int l= l1<l2 ? l1 : l2;
  while (--l>=0) if (*p++!=*q++) return *--p<*--q ? less : greater;
  return l1<l2 ? prefix : l1>l2 ? extension : equal;
}

local mod_pointer make_mod_node (char* name)
{ mod_pointer node=mod_ptr; /* allocate new node */
  if (mod_ptr++>=mod_table_end) overflow ("module name");
  name_begin(node)=name;
  node->llink=NULL; node->rlink=NULL;
  init_module_name(node); /* initialise new node */
  return node;
}

local mod_pointer mod_name_lookup (char* name, int l)
  /* finds complete module name */
{ mod_pointer p; /* current node of the search tree */
  mod_pointer* loc=&root; /* |p| will come from this location */
  while ((p=*loc)!=NULL)
  { int l0=p->key_length; char* key=name_begin(p);
    switch (mod_name_cmp(name,l,key,l0))
    { case less: loc=&p->llink; break;
      case greater: loc=&p->rlink; break;
      case equal: case extension:
	
	{ enum mod_comparison cmp=
	    mod_name_cmp(name+l0,l-l0,key+l0,(int)strlen(key+l0));
	  switch(cmp)
	  { case less: case greater:
	      err_print("! Incompatible module name"); 
	      print("\nName inconsistently extends <%.*s...>.\n",l0,key);
		     
	      return NULL;
	    case extension: case equal:
	      if (complete_name(p))
		if (cmp==equal) return p;
		else
		{ err_print("! Incompatible module name"); 
		  print("\nPrefix exists: <%s>.\n",key); return NULL;
			 
		}
	      name_begin(p)=store_string(name,l);
	        /* install |name| in place of |key| */
	      
	         free(key-1);
	      return p;
	  }
	}
      case prefix:
	err_print("! Incompatible module name"); 
	print("\nName is a prefix of <%s%s>.\n" 
	     ,key, complete_name(p) ? "" : "...");
      return NULL; /* dummy module name */
    }
  }
  
  { (p=make_mod_node(store_string(name,l)))->key_length=l; /* prepare new node */
    return *loc=p; /* install new node into tree */
  }
}

local mod_pointer prefix_lookup (char* name,int l)
  /* finds module name given a prefix */
{ mod_pointer p=root,* loc=&root; /* current node and where it comes from */
  mod_pointer match=NULL; /* the first matching node, if any */
  mod_pointer saved=NULL; /* another subtree that might have matches */
  while (p!=NULL)
  { int l0=p->key_length; char* key=name_begin(p);
    switch (mod_name_cmp(name,l,key,l0))
    { case less: p=*(loc=&p->llink); break;
      case greater: p=*(loc=&p->rlink); break;
      case equal: return p; /* a match, and no other matches are possible */
      case extension:
	
	{ enum mod_comparison cmp=
	    mod_name_cmp(name+l0,l-l0,key+l0,(int)strlen(key+l0));
	  switch(cmp)
	  { case less: case greater:
	      err_print("! Incompatible module name"); 
	      print("\nName inconsistently extends <%.*s...>.\n",l0,key);
		     
	      return NULL;
	    case prefix: case equal: return p;
	    case extension:
	      if (complete_name(p))
	      { err_print("! Incompatible module name"); 
		print("\nPrefix exists: <%s>.\n",key); return NULL; 
	      }
	      
	      { 
	           free(key-1);
	        if ((key=(char*)malloc(l+2))==NULL) fatal("Out of dynamic memory!");
	        *key++='\1'; /* ensure that |complete_name(p)| is false afterwards */
	        strncpy(key,name,l); key[l]='\0'; /* store the incomplete name */
	        name_begin(p)=key; /* install new name in node |p| */
	      }
	      return p;
	  }
	}
      case prefix:
	if (match!=NULL)
	{ err_print("! Ambiguous prefix"); return NULL; }
		       
	match=p; saved=p->rlink; p=p->llink; /* |loc| is irrelevant now */
   }
    if (p==NULL && match!=NULL)
      p=saved, saved=NULL; /* search other subtree */
  }
  if (match==NULL)
    
    { char* key=(char*)malloc(l+2);
      if (key==NULL) fatal("Out of dynamic memory!");
      *key++='\1'; /* ensure that |complete_name(p)| is false afterwards */
      strncpy(key,name,l); key[l]='\0'; /* store the incomplete name */
      (p=make_mod_node(key))->key_length=l; /* prepare new node */
      return *loc=p; /* install new node into tree */
    }
  match->key_length=l; /* |name| is a shorter prefix than used before */
  return match;
}

mod_pointer get_module_name (void)
{ 
  { eight_bits c; char* k=mod_text; /* points to last recorded character */
    do
    { if (!find_char())
      { err_print("! Input ended in module name"); break; }
  		   
      c=*loc++;
      
      if (c=='@')
      { if ((c=*loc++)=='>') break;
        if (isspace(c) || c=='*' || c=='~')
        { err_print("! Module name didn't end"); loc-=2; break; }
      		   
        if (k<mod_text_end-1) *++k='@';
          /* record the `\.{@}', now |c==loc[-1]| again */
      }
      if (isspace(c)) c=' '; /* convert tabs, newlines etc. */
      if (k<mod_text_end-1 && !(c==' ' && *k==' ') ) *++k=c;
    } while(true);
    id_first=&mod_text[1];
    if (k>=mod_text_end-1)
  { print("\n! Module name too long: "); 
      term_write(id_first,25); err_print("..");
    }
    id_loc= *k==' ' && k>mod_text ? k : k+1;
      /* point after last non-space character */
  }
  { int l=(int)(id_loc-id_first);
    return l>=3 && strncmp(id_loc-3,"...",3)==0
    	? prefix_lookup(id_first,l-3) : mod_name_lookup(id_first,l);
  }
}

boolean get_control_text(void)
{ char c,* k=id_first=&mod_text[1]; /* points after last recorded character */
  do
    if ((*k++=*loc++)=='@')
      if ((c=*loc++)!='@')
      { if (c!='>')
	  err_print("! Control codes are forbidden in control text");
		     
	return (id_loc=k-1)==id_first;
      }
  while(loc<=limit);
  err_print("! Control text didn't end"); 
  return (id_loc=k)==id_first;
}

void get_string(void)
{ char c, delim = loc[-1]; /* what started the string */
  id_loc=id_first = &mod_text[1]; copy_char(delim);
  if (delim=='L')
    *id_loc++=delim=*loc++; /* after `\.L' comes the real delimiter */
  else if (delim=='<') delim='>'; /* for file names in \&{\#include} lines */
  do
  { if (loc>=limit)
    { err_print("! String didn't end"); loc=limit; break; }
		   
    copy_char(c=*loc++);
    if (c=='\\')
      if (loc<limit) copy_char(*loc++);
        /* record `\.\\|c|' with |c!='\n'| literally */
      else if (get_line())
        if (program==cweave) --id_loc; /* |CWEAVE| erases backslash */
	else copy_char('\n'); /* but |CTANGLE| copies escaped newline */
      else
      { loc=buffer;
        err_print("! Input ended in middle of string");
		   
        break;
      }
    else if (!including_header_file && c=='@')
      if (*loc=='@') ++loc; /* undouble \:@ */
      else err_print("! Double @ required in strings");
		      
  }
  while (c!=delim);
  if (id_loc>=mod_text_end)
  { print("\n! String too long: "); 
    term_write(mod_text+1,25); err_print("..");
  }
}

void err_print (char *s) /* prints `\..' and location of error message */
{ print(*s=='!' ? "\n%s." : "%s.",s);
  if (web_file_open) 
                     { char *k, *l=(loc<limit) ? loc : limit; /* pointers into |buffer| */
                       if (changing) printf(" (l. %d of change file)\n", change_line);
                       else if (include_depth==0) printf(" (l. %d)\n", cur_line);
                       else printf(" (l. %d of include file %s)\n", cur_line, cur_file_name);
                       if (l>buffer)
                       { for (k=buffer; k<l; k++) putchar(*k=='\t' ? ' ': *k);
                         new_line();
                         for (k=buffer; k<l; k++) putchar(' '); /* space out the next line */
                       }
                       for (k=l; k<limit; k++) putchar(*k); /* print the part not yet read */
                     }
  update_terminal(); mark_error();
}

void wrap_up (void)
{
#ifdef STAT
  if (show_stats) print_stats(); /* print statistics about memory usage */
#endif
  
  { static char* mess[]=
    { "No errors were found.",
      "Did you see the warning message above?",
      "Pardon me, but I think I spotted something wrong",
      "That was a fatal error, my friend."
    };
    if (show_happiness || history>0) print("\n(%s)\n",mess[history]);
  }
  exit(history>harmless_message);
}

void fatal(char* s,...)
{ va_list p; va_start(p,s);
  vprintf(s,p); va_end(p); err_print("");
    /* print reason and location of fatal stop */
  history=fatal_message; wrap_up();
}

local void scan_args (int argc,char** argv)
{ char *dot_pos; /* position of rightmost |'.'| in the argument */
  int files_found=0, paths_found= at_h_path[0].name==NULL ? 0 : 1;
  while (--argc>0) /* first ``argument'' (program name) is irrelevant */
    if (((*++argv)[0]=='+' || (*argv)[0]=='-') && (*argv)[1]!='\0')
      
      { boolean flag_change=(**argv == '+');
        char* p=&(*argv)[1]; unsigned char c;
        while ((c=*p++)!='\0')
          if ((c=tolower(c))!='i') flags[c]=flag_change;
          else 
               { size_t l=strlen(p);
                 if (l==0) err_print("! Empty include path");
               		       
                 else if (l>max_path_length) err_print("! Include path too long");
                 					 
                 else if (paths_found>=max_include_paths)
                   err_print("! Too many include paths");
               	       
                 else
                 { at_h_path[paths_found].length=(int)l;
                   at_h_path[paths_found++].name=strcpy(byte_ptr,p);
                   byte_ptr+=l+1;
                 }
                 break;
               }
      }
    else
    { if (strlen(*argv)+5>max_file_name_length)
	/* we need room to add things like `\.{.web}' */
	fatal("! Filename too long:\n%s", *argv);
      dot_pos=strrchr(*argv,'.');
      switch (++files_found)
      { case 1:
	
	#ifndef CPPEXT
	#define CPPEXT "cpp"
	  /* extension for \Cpp\ file names; should not exceed 3 characters */
	#endif
	{ if (dot_pos==NULL) sprintf(web_file_name,"%s.w",*argv);
	  else
	  { sprintf(web_file_name,"%s",*argv); /* use file name and extension */
	    *dot_pos='\0'; /* truncate the name before the dot */
	  }
	  sprintf(alt_web_file_name,"%s.web",*argv);
	  sprintf(change_file_name,"%s.ch",*argv);
	  if (program==ctangle)
	    sprintf(C_file_name,"%s.%s",*argv, C_plus_plus ? CPPEXT : "c");
	  else
	  { sprintf(tex_file_name,"%s.tex",*argv);
	    sprintf(idx_file_name,"%s.idx",*argv);
	    sprintf(scn_file_name,"%s.scn",*argv);
	  }
	}
 
      break; case 2: 
                     if ((*argv)[0]=='-') change_file_name[0]='\0';
                     else if ((*argv)[0]!='+')
                     { change_file_explicit=true;
                       sprintf(change_file_name,dot_pos==NULL ? "%s.ch" : "%s", *argv);
                     }
 
      break; case 3: 
                     if (program==ctangle)
                       if (dot_pos!=NULL) sprintf(C_file_name, "%s", *argv);
                       else sprintf(C_file_name,"%s.%s", *argv, C_plus_plus ? CPPEXT : "c");
                     else
                     { if (dot_pos!=NULL)
                         { sprintf(tex_file_name, "%s", *argv); *dot_pos='\0'; }
                       else sprintf(tex_file_name,"%s.tex", *argv);
                       sprintf(idx_file_name,"%s.idx",*argv);
                       sprintf(scn_file_name,"%s.scn",*argv);
                     }
 
      break; default: 
                      fatal("! Usage:\n" 
                      "c%se [(+|-)options] cwebfile[.w] [(changefile[.ch]|+|-) [outputfile[.%s]]]"
                            , program==ctangle ? "tangl" : "weav"
                            , program==ctangle ? "c" : "tex");
      }
    }
  if (files_found==0) 
                    fatal("! Usage:\n" 
                    "c%se [(+|-)options] cwebfile[.w] [(changefile[.ch]|+|-) [outputfile[.%s]]]"
                          , program==ctangle ? "tangl" : "weav"
                          , program==ctangle ? "c" : "tex");
  if (paths_found<max_include_paths)
    at_h_path[paths_found].name=NULL; /* mark end of list */
}

void open_output_file(void)
{ char* name; FILE** file;
  if (program==ctangle) { name=C_file_name; file=&C_file; }
  else  { name=tex_file_name; file=&tex_file; }
  if ((*file=fopen(name,"w"))==NULL)
    fatal("! Cannot open \"%s\" as output file",name);
	   
}

void print(char* s,...)
{ va_list p; va_start(p,s);
  if (term_line_empty && *s=='\n') ++s; /* avoid printing empty line */
  vprintf(s,p); va_end(p); /* print formatted value */
  term_line_empty= s[strlen(s)-1]=='\n'; update_terminal();
}

void print_progress (char* s) { if (show_progress) print(s); }

void print_section_progress (void)
{ if (show_progress) print("*%u",section_count); }

