#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "common.h"
#define max_bytes 50000L
#define max_modules 1000
#define max_idents 5000
#define max_sections 4000
#define hash_size 353
#define buf_size 200
#define longest_name 1000
#define long_buf_size (buf_size+longest_name)
#define local static
#define array_size(a)((int)(sizeof(a)/sizeof(a[0])))
#define false (boolean)0
#define true (boolean)1
#define ctangle 0
#define cweave 1
#define and_and 04
#define lt_lt 020
#define gt_gt 021
#define plus_plus 013
#define minus_minus 01
#define minus_gt 031
#define not_eq 032
#define lt_eq 034
#define gt_eq 035
#define eq_eq 036
#define or_or 037
#define find_char()(loc<=limit||get_line())
#define id_index(p)((sixteen_bits)((p)-id_table))
#define id_at(i)(&id_table[i])
#define mod_index(p)((sixteen_bits)((p)-mod_table))
#define mod_at(i)(&mod_table[i])
#define name_begin(p)((p)->byte_start)
#define length(p)((int)(strlen(name_begin(p))))
#define name_end(p)(name_begin(p)+length(p))
#define complete_name(p)((p)->byte_start[-1]=='\0')
#define print_mod(p) \
printf(": <%s%s>",name_begin(p),complete_name(p)?"":"...")
#define spotless 0
#define harmless_message 1
#define error_message 2
#define fatal_message 3
#define mark_harmless() \
if(history==spotless)history=harmless_message;else
#define mark_error()(history=error_message)
#define overflow(t)fatal("\n! Sorry, %s capacity exceeded",t)
#define confusion(s)fatal("\n! This can't happen: %s",s)
#define show_banner flags['b']
#define show_happiness flags['h']
#define show_progress flags['p']
#define show_stats flags['s']
#define C_plus_plus flags['+']
#define compatibility_mode flags['c']
#define update_terminal()fflush(stdout)
#define new_line()putchar('\n')
#define term_write(string,leng)printf("%.*s",(int)(leng),string)
#define buffer_end (&buffer[buf_size-2])
#define max_include_depth 32
#define max_include_paths 16
#define max_path_length 200
#define lines_match() \
(change_limit-change_buffer==limit-buffer \
&&strncmp(buffer,change_buffer,limit-buffer)==0)
#define byte_mem_end (&byte_mem[max_bytes])
#define id_table_end (&id_table[max_idents])
#define mod_table_end (&mod_table[max_modules])
#define copy_char(c)if(id_loc<mod_text_end)*id_loc++=c;else(void)(c)
#include <stdarg.h>
/*13:*//*8:*/
#line 122 "common.inc"
boolean names_match(id_pointer,char*,int,int);
void init_id_name(id_pointer,int);
void init_module_name(mod_pointer);/*:8*//*17:*/
#line 111 "common.w"
int program,phase;/*:17*//*22:*/
#line 158 "common.w"
char buffer[long_buf_size];
char*loc=buffer;

char*limit=buffer;/*:22*//*25:*/
#line 244 "common.w"
struct f file[max_include_depth];
struct f change;
local char web_file_name[max_file_name_length]
,change_file_name[max_file_name_length]
,alt_web_file_name[max_file_name_length];
int include_depth;
boolean input_has_ended;
boolean changing;
boolean web_file_open=false;
boolean print_where=false;
local struct{char*name;int length;}
at_h_path[max_include_paths],at_i_path;

boolean including_header_file=false;/*:25*//*31:*/
#line 420 "common.w"
local boolean saved_changing;
local char*saved_change_limit;
local int saved_include_depth=0;/*:31*//*34:*/
#line 492 "common.w"
local char change_buffer[buf_size];
local char*change_limit;/*:34*//*46:*/
#line 672 "common.w"
sixteen_bits section_count;
eight_bits changed_section[(max_sections+7)/8];/*:46*//*55:*/
#line 853 "common.w"
char byte_mem[max_bytes];
char*byte_ptr= &byte_mem[0];
id_info id_table[max_idents];
id_pointer id_ptr= &id_table[0];
mod_info mod_table[max_modules];
mod_pointer mod_ptr= &mod_table[0];/*:55*//*61:*/
#line 935 "common.w"
id_pointer hash[hash_size];/*:61*//*70:*/
#line 1051 "common.w"
mod_pointer root=NULL;/*:70*//*83:*/
#line 1318 "common.w"
char mod_text[longest_name+1];
char*id_first;
char*id_loc;/*:83*//*91:*/
#line 1481 "common.w"
int history=spotless;/*:91*//*98:*/
#line 1579 "common.w"
boolean flags[UCHAR_MAX+1];
char C_file_name[max_file_name_length];
local char tex_file_name[max_file_name_length];
char idx_file_name[max_file_name_length];
char scn_file_name[max_file_name_length];
local boolean change_file_explicit=false;/*:98*//*109:*/
#line 1748 "common.w"
FILE*C_file;
FILE*tex_file;/*:109*//*112:*/
#line 1774 "common.w"
local boolean term_line_empty=true;/*:112*//*71:*/
#line 1063 "common.w"
enum mod_comparison
{less,

equal,
greater,

prefix,
extension
};/*:71*//*99:*/
#line 1591 "common.w"
local void scan_args(int argc,char* *argv);/*:99*//*:13*//*19:*/
#line 121 "common.w"
void common_init(int argc,char* *argv,char*version)
{/*20:*/
#line 133 "common.w"
if(argc==2&&strcmp(argv[1],"--version")==0)
{print("%s\n",version);exit(0);}/*:20*//*30:*/
#line 390 "common.w"
{char*cwebinputs=getenv("CWEBINPUTS");
at_h_path[0].name=at_i_path.name=NULL;
#ifdef CWEBHEADERS
at_h_path[0].name=CWEBHEADERS;
at_h_path[0].length=(int)strlen(CWEBHEADERS);
#endif
if(cwebinputs!=NULL)
{at_i_path.length=(int)strlen(cwebinputs);
at_i_path.name=strcpy(byte_ptr,cwebinputs);
byte_ptr+=at_i_path.length+1;
}
else
{
#ifdef CWEBINPUTS
at_i_path.name=CWEBINPUTS;at_i_path.length=(int)strlen(CWEBINPUTS);
#endif
}
}/*:30*//*62:*/
#line 943 "common.w"
{int i=hash_size;do hash[--i]=NULL;while(i>0);}/*:62*//*68:*/
#line 1042 "common.w"
*byte_ptr++='\0';/*:68*//*85:*/
#line 1349 "common.w"
mod_text[0]=' ';/*:85*//*104:*/
#line 1693 "common.w"
show_banner=show_happiness=show_progress=true;/*:104*/
#line 125 "common.w"
scan_args(argc,argv);
}/*:19*//*23:*/
#line 175 "common.w"
local boolean input_ln(FILE*f)

{register int c;
register char*k=limit=buffer;
while((c=getc(f))!='\n'&&c!=EOF)
if(k<=buffer_end){*k++=c;if(!isspace(c))limit=k;}
if(k>buffer_end)
{loc= &buffer[0];
err_print("! Input line too long");
if(limit>buffer_end)limit=buffer_end;
}
if(buffer[0]=='@'&&limit> &buffer[1]&&strchr("IXYZ",buffer[1])!=NULL)
buffer[1]=tolower(buffer[1]);
return c!=EOF||limit>buffer;
}/*:23*//*26:*/
#line 272 "common.w"
boolean locate_file_name()
{char delim=' ';
while(loc<limit&&(isspace((eight_bits)*loc)))++loc;
if(*loc=='"')delim= *loc++;
else if(*loc=='<')delim='>',++loc;
if(loc>=limit)
{err_print("! Include file name not given");

return false;
}
id_first=loc;
while(loc<limit&&(delim==' '?!isspace((eight_bits)*loc):*loc!=delim))
++loc;
id_loc=loc;
loc= &limit[1];
return true;
}/*:26*//*27:*/
#line 305 "common.w"
boolean push_input_file(boolean header,boolean suspend)
{boolean success=false;
boolean non_system= *id_loc!='>';

if(++include_depth>=max_include_depth)
{--include_depth;
err_print("! Too many nested includes");

print(" (%d)",max_include_depth);
}
else
{/*28:*/
#line 337 "common.w"
{char*k=cur_file_name;
while(id_first<id_loc)
if(k== &cur_file_name[max_file_name_length-1])
{err_print("! Include file name truncated");break;}

else*k++= *id_first++;
*k='\0';
}/*:28*/
#line 317 "common.w"
if(non_system&&(cur_file=fopen(cur_file_name,"r"))!=NULL)
success=true;
else/*29:*/
#line 358 "common.w"
{char name_buf[max_path_length+max_file_name_length];int i;
if(header)
for(i=0;i<max_include_paths;++i)
if(at_h_path[i].name==NULL)break;
else
{strcpy(name_buf,at_h_path[i].name);
strcpy(&name_buf[at_h_path[i].length],cur_file_name);
if((cur_file=fopen(name_buf,"r"))!=NULL){success=true;break;}
}
else if(at_i_path.name!=NULL)
{strcpy(name_buf,at_i_path.name);
strcpy(&name_buf[at_i_path.length],cur_file_name);
success=(cur_file=fopen(name_buf,"r"))!=NULL;
}
}/*:29*/
#line 320 "common.w"
if(success)
{cur_line=0;print_where=true;/*32:*/
#line 430 "common.w"
if(suspend)
{saved_changing=changing;changing=false;
saved_change_limit=change_limit;change_limit=change_buffer;
saved_include_depth=include_depth;
}/*:32*/
#line 323 "common.w"
}
else
{--include_depth;
if(non_system)
err_print("! Cannot open include file");

}
}
return success;
}/*:27*//*33:*/
#line 456 "common.w"
local boolean get_web_line(void)
{do
if(++cur_line,input_ln(cur_file))
if(!compatibility_mode
&&limit> &buffer[1]&&buffer[0]=='@'&&buffer[1]=='i')
{loc= &buffer[2];print_where=true;
if(locate_file_name())
push_input_file(false,false);
}
else return true;
else if(include_depth==0)
{input_has_ended=true;web_file_open=false;return false;}
else
{fclose(cur_file);print_where=true;
if(include_depth-- ==saved_include_depth)
{changing=saved_changing;change_limit=saved_change_limit;
saved_include_depth=0;including_header_file=false;
if(changing)return false;
}
}
while(true);
}/*:33*//*35:*/
#line 504 "common.w"
local void prime_the_change_buffer(void)
{change_limit=change_buffer;/*36:*/
#line 517 "common.w"
do
{if(++change_line,!input_ln(change_file))return;
if(limit> &buffer[1]&&buffer[0]=='@')
if(buffer[1]=='x')break;
else if(buffer[1]=='y'||buffer[1]=='z')
{loc= &buffer[2];
err_print("! Where is the matching @x?");
}
else/*37:*/
#line 537 "common.w"
{if(buffer[1]=='i'&&!compatibility_mode)
{loc= &buffer[2];err_print("! No includes allowed in change file");}
}/*:37*/
#line 526 "common.w"
}while(true);/*:36*//*38:*/
#line 546 "common.w"
do
if(++change_line,!input_ln(change_file))
{loc= &buffer[0];err_print("! Change file ended after @x");return;}
while(limit==buffer);/*:38*//*39:*/
#line 552 "common.w"
{int n=(int)(limit-buffer);change_limit=change_buffer+n;
strncpy(change_buffer,buffer,n);
}/*:39*/
#line 510 "common.w"
}/*:35*//*40:*/
#line 569 "common.w"
local void check_change(void)

{int n=0;
if(!lines_match())return;
print_where=true;
do
{changing=true;/*41:*/
#line 600 "common.w"
{if(++change_line,!input_ln(change_file))
{loc= &buffer[0];err_print("! Change file ended before @y");

change_limit=change_buffer;changing=false;return;
}
if(limit> &buffer[1]&&buffer[0]=='@')
if(buffer[1]=='y')break;
else if(buffer[1]=='x'||buffer[1]=='z')
{loc= &buffer[2];err_print("! Where is the matching @y?");}

else/*37:*/
#line 537 "common.w"
{if(buffer[1]=='i'&&!compatibility_mode)
{loc= &buffer[2];err_print("! No includes allowed in change file");}
}/*:37*//*39:*/
#line 552 "common.w"
{int n=(int)(limit-buffer);change_limit=change_buffer+n;
strncpy(change_buffer,buffer,n);
}/*:39*/
#line 612 "common.w"
}/*:41*/
#line 578 "common.w"
changing=false;
if(!get_web_line())
{loc= &buffer[0];
err_print("! CWEB file ended during a change");return;

}
if(!lines_match())++n;
}while(true);
if(n>0)
{loc= &buffer[2];
print("\n! Hmm... %d of the preceding lines failed to match",n);

err_print("");
}
}/*:40*//*43:*/
#line 625 "common.w"
void reset_input(void)
{boolean use_change_file=change_file_name[0]!='\0';/*44:*/
#line 641 "common.w"
{if((web_file=fopen(web_file_name,"r"))!=NULL)
strcpy(file[0].name,web_file_name);
else if((web_file=fopen(alt_web_file_name,"r"))!=NULL)
strcpy(file[0].name,alt_web_file_name);
else fatal("! Cannot open \"%s\" as input file",web_file_name);

web_file_open=true;
if(use_change_file)
if((change_file=fopen(change_file_name,"r"))!=NULL)
strcpy(change.name,change_file_name);
else if(!change_file_explicit)
use_change_file=false;
else fatal("! Cannot open \"%s\" as change file",change_file_name);

}/*:44*/
#line 628 "common.w"
cur_line=0;change_line=0;include_depth=0;
if(use_change_file){changing=true;prime_the_change_buffer();}

else change_limit=change_buffer;

limit=buffer;loc= &buffer[1];
changing=false;input_has_ended=false;
}/*:43*//*47:*/
#line 700 "common.w"
boolean get_line(void)
{
restart:
if(changing)mark_section_as_changed(section_count);
else/*49:*/
#line 757 "common.w"
{if(get_web_line()
&&change_limit>change_buffer
&&limit-buffer==change_limit-change_buffer
&&buffer[0]==change_buffer[0]
)check_change();
}/*:49*/
#line 705 "common.w"
if(changing)
{/*50:*/
#line 770 "common.w"
{if(++change_line,!input_ln(change_file))
{err_print("! Change file ended without @z");
buffer[0]='@';buffer[1]='z';limit= &buffer[2];
}
if(limit> &buffer[1]&&buffer[0]=='@')
if(buffer[1]=='z')
{prime_the_change_buffer();changing=false;print_where=true;}
else if(buffer[1]=='x'||buffer[1]=='y')
{loc= &buffer[2];err_print("! Where is the matching @z?");}

else/*37:*/
#line 537 "common.w"
{if(buffer[1]=='i'&&!compatibility_mode)
{loc= &buffer[2];err_print("! No includes allowed in change file");}
}/*:37*/
#line 781 "common.w"
}/*:50*/
#line 707 "common.w"
if(!changing)
{mark_section_as_changed(section_count);goto restart;}
}
loc= &buffer[0];*limit=' ';
if(compatibility_mode&&buffer[0]=='@'&&buffer[1]=='i')
{loc+=2;print_where=true;
if(locate_file_name())
push_input_file(false,changing);
goto restart;
}
if(limit-buffer>5
&&strncmp(buffer,"#line",5)==0&&isspace((eight_bits)buffer[5]))/*48:*/
#line 728 "common.w"
{sixteen_bits line=0;
print_where=true;
loc= &buffer[6];while(loc<limit&&isspace((eight_bits)*loc))++loc;
if(isdigit((eight_bits)*loc))
{do line=10*line+ *loc++ -'0';while(isdigit((eight_bits)*loc));
while(loc<limit&&isspace((eight_bits)*loc))++loc;
if(*loc++ =='"')
{int i=0;while(&loc[i]<limit&&loc[i]!='"')++i;
if(loc[i]=='"'&&i<max_file_name_length)
{struct f*cur_f=changing? &change:&file[include_depth];
cur_f->line=line-1;
strncpy(cur_f->name,loc,i);cur_f->name[i]='\0';
goto restart;
}
}
}
err_print("! Improper #line directive");goto restart;

}/*:48*/
#line 721 "common.w"
return!input_has_ended;
}/*:47*//*52:*/
#line 797 "common.w"
void check_complete(void)
{if(change_limit!=change_buffer)
{int l=(int)(change_limit-change_buffer);
strncpy(buffer,change_buffer,l);limit= &buffer[l];
changing=true;loc=buffer;web_file_open=true;

err_print("! Change file entry did not match");

}
}/*:52*//*56:*/
#line 865 "common.w"
char*store_string(char*s,int l)
{char*dest=byte_ptr;
if(byte_mem_end-byte_ptr<=l)overflow("byte memory");
byte_ptr+=l;*byte_ptr++='\0';return strncpy(dest,s,l);
}/*:56*//*63:*/
#line 956 "common.w"
id_pointer id_lookup(char*first,char*last,int ilk)

{int l,h;
if(last==NULL)last=first+(l=(int)strlen(first));
else l=(int)(last-first);/*64:*/
#line 974 "common.w"
{char*p=first;
h= *p;while(++p<last)h=((h<<1)+ *p)%hash_size;
}/*:64*//*65:*/
#line 984 "common.w"
{id_pointer p=hash[h];
while(p!=NULL&&!names_match(p,first,l,ilk))p=p->hash_link;
if(p==NULL)/*66:*/
#line 997 "common.w"
{p=id_ptr;
if(id_ptr++ >=id_table_end)overflow("identifier");
name_begin(p)=store_string(first,l);
if(program==cweave)init_id_name(p,ilk);
p->hash_link=hash[h];hash[h]=p;
}/*:66*/
#line 989 "common.w"
return p;
}/*:65*/
#line 964 "common.w"
}/*:63*//*72:*/
#line 1078 "common.w"
local enum mod_comparison mod_name_cmp
(char*p,int l1,char*q,int l2)
{int l=l1<l2?l1:l2;
while(--l>=0)if(*p++ != *q++)return* --p< * --q?less:greater;
return l1<l2?prefix:l1>l2?extension:equal;
}/*:72*//*73:*/
#line 1093 "common.w"
local mod_pointer make_mod_node(char*name)
{mod_pointer node=mod_ptr;
if(mod_ptr++ >=mod_table_end)overflow("module name");
name_begin(node)=name;
node->llink=NULL;node->rlink=NULL;
init_module_name(node);
return node;
}/*:73*//*74:*/
#line 1109 "common.w"
local mod_pointer mod_name_lookup(char*name,int l)

{mod_pointer p;
mod_pointer*loc= &root;
while((p= *loc)!=NULL)
{int l0=p->key_length;char*key=name_begin(p);
switch(mod_name_cmp(name,l,key,l0))
{case less:loc= &p->llink;break;
case greater:loc= &p->rlink;break;
case equal:case extension:/*76:*/
#line 1151 "common.w"
{enum mod_comparison cmp=
mod_name_cmp(name+l0,l-l0,key+l0,(int)strlen(key+l0));
switch(cmp)
{case less:case greater:
err_print("! Incompatible module name");
print("\nName inconsistently extends <%.*s...>.\n",l0,key);

return NULL;
case extension:case equal:
if(complete_name(p))
if(cmp==equal)return p;
else
{err_print("! Incompatible module name");
print("\nPrefix exists: <%s>.\n",key);return NULL;

}
name_begin(p)=store_string(name,l);/*79:*/
#line 1242 "common.w"
free(key-1);/*:79*/
#line 1170 "common.w"
return p;
}
}/*:76*/
#line 1123 "common.w"
case prefix:
err_print("! Incompatible module name");
print("\nName is a prefix of <%s%s>.\n"
,key,complete_name(p)?"":"...");
return NULL;
}
}/*75:*/
#line 1140 "common.w"
{(p=make_mod_node(store_string(name,l)))->key_length=l;
return*loc=p;
}/*:75*/
#line 1132 "common.w"
}/*:74*//*77:*/
#line 1194 "common.w"
local mod_pointer prefix_lookup(char*name,int l)

{mod_pointer p=root,*loc= &root;
mod_pointer match=NULL;
mod_pointer saved=NULL;
while(p!=NULL)
{int l0=p->key_length;char*key=name_begin(p);
switch(mod_name_cmp(name,l,key,l0))
{case less:p= *(loc= &p->llink);break;
case greater:p= *(loc= &p->rlink);break;
case equal:return p;
case extension:/*80:*/
#line 1251 "common.w"
{enum mod_comparison cmp=
mod_name_cmp(name+l0,l-l0,key+l0,(int)strlen(key+l0));
switch(cmp)
{case less:case greater:
err_print("! Incompatible module name");
print("\nName inconsistently extends <%.*s...>.\n",l0,key);

return NULL;
case prefix:case equal:return p;
case extension:
if(complete_name(p))
{err_print("! Incompatible module name");
print("\nPrefix exists: <%s>.\n",key);return NULL;
}/*81:*/
#line 1278 "common.w"
{/*79:*/
#line 1242 "common.w"
free(key-1);/*:79*/
#line 1279 "common.w"
if((key=(char*)malloc(l+2))==NULL)fatal("Out of dynamic memory!");
*key++='\1';
strncpy(key,name,l);key[l]='\0';
name_begin(p)=key;
}/*:81*/
#line 1266 "common.w"
return p;
}
}/*:80*/
#line 1209 "common.w"
case prefix:
if(match!=NULL)
{err_print("! Ambiguous prefix");return NULL;}

match=p;saved=p->rlink;p=p->llink;
}
if(p==NULL&&match!=NULL)
p=saved,saved=NULL;
}
if(match==NULL)/*78:*/
#line 1232 "common.w"
{char*key=(char*)malloc(l+2);
if(key==NULL)fatal("Out of dynamic memory!");
*key++='\1';
strncpy(key,name,l);key[l]='\0';
(p=make_mod_node(key))->key_length=l;
return*loc=p;
}/*:78*/
#line 1221 "common.w"
match->key_length=l;
return match;
}/*:77*//*84:*/
#line 1333 "common.w"
mod_pointer get_module_name(void)
{/*86:*/
#line 1357 "common.w"
{eight_bits c;char*k=mod_text;
do
{if(!find_char())
{err_print("! Input ended in module name");break;}

c= *loc++;/*87:*/
#line 1384 "common.w"
if(c=='@')
{if((c= *loc++)=='>')break;
if(isspace(c)||c=='*'||c=='~')
{err_print("! Module name didn't end");loc-=2;break;}

if(k<mod_text_end-1)* ++k='@';

}/*:87*/
#line 1364 "common.w"
if(isspace(c))c=' ';
if(k<mod_text_end-1&&!(c==' '&& *k==' '))* ++k=c;
}while(true);
id_first= &mod_text[1];
if(k>=mod_text_end-1)
{print("\n! Module name too long: ");
term_write(id_first,25);err_print("..");
}
id_loc= *k==' '&&k>mod_text?k:k+1;

}/*:86*/
#line 1336 "common.w"
{int l=(int)(id_loc-id_first);
return l>=3&&strncmp(id_loc-3,"...",3)==0
?prefix_lookup(id_first,l-3):mod_name_lookup(id_first,l);
}
}/*:84*//*88:*/
#line 1402 "common.w"
boolean get_control_text(void)
{char c,*k=id_first= &mod_text[1];
do
if((*k++= *loc++)=='@')
if((c= *loc++)!='@')
{if(c!='>')
err_print("! Control codes are forbidden in control text");

return(id_loc=k-1)==id_first;
}
while(loc<=limit);
err_print("! Control text didn't end");
return(id_loc=k)==id_first;
}/*:88*//*89:*/
#line 1428 "common.w"
void get_string(void)
{char c,delim=loc[-1];
id_loc=id_first= &mod_text[1];copy_char(delim);
if(delim=='L')
*id_loc++=delim= *loc++;
else if(delim=='<')delim='>';
do
{if(loc>=limit)
{err_print("! String didn't end");loc=limit;break;}

copy_char(c= *loc++);
if(c=='\\')
if(loc<limit)copy_char(*loc++);

else if(get_line())
if(program==cweave)--id_loc;
else copy_char('\n');
else
{loc=buffer;
err_print("! Input ended in middle of string");

break;
}
else if(!including_header_file&&c=='@')
if(*loc=='@')++loc;
else err_print("! Double @ required in strings");

}
while(c!=delim);
if(id_loc>=mod_text_end)
{print("\n! String too long: ");
term_write(mod_text+1,25);err_print("..");
}
}/*:89*//*92:*/
#line 1491 "common.w"
void err_print(char*s)
{print(*s=='!'?"\n%s.":"%s.",s);
if(web_file_open)/*93:*/
#line 1504 "common.w"
{char*k,*l=(loc<limit)?loc:limit;
if(changing)printf(" (l. %d of change file)\n",change_line);
else if(include_depth==0)printf(" (l. %d)\n",cur_line);
else printf(" (l. %d of include file %s)\n",cur_line,cur_file_name);
if(l>buffer)
{for(k=buffer;k<l;k++)putchar(*k=='\t'?' ':*k);
new_line();
for(k=buffer;k<l;k++)putchar(' ');
}
for(k=l;k<limit;k++)putchar(*k);
}/*:93*/
#line 1494 "common.w"
update_terminal();mark_error();
}/*:92*//*94:*/
#line 1525 "common.w"
void wrap_up(void)
{
#ifdef STAT
if(show_stats)print_stats();
#endif/*95:*/
#line 1538 "common.w"
{static char*mess[]=
{"No errors were found.",
"Did you see the warning message above?",
"Pardon me, but I think I spotted something wrong",
"That was a fatal error, my friend."
};
if(show_happiness||history>0)print("\n(%s)\n",mess[history]);
}/*:95*/
#line 1531 "common.w"
exit(history>harmless_message);
}/*:94*//*96:*/
#line 1554 "common.w"
void fatal(char*s,...)
{va_list p;va_start(p,s);
vprintf(s,p);va_end(p);err_print("");

history=fatal_message;wrap_up();
}/*:96*//*100:*/
#line 1607 "common.w"
local void scan_args(int argc,char* *argv)
{char*dot_pos;
int files_found=0,paths_found=at_h_path[0].name==NULL?0:1;
while(--argc>0)
if(((* ++argv)[0]=='+'||(*argv)[0]=='-')&&(*argv)[1]!='\0')/*105:*/
#line 1700 "common.w"
{boolean flag_change=(* *argv=='+');
char*p= &(*argv)[1];unsigned char c;
while((c= *p++)!='\0')
if((c=tolower(c))!='i')flags[c]=flag_change;
else/*106:*/
#line 1712 "common.w"
{size_t l=strlen(p);
if(l==0)err_print("! Empty include path");

else if(l>max_path_length)err_print("! Include path too long");

else if(paths_found>=max_include_paths)
err_print("! Too many include paths");

else
{at_h_path[paths_found].length=(int)l;
at_h_path[paths_found++].name=strcpy(byte_ptr,p);
byte_ptr+=l+1;
}
break;
}/*:106*/
#line 1705 "common.w"
}/*:105*/
#line 1613 "common.w"
else
{if(strlen(*argv)+5>max_file_name_length)

fatal("! Filename too long:\n%s",*argv);
dot_pos=strrchr(*argv,'.');
switch(++files_found)
{case 1:/*101:*/
#line 1641 "common.w"
#ifndef CPPEXT
#define CPPEXT "cpp"

#endif
{if(dot_pos==NULL)sprintf(web_file_name,"%s.w",*argv);
else
{sprintf(web_file_name,"%s",*argv);
*dot_pos='\0';
}
sprintf(alt_web_file_name,"%s.web",*argv);
sprintf(change_file_name,"%s.ch",*argv);
if(program==ctangle)
sprintf(C_file_name,"%s.%s",*argv,C_plus_plus?CPPEXT:"c");
else
{sprintf(tex_file_name,"%s.tex",*argv);
sprintf(idx_file_name,"%s.idx",*argv);
sprintf(scn_file_name,"%s.scn",*argv);
}
}/*:101*/
#line 1621 "common.w"

break;case 2:/*102:*/
#line 1668 "common.w"
if((*argv)[0]=='-')change_file_name[0]='\0';
else if((*argv)[0]!='+')
{change_file_explicit=true;
sprintf(change_file_name,dot_pos==NULL?"%s.ch":"%s",*argv);
}/*:102*/
#line 1622 "common.w"

break;case 3:/*103:*/
#line 1678 "common.w"
if(program==ctangle)
if(dot_pos!=NULL)sprintf(C_file_name,"%s",*argv);
else sprintf(C_file_name,"%s.%s",*argv,C_plus_plus?CPPEXT:"c");
else
{if(dot_pos!=NULL)
{sprintf(tex_file_name,"%s",*argv);*dot_pos='\0';}
else sprintf(tex_file_name,"%s.tex",*argv);
sprintf(idx_file_name,"%s.idx",*argv);
sprintf(scn_file_name,"%s.scn",*argv);
}/*:103*/
#line 1623 "common.w"

break;default:/*107:*/
#line 1732 "common.w"
fatal("! Usage:\n"
"c%se [(+|-)options] cwebfile[.w] [(changefile[.ch]|+|-) [outputfile[.%s]]]"
,program==ctangle?"tangl":"weav"
,program==ctangle?"c":"tex");/*:107*/
#line 1625 "common.w"
}
}
if(files_found==0)/*107:*/
#line 1732 "common.w"
fatal("! Usage:\n"
"c%se [(+|-)options] cwebfile[.w] [(changefile[.ch]|+|-) [outputfile[.%s]]]"
,program==ctangle?"tangl":"weav"
,program==ctangle?"c":"tex");/*:107*/
#line 1628 "common.w"
if(paths_found<max_include_paths)
at_h_path[paths_found].name=NULL;
}/*:100*//*110:*/
#line 1752 "common.w"
void open_output_file(void)
{char*name;FILE* *file;
if(program==ctangle){name=C_file_name;file= &C_file;}
else{name=tex_file_name;file= &tex_file;}
if((*file=fopen(name,"w"))==NULL)
fatal("! Cannot open \"%s\" as output file",name);

}/*:110*//*113:*/
#line 1788 "common.w"
void print(char*s,...)
{va_list p;va_start(p,s);
if(term_line_empty&& *s=='\n')++s;
vprintf(s,p);va_end(p);
term_line_empty=s[strlen(s)-1]=='\n';update_terminal();
}

void print_progress(char*s){if(show_progress)print(s);}

void print_section_progress(void)
{if(show_progress)print("*%u",section_count);}/*:113*/
