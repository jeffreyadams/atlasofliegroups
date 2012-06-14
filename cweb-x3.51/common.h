/*14:*//*15:*/
#line 95 "common.w"
typedef char boolean;
typedef unsigned char eight_bits;
typedef unsigned short sixteen_bits;/*:15*//*53:*/
#line 820 "common.w"
typedef struct id_info
{char*byte_start;/*57:*/
#line 902 "common.w"
struct variant*equiv_or_xref;/*:57*//*59:*/
#line 913 "common.w"
struct id_info*hash_link;
int ilk;/*:59*/
#line 822 "common.w"

}id_info,*id_pointer;

typedef struct mod_info
{char*byte_start;/*58:*/
#line 905 "common.w"
struct variant*equiv_or_xref;/*:58*//*67:*/
#line 1021 "common.w"
struct mod_info*llink,*rlink;

int key_length;/*:67*/
#line 827 "common.w"

}mod_info,*mod_pointer;/*:53*//*16:*/
#line 108 "common.w"
extern int program,phase;/*:16*//*18:*/
#line 117 "common.w"
void common_init(int argc,char* *argv,char*version);/*:18*//*21:*/
#line 146 "common.w"
extern char buffer[],*loc,*limit;/*:21*//*24:*/
#line 214 "common.w"
#define max_file_name_length 256
extern struct f
{FILE*file;char name[max_file_name_length];sixteen_bits line;}
file[],change;
extern int include_depth;
extern boolean input_has_ended,changing,web_file_open,print_where
,including_header_file;

boolean locate_file_name();
boolean push_input_file(boolean,boolean);
boolean get_line(void);

#define cur_file file[include_depth].file 
#define cur_file_name file[include_depth].name 
#define cur_line file[include_depth].line

#define web_file file[0].file
#define change_file change.file
#define change_line change.line/*:24*//*42:*/
#line 619 "common.w"
void reset_input(void);/*:42*//*45:*/
#line 665 "common.w"
extern sixteen_bits section_count;
extern eight_bits changed_section[];
#define mark_section_as_changed(n) (changed_section[(n)>>3]|=1<<((n)&7))
#define section_changed(n) ((changed_section[(n)>>3]&(1<<((n)&7)))!=0)/*:45*//*51:*/
#line 787 "common.w"
extern void check_complete(void);/*:51*//*54:*/
#line 838 "common.w"
extern char byte_mem[],*byte_ptr;
extern id_info id_table[],*id_ptr;
extern mod_info mod_table[],*mod_ptr;/*:54*//*60:*/
#line 928 "common.w"
extern id_pointer hash[];
#define hash_end  (&hash[hash_size]) 
id_pointer id_lookup(char*,char*,int);/*:60*//*69:*/
#line 1047 "common.w"
extern mod_pointer root;/*:69*//*82:*/
#line 1308 "common.w"
extern char mod_text[],*id_first,*id_loc;
#define mod_text_end (&mod_text[longest_name+1]) 
mod_pointer get_module_name(void);
boolean get_control_text(void);
void get_string(void);/*:82*//*90:*/
#line 1475 "common.w"
extern history;
extern void err_print(char*),wrap_up(void),print_stats(void),
fatal(char*,...);/*:90*//*97:*/
#line 1575 "common.w"
extern boolean flags[];
extern char C_file_name[],idx_file_name[],scn_file_name[];/*:97*//*108:*/
#line 1741 "common.w"
extern FILE*C_file,*tex_file;
void open_output_file(void);/*:108*//*111:*/
#line 1767 "common.w"
void print(char*,...),print_progress(char*),print_section_progress(void);/*:111*//*:14*/
