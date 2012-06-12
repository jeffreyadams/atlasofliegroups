

typedef char boolean;
typedef unsigned char eight_bits;
typedef unsigned short sixteen_bits;

typedef struct id_info
{ char *byte_start; /* beginning of the name in |byte_mem| */
  
  struct variant* equiv_or_xref; 
  
  struct id_info *hash_link; /* links identifiers with same hash code */
  int ilk; 

} id_info, *id_pointer;

typedef struct mod_info
{ char *byte_start; /* beginning of the name in |byte_mem| */
  
  struct variant* equiv_or_xref; 
  
  struct mod_info *llink,*rlink;
    /* left and right links in binary search tree */
  int key_length; 

} mod_info, *mod_pointer;


extern int program, phase;

void common_init (int argc,char** argv,char* version);

extern char buffer[], *loc, *limit;

#define max_file_name_length 60
extern struct f
{ FILE *file; char name[max_file_name_length]; sixteen_bits line; }
file[], change;
extern int include_depth;
extern boolean input_has_ended, changing, web_file_open, print_where
	     , including_header_file;

boolean locate_file_name();
boolean push_input_file(boolean,boolean); /* start a new level of input */
boolean get_line (void); /* get the next line of merged input */

#define cur_file file[include_depth].file /* current file */
#define cur_file_name file[include_depth].name /* current file name */
#define cur_line file[include_depth].line
  /* number of current line in current file */
#define web_file file[0].file
#define change_file change.file
#define change_line change.line

void reset_input (void);

extern sixteen_bits section_count;
extern eight_bits changed_section[];
#define mark_section_as_changed(n) (changed_section[(n)>>3]|=1<<((n)&7))
#define section_changed(n) ((changed_section[(n)>>3]&(1<<((n)&7)))!=0)

extern void check_complete (void);

extern char byte_mem[], *byte_ptr;
extern id_info id_table[], *id_ptr;
extern mod_info mod_table[], *mod_ptr;

extern id_pointer hash[];
#define hash_end  (&hash[hash_size]) /* end of |hash| */
id_pointer id_lookup(char*,char*,int);

extern mod_pointer root;

extern char mod_text[], *id_first, *id_loc;
#define mod_text_end (&mod_text[longest_name+1]) /* end of |mod_text| */
mod_pointer get_module_name (void);
boolean get_control_text(void);
void get_string(void);

extern history; /* indicates how bad this run was */
extern void err_print (char*), wrap_up (void), print_stats (void),
       fatal (char*,...);

extern boolean flags[];
extern char C_file_name[],idx_file_name[],scn_file_name[];

extern FILE *C_file, *tex_file;
void open_output_file(void);

void print(char*,...), print_progress(char*), print_section_progress(void);


