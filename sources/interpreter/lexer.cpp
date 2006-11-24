/*3:*/
#line 35 "lexer.w"
#include "lexer.h"
using namespace std;
namespace atlas
{namespace interpreter
{/*5:*/
#line 89 "lexer.w"
String_pool::String_pool(size_t b)
:block_size(b),start(NULL),point(NULL),chars_left(0),prev(NULL)
{}/*:5*//*6:*/
#line 99 "lexer.w"
char*String_pool::store(const char*s,size_t len)
{size_t n=len>=block_size?len+1:block_size;

if(len>=chars_left)
{if(start!=NULL)
{String_pool*p=new String_pool(block_size);
p->prev=prev;p->start=start;
prev=p;
}
start=point=new char[n];
chars_left=n;
}
strncpy(point,s,len);char*result=point;
chars_left-=len+1;
point+=len;*point++='\0';
return result;
}/*:6*//*7:*/
#line 122 "lexer.w"
String_pool::~String_pool(void)
{delete prev;delete[]start;}/*:7*//*9:*/
#line 177 "lexer.w"
Hash_table::Hash_table(size_t init,size_t block_size)
:pool(block_size)
,mod(init<2?2:init),hash_tab(mod),name_tab(0)
{name_tab.reserve(max_fill());for(int i=0;i<mod;++i)hash_tab[i]= -1;
}/*:9*//*10:*/
#line 207 "lexer.w"
Hash_table::id_type Hash_table::hash
(const char*s,size_t l)const
{long r=static_cast<unsigned char>(*s++)%mod;
for(--l;l>0;--l)r=(static_cast<unsigned char>(*s++)+(r<<8))%mod;
return static_cast<id_type>(r);
}/*:10*//*11:*/
#line 237 "lexer.w"
Hash_table::id_type Hash_table::do_match
(const char*str,size_t l,bool copy_string)
{id_type i,h=hash(str,l);
while((i=hash_tab[h])>=0)
if(strncmp(name_tab[i],str,l)==0&&name_tab[i][l]=='\0')return i;

else if(++h==mod)h=0;
if(name_tab.size()>=max_fill())
{/*12:*/
#line 266 "lexer.w"
{hash_tab=vector<id_type>(mod=2*mod-1);
for(int h=0;h<mod;++h)hash_tab[h]= -1;
for(size_t i=0;i<name_tab.size();++i)
{int h=hash(name_tab[i],strlen(name_tab[i]));
while(hash_tab[h]>=0)if(++h==mod)h=0;
hash_tab[h]=i;
}
name_tab.reserve(max_fill());
}/*:12*/
#line 246 "lexer.w"
h=hash(str,l);
while(hash_tab[h]>=0)
if(++h==mod)h=0;
}
hash_tab[h]=name_tab.size();
name_tab.push_back(copy_string?pool.store(str,l):str);
return hash_tab[h];
}/*:11*//*20:*/
#line 411 "lexer.w"
Lexical_analyser::Lexical_analyser
(BufferedInput&source,Hash_table&hash,const char* *keywords)
:input(source),id_table(hash),nesting(0)
,prevent_termination('\0'),state(initial)
{/*21:*/
#line 427 "lexer.w"
{for(size_t i=0;keywords[i]!=0;++i)
id_table.match_literal(keywords[i]);
keyword_limit=id_table.nr_entries();
}/*:21*/
#line 416 "lexer.w"
comment_start=comment_end=256;
}/*:20*//*22:*/
#line 440 "lexer.w"
bool Lexical_analyser::reset()
{nesting=0;state=initial;input.reset();return input.getline();}/*:22*//*23:*/
#line 466 "lexer.w"
void Lexical_analyser::skip_space(void)
{if(prevent_termination!='\0')input.push_prompt(prevent_termination);
do
{char c=input.shift();
if(isspace(c))
{if(c=='\n'&&prevent_termination=='\0'&&nesting==0)break;

continue;
}
if(c==comment_start)
{input.push_prompt(comment_start);
do c=input.shift();while(c!=comment_end&&c!='\0');
input.pop_prompt();
if(c=='\n')input.unshift();
}
else break;
}while(true);
if(prevent_termination!='\0')input.pop_prompt();
input.unshift();
}/*:23*//*24:*/
#line 497 "lexer.w"
char*Lexical_analyser::scan_quoted_string()
{const char*start=input.point(),*end;
int nr_quotes=0;
do
{char c;
do c=input.shift();while(c!='"'&&c!='\n'&&c!='\0');
if(c!='"')
{input.unshift();end=input.point();
int l0,c0,l1,c1;
input.locate(start,l0,c0);input.locate(end,l1,c1);
input.show_range(cerr,l0,c0,l1,c1);
cerr<<"Closing string denotation.\n";
break;
}
else if((c=input.shift())!='"')
{input.unshift();end=input.point()-1;break;}
else++nr_quotes;
}while(true);
size_t len=end-start-nr_quotes;
char*s=new char[len+1];
while(start<end)
if((*s++= *start++)=='"')++start;
*s='\0';return s-len;
}/*:24*//*25:*/
#line 588 "lexer.w"
int Lexical_analyser::get_token(YYSTYPE*valp,YYLTYPE*locp)
{if(state==ended){state=initial;return 0;}
skip_space();prevent_termination='\0';
input.locate(input.point(),locp->first_line,locp->first_column);
int code;char c=input.shift();
if(isalpha(c))/*26:*/
#line 610 "lexer.w"
{const char*p=input.point()-1;
do c=input.shift();while(isalpha(c)||isdigit(c)||c=='_');
input.unshift();
Hash_table::id_type id_code=id_table.match(p,input.point()-p);
if(id_code>=keyword_limit){valp->id_code=id_code;code=IDENT;}
else code=QUIT+id_code;
}/*:26*/
#line 594 "lexer.w"
else if(isdigit(c))/*27:*/
#line 622 "lexer.w"
{const char*p=input.point()-1;
do c=input.shift();while(isdigit(c));
input.unshift();
const char*end=input.point();
unsigned long val= *p++ -'0';while(p<end)val=10*val+(*p++ -'0');
valp->val=val;code=INT;
}/*:27*/
#line 595 "lexer.w"
else/*28:*/
#line 638 "lexer.w"
{switch(c)
{case '"':/*29:*/
#line 679 "lexer.w"
{valp->expression=make_string_denotation(scan_quoted_string());
code=STRING;
}/*:29*/
#line 640 "lexer.w"
break;case '(':
case '{':
case '[':++nesting;input.push_prompt(c);code=c;
break;case ')':
case '}':
case ']':--nesting;input.pop_prompt();code=c;
break;case '<':
case '>':
if(state==initial)
{code=c=='<'?FROMFILE:
input.shift()=='>'?ADDTOFILE:
(input.unshift(),TOFILE);/*30:*/
#line 690 "lexer.w"
if((skip_space(),c=input.shift())=='"')
{char*s=scan_quoted_string();file_name=s;delete[]s;}
else
{file_name="";
while(!isspace(c)){file_name+=c;c=input.shift();}
input.unshift();
}/*:30*/
#line 652 "lexer.w"
break;
}

case '=':
case '+':
case '-':
case '*':
case '%':
case '^':
case ':':prevent_termination=c;code=c;
break;case '/':prevent_termination=c;
code=input.shift()=='%'?DIVMOD:(input.unshift(),'/');
break;case '\n':state=ended;
default:code=c;
}
}/*:28*/
#line 596 "lexer.w"
input.locate(input.point(),locp->last_line,locp->last_column);
if(state==initial)state=normal;return code;
}/*:25*//*16:*/
#line 322 "lexer.w"
extern "C"
char*id_completion_func(const char*text,int state)
{static size_t l;static Hash_table::id_type i,n;
if(state==0)
{i=0;n=main_hash_table->nr_entries();
l=std::strlen(text);
}
while(i<n)

{std::string s=main_hash_table->name_of(i++);

if(s.compare(0,l,text,l)==0)
{char*res=static_cast<char* >(std::malloc(s.length()+1));

if(res==NULL)
return NULL;
s.copy(res,s.length());res[s.length()]='\0';
return res;
}
}
return NULL;
}/*:16*//*14:*/
#line 287 "lexer.w"
Hash_table*main_hash_table=NULL;/*:14*//*19:*/
#line 391 "lexer.w"
Lexical_analyser*lex=NULL;/*:19*/
#line 42 "lexer.w"

}
}/*:3*/
