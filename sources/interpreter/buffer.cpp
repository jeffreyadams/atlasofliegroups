#include <iostream>
#include "buffer.h"
#include <cctype>
#include <cstring>
#include <stdexcept>
/*2:*/
#line 50 "buffer.w"
namespace atlas
{namespace interpreter
{/*9:*/
#line 172 "buffer.w"
BufferedInput::BufferedInput(std::istream&s)
:
base_stream(s),line_buffer(),p(NULL),prompt(""),prompt2(""),temp_prompt(""),
readline(NULL),add_hist(NULL),line_no(1),cur_lines(0),input_stack()
{stream= &base_stream;}

BufferedInput::BufferedInput
(const char*pr,rl_type rl,add_hist_type ah,const char*pr2)
:
base_stream(std::cin),line_buffer(),p(NULL),
prompt(pr),prompt2(pr2),temp_prompt(""),
readline(rl),add_hist(ah),line_no(1),cur_lines(0),input_stack()
{stream= &base_stream;}/*:9*//*10:*/
#line 192 "buffer.w"
BufferedInput::~BufferedInput()
{while(not input_stack.empty())
{delete input_stack.top().first;input_stack.pop();}
}/*:10*//*11:*/
#line 205 "buffer.w"
void BufferedInput::pop_file()
{delete input_stack.top().first;input_stack.pop();
stream=input_stack.empty()? &base_stream:input_stack.top().first;
}

void BufferedInput::close_includes()
{while(not input_stack.empty())
{delete input_stack.top().first;input_stack.pop();}
stream= &base_stream;
}/*:11*//*12:*/
#line 222 "buffer.w"
void BufferedInput::push_file(const char*name)
{std::ifstream*new_file=new std::ifstream(name);
if(new_file->good())
{input_stack.push(std::make_pair(new_file,std::string(name)));
stream=new_file;
}
else std::cerr<<"failed to open input file '"<<name<<"'.\n";
}/*:12*//*14:*/
#line 276 "buffer.w"
bool BufferedInput::getline()
{std::string line;line_buffer="";
line_no+=cur_lines;cur_lines=0;
const char*pr=/*15:*/
#line 302 "buffer.w"
(temp_prompt.empty()?prompt:(temp_prompt+" "+prompt2).c_str())/*:15*/
#line 279 "buffer.w"
;
while(stream->good())
{/*16:*/
#line 309 "buffer.w"
{do/*17:*/
#line 337 "buffer.w"
if(stream== &std::cin)
{prompt_length=std::strlen(pr);
if(readline!=NULL)
{char*l=readline(pr);
if(l==NULL)
{line="";stream->setstate(std::ios_base::eofbit);
std::cout<<"^D\n";
}
else
{line=l;
if(add_hist!=NULL and*l!='\0')add_hist(l);
else free(l);
}
}
else{std::cout<<pr;std::getline(std::cin,line,'\n');}
}
else std::getline(*stream,line,'\n');/*:17*/
#line 310 "buffer.w"
while(line.empty()and not stream->good()and not input_stack.empty()
and(pop_file(),true));
if(line.empty()and not stream->good())break;

}/*:16*/
#line 283 "buffer.w"
++cur_lines;
std::string::size_type l=line.length();
while(l>0 and std::isspace(line[l-1]))--l;
if(l<line.length())line.erase(l);
if(l==0 or line[l-1]!='\\'){line_buffer+=line;break;}
line.erase(l-1);line_buffer+=line;pr="\\ ";
}
line_buffer.push_back('\n');
p=line_buffer.data();
return stream->good()or cur_lines>0;

}/*:14*//*20:*/
#line 404 "buffer.w"
void BufferedInput::locate(const char*p,int&line,int&column)const
{line=line_no;column=p-line_buffer.data();}/*:20*//*22:*/
#line 421 "buffer.w"
void BufferedInput::show_range
(std::ostream&out,unsigned long l0,int c0,unsigned long l1,int c1)
const
{if(l1==line_no)
{int pl=prompt_length;
if(stream!= &std::cin or cur_lines>1)
{if(not input_stack.empty())
out<<"In input file '"<<input_stack.top().second<<"', ";
if(stream!= &std::cin)out<<"line "<<l0<<":\n";
out<<line_buffer;pl=0;
}
if(l0<l1)c0=0;
for(int i=pl+c0;i>0;--i)out<<' ';
for(int i=c1;i>c0;--i)out<<'^';
out<<std::endl;
if(l0<l1)
if(l1-l0==1)out<<"Range started in previous line\n";
else out<<"Range started "<<(l1-l0)<<"lines above\n";
}
else
{out<<"Range from line "<<l0<<" column "<<c0
<<" to line "<<l1<<" column "<<c1<<".\n";
}
}/*:22*//*24:*/
#line 462 "buffer.w"
void BufferedInput::push_prompt(char c){temp_prompt.push_back(c);}
void BufferedInput::pop_prompt()
{size_t l=temp_prompt.length();
if(l>0)temp_prompt.erase(l-1);
}
void BufferedInput::reset(){temp_prompt="";}/*:24*//*6:*/
#line 117 "buffer.w"
BufferedInput*main_input_buffer=NULL;/*:6*/
#line 54 "buffer.w"

}
}/*:2*/
