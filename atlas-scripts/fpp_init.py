#!/bin/python3

import sys, os, getopt, subprocess,re
from subprocess import Popen, PIPE, STDOUT

#usage:  fpp_int.py -g "Sp(6,R)" -G sp6R 
#format atlas command, given as a text string, for passing to atlas via stdin.write
def format_cmd(atlas_command):
   return('{}'.format(atlas_command).encode('utf-8'))

#execute atlas command, read/write output to log file, until encountering final_string
def execute_atlas_command(atlas_cmd,final_string,my_proc):  #no my_log: stdout instead
   print("executing atlas_cmd: " + atlas_cmd + "\n")
   my_proc.stdin.write(format_cmd(atlas_cmd))
   print("wrote")
   my_proc.stdin.flush()
   print("flushed")
   while True:
      print("TRUE LOOP")
      print("proc: ", my_proc)
      #line=my_proc.stdout.readline().decode('ascii').strip()
      line=my_proc.stdout.readline().decode('ascii').strip()
      proc.stdout.flush()
      print("LINE: " + line)
      if not line:
         print("[empty line]\n")
      if line.find(final_string)>=0:
         print("got termination line with \"" + final_string + "\": " + line + "\n")
         print("finished executing: " + atlas_cmd + "\n")
         return(line)
      else:
         print(line + "\n")
         
#MAIN
def main(argv):
    executable_dir="/.ccs/u02/jdada11/atlasSoftware/to_ht_branch_jeff_2/"
    atlas_executable=executable_dir + "atlas"
    opts, args = getopt.getopt(argv, "g:G:d:o:r:")
    #g: group; G: group_name; d: data_directory (default "./data")
    #o: output_file (default  group_name_init.at
    data_directory="./data"
    group=""
    group_name=""
    output_file=""
    round=-1
    for opt, arg in opts:
        if opt in ('-g'):
            group=arg
            group_name=arg.strip()
            print("group= " + group)
        elif opt in ('-r'):
           round=arg
        elif opt in ('-o'):
           output_file=arg;
           print("output file: ", output_file, "\n")
        elif opt in ('-G'):
            group_name=arg
            print("set group_name: ", group_name)
        elif opt in ('-d'):
            data_directory=arg
            print("data directory: " + data_directory)
    print("group: ", group, "\ngroup_name: ", group_name, "\ndata directory: ", data_directory,"\nround: ", round)
    if group=="":
       print("group must be defined with -g");
       exit()

    #write data/group_name_definition.at
    arg=[atlas_executable,"all.at"]
    #print("arg: ", arg)
    proc=subprocess.Popen(arg,stdin=PIPE,stdout=PIPE)
    atlas_cmd="write_real_form_plus(" + group + ",\"G_temp\")" + "\n"  #plus: in FPP.at: includes j line
    #print("atlas_cmd: ",  atlas_cmd)
    proc.stdin.write(format_cmd(atlas_cmd))
    proc.stdin.flush()
    group_definition_file=data_directory + "/" + group_name + "_definition.at"
    #print("group_definition_file: ", group_definition_file)
    group_definition=open(group_definition_file,"w",buffering=1)
    #print("Opened" , group_definition_file)
    group_definition.write("<groups.at\n")
    while True:
       line = proc.stdout.readline().decode('ascii').strip()
       #print("line in definition process: ", line)
       group_definition.write(line + "\n")
       if line.find("rf_number")>=0:
          break
          group_definition.close()
    #git list of directories: data_directory/e8q_j j=1,2,3... (for example)
    #also read G2_s_init.at if it exists
    #dirs needs to be a string for passing as an argument to grep
    dirs=""
    for entry in os.listdir(data_directory):
       if os.path.isdir(os.path.join(data_directory,entry)) and entry.startswith(group_name):
          dirs= dirs + os.path.join(data_directory,entry + "/*at ")
    #use grep to capture all relevant lines from #.at files in dirs
    #relevant: big_unitary_hash.finish_num(G_temp,158,131072)
    #          void: big_unitary_hash.long_match(parameter(G_temp,264,[0,1,1,0,1,1,1,1]/1,[0,0,0,0,0,0,0,0]/1),j)
    #not relevant: set j = big_unitary_hash.rf_number(G_temp)
    temp_file=group_name + "_grep.at"
    grep_arg= "grep -E -h \"big_unitary_hash.(finish|long)\" " + dirs + "> " + temp_file
    #print("grep arg: ", grep_arg)
    data=subprocess.run(grep_arg,shell=True,stdout=subprocess.PIPE)
    files_to_read=["all.at",group_definition_file,temp_file]
    init_file=group_name + "_init.at"
    if os.path.exists(init_file):
        files_to_read.append(init_file)
        #print("adding " + init_file + " to list of files\n")
    arg=[atlas_executable] + files_to_read
    #print("loading atlas, definition file, init file")
    #print("files_to_read: ", files_to_read)
    proc=subprocess.Popen(arg,stdin=PIPE,stdout=PIPE)
    main_dir=data_directory + "/" + group_name + "_" + str(round)
    atlas_cmd="prints(G_temp)\n"
    #print("atlas_cmd1: ", atlas_cmd)
    proc.stdin.write(format_cmd(atlas_cmd + "\n"))
    proc.stdin.flush()
    line = proc.stdout.readline().decode('ascii').strip()
    #print("gtemp line: ", line)
    #print("main_dir: ", main_dir)
    if len(output_file)==0:
        output_file=group_name + "_init.at"
        #print("output_file is now: ", output_file)
    #print("sending output to ", output_file, "\n")
    atlas_cmd="> " + output_file + " big_unitary_hash.write()"
    atlas_cmd="prints(\"result:\", big_unitary_hash.uhash_sizes())"
    proc.stdin.write(format_cmd(atlas_cmd + "\n"))
    proc.stdin.flush()
    line = proc.stdout.readline().decode('ascii').strip()
    atlas_cmd="> " + output_file + " big_unitary_hash.write()"
    proc.stdin.write(format_cmd(atlas_cmd + "\n"))
    proc.stdin.flush()

if __name__ == "__main__":
   main(sys.argv[1:])




