#!/bin/python3

import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math,psutil,signal
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
    executable_dir="/.ccs/u02/jdada11/atlasSoftware/david_facets/"
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
    print("arg: ", arg)
    proc=subprocess.Popen(arg,stdin=PIPE,stdout=PIPE)
    print("proc defined: ", proc)
    atlas_cmd="write_real_form_plus(" + group + ",\"G_temp\")" + "\n"  #plus: in FPP.at: includes j line
    print("atlas_cmd: ",  atlas_cmd)
    proc.stdin.write(format_cmd(atlas_cmd))
    proc.stdin.flush()
    group_definition_file=data_directory + "/" + group_name + "_definition.at"
    print("gdf: ", group_definition_file)
    group_definition=open(group_definition_file,"w",buffering=1)
    group_definition.write("<groups.at\n")
    while True:
        line = proc.stdout.readline().decode('ascii').strip()
        print("line in definition process: ", line)
        group_definition.write(line + "\n")
        if line.find("rf_number")>=0:
            break
            group_definition.close()
    
    #read the directories ./data/G2_s_j j=0,1,2,3 (for example)
    #also read G2_s_init.at if it exists
    dirs=[]
    for entry in os.listdir(data_directory):
       print("entry: ", entry)
       if os.path.isdir(os.path.join(data_directory,entry)) and entry.startswith(group_name):
          print("adding", entry)
          dirs.append(os.path.join(data_directory,entry))
    print("dirs: ", dirs)
    files_to_read=["all.at",group_definition_file]
    init_file=group_name + "_init.at"
    print("look for: " + init_file)
    if os.path.exists(init_file):
        files_to_read.append(init_file)
        print("adding " + init_file + " to list of files\n")
    arg=[atlas_executable] + files_to_read
    print("loading atlas, definition file, init file")
    print("files_to_read: ", files_to_read)
    proc=subprocess.Popen(arg,stdin=PIPE,stdout=PIPE)
    main_dir=data_directory + "/" + group_name + "_" + str(round)
    atlas_cmd="prints(G_temp)\n"
    print("atlas_cmd1: ", atlas_cmd)
    proc.stdin.write(format_cmd(atlas_cmd + "\n"))
    proc.stdin.flush()
    line = proc.stdout.readline().decode('ascii').strip()
    print("gtemp line: ", line)
    print("main_dir: ", main_dir)
    dirs=[main_dir]
    for dir in dirs:
        print("dir: ", dir)
        files=os.listdir(dir)
        print("files: ", files)
        for file in files:
            if file.endswith("at"):
                files_to_read.append(dir + "/" + file)
    for file in files_to_read:
       #print("loading file: ", file)
       atlas_cmd="<\"" + file + "\"\n"
       proc.stdin.write(format_cmd(atlas_cmd))
       proc.stdin.flush()
       execute_atlas_command(atlas_cmd,"Completely",proc)
    if len(output_file)==0:
        output_file=group_name + "_init.at"
        print("output_file is now: ", output_file)
    print("sending output to ", output_file, "\n")
#    atlas_cmd="> " + output_file + " big_unitary_hash.writeG("+ group + ")"
    atlas_cmd="> " + output_file + " big_unitary_hash.write()"
#    atlas_cmd="> " + output_file + " prints(111) "
    atlas_cmd="prints(\"result:\", big_unitary_hash.uhash_sizes())"
    print("atlas_cmd: " + atlas_cmd + "\n")
    proc.stdin.write(format_cmd(atlas_cmd + "\n"))
    proc.stdin.flush()
    line = proc.stdout.readline().decode('ascii').strip()
    print("LINE:", line)
    atlas_cmd="> " + output_file + " big_unitary_hash.write()"
    #atlas_cmd="> " + output_file + " prints(11111111)"
    print("atlas_cmd: " + atlas_cmd + "\n")
    proc.stdin.write(format_cmd(atlas_cmd + "\n"))
    proc.stdin.flush()
    print("done")

if __name__ == "__main__":
   main(sys.argv[1:])




