#!/bin/python3

import getopt,os,sys, re

def main(argv):
    directory=""
    outputfile=""
    opts,args=getopt.getopt(argv,"d:o:")
    for opt, arg in opts:
        if opt in ('-d'):
            directory=arg
        elif opt in ('-o'):
            outputfile=arg
    if directory=="" or outputfile=="":
        print("Usage: \n-d: directory\n-o: output file")
        exit()
    print("reading directory: ", directory,"\nwriting to file: ", outputfile)
    all_files=os.listdir(directory)
    params=[]
    for file in all_files:
         if re.match(r'[0-9]*\.at',file):
            at_file = directory + "/" + file;
            print("opening at file: ",at_file)
            at_data=open(at_file,"r")
            while True:
               line=at_data.readline().strip()
               line=re.sub(".*:=","",line) 
               line=re.sub("\)[^,].*",")",line)  #get rid of {dual}
               if len(line)==0:
                  break
               elif line.find("parameter")==-1:
                   ()
               else:
                  params.append(line)
    print("Number of params: ", len(params))
    out=open(outputfile,"w")
    group_data=open(directory + "/the_group.at")
    group_lines=group_data.read()
    out.write(group_lines)
    out.write("set params=[")
    for i in range(len(params)-1):
        out.write("\n" + params[i] + ",")
    out.write("\n" + params[-1] + "]\n")
    out.close()

if __name__ == "__main__":
   main(sys.argv[1:])
