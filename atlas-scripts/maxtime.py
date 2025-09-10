#!/bin/python3

import getopt,os,sys, re,subprocess,decimal,time
from subprocess import Popen, PIPE, STDOUT
from decimal import Decimal

def main(argv):

    base_directory="latest" #default; input directory will be directory/logs
    outputfile=""
    opts,args=getopt.getopt(argv,"d:o:h")
    
    for opt, arg in opts:
        if opt in ('-d'):
            base_directory=arg
        elif opt in ('-o'):
            outputfile=arg
        elif opt in ('-h'):
            print("Usage: \n-d: directory\n-o: output file")
            print("reads data from directory/logs, writes to directory/outputfile")
            print("Default: read from latest/logs, output to latest/(actual name of directory).times.txt")
            exit()
    directory=base_directory + "/logs"
    if outputfile=="":
        data=subprocess.run("ls -l " + base_directory,shell=True,stdout=subprocess.PIPE).stdout.splitlines()[0].decode('ascii')
        name=data
        print("name: ", data)
        name=name.strip("/")
        name=re.sub(".*/",'',name)
        outputfile=directory   + "/" + name + ".times.txt"
    print("reading directory: ", directory,"\nwriting to file: ", outputfile)

    tempfile="maxtime_tmp"
    output_file_tmp=open(tempfile,'a')
    #grep_arg= "grep -n \"max time\" " + directory + "/*txt"
    grep_arg= "grep -n \"^|\" " + directory + "/*txt"
    print("grep_arg: ", grep_arg)
    data=subprocess.run(grep_arg,shell=True,stdout=subprocess.PIPE)
    d=data.stdout.splitlines()
    #    line without kgb_number: Total time = 116.073sec; max time = 115.100sec for (x,lambda) = ((KGB element #20895,[  2,  0, -1,  2,  1,  1,  1 ]/1))
    total_time=0
    #print("DATA: ", data)
    for line in d:
        line=line.decode('ascii')
        #print("line: ", line)
        if "job" in line:
            #print("skipping")
            junk=0
        else:
            #output_file_tmp.write("LINE: " + line + "\n")
            (file_name_and_line_number,job_number,round,pair_number,x,lambda_,time)=line.split('|')
            (file_name,line_number,nothing)=file_name_and_line_number.split(':')
            file_name=file_name.replace(directory + "/","")
            #output_file_tmp.write("file_name: " + file_name+"\n")
            #output_file_tmp.write("line_number: " + line_number+"\n")
            # output_file_tmp.write("round: " + round+"\n")
            # output_file_tmp.write("pair_number: "+ pair_number+"\n")
            # output_file_tmp.write("x: "+ x+"\n")
            # output_file_tmp.write("lambda_: " + lambda_+"\n")
            # output_file_tmp.write("time: " + time+"\n")
            line =re.sub(".*element ",'',line)
            output_file_tmp.write(time + "|" + file_name + "|" + line_number + "|" + job_number + "|" + round + "|" + pair_number + "|" + x + "|" + lambda_ +   "\n")
            #total_time+=Decimal(time)
    output_file_tmp.close()
    print("total time: ", str(total_time))
    header_arg="echo \"time|file|line number|job|round|pair_number|x|lambda|\" + \"\ndata directory:\" " + directory + "\"\n\" > " + outputfile
    print("header_arg: ", header_arg)
    data=subprocess.run(header_arg,shell=True,stdout=subprocess.PIPE)
    sort_arg= "sort -r -n " + tempfile + " >> " + outputfile
    print("sort_arg: ", sort_arg)
    data=subprocess.run(sort_arg,shell=True,stdout=subprocess.PIPE)
    os.remove(tempfile)
    output_file=open(outputfile,'a')
    output_file.write("Total time: " + str(total_time))
    output_file.close()
    print("done")

if __name__ == "__main__":
    main(sys.argv[1:])
