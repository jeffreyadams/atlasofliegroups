#!/bin/python3

import getopt,os,sys, re,subprocess,time, datetime
from subprocess import Popen, PIPE, STDOUT

def help(base_directory):
    print("Usage: \n-d: directory (default: ", base_directory,")")
    print("-G: group_name - reads from directory/group_name_j/logs/[number].txt")
    print("-o: output file (default: ./group_name.times.txt)")
    print("-c: cutoff (only keep lines with time > cutoff, default=-1)")
    exit()

def main(argv):
    base_directory="data" #default; input directory will be directory/logs
    outputfile=""
    group_name=""
    number_processors=1
    cutoff=-1
    excluded=0
    number_of_pairs=0
    line_start_string=".|"  #look for lines starting with this
                           #default: |  to read single line running reports use ".|"
    opts,args=getopt.getopt(argv,"o:d:n:h:l:G:c:")
    for opt, arg in opts:
        if opt in ('-d'):
            base_directory=arg
        elif opt in ('-n'):
            number_processors=arg
        elif opt in ('-o'):
            outputfile=arg
        elif opt in ('-G'):
            group_name=arg
        elif opt in ('-l'):
            line_start_string=arg
        elif opt in ('-c'):
            cutoff=arg;
        elif opt in ('-h'):
            help(base_directory)
    if group_name=="":
        help(base_directory)
    directories=base_directory + "/" + group_name + "_[0-9]*/logs"
    files=directories + "/[0-9]*.txt"
    if outputfile=="":
        outputfile="./" + group_name + ".times.txt"
    print("reading files: ", files,"\nwriting to file: ", outputfile)
    tempfile="time_report.tmp"
    output_file_tmp=open(tempfile,'a')
    #grep_arg= "grep -n \"^|\" " + directories + "/[0-9]*txt"
    grep_arg= "grep -n \"^" + line_start_string + "\" " + directories + "/[0-9]*txt"
    print("running: ", grep_arg)
    data=subprocess.run(grep_arg,shell=True,stdout=subprocess.PIPE)
    print("finished running ", grep_arg)
    d=data.stdout.splitlines()
    total_time=0
    #|job_number|round|restart|pair number|x|lambda|time
    #  |0|1||1|26170|978|[ 1, 3, 2, 0, 2, 1 ]/1|0:00:01
    print("processing lines")
    for line in d:
        line=line.decode('ascii')
        #print("line: ", line)
        if "job" in line:
            #print("skipping")
            junk=0
        else:
            #output_file_tmp.write("LINE: " + line + "\n")
            #print("line: ", line)
            #data/e6s_1/logs/9.txt:24317:|9|1||1|26275|981|[ 1, 1, 2, 2, 2, 1 ]/1|0:00:04
            (file_name_and_line_number,job_number,round,junk,restart,pair_number,x,lambda_,days_and_time)=line.split('|')
            (file_name,line_number,nothing)=file_name_and_line_number.split(':')
            days=0
            if "day" in days_and_time:
                (days,time)=days_and_time.split(",")
                days=re.sub("s","",days)
                days=int(re.sub(" day","",days))
            else:
                time=days_and_time
            file_name=file_name.split("/")[-1]
            verbose=False
            if verbose:
                print("file_name: ", file_name)
                print("job_number: ", job_number)
                print("round: ", round)
                print("restart: ", restart)
                print("pair_number: ", pair_number)
                print("x: ", x)
                print("lambda: ", lambda_)
                print("days_and_time: ", days_and_time)
                print("days: ", days)
                print("time: ", time)
                print("junk: ", junk)
            
            #output_file_tmp.write("file_name: " + file_name+"\n")
            #output_file_tmp.write("line_number: " + line_number+"\n")
            # output_file_tmp.write("round: " + round+"\n")
            # output_file_tmp.write("pair_number: "+ pair_number+"\n")
            # output_file_tmp.write("x: "+ x+"\n")
            # output_file_tmp.write("lambda_: " + lambda_+"\n")
            # output_file_tmp.write("time: " + time+"\n")
            line =re.sub(".*element ",'',line)
            h, m, s = map(int, time.split(':'))
            if days>0:
                h=h+24*days
                print("added days, h=", h)
            total_seconds = h * 3600 + m * 60 + s
            total_time+=total_seconds
            if total_seconds>int(cutoff):
                number_of_pairs+=1
                output_file_tmp.write(days_and_time + "|" + file_name + "|" + line_number + "|" + job_number + "|" + restart + "|" + round + "|" + pair_number + "|" + x + "|" + lambda_ +   "\n")
            else:
                excluded+=1
    output_file_tmp.close()
    print("finished processing lines")
    average=total_time/int(number_processors)
    print("total time: ", total_time)
    print("average: ", average)
    print("excluded (<=", cutoff, "seconds): ", excluded)
    total_time_string=str(datetime.timedelta(seconds=total_time))
    time_per_processor_string=str(datetime.timedelta(seconds=average))
    print("raw: ", datetime.timedelta(seconds=average))
    print("total time: ", total_time_string)
    print("per: " , time_per_processor_string)
    #header_arg="echo \"time|file|line number|job|round|pair_number|x|lambda|\" + \"\ndata directory:\" " + directories + "\"\n\" > " + outputfile

    header_arg="echo \"data directory:\" " + directories + "\"\ntime|file|line number|job|round|pair_number|x|lambda|\n\" > " + outputfile
    print("header_arg: ", header_arg)
    data=subprocess.run(header_arg,shell=True,stdout=subprocess.PIPE)
    sort_arg= "sort -r -n " + tempfile + " >> " + outputfile
    print("sort_arg: ", sort_arg)
    data=subprocess.run(sort_arg,shell=True,stdout=subprocess.PIPE)
    os.remove(tempfile)
    output_file=open(outputfile,'a')
    excluded_string=f"{excluded:,}"
    output_file.write("\nnumber of pairs: " + str(number_of_pairs) + "\n")
    output_file.write("excluded (<=" + str(cutoff) + "seconds): " + excluded_string+ "\n")
    output_file.write("Total time in seconds: " + str(total_time) + "\n")
    output_file.write("Total time: " + total_time_string + "\n")
    print("outputfile: ", outputfile)
    output_file.write("Average time per processor: " +  time_per_processor_string + "\n")
    output_file.close()
    print("done")

if __name__ == "__main__":
    main(sys.argv[1:])


