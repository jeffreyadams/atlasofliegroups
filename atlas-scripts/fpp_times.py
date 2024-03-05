#!/bin/python3

import getopt,os,sys, re
import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math
from subprocess import Popen, PIPE, STDOUT

def main(argv):
    directory=""
    outputfile="fpp_times.log"
    force=False
    opts,args=getopt.getopt(argv,"d:o:f")
    for opt, arg in opts:
        if opt in ('-d'):
            directory=arg.strip()
        elif opt in ('-o'):
            outputfile=arg
        elif opt in ('-f'):
            force=True
    if directory=="":
        print("Usage: \n-d: directory\n-o: output file (optional: default is directory/fpp.log)\n-f: force (overwrite output file)")
        exit()
    outputfile=directory + "/" + outputfile
    if os.path.isfile(outputfile) and not force:
        print("File ", outputfile , " exists, (abort)")
        exit()
#    print("out:", outputfile)
    print("reading directory: ", directory,"\nwriting to file: ",  outputfile)
    all_files=os.listdir(directory)
    kgb=[]
    for file in all_files:
         if re.match(r'[0-9]*\.txt',file):
            txt_file = directory + "/" + file;
#            print("opening txt file: ",txt_file)
            txt_data=open(txt_file,"r")
            while True:
               line=txt_data.readline().strip()
               if len(line)==0:
                   break
#               print(line)
               if re.search('list_number:',line):
                   list_number=re.sub('list_number:','',line)
                   list_number=re.sub(' .*','',list_number)
                   time=re.sub('list_number:[0-9]* ','',line)   #1 day, 2:04:34
                   if re.search('day',time):
                       (day,hms)=time.split(', ')
                       days=re.sub("[^0-9]*",'',day)
                       (hours,minutes,seconds)=hms.split(':')
                       time_seconds=int(days)* 24 * 60 * 60 + int(hours) * 60 * 60 + int(minutes) * 60 + int(seconds)
                       kgb.append((list_number,time_seconds,time))
                   else:
                       (hours,minutes,seconds)=time.split(':')
#                       print("hours: ", hours, "minutes", minutes, "seconds", seconds)
                       time_seconds=int(hours) * 60 * 60 + int(minutes) * 60 + int(seconds)
#                       print("ts: ", int(time_seconds))
                       kgb.append((list_number,time_seconds,time))
    print("Number of KGB elements: ", len(kgb))
    kgb=sorted(kgb,  key=lambda a:a[1], reverse=True)
    out=open(outputfile,"w")
    out.write("list number/time in seconds/time in (days) hours:minutes:seconds\n\n")
    for i in range(len(kgb)):
        (x,time_seconds,time)=kgb[i]
        out.write("list number=" + x + "  " + str(time_seconds) + "  " + time + "\n")
    out.close()

if __name__ == "__main__":
   main(sys.argv[1:])
