#!/bin/python3

import getopt,os,sys, re
import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math
from subprocess import Popen, PIPE, STDOUT

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
               if re.search('x:',line):
                   x=re.sub('x:','',line)
                   x=re.sub(' .*','',x)
                   time=re.sub('x:[0-9]* ','',line)   #1 day, 2:04:34
                   if re.search('day',time):
                       (day,hms)=time.split(', ')
                       days=re.sub("[^0-9]*",'',day)
                       (hours,minutes,seconds)=hms.split(':')
                       time_seconds=str(int(days)* 24 * 60 * 60 + int(hours) * 60 * 60 + int(minutes) * 60 + int(seconds))
                       kgb.append((x,time_seconds,time))
                   else:
                       (hours,minutes,seconds)=time.split(':')
#                       print("hours: ", hours, "minutes", minutes, "seconds", seconds)
                       time_seconds=str( int(hours) * 60 * 60 + int(minutes) * 60 + int(seconds))
#                       print("ts: ", time_seconds)
                       kgb.append((x,time_seconds,time))
    print("Number of KGB elements: ", len(kgb))
    out=open(outputfile,"w")
    out.write("set kgb_times=[")
    for i in range(len(kgb)-1):
        (x,time_seconds,time)=kgb[i]
        out.write("(" + x + "," + time_seconds + ",\"" + time + "\"),\n")
    (x,time_seconds,time)=kgb[-1]
    out.write("(" + x + "," + time_seconds + ",\"" + time + "\")]\n")
    out.close()

if __name__ == "__main__":
   main(sys.argv[1:])
