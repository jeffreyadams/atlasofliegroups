#!/bin/python3.6

import os, time, getopt,sys, subprocess, glob
import time
from datetime import datetime

interval=600              #number of seconds between reports (default, reset with -s)
filename="freereport.txt" #default, reset with -f
directory=""              #set with -d

argv=sys.argv[1:]
opts, args = getopt.getopt(argv,"s:d:f:t")
testonly=False

for opt, arg in opts:
    if opt in ('-s'):
        interval=int(arg)
    elif opt in ('-d'):
        directory=arg
    elif opt in ('-f'):
        filename=arg
    elif opt in ('-t'):
        testonly=True

if directory=="":
      print("Usage: \n-d: directory (required)\n-s: interval in seconds (default 600)\n-f: filename (default freereport.txt)\n-t: test only")
      exit()

outputfile = directory + "/" + filename
if testonly:
    print("Testing only\nOutput Directory: " + directory + "\nOutput file: " + outputfile + "\nInterval in seconds: " + str(interval) )
    exit()

#nice format of datetime: 2022-05-23 16:13:19
my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')

file=open(outputfile,'w')
file.write("free report started at " + my_time + "\n")
file.write("interval in seconds: " + str(interval) + "\n")
file.close()
#run top
#-b: batch mode
#-o %CPU sort by CPU usage
#head -n 20: only keep 20 lines

while (1):
    finished=0
    my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    top=os.popen('free -h')
    top_text = top.read()
    top.close()
    file=open(outputfile,'a')
    file.write(my_time + "\n")
    file.write(top_text)
    files=os.listdir(directory)
    for filename in files:
        if filename.endswith("txt"):
            finish= subprocess.run(['grep','Finish',str(directory + "/" + filename)],stdout=subprocess.PIPE)
            number=len(finish.stdout.splitlines())
            finished += number
    file.write("#KGB elements: " + str(finished) + "\n\n")
    file.close()
    time.sleep(interval)


