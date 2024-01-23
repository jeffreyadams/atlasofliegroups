#!/bin/python3.6

import os, time, getopt,sys, subprocess, glob
import time,re
import datetime

interval=600              #number of seconds between reports (default, reset with -s)
filename="freereport.log" #default, reset with -f
directory=""              #set with -d

argv=sys.argv[1:]
opts, args = getopt.getopt(argv,"s:d:f:t:n:")
testonly=False
number_of_kgb_elements=-1

for opt, arg in opts:
    if opt in ('-s'):
        interval=int(arg)
    elif opt in ('-d'):
        directory=arg
    elif opt in ('-f'):
        filename=arg
    elif opt in ('-t'):
        testonly=True
    elif opt in ('-n'):
        number_of_kgb_elements=int(arg)



if directory=="":
      print("Usage: \n-d: directory (required)\n-s: interval in seconds (default 600)\n-f: filename (default freereport.txt)\n-n: number of kgb elements\n-t: test only")
      exit()

outputfile = directory + "/" + filename
if testonly:
    print("Testing only\nOutput Directory: " + directory + "\nOutput file: " + outputfile + "\nInterval in seconds: " + str(interval) )
    exit()

#nice format of datetime: 2022-05-23 16:13:19
start_time_raw=time.time()
start_time=time.ctime()

file=open(outputfile,'w')
file.write("free report started at " + start_time + "\n")
file.write("interval in seconds: " + str(interval) + "\n")
file.close()
#run top
#-b: batch mode
#-o %CPU sort by CPU usage
#head -n 20: only keep 20 lines

def nice_time(t):
   return(re.sub("\..*","",str(datetime.timedelta(seconds=t))))

previous_time_raw=0
previous_finished=0

while (1):
    finished=0
    my_time_raw=time.time()
    my_time=time.ctime()
    elapsed_time_in_seconds = my_time_raw-start_time_raw
    local_elapsed_time_in_seconds=my_time_raw-previous_time_raw
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
    file.write("#KGB elements: " + str(finished))
    if number_of_kgb_elements !=-1:
        file.write("/" + str(number_of_kgb_elements) + " (remaining: " +  str(number_of_kgb_elements-finished) + ")")
    file.write("\n")
    avg_time=finished/(elapsed_time_in_seconds/60)
    file.write("x finished this step: " + str(finished-previous_finished) + "\n")
    local_avg_time=(finished-previous_finished)/(elapsed_time_in_seconds/60)
    file.write("total time: " + "{:.2f}".format(elapsed_time_in_seconds) + " seconds (" + nice_time(elapsed_time_in_seconds) + ")\n")
    file.write("time this step: " + nice_time(local_elapsed_time_in_seconds) + "\n")
    file.write("average number of x per minute from start: " +  "{:.2f}".format(avg_time) + "\n")
    file.write("average number of x per minute this step: " +  "{:.2f}".format(local_avg_time) + "\n")
    file.close()
    time.sleep(interval)
    previous_time_raw=my_time_raw
    previous_finished=finished


