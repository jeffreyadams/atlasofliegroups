#!/bin/python3.6
import os
import time
from datetime import datetime

filename="topreport2.txt"  #output to this file
interval=600               #number of seconds between reports

#nice format of datetime: 2022-05-23 16:13:19
my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')


file=open(filename,'w')
file.write("top report started at " + my_time + "\n")
file.write("interval in seconds: " + str(interval) + "\n")
file.close()
#run top
#-b: batch mode
#-o %CPU sort by CPU usage
#head -n 20: only keep 20 lines

while (1):
    top=os.popen('top -u jdada11 -o %CPU | head -1200')
    top_text = top.read()
    print(top_text)
    top.close()
    file=open(filename,'a')
    file.write(top_text)
    file.close()
    time.sleep(interval)


