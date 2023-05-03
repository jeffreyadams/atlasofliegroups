#!/bin/python3

#python script to compute the fpp
#usage: facets.py -g group -d dimension -o output_directory
#example: facets.py -g E8_s -d 4 -c 20 -o facets_E8
#compute number #ff of fundamental facets of dimension d
#write files facets_E8/E8_s_dim_3_i (0\le i\le #ff-1)

import sys, time, os, getopt, subprocess, queue,re
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing   #only for cpu count
from datetime import datetime

progress_step_size=100 #how often to report progress

#def atlas_compute(group,output_dir,procs,dim,i):
def atlas_compute(i):
   my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
   print("starting process #",i, " at ", my_time)
   proc=procs[i]
   atlas_arg=b"".join([b"\n  facets(",bytes(group,'utf-8'),b",",bytes(str(dim),'utf-8'),b",",bytes(str(i),'utf-8'),b")\n"])
#   atlas_arg='{}'.format("\n  facets(" + group + "," + str(dim) + "," + str(i) + ")\n").encode('utf-8')
   outputfile=  output_dir + "/facets_" + group + "_dim_" + str(dim) + "_" + str(i)
#   print("atlas_arg; ", atlas_arg)
   f=open(outputfile,"wb", buffering=4000000)
   my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
   print("running atlas function, process ", i, " at ", my_time)
   proc.stdin.write(atlas_arg)
   my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
   print("done running atlas function, process ", i, " at ", my_time)
   proc.stdin.flush()
   print("first line of process ", i, " written at ", time.ctime())
   while True:   #read until line with "Value"
      line = proc.stdout.readline().decode('ascii')
      if "Value" in line:
         break
      line=re.sub(r"[^\S\r\n]","",line)  #remove spaces but not newlines
      f.write(line.encode())
   print("Done with dimension ", dim, " i=", i)

def main(argv):
   
   global group
   global dim
   global output_dir
   global procs
   opts, args = getopt.getopt(argv, "g:c:d:o:")
   if len(opts)==0:
      print("Usage: \n-g: group\n-d: dimension\n-o: output directory")
      exit()
   for opt, arg in opts:
       if opt in ('-g'):
          group=arg
       elif opt in ('-d'):
          dim=int(arg)
       elif opt in ('-o'):
          output_dir=arg
   if not os.path.exists(output_dir):
      os.makedirs(output_dir)
   my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
   print("starting facet computation: ", my_time)
   print("group: ", group)
   print("dimension: ", dim)
   print("output directory: ", output_dir)
   #run atlas once to determine number of fundamental facets
   proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE,stdout=PIPE)
   atlas_arg='{}'.format("\n prints(#facets_fundamental(" + group + "," + str(dim) + "))\n").encode('utf-8')
   print("atlas_arg to count fundamental facets: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   number_fundamental_facets =int(proc.stdout.readline().splitlines()[0])
   print("number of fundamental facets: ", number_fundamental_facets)
   proc.kill()
   procs=[]
   for i in range(number_fundamental_facets):
      proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE,stdout=PIPE)
      procs.append(proc)
   print("number of processes started: ", len(procs))
   with concurrent.futures.ProcessPoolExecutor(number_fundamental_facets) as P:
      for i in range(number_fundamental_facets):
         proc=procs[i]
         Q=P.submit(atlas_compute,i)
   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



