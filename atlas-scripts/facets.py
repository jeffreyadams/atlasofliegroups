#!/bin/python3

#python script to compute the fpp
#usage: facets.py -g group -d dimension -o output_directory
#example: facets.py -g E8_s -d 4 -c 20 -o facets_E8
#compute number #ff of fundamental facets of dimension d
#write files facets_E8/E8_s_dim_3_i (0\le i\le #ff-1)

import sys, time, os, getopt, subprocess, queue
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing   #only for cpu count

progress_step_size=100 #how often to report progress

def atlas_compute(group,output_dir,procs,dim,i):
   atlas_arg='{}'.format("\n  facets(" + group + "," + str(dim) + "," + str(i) + ")\n").encode('utf-8')
   outputfile=  output_dir + "/facets_" + group + "_dim_" + str(dim) + "_" + str(i)
#   print("atlas_arg; ", atlas_arg)
   f=open(outputfile,"w")
   proc=procs[i]
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   while True:
      line = proc.stdout.readline().decode('ascii')
      if "Value" in line:
         break
      z=f.write(line)
   print("Done with dimension ", dim, " i=", i)

def main(argv):
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
   T=[]   #array of results from the atlas processes
   for i in range(number_fundamental_facets):
      proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE,stdout=PIPE)
      procs.append(proc)
   print("number of processes started: ", len(procs))
   P=concurrent.futures.ThreadPoolExecutor()
   for i in range(number_fundamental_facets):
         proc=procs[i]
         T.append(P.submit(atlas_compute,group,output_dir,procs,dim,i))
   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



