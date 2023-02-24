#!/bin/python3
#python script to compute K_characters in parallel
#usage:  Kpols.py group input_file output_file
#example: Kpols.py simply_connected(G2)  G2_facets.txt G2_K_characters.txt

import sys
import concurrent.futures
import subprocess, queue
from subprocess import Popen, PIPE, STDOUT
import multiprocessing   #only for cpu count

cpu_count=multiprocessing.cpu_count()
print("Number of cores: ", cpu_count)
max_cores=128   #set this to the number of cores to use
print("Using at most ", max_cores, " cores")


#call the atlas process running on proc, with argument group, output_file (strings) and q (a queue)
#q is an array of strings, each string is of the form [ratvec,vec,vec]=[facet, bup, bdown]
#for example (in F4):
#[[ 2, 2, 1, 2 ]/3,[ 12, -7,  8, -7 ],[  1, -1, -5, -1 ]]
#facet is a point in a facet, bup and bdown are shifts
def atlas_compute(group,output_file,q,proc,i):
   count=0
   while not q.empty():
      print("i: ", i, " q: ", q.qsize())
      count+=1
      facet=q.get()
      atlas_arg='{}'.format("\n>> " + output_file +" prints(\"(\"," + facet+ "[0],\",\",K_data(K_char(" + group + "," + facet + "))," "\")\")" ).encode('utf-8')
      #sends this input line to running atlas process (for example):
      # >>       G2_Kchar.at prints("(",[[ 15, 16 ]/18,[ -3, -3 ],[ -3, -3 ]][0],",",K_data(K_char(G2_s,[[ 15, 16 ]/18,[ -3, -3 ],[ -3, -3 ]])),")" )
      #print("atlas_arg: ", atlas_arg)
      z=proc.stdin.write(atlas_arg)
#   proc.kill()
#   return(i,count)

def main(argv):
   args=sys.argv
   print("args: ",args)
   group=args[1]
   facet_file=args[2]
   output_file=args[3]
   print("group: ", group)
   print("facet_file: ", facet_file)
   print("output_file: ", output_file)
   q=queue.Queue()
   file=open(facet_file,"r")
   #read facet_file and put each entry on q
   data=file.read().splitlines()
#   print("data: ", data)
   for d in data:
      q.put(d)
   print("length of q: ", q.qsize())
   #initialize array of max_cores atlas processes
   procs=[]
   for i in range(max_cores):
      proc=subprocess.Popen(["../atlas","polsSMALLEST.at"], stdin=PIPE)
      procs.append(proc)
   print("number of procs: ", len(procs))
   #This is the object which manages the threads
   P=concurrent.futures.ThreadPoolExecutor()
   T=[]   #array of results from the atlas processes
   for i in range(max_cores):
         print("submitting job i=",i)
         T.append(P.submit(atlas_compute,group,output_file,q,procs[i],i))
   print("T: ", len(T))
#   for t in T:
#      print(t)
#      try:
#         data = t.result()
#      except Exception as exc:
#         print('%r generated an exception: %s' % ("x",exc))
#      else:
#         print("number facets: ", data)

if __name__ == "__main__":
   main(sys.argv[1:])



