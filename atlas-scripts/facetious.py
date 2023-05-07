#!/bin/python3

#python script to compute the facets and K_characters
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
   atlas_arg_1=b"".join([b"set vd=vertex_data(",bytes(group,'utf-8'),b")\n"])
   atlas_arg_2=b"".join([b" void:facets_big(",bytes(group,'utf-8'),b",",bytes(str(dim),'utf-8'),b",",bytes(str(i),'utf-8'),b")\n"])
   atlas_arg_3=b"".join([b" prints(K_data(K_char(",bytes(group,'utf-8'),b",",bytes(facet,'utf-8'),b"\n"])
   atlas_arg_4=b"".join([b" prints(\"end\")\n"])
   atlas_arg=atlas_arg_1+atlas_arg_2
   print("atlas_arg: ", repr(atlas_arg))
   z=proc.stdin.write(atlas_arg)
   print("z=",z)
   proc.stdin.flush()
#   z=proc.stdin.write(atlas_arg)
#   proc.stdin.flush()
   
#   print("z=",z)
#   z1=proc.stdin.write(atlas_vd_arg_1)
#   z2=proc.stdin.write(atlas_vd_arg_2)
   outputfile= data_dir + "/facets_" + group + "_dim_" + str(dim) + "_" + str(i)
   print("outputfile: ", outputfile)
   f=open(outputfile,"wb", buffering=4000000)
   firstline=True
   while True:   #read until line with "Value"
      line = proc.stdout.readline().decode('ascii')
      if not "Variable" in line:
         if "end" in line:
            break
         if firstline:
            print("first line of process ", i, " written at ", time.ctime())
            f.write("set v=[\n".encode())
            line=re.sub(r"[^\S\r\n]","",line)  #remove spaces but not newlines
            f.write(line.encode())
            firstline=False
         else:
            line=","+line
            f.write(line.encode())
   f.write("]".encode())
   f.close()
   print("closed")


   my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
   print("running atlas_compute/process#", i, " at ", my_time)
   proc.stdin.write(atlas_arg)
   my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
   print("done running atlas_compute/process#", i, " at ", my_time)
   proc.stdin.flush()
   print("first line of process ", i, " written at ", time.ctime())
   while True:   #read until line with "Value"
      line = proc.stdout.readline().decode('ascii')
      if "int_part" in line:
         break
      line=re.sub(r"[^\S\r\n]","",line)  #remove spaces but not newlines
      f.write(line.encode())
   print("Done with dimension ", dim, " i=", i)

def main(argv):
   global group
   global dim
   global data_dir
   global procs
   opts, args = getopt.getopt(argv, "g:c:d:D:")
   if len(opts)==0:
      print("Usage: \n-g: group\n-d: dimension\n-D: data directory")
      exit()
   for opt, arg in opts:
       if opt in ('-g'):
          group=arg
          data_dir=group #default, overridden by -D
       elif opt in ('-d'):
          dim=int(arg)
       elif opt in ('-D'):
          data_dir=arg
   if not os.path.exists(data_dir):
      print("making data directory: " + data_dir)
      os.makedirs(data_dir)
   my_time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
   print("starting facet computation: ", my_time)
   print("group: ", group)
   print("dimension: ", dim)
   print("data directory: ", data_dir)
   #run atlas once to determine number of fundamental facets
   proc=subprocess.Popen(["../atlas","facetious_jeff.at"], stdin=PIPE,stdout=PIPE)
   atlas_arg='{}'.format("\n prints(#facets_fundamental(" + group + "," + str(dim) + "))\n").encode('utf-8')
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   number_fundamental_facets =int(proc.stdout.readline().splitlines()[0])
   print("number of fundamental facets: ", number_fundamental_facets)
   proc.kill()
   procs=[]
   for i in range(number_fundamental_facets):
      proc=subprocess.Popen(["../atlas","facetious_jeff.at"], stdin=PIPE,stdout=PIPE)
      print("opened atlas process #", i)
      procs.append(proc)
   print("number of processes started: ", len(procs))
   with concurrent.futures.ProcessPoolExecutor(number_fundamental_facets) as P:
      for i in range(number_fundamental_facets):
         proc=procs[i]
         print("submiting process #",i, " to atlas_compute")
         Q=P.submit(atlas_compute,i)
   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



