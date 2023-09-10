#!/bin/python3

#test file of ostensibly non-unitary representations
#mainly intended for checking 67M non-unitary principal series for E8 calculated by Steve Miller

import sys, time, os, getopt, subprocess, gc,re
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp   #only for cpu count
progress_step_size=10 #how often to report progress

#simple reporting function
def report(q,i,proc):
   return("(i)[size of q]<status_i>: (" + str(i) + ")[:" + str(q.qsize()) +  "]<" + str(proc.poll()) + ">")

#call the atlas process proc, with arguments:
#q: queue of facet data, obtained from facet_file
#procs: array of processes
#i: number of process
def atlas_compute(i):
   print("starting atlas_compute process #", i)
   proc=procs[i]

   #first read in file, get N=total number of parameters
   atlas_arg='{}'.format("< " + input_file+ "\n").encode('utf-8')
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   while True:   #read until line with "Variable"
      line = proc.stdout.readline().decode('ascii')
      if "Variable" in line:
         break
      f.write(line)
   print("Done reading in file " + str(input_file) + "(job #" + str(i) + ")")
   
   output_file=directory + "/" + str(i)
   f=open(output_file,"w", buffering=1)
   f.write("Starting computation at " + str(time.ctime()) + "\n")
   f.write("Computing FPP_unitary_hash_one_level(G,hash) for G=" + group)
   line= "Job number: "  + str(i) + "\n" + "Number of jobs: " + str(number_jobs) + "\n" + "Number of parameters: " + str(number_parameters) + "\n" + "Number of parameters this job: " + str(N) + "\n"
   f.write(line)
   starttime=time.time()
   print("Reading in file  (job #" + str(i) + "): " + str(input_file))
   #get number of parameters in file
   atlas_arg='{}'.format("#parameters" + "\n").encode('utf-8')
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   N=int(line)
   print("N=",N)
   f.write("Number of parameters; " + str(N))
   atlas_arg ='{}'.format("set x=is_unitary_to_level_one(parameters,"  + str(i*N) + "," + "N" + ")" + "\n").encode('utf-8')
   print("Atlas arg: ", atlas_arg)
   
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()

   while True:   #read until line with "Variable"
      line = proc.stdout.readline().decode('ascii')
      if "Variable" in line:
         break
      line = line.strip() + "  " +  time.strftime("%H:%M:%S", time.gmtime(time.time()-starttime)) + "\n"
      f.write(line)
   print("Done computing nonunitary of parameters (job #" + str(i) + ")")

   stoptime=time.time()
   f.write("Finished computation at " + str(time.ctime()) + "\n")
   elapsed = time.strftime("%H:%M:%S", time.gmtime(stoptime-starttime))
   print("Finished ",group, "  time: ", elapsed)
   f.write("elapsed time: "+ elapsed)
   f.close()
   return()

def main(argv):
   global directory, number_jobs,input_file,group
   group="E8_s"
   opts, args = getopt.getopt(argv, "d:i:n:")
   if len(opts)<3:
      print("Usage: \n-d: directory\n-i: input_file\n-n: number jobs")
      exit()
   print("----------------------------")
   print("Starting at ", time.ctime())
   print("Computing nonunitary of parameters, on david_facets branch")
   cpu_count=mp.cpu_count()
   print("Number of cores: ", cpu_count)
   for opt, arg in opts:
       if opt in ('-d'):
          directory=arg
          print("directory: ", directory)
          if not os.path.exists(directory):
             os.makedirs(directory)
       elif opt in ('-i'):
          input_file=arg
          print("input_file: ", input_file)
       elif opt in ('-n'):
          number_jobs=int(arg)
          print("number_jobs: ", number_jobs)
   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   for i in range(number_jobs):
      proc=subprocess.Popen(["../atlas","jeff.at"], stdin=PIPE,stdout=PIPE)
      procs.append(proc)
#   print("number of processes: ", len(procs))
   with concurrent.futures.ProcessPoolExecutor(number_jobs) as P:
      for i in range(number_jobs):
         proc=procs[i]
         Q=P.submit(atlas_compute,i)
         T.append(Q)
#   for t in T:
#      try:
#         data = t.result()
#      except Exception as exc:
#         print("exception: ", exc)
#         print("data: ", data)
#         print('%r generated an exception: %s' % ("x",exc))
#         print("got data: ", data)
#      else:
#         print("OK")

if __name__ == "__main__":
   main(sys.argv[1:])



