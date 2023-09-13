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

   output_file=directory + "/" + str(i)
   f=open(output_file,"w", buffering=1)
   f.write("{Starting computation at " + str(time.ctime()) + "\n")
   f.write("Computing is_unitary_to_level_one(G,[Param])  for G=" + group + "\n")

    #first read in file, get N=total number of parameters
   atlas_arg='{}'.format("< " + input_file+ "\n").encode('utf-8')
#   print("atlas_arg 1: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   count=0
   while True:   #read until line with "Completely"
      line = proc.stdout.readline().decode('ascii')
      if "Completely" in line:
         break
#     f.write(line)
#      print(line)
      count += 1
      if count%1000==0:
         f.write(str(count) + "  " + line)
   print("Done reading in file " + str(input_file) + "(job #" + str(i) + ")")
   atlas_arg='{}'.format("#parameters" + "\n").encode('utf-8')
#   print("atlas_arg 2: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   N=int(re.sub('\D','',line))
   print("Total number of parameters: ", N)
   n=N//number_jobs
   f.write("Total number of parameters: " + str(N) + "\n")
   f.write("Number of parameters per job: " + str(n) + "\n")
   start_index=i*n
   number_parameters_this_job=min(n,N-(i)*n)
   end_index=start_index+number_parameters_this_job-1
   print("Number of parameters this job: ", number_parameters_this_job)
   f.write("Number of parameters this job: " + str(number_parameters_this_job) + "\n")
   f.write("Starting index: " +  str(start_index) + "\n")
   f.write("Ending index: " + str(end_index) + "\n")
   line= "Job number: "  + str(i)  + "\n" + "Number of jobs: " + str(number_jobs) + "\n" + "Number of parameters: " + str(N) + "\n" + "Number of parameters this job: " + str(number_parameters_this_job) + "\n" + "Starting index: " + str(start_index) + "\n" + "Ending index: " + str(end_index) + "\n"
   print("line=", line)
   f.write(line + "}\n")
   f.write("passed_" + str(i) + "=[\n")
   starttime=time.time()
   atlas_arg_0='{}'.format("set comma=\"\"" + "\n").encode('utf-8')
   proc.stdin.write(atlas_arg_0)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #get rid of one line
   atlas_arg_1="for j:" + str(number_parameters_this_job) + " from " + str(start_index)
   atlas_arg_2=" do let test=is_unitary_to_level_one(parameters[j]) in "
   atlas_arg_3=  "if  true then prints(comma + \"(\"+ j+ \",\"+ test.to_string + \")\");comma:=\",\" fi od" + "\n"
   atlas_arg=    atlas_arg_1 + atlas_arg_2 + atlas_arg_3 + "\n"
   atlas_arg='{}'.format(atlas_arg).encode('utf-8')
 #  print(atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   while True:   #read until line with "Value"
      line = proc.stdout.readline().decode('ascii')
#      print("MY LINE: ", line)
      if "Value" in line:
         break
      line = line.strip() + "  " +  time.strftime("%H:%M:%S", time.gmtime(time.time()-starttime)) + "\n"
      f.write(line)
#      print("line=", line)
#   print("NOW line: ", line)
   f.write("]\n")
   stoptime=time.time()
   f.write("{Finished computation at " + str(time.ctime()) + "}\n")
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
   exit()
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



