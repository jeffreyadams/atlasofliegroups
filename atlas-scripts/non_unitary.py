#!/bin/python3

#output to a valid atlas file, of the form
#K_chars=[...]
#sample usage: Kpols.py -g G2_s -d 2 -f "facets_G2/facets_G2_s_dim
#compute number of fundamental facets of dimension d, then look for files
#facets_G2/facets_G2_s_dim2_i (0\le i\le d-1)
#command line arguments:
#-g: group
#-d: dimension (for naming output files)
#-f: facet_file_stem
#-o: output_file_stem (optional: override default naming convention)
#default output files: directory "group"/"group_dim" + dim + "_" + j + ".at"
#example: G2_s/G2_s_dim1_3.at

#to create a facet file (example):
#in atlas:
#>facetsE6dim4 facets(E6_s,4)   {4 dim facets}
#>facetsE6dim4 facets(E6_s)     {all facets}
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
   N=number_parameters/number_jobs
   proc=procs[i]
   output_file=directory + "/" + str(i)
#   f=open(output_file,"w", buffering=4000000)
   f=open(output_file,"w", buffering=1)
   f.write("Starting computation at " + str(time.ctime()) + "\n")
   f.write("Computing FPP_unitary_hash_one_level(G,hash) for G=E8_s\n")
   line= "Job number: "  + str(i) + "\n" + "Number of jobs: " + str(number_jobs) + "\n" + "Number of parameters: " + str(number_parameters) + "\n" + "Number of parameters this job: " + str(N) + "\n"
   f.write(line)
   starttime=time.time()
   print("Reading in file  (job #" + str(i) + "): " + str(input_file))

   atlas_arg='{}'.format("< " + input_file+ "\n").encode('utf-8')
#   print("atlas arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   while True:   #read until line with "Variable"
      line = proc.stdout.readline().decode('ascii')
      if "Variable" in line:
         break
      f.write(line)
   print("Done reading in file (job #" + str(i) + ")")

   atlas_arg='{}'.format("prints(\"number of parameters in file: \", #parameters)" + "\n").encode('utf-8')
#   print("atlas arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
#   print("line="+line)   
   f.write(line)
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
   global directory, number_jobs,number_parameters,input_file
   opts, args = getopt.getopt(argv, "d:i:n:N:")
   if len(opts)<3:
      print("Usage: \n-d: directory\n-i: input_file\n-n: number jobs\n-N:number of parameters to test")
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
       elif opt in ('-N'):
          number_parameters=int(arg)
          print("number_parameters: ", number_parameters)
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



