#!/bin/python3

#test file of ostensibly non-unitary representations
#mainly intended for checking 67M non-unitary principal series for E8 calculated by Steve Miller

import sys, time, datetime, os, getopt, subprocess, gc,re
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
   f.write("Computing is_unitary_to_level(" + level + ",G,p)  for G=" + group + "\n")
   f.write("Level <0 <-> is_unitary(p) \n")

   #get N=total number of parameters from number of lines in input file
   N=int(re.search(r'\d+',(os.popen("wc -l " + input_file).read())).group())
   print("Number of all of parameters: ", N)
   #set n=number_jobs
   #define m,r: N=nm + r   (r<n)
   #number of jobs each process:  (m+1,m+1,...,m+1,m,...,m)  (r ones)
   #j   start_index(j)   end_index(j)
   #0   0                m
   #1   m+1              2m+1
   #j   j(m+1)           (j+1)m+j   (j\le r-1)
   #r-1 (r-1)(m+1)       rm+r-1    so far: rm+r-1+1=rm+r parameters
   #r   rm+r             (r+1)m+(r-1)
   #r+1 (r+1)m+r         (r+2)m+(r-1)
   #n-1 (n-1)m+r         nm+r-1     total: nm+r=N
   #formula:
   #start_index(i)=im+min(i,r)
   #end_index(i)  =(i+1)m+min(i,r-1)
   n=number_jobs
   m=N//n  #all jobs approximately same size
   r=N-n*m
   start_index=i*m + min(i,r)
   end_index=(i+1)*m + min(i,r-1)
   print(i, " ", start_index, " ", end_index)
   f.write("Number of all parameters: " + str(N) + "\n")
   f.write("Number of parameters per job: " + str(m) + "/" + str(m+1)  + "\n")
   number_parameters_this_job=end_index-start_index+1

   f.write("Number of parameters this job: " + str(number_parameters_this_job) + "\n")
   f.write("Starting index: " +  str(start_index) + "\n")
   f.write("Ending index: " + str(end_index) + "\n}\n")
   f.write("passed_" + str(i) + "=[\n")
   starttime=time.time()
   #some initialization
   atlas_arg='{}'.format("set G=" + group + "\n").encode('utf-8')
#   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
   atlas_arg_0='{}'.format("set comma=\"\"" + "\n").encode('utf-8')
   proc.stdin.write(atlas_arg_0)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
#   print("discard: ", line)
   sed_arg="sed -n \'" + str(start_index+1) + "," + str(end_index+1) + " p\' " + input_file
#   print("sed_arg: ", sed_arg)
   parameters_from_file =os.popen(sed_arg).read()
   parameters=parameters_from_file.splitlines()
   comma=""
   for j in range(len(parameters)):
      atlas_arg='{}'.format("let test=is_unitary_to_level(" + level + "," + parameters[j] + ") in prints(test, \" \", deformLookupCounter.use_count(), \" \", deformCalcCounter.use_count()) \n").encode('utf-8')
#      atlas_arg='{}'.format("let test=is_unitary_to_level(" + level + "," + parameters[j] + ") in prints(test) \n").encode('utf-8')
      print("atlas_arg: ", atlas_arg)
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()

      line = proc.stdout.readline().decode('ascii').strip()
      print("line: ", line)
      line = comma + "(" + str(start_index + j) + "," + parameters[j] + "," + line + ")"
#      line = line + "  {" +  time.strftime("%H:%M:%S", time.gmtime(time.time()-starttime)) + "}\n"
      line = line + "  {"  +str(datetime.timedelta(seconds=384374974+time.time()-starttime)) + "}\n"
      f.write(line)
      comma=","
   f.write("]\n")

   stoptime=time.time()
   f.write("{Finished computation at " + str(time.ctime()) + "}\n")
   elapsed = time.strftime("%H:%M:%S", time.gmtime(stoptime-starttime))
   print("Finished ",group, "  time: ", elapsed)
   f.write("elapsed time: "+ elapsed)
   f.close()
   return()

def main(argv):
   global directory, number_jobs,input_file,group, level
   group="E8_s"
   opts, args = getopt.getopt(argv, "d:i:n:g:l:")
   if len(opts)<3:
      print("Usage: \n-d: directory\n-i: input_file\n-n: number jobs\n-g: group\n-l: level  (level <0 <-> is_unitary(p))")
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
       elif opt in ('-g'):
          group=arg
       elif opt in ('-l'):
          level=arg
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



