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

   data_file=directory + "/" + str(i) + ".at"
   log_file=directory + "/logs/" + str(i)
   data=open(data_file,"w", buffering=1)
   log=open(log_file,"w", buffering=1)
   log.write("{Starting computation at " + str(time.ctime()) + "\n")
   log.write("Computing FPP for " + group + "\n")
   atlas_arg='{}'.format("KGB_size(" + group + ")\n").encode('utf-8')
#   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   N=int(re.sub(".* ","",proc.stdout.readline().decode('ascii')))
#   print("N", str(N))
   
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
#   print("r=",r,"  m=", m, "i= ", i)
   start_index=int(start_KGB)+i*m + min(i,r)
   end_index=int(start_KGB)+(i+1)*m + min(i,r-1)
#   print(i, " ", start_index, " ", end_index)
   log.write("Number of KGB elements: " + str(N) + "\n")
   log.write("Number of KGB elements per job: " + str(m) + "/" + str(m+1)  + "\n")
   number_kgb_elements_this_job=int(end_index)-int(start_index)+1
   log.write("Number of KGB Elements this job: " + str(number_kgb_elements_this_job) + "\n")
   log.write("Starting index: " +  str(start_index) + "\n")
   log.write("Ending index: " + str(end_index) + "\n}\n")

   starttime=time.time()
   #some initialization

   atlas_arg='{}'.format("set G=" + group + "\n").encode('utf-8')
#   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
#   print("LINE: ", line)
   
#   print("discard: ", line)
   atlas_arg='{}'.format("set list=FPP(" + group + "," + str(start_index) + "," + str(number_kgb_elements_this_job) + ")"+ "\n").encode('utf-8')
#   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   newtime=starttime
   while True:
      oldtime=newtime
      newtime=time.time()
      timediff=newtime-oldtime
      line = proc.stdout.readline().decode('ascii').strip()
      line =line + "  " + time.strftime("%H:%M:%S", time.gmtime(timediff)) + "  " + time.strftime("%H:%M:%S", time.gmtime(newtime-starttime))
      if line.find("Variable")>=0:
         print("BREAK")
         break
      log.write(line + "\n")
#   f.write("]\n")
   atlas_arg='{}'.format("write_param_list_jda(##list,\"p" + str(i) + "\");print(\"done\")"+ "\n").encode('utf-8')
   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   while True:
      line = proc.stdout.readline().decode('ascii').strip()
#      print("my line: ", line)
      if line.find("done")>=0:
         break
      data.write(line + "\n")
   print("DONE")

   stoptime=time.time()
   log.write("{Finished computation at " + str(time.ctime()) + "}\n")
   elapsed = time.strftime("%H:%M:%S", time.gmtime(stoptime-starttime))
   print("Finished ",group, "  time: ", elapsed)
   log.write("elapsed time: "+ elapsed)
   log.close()
   data.close()
   return()

def main(argv):
   global directory, number_jobs,group, start_KGB
   group="E8_s"
   start_KGB=0
   opts, args = getopt.getopt(argv, "d:n:g:k:i:")
   if len(opts)<4:
      print("Usage: \n-d: directory\n-n: number jobs\n-g: group\n-k: start KGB\n")
      exit()
   print("----------------------------")
   print("Starting at ", time.ctime())
   cpu_count=mp.cpu_count()
   print("Number of cores: ", cpu_count)
   for opt, arg in opts:
       if opt in ('-d'):
          directory=arg
          print("directory: ", directory)
          if not os.path.exists(directory):
             os.makedirs(directory)
          if not os.path.exists(directory + "/logs"):
             os.makedirs(directory + "/logs")
       elif opt in ('-g'):
          group=arg
       elif opt in ('-k'):
          start_KGB=arg
          print("starting KGB: ", start_KGB)
       elif opt in ('-n'):
          number_jobs=int(arg)
          print("number_jobs: ", number_jobs)
   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   for i in range(number_jobs):
      proc=subprocess.Popen(["../atlas","FPP.at"], stdin=PIPE,stdout=PIPE)
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



