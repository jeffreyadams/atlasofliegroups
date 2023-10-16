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
#procs: array of processes
#i: number of process
def atlas_compute(i,pid):
   print("starting atlas_compute process #", i, "pid: ", str(pid))
   proc=procs[i]

   data_file=directory + "/" + str(i) + ".at"
   log_file=directory + "/logs/" + str(i)
   data=open(data_file,"w", buffering=1)
   log=open(log_file,"w", buffering=1)
   log.write("{Starting computation at " + str(time.ctime()) + "\n")
   log.write("Computing FPP for " + group + "\n")
   log.write("Job number: " + str(i))

   starttime=time.time()
   #some initialization

   atlas_arg='{}'.format("set G=" + group + "\n").encode('utf-8')
#   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
   time.sleep(1)
   while True:
      print("TRUE")
      head_arg="head -1 " + input_file 
      print("head_arg:", head_arg)
      kgb_line=os.popen(head_arg).read()
      print("kgb_line: ", kgb_line + "pid=" + str(pid))
      if len(kgb_line)==0:
         print("no more kgb elements")
         break
      sed_arg="sed  1d -i " + input_file
      print("sed_arg: ", sed_arg)
      delete_line=os.popen(sed_arg)
      print("lines:")
#      atlas_arg='{}'.format("set list=FPP(" + group + "," + str(KGB_number) + ")" + "\n").encode('uft-8')
      aa="set list=FPP(" + str(group) +  str(KGB_number)
      print("aa=",aa)
      atlas_arg='{}'.format(aa).encode('uft-8')
      print("OKOKOK")
      print("atlas_arg: ", atlas_arg)
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      newtime=starttime
   print("NOW")
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
   global directory, number_jobs,group, start_KGB,end_KGB,input_file
   group="E8_s"
   start_KGB=0
   input_file=""
   end_KGB=1000  #should probably be KGB_size(G)
   opts, args = getopt.getopt(argv, "d:n:g:k:e:i:")
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
       elif opt in ('-i'):
          input_file=arg
       elif opt in ('-k'):
          start_KGB=arg
          print("starting KGB: ", start_KGB)
       elif opt in ('-e'):
          end_KGB=arg
          print("ending KGB: ", end_KGB)
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
      pid=proc.pid
      print("pid: ", pid)
#   print("number of processes: ", len(procs))
   with concurrent.futures.ProcessPoolExecutor(number_jobs) as P:
      for i in range(number_jobs):
         proc=procs[i]
         Q=P.submit(atlas_compute,i,proc.pid)
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



