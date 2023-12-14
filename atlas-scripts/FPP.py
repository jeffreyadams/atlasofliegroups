#!/bin/python3

import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp   #only for cpu count
from multiprocessing import Process, Queue
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
   test_function="FPP_unitary_hash"
   data_file=directory + "/" + str(i) + ".at"
   log_file=directory + "/logs/" + str(i) + ".txt"
   data=open(data_file,"w", buffering=1)
   log=open(log_file,"w", buffering=1)
   log.write("Starting computation at " + str(time.ctime()) + "\n")
   log.write("Computing FPP for " + group + "\n")
   log.write("Job number: " + str(i) + "\n")
   log.write("Job pid: " +  str(pid) + "\n")
   log.write("function: " +  test_function + "\n")
   vars=['coh_ind_flag','revert_flag','deform_flag','every_KGB_flag','every_lambda_flag','every_lambda_deets_flag','test_slightly_verbose','test_verbose']

   #for var in vars:
#      print("var=", var)
#      atlas_arg='{}'.format("prints(\"" + var + "\"," +  var + ")").encode('utf-8')
#      print("atlas_arg: ", atlas_arg)
#      proc.stdin.write(atlas_arg)
#      proc.stdin.flush()
#      line = proc.stdout.readline().decode('ascii')
#      log.write("my line=" + line)
#      print("wrote line: ", line)

   starttime=time.time()
   #some initialization

   atlas_arg='{}'.format("set G=" + group + "\n").encode('utf-8')
#   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
#   time.sleep(.1)
   include_header=True
   while kgb_queue.qsize()>0:
#      print("q size: ", kgb_queue.qsize())
      log.write("==================================================================\n")
      log.write("size of kgb_queue: " + str(kgb_queue.qsize()) + "\n")
      x=kgb_queue.get()
      log.write("x=" + str(x) + "\n")
      atlas_arg_txt="set list=" + test_function + "(" + group + "," + str(x) + ")\n"
#      print("atlas_arg: ", atlas_arg_txt)
      atlas_arg='{}'.format(atlas_arg_txt).encode('utf-8')
      proc.stdin.write(atlas_arg)
#      print("Done executing atlas_arg")
      proc.stdin.flush()
      newtime=starttime
      while True:
         oldtime=newtime
         newtime=time.time()
         timediff=newtime-oldtime
         line = proc.stdout.readline().decode('ascii').strip()
#         print("THE LINE: ", line, " x=", x)
         if line.find("Variable")>=0:
#            print("BREAK")
            break
         line =line + "  [" + time.strftime("%H:%M:%S", time.gmtime(timediff)) + "  " + time.strftime("%H:%M:%S", time.gmtime(newtime-starttime)) + "]"
         log.write(line + "\n")
#      print("DONE WITH FIRST LOOP")
      proc.stdout.flush()
      if include_header:
         header_arg="true"
         include_header=False
      else:
         header_arg="false"
      atlas_arg_txt="write_param_list_jda(list,\"p" + str(x) + "\"," + header_arg + ");print(\"done\")"+ "\n"
#      print("atlas_arg: ", atlas_arg_txt)
      atlas_arg='{}'.format(atlas_arg_txt).encode('utf-8')
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      #  exit()
      while True:
         line = proc.stdout.readline().decode('ascii').strip()
         if line.find("done")>=0:
            break
         data.write(line + "\n")
      log.write("Finished KGB element " + str(x)   +  "  [" + time.strftime("%H:%M:%S", time.gmtime(timediff)) + "  " + time.strftime("%H:%M:%S", time.gmtime(newtime-starttime)) + "]\n")
      
   stoptime=time.time()
   log.write("Finished computation at " + str(time.ctime()) + "\n")
   elapsed = time.strftime("%H:%M:%S", time.gmtime(stoptime-starttime))
   log.write("Finished " + group)
   log.write(" elapsed time: "+ elapsed + "\n")
   log.close()
   data.close()
   return()

def main(argv):
   global directory, number_jobs,group, start_KGB,end_KGB, kgb_queue, kgb_file, step_size
   executable_dir="~/atlasSoftware/FPP_jeff/"
   group="E8_s"
   directory=""
   number_jobs=0
   start_KGB=0
   end_KGB=1000  #should probably be KGB_size(G)
   step_size=1
   kgb_file=""
   kgb_list=[]
   opts, args = getopt.getopt(argv, "d:n:g:k:e:s:S:")
   print("Starting at ", time.ctime())
   for opt, arg in opts:
      if opt in ('-g'):
          group=arg
          print("G=",arg)
      elif opt in ('-d'):
         directory=arg
         print("directory: ", directory)
         if not os.path.exists(directory):
            os.makedirs(directory)
         if not os.path.exists(directory + "/logs"):
            os.makedirs(directory + "/logs")
         #if symlinks dir exists, remove it, and make it again empty
         if os.path.exists(directory + "/symlinks"):
            shutil.rmtree(directory + "/symlinks")
         os.makedirs(directory + "/symlinks")
      elif opt in ('-S'):
          step_size=arg
      elif opt in ('-s'):
          start_KGB=arg
          print("starting KGB: ", start_KGB)
      elif opt in ('-e'):
         end_KGB=arg
         print("ending KGB: ", end_KGB)
      elif opt in ('-k'):
         kgb_file=arg
         print("loading KGB file ", kgb_file)
         data=open("kgb_file","r")
         while True:
            line=data.readline()
            if not line:
               break
            else:
               kgb_list.append(int(line.strip()))
         print("loaded kgb elements from file", kgb_file)
      elif opt in ('-n'):
          number_jobs=int(arg)
          print("number_jobs: ", number_jobs)
   if kgb_file=="":
      for i in range(int(start_KGB),int(end_KGB)):
         kgb_list.append(i)
      print("made kgb queue from ", start_KGB, " to ", end_KGB)
   print("Number of cores being used: ", number_jobs)
   kgb_number=len(kgb_list)  #12833 # size of kgb array
   print("number of kgb elements: ", kgb_number)
   if group=="" or directory=="" or len(kgb_list)==0 or number_jobs==0:
#      print(group, " ", directory, " ", len(kgb_list), " ", number_jobs)
      print("Usage: \n-d: directory\n-n: number jobs\n-g: group\n-s: start KGB\n-e: end KGB\n-k: file of KGB elements\n-S: step size (each job does x=s,s+S,s+2S...(mod n)), default is 1")
      print("\n-g group, -d directory, -n number jobs are required")
      print("-s/-e (start/end kgb elements) or -k (kgb file) are required")
      exit()
   print("----------------------------")
   cpu_count=mp.cpu_count()
   print("Number of cores available: ", cpu_count)
   shutil.copy("FPP.at",directory + "/logs")
   print("Copied FPP.at to logs directory")

#make kgb_queue from kgb_list
   kgb_queue=Queue()
   print("step size: ", step_size)
   gcd=math.gcd(kgb_number,int(step_size))
   print("gcd=",gcd)
   if gcd>1:
      print("Error: gcd(kgb_number,step_size)>1 (=",gcd,")")
      exit()
   kgb_list_reordered=[]  #just for debugging
   for i in range(len(kgb_list)):
      entry=kgb_list[i*int(step_size)%kgb_number]
      kgb_list_reordered.append(entry)
      kgb_queue.put(entry)
   print("size of queue: ", kgb_queue.qsize())
   print("kgb list (original order):")
   print(kgb_list)

   print("order of kgb elements in queue: ")
   print(kgb_list_reordered)
   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   print("starting ", number_jobs, " atlas processes")
   symlinks_dir=directory + "/symlinks"
   for i in range(number_jobs):
      atlas_cmd=symlinks_dir + "/atlas_" + str(i)
      symlink_cmd="ln -s " + executable_dir + "atlas " + atlas_cmd
#      print("symlink cmd: ", symlink_cmd)
      os.system(symlink_cmd)
      proc=subprocess.Popen([atlas_cmd,"FPP.at"], stdin=PIPE,stdout=PIPE)
      procs.append(proc)
      pid=proc.pid
   print("executing ", number_jobs, " atlas functions")
   with concurrent.futures.ProcessPoolExecutor(number_jobs) as P:
      for i in range(number_jobs):
         proc=procs[i]
         Q=P.submit(atlas_compute,i,proc.pid)
         T.append(Q)
   print("removing: ", symlinks_dir)
   shutil.rmtree(symlinks_dir)

   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



