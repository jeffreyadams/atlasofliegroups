#!/bin/python3

#test file of ostensibly non-unitary representations
#mainly intended for checking 67M non-unitary principal series for E8 calculated by Steve Miller

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
   if i%2==0:
      test_function="local_test_GEO_hash@(KGBElt, ratvec, Param_hash)"
   else:
      test_function="local_test_GEO_hash_one_level@(KGBElt, ratvec, Param_hash)"
   data_file=directory + "/" + str(i) + ".at"
   log_file=directory + "/logs/" + str(i) + ".txt"
   data=open(data_file,"w", buffering=1)
   log=open(log_file,"w", buffering=1)
   log.write("Starting computation at " + str(time.ctime()) + "\n")
   log.write("Computing FPP2 for " + group + "\n")
   log.write("Job number: " + str(i) + "\n")
   log.write("Job pid: " +  str(pid) + "\n")
   log.write("function: " +  test_function + "\n")
   starttime=time.time()
   #some initialization

   atlas_arg='{}'.format("set G=" + group + "\n").encode('utf-8')
   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
#   time.sleep(.1)
   include_header=True
   while kgb_queue.qsize()>0:
      print("q size: ", kgb_queue.qsize())
      log.write("size of kgb_queue: " + str(kgb_queue.qsize()) + "\n")
      x=kgb_queue.get()
      log.write("x=" + str(x) + "\n")
      log.write("KGB loop: "  + str(x) +  " size of queue: " + str(kgb_queue.qsize()))
#      atlas_arg_txt="set list=FPP2(" + group + "," + str(x) + "," + test_function + ")" + "\n"
      atlas_arg_txt="set list=FPP2(" + group + "," + str(x) + "," + test_function + ",unitary_hash)" + "\n"
#      print("atlas_arg: ", atlas_arg_txt)
      atlas_arg='{}'.format(atlas_arg_txt).encode('utf-8')
      proc.stdin.write(atlas_arg)
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
         line =line + "  " + time.strftime("%H:%M:%S", time.gmtime(timediff)) + "  " + time.strftime("%H:%M:%S", time.gmtime(newtime-starttime))
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
   print("DONE")

   stoptime=time.time()
   log.write("Finished computation at " + str(time.ctime()) + "\n")
   elapsed = time.strftime("%H:%M:%S", time.gmtime(stoptime-starttime))
   print("Finished ",group, "  time: ", elapsed)
   log.write("elapsed time: "+ elapsed)
   log.close()
   data.close()
   return()

def main(argv):
   global directory, number_jobs,group, start_KGB,end_KGB, kgb_queue, step_size
   group="E8_s"
   start_KGB=0
   end_KGB=1000  #should probably be KGB_size(G)
   step_size=1
   opts, args = getopt.getopt(argv, "d:n:g:k:e:s:")
   if len(opts)<4:
      print("Usage: \n-d: directory\n-n: number jobs\n-g: group\n-k: start KGB\n-e: end KGB\n-s: step size (each job does x=k,k+s,k+2s...)")
      exit()
   print("----------------------------")
   print("Starting at ", time.ctime())
   cpu_count=mp.cpu_count()
   print("Number of cores available: ", cpu_count)
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
          if not os.path.exists(directory + "/symlinks"):
             os.makedirs(directory + "/symlinks")
      elif opt in ('-s'):
          step_size=arg
      elif opt in ('-k'):
          start_KGB=arg
          print("starting KGB: ", start_KGB)
      elif opt in ('-e'):
          end_KGB=arg
          print("ending KGB: ", end_KGB)
      elif opt in ('-n'):
          number_jobs=int(arg)
          print("number_jobs: ", number_jobs)
   kgb_queue=Queue()
   print("Number of cores being used: ", number_jobs)
   kgb_list = [428,433,434,435,1089,3624,4439,4514,4547,4552,5445,5490,5491,6486,6526,
 6532,6598,7068,7106,7624,7636,7654,7655,7656,7817,8804,8806,8810,9962,9983,9993,
 9994,11154,11169,12330,12775,13189,13225,13409,13426,13464,14323,14436,14500,14518,
 14527,15410,15456,15505,16431,17267,17268,17269,17988,18011,18012,18656,18657,18658,
 19204,19658,20024,20027]

#   for i in range(int(start_KGB),int(end_KGB)):
#      kgb_list.append(i)
   kgb_number=len(kgb_list)  #12833 # size of kgb array
   print("kgb list (original order):")
   print(kgb_list)
   print("number of kgb elements: ", kgb_number)
   print("step size: ", step_size)
   gcd=math.gcd(kgb_number,int(step_size))
   print("gcd=",gcd)
   if gcd>1:
      print("Error: gcd(kgb_number,step_size)>1 (=",gcd,")")
      exit()
   kgb_list_reordered=[]  #just for debugging
   for i in range(len(kgb_list)):
#      entry=int(start_KGB)+i*int(step_size)%(int(end_KGB)-int(start_KGB))
      entry=kgb_list[i*int(step_size)%kgb_number]
      kgb_list_reordered.append(entry)
      kgb_queue.put(entry)
   print("size of queue: ", kgb_queue.qsize())
   print("order of kgb elements in queue: ")
   print(kgb_list_reordered)

#   for i in range(len(kgb)):
#      kgb_queue.put(kgb[i*23%kgbnumber])
#   print("queue size: ", kgb_queue.qsize())
#   print("Here is the kgb_queue:")
#   for i in range(kgb_queue.qsize()):
#      x=kgb_queue.get()
#      print(x)
   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   for i in range(number_jobs):
      symlinks_dir=directory + "/symlinks"
      atlas_cmd=symlinks_dir + "/atlas_" + str(i)
#      print("atlas_cmd: ",atlas_cmd)
      symlink_cmd="ln -s ../../../../atlas " + atlas_cmd
      print("symlink cmd: ", symlink_cmd)
      os.system(symlink_cmd)
      proc=subprocess.Popen([atlas_cmd,"FPP3.at"], stdin=PIPE,stdout=PIPE)
      procs.append(proc)
      pid=proc.pid
   with concurrent.futures.ProcessPoolExecutor(number_jobs) as P:
      for i in range(number_jobs):
         proc=procs[i]
         Q=P.submit(atlas_compute,i,proc.pid)
         T.append(Q)
   print("removing: ", symlinks_dir)
   shutil.rmtree(symlinks_dir)
   
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



