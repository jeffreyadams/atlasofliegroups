#!/bin/python3

import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp   #only for cpu count
from multiprocessing import Process, Queue

FPP_at_file="FPP.at" #default
FPP_py_file="FPP.py" #default
#unitary_hash_function="FPP_unitary_hash_bottom_layer"   get this from an atlas variable

def nice_time(t):
   return(re.sub("\..*","",str(datetime.timedelta(seconds=t))))

def elapsed_time(start_time):
   return("[" + nice_time(time.time()-start_time) + "]\n")

#format atlas command, given as a text string, for passing to atlas via stdin.write
def format_cmd(atlas_command):
   return('{}'.format(atlas_command).encode('utf-8'))

#call the atlas process proc, with arguments:
#procs: array of processes
#i: number of process
def atlas_compute(i,pid):
   print("starting atlas_compute process #:",i, "pid: ", str(pid), " at ", str(time.ctime()))
   proc=procs[i]
   log_file=directory + "/logs/" + str(i) + ".txt"
   log=open(log_file,"w",buffering=1)
   log.write("Starting computation at " + str(time.ctime()) + "\n")
   log.write("Computing FPP for " + group + "\n")
   log.write("(default) main file: " +  FPP_at_file + "\n")
   log.write("other at files: ")
   if len(extra_files)==0:
      log.write("none")
   else:
      for file in extra_files:
         log.write("\n  " + file)
   log.write("\nJob number: " + str(i) + "\n")
   log.write("Job pid: " +  str(pid) + "\n")
   atlas_cmd="prints(unitary_hash_function)"
#   log.write("atlas_cmd: " + atlas_cmd + "\n")
   proc.stdin.write(format_cmd(atlas_cmd + "\n"))
   proc.stdin.flush()
   line=proc.stdout.readline().decode('ascii')
   log.write("unitary hash function line: " + line + "\n")
   line=line[:-1]
   log.write("unitary has function line revised: " + line + "\n")
   unitary_hash_function = line
   vars=['test_verbose','low_frac','unip_flag','load_unipotents_flag','every_lambda_flag','every_lambda_deets_flag','every_KGB_flag','by_zero_flag','more_shift_level','more_flag','unitary_hash_function','interrupt_flag', 'bl_interrupt_flag','bl_step_size']

   log.write("some variables:\n")
   for var in vars:
      atlas_cmd="prints(\"" + var + ": \"," +  var + ")" + "\n"
#      log.write("atlas_cmd: " + atlas_cmd)
      proc.stdin.write(format_cmd(atlas_cmd))
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
      log.write(line)
   starttime=time.time()
   #some initialization
   log.write("some initialization\n")
   atlas_cmd="set G=" + group + "\n"
   #log.write("atlas_cmd: " + atlas_cmd)
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   log.write("defined group: " + line + "\n")
   proc.stdout.flush()
   proc.stdin.write(format_cmd(atlas_cmd))
   #log.write("atlas_cmd: " + atlas_cmd + "\n")
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   log.write("hash: " + line + "\n")
   reporting_data=[]
   atlas_cmd="prints(big_unitary_hash.uhash_sizes())\n"
   proc.stdin.write(format_cmd(atlas_cmd))
   log.write("atlas_cmd: " + atlas_cmd + "\n")
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   print("starting size of unitary_hash: (list number ",str(i), "): ", line)
   log.write("starting size of unitary_hash: (list number " + str(i) + "): " + line + "\n")

   log.write("loading file: " + xlambdas_todo_file + "\n")   
   atlas_cmd="< " + xlambdas_todo_file + "\n"
   log.write("atlas_cmd: " + atlas_cmd + "\n")
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   #log.write("discarding line 0: " + line)
   line = proc.stdout.readline().decode('ascii')
   #log.write("discarding line 00: " + line)
   line = proc.stdout.readline().decode('ascii')
   #log.write("discarding line 000: " + line)
   log.write("loaded file\n")
   log.write("kgb_lambda_queue_size: " + str(kgb_lambda_queue.qsize()) + "\n")
   if kgb_lambda_queue.qsize()==0:
      log.write("KGB lambda queue was empty from the start\n")
      log.write("Killing process at " + str(time.ctime()) + "\n")
      log.close()
      data.close()
      print("killing process ",i, " at " + str(time.ctime()) + "\n")
      proc.kill()
      return()
   else:
      log.write("Now processing kgb_lambda_queue of size: " + str(kgb_lambda_queue.qsize()))
      data_file=directory + "/" + str(i) + "pair.at"
      log.write("data file: " + data_file)
      atlas_cmd="prints(unitary_hash_pair_function)"
      log.write("atlas_cmd: " + atlas_cmd)
      proc.stdin.write(format_cmd(atlas_cmd + "\n"))
      proc.stdin.flush()
      line=proc.stdout.readline().decode('ascii')
      line=line[:-1]
      unitary_hash_pair_function = line
      kgb_lambda_queue_count=0
      reporting_data=[]
      log.write("starting kgb_lambda_queue size: " + str(kgb_lambda_queue.qsize()) + "\n")
      atlas_cmd="prints(unitary_hash.size())\n"
      log.write("atlas_cmd: " + atlas_cmd)
      proc.stdin.write(format_cmd(atlas_cmd))
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
      log.write("Line: " + line)
      print("starting size of unitary_hash: (list number ",str(i), "): ", line)
      proc.stdout.flush()
      log.write("starting size of unitary_hash: (list number " + str(i) + "): " + line + "\n")
      while kgb_lambda_queue.qsize()>0:
         log.write("==================================================================\n")
         log.write(elapsed_time(starttime))
         log.write("size of kgb_lambda_queue 1: " + str(kgb_lambda_queue.qsize()) + "\n")
         try:
            log.write("try to get kgb_lambda_number\n")
            kgb_lambda_number=kgb_lambda_queue.get()
         except Queue.Empty:
            log.write("Exception\n")
            log.write("KGB lambda queue is empty\n")
         else:
            log.write("succeeded in getting kgb_lambda_number\n")
            log.write("new kgb_lambda_number: " + str(kgb_lambda_number) + "\n")
            kgb_lambda_queue_count+=1
            x_lambda_start_time=time.time()
            #get (kgb_number,lambda)
            atlas_cmd="set (kgb_number,lambda)=xlambdas_todo[" + str(kgb_lambda_number) + "]\n" 
            log.write("atlas_cmd: " + atlas_cmd)
            log.write("x:" + str(kgb_lambda_number) + "\n")
            proc.stdin.write(format_cmd(atlas_cmd))
            proc.stdin.flush()
            atlas_cmd="prints(\"kgb_number=\", kgb_number)" + "\n"
            log.write("atlas_cmd(1): " + atlas_cmd)
            #log.write("finished atlas_cmd" + "\n")
            proc.stdin.write(format_cmd(atlas_cmd))
            proc.stdin.flush()
            #log.write("finished atlas_cmd2" + "\n")
            while True:
               #log.write("In True loop" + "\n")
               line = proc.stdout.readline().decode('ascii').strip()
               #log.write("my line: " + line + "\n")
               if line.find("kgb_number=")>=0:
                  log.write("got line with kgb number: " + line + "\n")
                  break;
               else:
                  log.write("(skipping line): " + line + "\n")
            log.write("line: " + line)
            kgb_number=re.sub(".*=","",line)
            log.write("kgb_number: " + kgb_number + "\n")
            proc.stdout.flush()
            atlas_cmd="prints(\"lambda=\", lambda)" + "\n"
            log.write("atlas_cmd: " + atlas_cmd)
            proc.stdin.write(format_cmd(atlas_cmd))
            proc.stdin.flush()
            while True:
               line = proc.stdout.readline().decode('ascii').strip()
               if line.find("lambda=")>=0:
                  log.write("got line with lambda=: " + line + "\n")
                  break;
               else:
                  log.write("(skipping line): " + line + "\n")
            log.write("line: " + line + "\n")
            my_lambda=re.sub(".*=","",line)
            log.write("my_lambda: " + my_lambda + "\n")
            proc.stdout.flush()
            atlas_cmd="set list=" + unitary_hash_pair_function + "(" + group + ','  + kgb_number + "," + my_lambda + ")\n"
            log.write("atlas_cmd: " + atlas_cmd)
            log.write("main atlas_pair_arg: " + atlas_cmd)
            proc.stdin.write(format_cmd(atlas_cmd))
            proc.stdin.flush()
            proc.stdout.flush()
            x_lambda_end_time=time.time()
            x_lambda_total_time=x_lambda_end_time-x_lambda_start_time

            #use write_one_x from FPP.at, first to temp file
            #file_name_temp=directory + "/" + str(kgb_lambda_number) + ".pair.temp.at"
            #file_name_final=directory + "/" + str(kgb_lambda_number) + "pair.at"
            #atlas_cmd=">\"" + file_name_final  + "\" write_one_pair(big_unitary_hash," + group + "," + kgb_number + ","  +  my_lambda  +")\n"
            atlas_cmd=">>\"" + data_file  + "\" write_one_pair(big_unitary_hash," + group + "," + kgb_number + ","  +  my_lambda  +")\n"
            log.write("atlas_cmd: " + atlas_cmd)
            proc.stdin.write(format_cmd(atlas_cmd))
            proc.stdin.flush()

            log.write("Finished x lambda pair " + kgb_number + "," + my_lambda    + "\nelapsed time: " + nice_time(x_lambda_total_time) +
                      "\ntime from start: " + nice_time(x_lambda_end_time-starttime) + "\n")
            reporting_data.append((kgb_number,my_lambda,x_lambda_total_time))
            log.write("end of loop\n")
            log.write("Get another pair from queue\n")
            proc.stdout.flush()
         
      log.write("No more (x,lambda) pairs to do; time: " + str(time.ctime()) + "\n")
      while True:
         line = proc.stdout.readline().decode('ascii').strip()
         if line.find("Variable")>=0:
            log.write("got line with Variabla=: " + line + "\n")
            break;
         else:
            log.write("(skipping line): " + line + "\n")
      log.write("report on times:\n")
      for (a,b,c) in reporting_data:
         log.write(str(a) + " " + str(b) + " " + str(nice_time(c)) + "\n")
      atlas_cmd="big_unitary_hash.uhash_size(" + group + ")\n"
      log.write("atlas_cmd: " + atlas_cmd)
      proc.stdin.write(format_cmd(atlas_cmd))
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii').strip()
      if line.find("void")>=0:
         line = proc.stdout.readline().decode('ascii').strip()
      proc.stdout.flush()
      log.write("size of unitary_hash after all pairs: (job number " + str(i) +  "): " +  line + "\n")
      stoptime=time.time()
      elapsed = nice_time(stoptime-starttime)
      log.write("size of kgb_lambda_queue 2: " + str(kgb_lambda_queue.qsize()) + "\n")
      #         for (x,t) in reporting_data:
      #            log.write("x:" + str(x) + " " + nice_time(t) + "\n")
      #         log.write("Total time for " + str(kgb_lambda_count) + " x/lambda pairs: "+ elapsed + "\n")
      log.write("FINISHED kgb_lambda queue\n")
      log.write("Killing process at " + str(time.ctime()) + "\n")
      log.close()
      data.close()
      print("killing process ",i, " at " + str(time.ctime()) + "\n")
      proc.kill()
      return()

#MAIN
def main(argv):
   global directory, number_jobs,group,extra_files,kgb_lambda_queue, known_kgb, generate_kgb,xlambdas_todo_file
   message=""
   directory=""
   group=""
   number_jobs=0
   extra_files=[]
   kgb_list=[]
   known_kgb=[]
   kgb_file=""
   xlambdas_todo=[]
   xlambdas_todo_file=""
   generate_kgb=True
   executable_dir="/.ccs/u02/jdada11/david_facets/"

   opts, args = getopt.getopt(argv, "d:n:g:x:m:k:p:")
   for opt, arg in opts:
      if opt in ('-d'):
         directory=arg
         if os.path.exists(directory):
            print("directory  exists: " + directory)
            exit()
         else:
            os.makedirs(directory)
         if not os.path.exists(directory + "/logs"):
            os.makedirs(directory + "/logs")
         if os.path.exists(directory + "/symlinks"):
            shutil.rmtree(directory + "/symlinks")
         os.makedirs(directory + "/symlinks")
      elif opt in ('-g'):
          group=arg
      elif opt in ('-p'):
         xlambdas_todo_file=arg
         print("loading (x,lambda) pairs from: ", xlambdas_todo_file)
         print("loaded xlambdas_todo_file=", xlambdas_todo_file)
      elif opt in ('-n'):
         number_jobs=int(arg)
         print("number_jobs: ", number_jobs)
      elif opt in ('-m'):
         message=arg
         print("message: ", message)
      elif opt in ('-x'):
         extra_file=arg
         extra_files.append(extra_file)
      elif opt in ('-k'):
         kgb_file=arg
         generate_kgb=False
         print("not computing KGB (getting it from file)")
         print("file: ", kgb_file)
   if group=="" or directory=="" or number_jobs==0:
      print("group -g, directory -d, number_jobs -n are required")
      print("-m message, -x extra files, -k are optional")
      exit()
   print("Starting at ", time.ctime())
   print("command line 2: ", " ".join(sys.argv))
   log_file=directory + "/logs/logfile.txt"
   print("redirecting output to ", log_file)
   sys.stdout = open(log_file, 'w')
   sys.stderr = sys.stdout
   print("Starting at ", time.ctime())
   if len(message)>0:
      print(message)
   print("command line: ", " ".join(sys.argv))
   print("Number of cores being used: ", number_jobs)
   kgb_lambda_queue=Queue()
   print("created kgb_lambda_queue")
   if len(xlambdas_todo_file)>1:
      xlambdas_todo=[]
      print("loading file of (x,lambda) pairs:", xlambdas_todo_file)
      atlas_cmd=executable_dir + "atlas"
      print("atlas_cmd: " + atlas_cmd)
      proc=subprocess.Popen([atlas_cmd,xlambdas_todo_file], stdin=PIPE,stdout=PIPE)
      print("DONE")
      atlas_cmd="prints(#xlambdas_todo)\n"
      proc.stdin.write(format_cmd(atlas_cmd))
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii').strip()
      print("line: ", line)
      for i in range(int(line)):
         kgb_lambda_queue.put(i)
   print("size of kgb_lambda_queue 3: ",str(kgb_lambda_queue.qsize()))

   freeproc=subprocess.Popen(["./free.py","-d" + directory + "/logs"], stdin=PIPE,stdout=PIPE)
   print("----------------------------")
   cpu_count=mp.cpu_count()
   print("Number of cores available: ", cpu_count)
   shutil.copy(FPP_at_file,directory + "/logs")
   print("Copied " + FPP_at_file + " to logs directory")
   shutil.copy(FPP_py_file,directory + "/logs")
   print("Copied " + FPP_py_file + " to logs directory")
   print("extra_files: ", extra_files)
   all_files=[FPP_at_file] + extra_files
   print("all_files: ", all_files)
   print("atlas will load: " + " ".join(all_files))
   #run git log -n 1 to get last log entry and print it out
   print("git log -n 1:")
   git_cmd=['git','log','-n','1']
   proc=subprocess.Popen(git_cmd, stdin=PIPE,stdout=PIPE)
   proc.stdin.flush()
   gitlog=proc.stdout.read().decode('ascii')
   print(gitlog)

   #   exit()
   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   #Write group definition file
   print("running one atlas job to write group definition file")
   atlas_cmd=executable_dir + "atlas"
   print("atlas_cmd: " + atlas_cmd)
   proc=subprocess.Popen([atlas_cmd,"writeFiles.at","FPP.at"], stdin=PIPE,stdout=PIPE)
   atlas_cmd="<FPP.at\n"
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()   
   atlas_cmd="write_real_form_plus(" + group + ",\"G_temp\")" + "\n"  #plus: in FPP.at: includes j line
   print("atlas_cmd: ",  atlas_cmd)
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   group_definition_file=directory + "/" + "the_group.at"
   group_definition=open(group_definition_file,"w",buffering=1)
   group_definition.write("<groups.at\n")
   while True:
      line = proc.stdout.readline().decode('ascii').strip()
      group_definition.write(line + "\n")
      if line.find("rf_number")>=0:
         break
         group_definition.close()

   #launch atlas processes
   print("starting ", number_jobs, " atlas processes")
   symlinks_dir=directory + "/symlinks"
   for i in range(number_jobs):
      atlas_cmd=symlinks_dir + "/atlas_" + str(i)
#      print("atlas_cmd: " + atlas_cmd)
      symlink_cmd="ln -s " + executable_dir + "atlas " + atlas_cmd
#      print("atlas_cmd: " + atlas_cmd)
      os.system(symlink_cmd)
      myarg=[atlas_cmd]+all_files
      print("myarg: ",  myarg)
      proc=subprocess.Popen(myarg, stdin=PIPE,stdout=PIPE)
      procs.append(proc)
      pid=proc.pid
   print("executing ", number_jobs, " atlas functions")
   with concurrent.futures.ProcessPoolExecutor(number_jobs) as P:
      for i in range(number_jobs):
         proc=procs[i]
         Q=P.submit(atlas_compute,i,proc.pid)
         print("submitted atlas job ", i, " ", proc.pid)
         T.append(Q)
   shutil.rmtree(symlinks_dir)
   freeproc.kill()
   print("Ending computation at " + str(time.ctime()) + "\n")
   exit()


if __name__ == "__main__":
   main(sys.argv[1:])



