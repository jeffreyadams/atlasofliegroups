#!/bin/python

import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp   #only for cpu count
from multiprocessing import Process, Queue

progress_step_size=10 #how often to report progress

FPP_at_file="FPP.at" #default
FPP_py_file="FPP.py" #default
coh_ind_flag=True  #False: load from the file/True: don't load from the file
def nice_time(t):
   return(re.sub("\..*","",str(datetime.timedelta(seconds=t))))
 
#simple reporting function
def report(q,i,proc):
   return("(i)[size of q]<status_i>: (" + str(i) + ")[:" + str(q.qsize()) +  "]<" + str(proc.poll()) + ">")

#format atlas command, given as a text string, for passing to atlas via stdin.write
def format_cmd(atlas_command):return('{}'.format(atlas_command).encode('utf-8'))

#call the atlas process proc, with arguments:
#procs: array of processes
#i: number of process
def atlas_compute(i,pid):
   print("starting atlas_compute process #:",i, "pid: ", str(pid), " at ", str(time.ctime()))
   print("xlambdalists_flag: ", xlambdalists_flag)
   print("more_shift: ", more_shift)
   print("max_time_human (in minutes): ", max_time_human)
   print("max_time (in milliseconds) ", max_time)
   if xlambdalists_flag:
      print("Using xlambdalists")
   else:
      print("Using list of KGB elements")
   x=i  #this is a hack: should change occurences of x to i
   proc=procs[i]
   data_file=directory + "/" + str(i) + ".at"
   log_file=directory + "/logs/" + str(i) + ".txt"
   data=open(data_file,"w", buffering=1)
   log=open(log_file,"w", buffering=1)
   log.write("Starting computation at " + str(time.ctime()) + "\n")
   log.write("Computing FPP for " + group + "\n")
   log.write("(default) main file: " +  FPP_at_file + "\n")
   log.write("other at files: ")
   log.write("\nJob number: " + str(i) + "\n")
   log.write("Job pid: " +  str(pid) + "\n")
   log.write("function: " +  unitary_hash_function + "\n")
   atlas_cmd=format_cmd("set more_shift=" + str(more_shift) + "\n")
   proc.stdin.write(atlas_cmd)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
   atlas_cmd=format_cmd("set max_time=" + str(max_time) + "\n")
#   print(atlas_cmd)
   proc.stdin.write(atlas_cmd)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
   vars=['coh_ind_flag','revert_flag','deform_flag','every_KGB_flag','every_lambda_flag','every_lambda_deets_flag','facet_verbose','test_slightly_verbose','test_verbose','global_top','low_KGB_frac','bottom_length_frac','global_facets_file_loaded','coh_ind_file_loaded', 'more_flag', 'old_proj_flag','to_ht_frac']
   for var in vars:
      atlas_arg='{}'.format("prints(\"" + var + ": \"," +  var + ")" + "\n").encode('utf-8')
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
      log.write(line)
      time.sleep(.1)
   starttime=time.time()
   #some initialization
   proc.stdout.flush()
   atlas_arg='{}'.format("set G=" + group + "\n").encode('utf-8')
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
   queue_count=0
   reporting_data=[]
   print("starting Q size: ", main_queue.qsize())
   #add all parameters from coh_ind_params to unitary_hash
   if not coh_ind_flag:
      #here means: coh_ind_flag False; loading from the file
      print("loading params from coh_ind_params")
      atlas_cmd=format_cmd("update(unitary_hash,coh_ind_params)\n")
      proc.stdin.write(atlas_cmd)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')  #discard
      atlas_cmd=format_cmd("unitary_hash.size()\n")
      proc.stdin.write(atlas_cmd)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
   else:
      #otherwise: coh_ind_flag=True: not loading from the file
      print("not loading params from coh_ind_params")
   print("starting size of unitary_hash: (list number",str(i), ")", line)
   while main_queue.qsize()>0:
#      print("MAIN LOOP")
      log.write("==================================================================\n")
      log.write("size of main_queue: " + str(main_queue.qsize()) + "\n")
      atlas_cmd=format_cmd("unitary_hash.size()\n")
      proc.stdin.write(atlas_cmd)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii').strip()
      log.write("current size of unitary_hash: " + line + "\n")

      atlas_cmd=format_cmd("surprise_hash.size()\n")
      proc.stdin.write(atlas_cmd)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii').strip()
      log.write("Current size of surprise_hash: " + line + "\n")

      next_queue_entry=main_queue.get()
      if xlambdalists_flag:
         list_number=next_queue_entry
         print("list_number: ", list_number)
         log.write("list_number=" + str(list_number) + "\n")
         #get size of list
         atlas_cmd=format_cmd("#xlambdalists[" + str(list_number) + "]" + "\n")
         print("atlas_cmd: ", atlas_cmd)
         log.write("in xlambdalists_flag loop: ")
         proc.stdin.write(atlas_cmd)
         proc.stdin.flush()
         line = proc.stdout.readline().decode('ascii').strip()
         line = re.sub(".*:",'',line)
         log.write("size of list: "+ line + "\n")
      else:
         x=next_queue_entry
         log.write("x=" + str(x) + "\n")
      queue_count+=1
      start_time=time.time()

      #primary atlas function
      #if xlambdalists_flag: pass 2nd argument = xlambdalists, then 3rd argument is the listnumber
      if xlambdalists_flag:
         atlas_cmd=format_cmd("set list=FPP_unitary_hash_timeout(" + group + ",xlambdalists,"  + str(list_number) + ",functions,function_params" + ")\n")
         print("ATLAS CMD: ", atlas_cmd)
      print("executing ATLAS_CMD")
      proc.stdin.write(atlas_cmd)
      print("EXECUTED ATLAS_CMD")
      proc.stdin.flush()
      newtime=starttime
      while True:
         oldtime=newtime
         newtime=time.time()
         timediff=newtime-oldtime
         line = proc.stdout.readline().decode('ascii').strip()
         if line.find("Variable")>=0:
            break
         line =line + "  [" + nice_time(timediff) + "  " + nice_time(newtime-starttime) + "]"
         log.write(line + "\n")
      proc.stdout.flush()
      atlas_cmd=format_cmd("write_param_list_jda(list,\"p" + str(list_number) + "\");prints(\"done\")" + "\n")
      proc.stdin.write(atlas_cmd)
      proc.stdin.flush()
      while True:
         line = proc.stdout.readline().decode('ascii').strip()
#         print(line)
         if line.find("done")>=0:
            break
         data.write(line + "\n")
      job_time=newtime-start_time
      log.write("Finished_KGB element " + str(i)   + ": elapsed time: " + nice_time(job_time) + 
                "  time from start: " + nice_time(newtime-starttime) + "\n")
      reporting_data.append((list_number,job_time))
#      print("end of loop, Qsize is ", main_queue.qsize())
#      print("Get another item from queue")
      log.write("Get another item from queue\n")
   log.write("No more KGB elements to do; time: " + str(time.ctime()) + "\n")
   stoptime=time.time()
   elapsed = nice_time(stoptime-starttime)
   log.write("Times:\n")
   for (list_number,t) in reporting_data:
      log.write("Time for list number:" + str(list_number) + " " + nice_time(t) + "\n")
   log.write("Total time for " + str(queue_count) + " lists: "+ elapsed + "\n")
   #report on unitary_hash
   atlas_cmd=format_cmd("unitary_hash.size()\n")
   proc.stdin.write(atlas_cmd)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii').strip()
#   print("Size of unitary_hash: (list number ",str(i), ")", line)
   log.write("Size of unitary_hash: " + line + "\n")
   #report on suprise_hash
   atlas_cmd=format_cmd("surprise_hash.size()\n")
   proc.stdin.write(atlas_cmd)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii').strip()
#   print("Size of surprise_hash: (list number ",str(i), ")", line)
   log.write("Size of surprise_hash: " + line + "\n")
   
   log.write("Killing process at " + str(time.ctime()) + "\n")
   log.close()
   data.close()

   print("killing process ",i, " at " + str(time.ctime()) + "\n")
   proc.kill()
   print("process ", i, "killed, waiting")
   proc.wait()   #see https://stackoverflow.com/questions/41961430/how-to-cleanly-kill-subprocesses-in-python
   print("done waiting ",i)
   return()

def main(argv):
   global directory, number_jobs,group, start_KGB,end_KGB, main_queue, kgb_file, step_size,all_kgb, read_at_files,log_file, unitary_hash_function, xlambdalists_flag, global_facets_file_loaded, xlambdalists_file,more_shift,max_time, max_time_human
   unitary_hash_function="FPP_unitary_hash_bottom_layer"
   executable_dir="/.ccs/u02/jdada11/atlasSoftware/FPP_jeff/"
   group=""
   FPP_at_file="FPP.at"
   directory=""
   number_jobs=0
   start_KGB=0
   end_KGB=0
   kgb_number=0
   more_shift=5
   max_time_human=10   #minutes
   max_time=max_time_human*60*1000  #millliseconds
   all_kgb=False
   step_size=1
   kgb_file=""
   xlambdalists_file=""
   xlambdalists_flag=False
   kgb_list=[]
   dry_run=False
   read_at_files=False
   other_file=""
   extra_files=[]
   log_file=""
   opts, args = getopt.getopt(argv, "d:n:g:S:f:l:x:t:m:o:D")
   for opt, arg in opts:
      if opt in ('-d'):
         directory=arg
         if not os.path.exists(directory):
            os.makedirs(directory)
         if not os.path.exists(directory + "/logs"):
            os.makedirs(directory + "/logs")
         if os.path.exists(directory + "/symlinks"):
            shutil.rmtree(directory + "/symlinks")
         os.makedirs(directory + "/symlinks")
      elif opt in ('-l'):
         log_file=directory + "/logs/" + arg
      elif opt in ('-g'):
          group=arg
      elif opt in ('-f'):
         FPP_at_file=arg
      elif opt in ('-D'):
         dry_run=True
      elif opt in ('-S'):
          step_size=arg
      elif opt in ('-m'):
         more_shift=arg
      elif opt in ('-t'):
         max_time_human=arg
         max_time=max_time_human*60*1000
      elif opt in ('-x'):
         xlambdalists_file=arg
         xlambdalists_flag=True
         print("xlambdalists file: ", xlambdalists_file)
         extra_files.append(xlambdalists_file)
      elif opt in ('-o'):
         other_file=arg   #one extra file: loads other at files and sets some variables
         extra_files.append(other_file)
         #example (e7other.at):
         #<global_facets_E7.at
         #<coh_ind_E7_to_hash.at
         #set functions = [FPP_unitary_hash_bottom_layer@(RealForm,[[(KGBElt,ratvec)]],int,int),
         #                 FPP_unitary_hash_bottom_layer@(RealForm,[[(KGBElt,ratvec)]],int,int),
         #                 FPP_unitary_hash2@(RealForm,[[(KGBElt,ratvec)]],int,int)]
         #set function_params=[(5,10),(15,20),(-1,30)]

      elif opt in ('-n'):
          number_jobs=int(arg)
          print("number_jobs: ", number_jobs)
   if log_file != "":
      print("redirecting output to ", log_file)
      sys.stdout = open(log_file, 'w')
      sys.stderr = sys.stdout
   else:
      print("output to terminal")
   print("Starting at ", time.ctime())
   print("command line: ", " ".join(sys.argv))
   print("unitary hash function: ", unitary_hash_function)
   print("extra files:")
   if len(extra_files)==0:
      print(" none")
   else:
      for file in extra_files:
         print("  ",file)
   print("running free.py")
   if group=="" or directory=="" or (xlambdalists_file=="" and len(kgb_list)==0) or number_jobs==0:
      print("Usage: \n-d: directory\n-n: number jobs\n-g: group\n-a: x: file containing xlambdalists\n-o: other at file to load\n-f: main atlas file (default FPP.at)\n-l: logfile (default: directory/logs/fpp_group.log)\n-S: step size (each job does x=s,s+S,s+2S...(mod n)), default is 1\n-D: dry run only\n-t: max_time (in minutes)\n-m: more_shift (default 5)")
      print("\n-g group, -d directory, -n number jobs are required")
      print("-s/-e (start/end kgb elements) or -k (kgb file) or -a (all kgb elements) or -x (xlambdalists_file) are required")
      exit()
   
 
   if all_kgb and kgb_file=="":
      print("Doing all KGB elements")
      #Run atlas twice to get [KGB_fixed] and [KGB_non_fixed]
      #KGB_fixed:  cross(w_0,x)=x
      #KGB_non_fixed:  one of each non-fixed pair
      atlas_cmd=executable_dir + "atlas"
      proc=subprocess.run([atlas_cmd,"FPP.at"], stdin=PIPE,stdout=PIPE)
      atlas_arg='{}'.format("kgb_fixed(" + group + ")\n").encode('utf-8')
#      print("atlas_arg: ", atlas_arg)
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii').strip()
      temp=re.sub("\].*","",re.sub(".*\[","",line));kgb_fixed=temp.split(',') if temp!= "" else []
      print("kgb_fixed: " +  line)
      print("number of fixed elements: ", len(kgb_fixed))
#      print("fixed elements: ", kgb_fixed)

      atlas_arg='{}'.format("kgb_non_fixed(" + group + ")\n").encode('utf-8')
#      print("atlas_arg: ", atlas_arg)
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii').strip()
      print("kgb_non_fixed: " +  line)
      temp=re.sub("\].*","",re.sub(".*\[","",line));kgb_non_fixed=temp.split(',') if temp!= "" else []
      print("number of non fixed pairs of elements: ", len(kgb_non_fixed))
      print("total: ",len(kgb_fixed) + 2*len(kgb_non_fixed))
      kgb_size=len(kgb_fixed) + 2*len(kgb_non_fixed)
      for a in kgb_fixed + kgb_non_fixed:
         kgb_list.append(int(a))
      print("kgb_list: ", kgb_list)
      print("KGB size: " + str(kgb_size))
      kgb_number=len(kgb_list)  #12833 # size of kgb array
   if (end_KGB>0):
      print("start KGB: ", start_KGB)
      print("end KGB: ", end_KGB)
      for i in range(int(start_KGB),int(end_KGB)):
         kgb_list.append(i)
      print("made kgb queue from ", start_KGB, " to ", end_KGB)
      print("kgb_list: ", kgb_list)
      print("KGB size: " + str(kgb_size))
      kgb_number=len(kgb_list)  #12833 # size of kgb array
   print("Number of cores being used: ", number_jobs)

   print("number of kgb elements: ", kgb_number)

   freeproc=subprocess.Popen(["./free.py","-d" + directory + "/logs"], stdin=PIPE,stdout=PIPE)   
   print("----------------------------")
   cpu_count=mp.cpu_count()
   print("Number of cores available: ", cpu_count)
   shutil.copy(FPP_at_file,directory + "/logs")
   print("Copied " + FPP_at_file + " to logs directory")
   shutil.copy(FPP_py_file,directory + "/logs")
   if other_file != "":
      shutil.copy(other_file,directory + "/logs")
      print("Copied " + other_file + " to logs directory")
   print("Copied " + FPP_py_file + " to logs directory")
#   all_files=" ".join([FPP_at_file] + extra_files)
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

#make main _queue from kgb_list OR xlambdalists
   main_queue=Queue()
   xlambdalists_queue=Queue()
   if xlambdalists_flag:
      print ("running atlas once to get xlambdalists_size from file: ", xlambdalists_file)
      atlas_cmd=executable_dir + "atlas" 
      proc=subprocess.Popen([atlas_cmd,xlambdalists_file], stdin=PIPE,stdout=PIPE)
      atlas_arg=format_cmd("prints(#xlambdalists)" + "\n")
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii').strip()
      xlambdalists_size=int(line)
      print("Size of xlambdalists: ", xlambdalists_size)
      for i in range(xlambdalists_size):
         main_queue.put(i)
      print("size of xlambdalists_queue: ", xlambdalists_queue.qsize())
      
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
      main_queue.put(entry)
   print("size of queue: ", main_queue.qsize())
   print("kgb list (original order):")
   print(kgb_list)
   print("order of kgb elements in queue: ")
   print(kgb_list_reordered)
   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   #Run atlas once to write the file defining the group
   print("running one atlas job to write group definition file")
   atlas_cmd=executable_dir + "atlas" 
   #   proc=subprocess.Popen([atlas_cmd,"writeFiles.at"], stdin=PIPE,stdout=PIPE)
   proc=subprocess.Popen([atlas_cmd,"writeFiles.at"], stdin=PIPE,stdout=PIPE)
   atlas_arg='{}'.format("write_real_form(" + group + ",\"G\")" + "\n").encode('utf-8')
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   group_definition_file=directory + "/" + "the_group.at"
   group_definition=open(group_definition_file,"w",buffering=1)
   group_definition.write("<groups.at\n")
   while True:
      line = proc.stdout.readline().decode('ascii').strip()
      group_definition.write(line + "\n")
      if line.find("real_form")>=0:
         break
   group_definition.close()
   #      proc.kill()
   if dry_run:
      print("Dry run only: exiting")
      exit()
   print("starting ", number_jobs, " atlas processes")
   symlinks_dir=directory + "/symlinks"
   for i in range(number_jobs):
      atlas_cmd=symlinks_dir + "/atlas_" + str(i)
      #      print("atlas_cmd: ", atlas_cmd)
      symlink_cmd="ln -s " + executable_dir + "atlas " + atlas_cmd
      #      print("symlink cmd: ", symlink_cmd)
      os.system(symlink_cmd)

      print("all_files: ", all_files)
      myarg=[atlas_cmd]+all_files
      #     proc=subprocess.Popen([atlas_cmd,"writeFiles.at"], stdin=PIPE,stdout=PIPE)
      #      proc=subprocess.Popen(myarg stdin=PIPE,stdout=PIPE)
      #      proc=subprocess.Popen([atlas_cmd,"all.at"], stdin=PIPE,stdout=PIPE)
      proc=subprocess.Popen(myarg, stdin=PIPE,stdout=PIPE)
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
   freeproc.kill()
   print("Ending computation at " + str(time.ctime()) + "\n")
   exit()

if __name__ == "__main__":
   main(sys.argv[1:])


