#!/bin/python3

import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp   #only for cpu count
from multiprocessing import Process, Queue
progress_step_size=10 #how often to report progress

FPP_at_file="FPP.at" #default
FPP_py_file="FPP.py" #default
#extra_files=["global_facets_E7.at","coh_ind_E7_to_hash.at","E7except589to15851.at"]
extra_files=["global_facets_E7.at","coh_ind_E7_to_hash.at"]
extra_files=[]

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
   print("starting atlas_compute process #:",i, "pid: ", str(pid))
   proc=procs[i]
   test_function="FPP_unitary_hash_one_level"
   data_file=directory + "/" + str(i) + ".at"
   log_file=directory + "/logs/" + str(i) + ".txt"
   data=open(data_file,"w", buffering=1)
   log=open(log_file,"w", buffering=1)
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
   log.write("function: " +  test_function + "\n")
   vars=['coh_ind_flag','revert_flag','deform_flag','every_KGB_flag','every_lambda_flag','every_lambda_deets_flag','facet_verbose','test_slightly_verbose','test_verbose','global_top','global_facets_file_loaded','coh_ind_file_loaded','low_KGB_frac']

   for var in vars:
      atlas_arg='{}'.format("prints(\"" + var + ": \"," +  var + ")" + "\n").encode('utf-8')
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
      log.write(line)
      time.sleep(.1)
#   atlas_arg='{}'.format("prints(\"global_facets"  + ": \",#global_facets"  + ")" + "\n").encode('utf-8')
#   proc.stdin.write(atlas_arg)
#   proc.stdin.flush()
#   line = proc.stdout.readline().decode('ascii')
#   log.write(line)

   starttime=time.time()
   #some initialization

   atlas_arg='{}'.format("set G=" + group + "\n").encode('utf-8')
#   print("atlas_arg: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')  #discard this line
#   time.sleep(.1)
   kgb_count=0
   reporting_data=[]
   while kgb_queue.qsize()>0:
      log.write("==================================================================\n")
      log.write("size of kgb_queue: " + str(kgb_queue.qsize()) + "\n")
      x=kgb_queue.get()
      kgb_count+=1
      x_start_time=time.time()
      log.write("x=" + str(x) + "\n")
      # atlas_cmd=format_cmd("<\"" + directory + "/the_group.at\"" + "\n")
      # print("atlas_cmd: ", atlas_cmd)
      # proc.stdin.write(atlas_cmd)
      # proc.stdin.flush()
      # print("okok")
      # while True:
      #    print("in true loop")
      #    line = proc.stdout.readline().decode('ascii').strip()
      #    print("line: ",line)
      #    print("flag: ", line.find("Completely"))
      #    if line.find("Completely")>-1:
      #       print("DONE with loop")
      #       break
      # print("moving on")
      # #read all at files
      # print("IN LOOP")
      all_files=os.listdir(directory)
      atlas_cmd=format_cmd("prints(unitary_hash.size())" + "\n")
      proc.stdin.write(atlas_cmd)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
      print("TOP process id: ", i, "Size of unitary hash: ",line)
      log.write("TOP process id: " + str(i) +  "\nSize of unitary hash: " + line)
      print("files: ", all_files)
      params_to_add=[]
      
      for file in all_files:
         if re.match(r'[0-9]*\.at',file):
#            print("doing",file)
            at_file = directory + "/" + file;
            print("opening at file: (",i,") ",at_file)
            log.write("opening at file: " + at_file + "\n")
            at_data=open(at_file,"r")
#            print("data: ", at_data)
            while True:
#               print("IN LOOP")
               line=at_data.readline().strip()
#               print("line is: ", line)
               line=re.sub(".*:=","",line) 
               line=re.sub("\)[^,].*",")",line)  #get rid of {dual}
#               print("line is now: ", line)
               if len(line)==0:
#                  print("no lines, breaking")
                  break
               elif line.find("parameter")==-1:
                   ()
#                  print("no parameter in: ", line)
               else:
#                  print("got something:", line)
                  params_to_add.append(line)
      print("Done reading at files")
      print("number of new parameters to possibly add to unitary_hash:", len(params_to_add))
      log.write("number of new parameters to possibly add to unitary_hash:" + str(len(params_to_add)) + "\n")
      if len(params_to_add)>0:
#         print("adding params to unitary_hash")
         params_string="[" + ','.join(params_to_add) + "]"
#         print("params_string: ", params_string)
         atlas_cmd=format_cmd("update(unitary_hash," + params_string + ")" + "\n")
#         print("atlas_cmd: ", atlas_cmd)
         proc.stdin.write(atlas_cmd)
         proc.stdin.flush()
         line = proc.stdout.readline().decode('ascii')  #(number of new params,[Param])
         line=re.sub(".*Value: *\(","",line)
         line=re.sub(",.*","",line)
         print("number of parameters added to unitary_hash: ", line.strip())
         log.write("number of parameters added to unitary_hash: " + line)
         proc.stdout.flush()
#         print("atlas result: ", line)
         atlas_cmd=format_cmd("prints(unitary_hash.size())" + "\n")
#         print("atlas_cmd: ", atlas_cmd)
         proc.stdin.write(atlas_cmd)
         proc.stdin.flush()
         line = proc.stdout.readline().decode('ascii')
         print("[",i,"] Size of unitary hash: ",line)
         log.write("Size of unitary hash: " + line)

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
#         line =line + "  [" + time.strftime("%H:%M:%S", time.gmtime(timediff)) + "  " + time.strftime("%H:%M:%S", time.gmtime(newtime-starttime)) + "]"
         line =line + "  [" + nice_time(timediff) + "  " + nice_time(newtime-starttime) + "]"
         log.write(line + "\n")
#      print("DONE WITH FIRST LOOP")
      proc.stdout.flush()
#      print("OK")
#      atlas_arg_txt="write_param_list_jda(list,\"p" + str(x) + "\"," + str(x) +  "," + str(x) ;
      atlas_arg_txt="write_param_list_jda(list,\"p"  + str(x) + "\"," + str(x) + ");prints(\"done\")" + "\n"
#      print("atlas_arg: ", atlas_arg_txt)
      atlas_arg='{}'.format(atlas_arg_txt).encode('utf-8')
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      #  exit()
      while True:
         line = proc.stdout.readline().decode('ascii').strip()
#         print("line: ", line)
         if line.find("done")>=0:
            break
         data.write(line + "\n")
      x_time=newtime-x_start_time
      log.write("Finished_KGB element " + str(x)   + ": elapsed time: " + nice_time(x_time) + 
#                time.strftime("%H:%M:%S", time.gmtime(x_time)) +
                "  time from start: " + nice_time(newtime-starttime+10000000) + "\n")
#                time.strftime("%H:%M:%S", time.gmtime(newtime-starttime)) + "\n")
      reporting_data.append((x,x_time))
   log.write("No more KGB elements to do; time: " + str(time.ctime()) + "\n")
   stoptime=time.time()
#   elapsed = time.strftime("%H:%M:%S", time.gmtime(stoptime-starttime))
   elapsed = nice_time(stoptime-starttime)
   log.write("Times:\n")
#   print("data: ", reporting_data)
   for (x,t) in reporting_data:
#      log.write("x:" + str(x) + " " + time.strftime("%H:%M:%S", time.gmtime(t)) + "\n")
      log.write("x:" + str(x) + " " + nice_time(t) + "\n")
   log.write("Total time for " + str(kgb_count) + " KGB elements: "+ elapsed + "\n")
   log.close()
   data.close()
   return()

def main(argv):
   global directory, number_jobs,group, start_KGB,end_KGB, kgb_queue, kgb_file, step_size,all_kgb
   #   executable_dir="~/atlasSoftware/FPP_jeff/"
   executable_dir="/.ccs/u02/jdada11/atlasSoftware/FPP_jeff/"
   group=""
   FPP_at_file="FPP.at"
   directory=""
   number_jobs=0
   start_KGB=0
   end_KGB=0
   all_kgb=False
   step_size=1
   kgb_file=""
   kgb_list=[]
   dry_run=False
   opts, args = getopt.getopt(argv, "d:n:g:k:e:s:S:f:aD")
   print("Starting at ", time.ctime())
   print("command line: ", " ".join(sys.argv))
   print("extra files:")
   if len(extra_files)==0:
      print(" none")
   else:
      for file in extra_files:
         print("  ",file)
   for opt, arg in opts:
      if opt in ('-g'):
          group=arg
          print("G=",arg)
      if opt in ('-f'):
          FPP_at_file=arg
          print("main FPP file: " + FPP_at_file)
      elif opt in ('-D'):
         dry_run=True
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
      elif opt in ('-e'):
         end_KGB=arg
      elif opt in ('-a'):
         all_kgb=True
      elif opt in ('-k'):
         kgb_file=arg
         print("loading KGB file ", kgb_file)
         data=open("kgb_file","r")
         line=data.readline().decode('ascii')
         print("MY LINE: ", line)
         while True:
            line=data.readline()
            if not line:
               break
            else:
               kgb_list.append(int(line.strip()))
         print("loaded kgb elements from file", kgb_file)
         print("number of kgb elements: ", len(kgb_list))
      elif opt in ('-n'):
          number_jobs=int(arg)
          print("number_jobs: ", number_jobs)
   if all_kgb and kgb_file=="":
      print("Doing all KGB elements")
      #Run atlas twice to get [KGB_fixed] and [KGB_non_fixed]
      #KGB_fixed:  cross(w_0,x)=x
      #KGB_non_fixed:  one of each non-fixed pair
      atlas_cmd=executable_dir + "atlas"
      proc=subprocess.Popen([atlas_cmd,"FPP.at"], stdin=PIPE,stdout=PIPE)
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
   if kgb_file=="":
      print("start KGB: ", start_KGB)
      print("end KGB: ", end_KGB)
      for i in range(int(start_KGB),int(end_KGB)):
         kgb_list.append(i)
      print("made kgb queue from ", start_KGB, " to ", end_KGB)
   print("Number of cores being used: ", number_jobs)
   kgb_number=len(kgb_list)  #12833 # size of kgb array
   print("number of kgb elements: ", kgb_number)
   if group=="" or directory=="" or len(kgb_list)==0 or number_jobs==0:
#      print(group, " ", directory, " ", len(kgb_list), " ", number_jobs)
      print("Usage: \n-d: directory\n-n: number jobs\n-g: group\n-a: all kgb elements\n-s: start KGB\n-e: end KGB\n-k: file of KGB elements\n-S: step size (each job does x=s,s+S,s+2S...(mod n)), default is 1\n-D: dry run only")
      print("\n-g group, -d directory, -n number jobs are required")
      print("-s/-e (start/end kgb elements) or -k (kgb file) or -a (all kgb elements) are required")
      exit()
   print("----------------------------")
   cpu_count=mp.cpu_count()
   print("Number of cores available: ", cpu_count)
   shutil.copy(FPP_at_file,directory + "/logs")
   print("Copied " + FPP_at_file + " to logs directory")
   shutil.copy(FPP_py_file,directory + "/logs")
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
#      print("myarg: ", myarg)
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

   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



