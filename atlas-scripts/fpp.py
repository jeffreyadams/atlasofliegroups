#!/bin/python3

import sys, time, datetime, os, getopt, subprocess, gc,re,shutil, math,psutil,signal
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp   #only for cpu count
from multiprocessing import Process, Queue

FPP_at_file="FPP.at" #default
FPP_settings_file=""
files_to_copy=[FPP_at_file,FPP_settings_file]

#unitary_hash_function="FPP_unitary_hash_bottom_layer"   #get this from an atlas variable

#max_memory=20000  #in megabytes: 20,000 = 20 gigabytes
#max_memory=50000   #50,000 megabyte = 50 gigabytes, x 250 jobs=12.5 terabytes
#max_memory=60000   #60,000 megabyte = 60 gigabytes, x 250 jobs=15 terabytes
max_memory=10000   #10,000 megabyte = 10 gigabytes, x 1000 jobs=10 terabytes

def nice_time(t):
   return(re.sub("\..*","",str(datetime.timedelta(seconds=t))))

def elapsed_time(start_time):
   return("[" + nice_time(time.time()-start_time) + "]\n")

#format atlas command, given as a text string, for passing to atlas via stdin.write
def format_cmd(atlas_command):
   return('{}'.format(atlas_command).encode('utf-8'))


#execute atlas command, read/write output to log file, until encountering final_string
def execute_atlas_command(atlas_cmd,final_string,my_log,my_proc):
   my_log.write("executing atlas_cmd: " + atlas_cmd + "\n")
   my_proc.stdin.write(format_cmd(atlas_cmd))
   my_proc.stdin.flush()
   while True:
      line=my_proc.stdout.readline().decode('ascii').strip()
      #my_log.write("LINE: " + line)
      if not line:
         my_log.write("[empty line]\n")
      if line.find(final_string)>=0:
         my_log.write("got termination line with " + final_string + ": " + line + "\n")
         my_log.write("finished executing: " + atlas_cmd + "\n")
         return(line)
      else:
         my_log.write(line + "\n")

#call the atlas process proc, with arguments:
#procs: array of processes
#i: number of process
def atlas_compute(job_number,pid):
   print("starting atlas_compute process #:",job_number, "pid: ", str(pid), " at ", str(time.ctime()))
   proc=procs[job_number]
   #print("output_dir: ", output_dir, flush=True)
   log_file=output_dir + "/logs/" + str(job_number) + ".txt"
   log=open(log_file,"w",buffering=1)
   log.write("test_var:" + test_var)
   round=test_var
   log.write("round:" + round)   
   log.write("Starting computation at " + str(time.ctime()) + "\n")
   log.write("Computing FPP for " + group + "\n")
   log.write("\nJob number: " + str(job_number) + "\n")
   log.write("Job pid: " +  str(pid) + "\n")
   vars=['algorithm_flag',
         'big_unitary_hash_flag',
         'bl_interrupt_flag',
         'blob_time_flag',
         'bl_step_flag',
         'branch_hash_flag',
         'by_zero_flag',
         'check_pols_flag',
         'coh_ind_all_flag',
         'coh_ind_most_flag',
         'coh_ind_report_flag',
         'def_flag',
         'deformed_hash_flag',
         'deform_flag',
         'Dirac_best_flag',
         'Dirac_best_local_flag',
         'Dirac_easy_flag',
         'Dirac_flag',
         'Dirac_max_flag',
         'Dirac_pro_flag',
         'drain_plug_flag',
         'dry_flag',
         'DV_K_type_flag',
         'easy_up_flag',
         'edge_report_flag',
         'every_KGB_flag',
         'every_lambda_deets_flag',
         'every_lambda_flag',
         'every_MBEF_flag',
         'fewer_reducible_unitary_flag',
         'FPP_lambda_table_flag',
         'FPP_report_flag',
         'ICPNU_flag',
         'ICPU_flag',
         'interrupt_flag',
         'is_FPP_unitary_flag',
         'KGB_frac_flag',
         'KGB_path_hash_flag',
         'kill_for_big_uTime_flag',
         'KNU_flag',
         'KU_flag',
         'KUKU_flag',
         'long_out_flag',
         'Lucas_fast_flag',
         'max_loc_flag',
         'min_after_flag',
         'mix_flag',
         'more_after_flag',
         'more_flag',
         'MvL_flag',
         'my_branch_flag',
         'my_formula_flag',
         'my_unitary_hash_flag',
         'next_taus_flag',
         'nuhash_flag',
         'old_K_type_flag',
         'old_local_test_flag',
         'old_local_testK_flag',
         'old_proj_flag',
         'old_test_interrupt_flag',
         'old_to_ht_flag',
         'one_level_revert_flag',
         'one_two_reverse_flag',
         'position_flag',
         'pre_bottom_flag',
         'pre_def_dumb_flag',
         'prefer_diff_flag',
         'prefer_gamma_flag',
         'print_big_B_flag',
         'print_int_flag',
         'proj_hash_flag',
         'Qs2B_flag',
         'Qs2_flag',
         'quick_flag',
         'quitNIH_flag',
         'real_flag',
         'recursive_tester_flag',
         'red_count_flag',
         'red_count_old_flag',
         'revert_flag',
         'R_packet_flag',
         'short_hts_flag',
         'short_mu_flag',
         'small_test_flag',
         'sort_LFD_flag',
         'step_flag',
         'test_bl_flag',
         'test_interrupt_flag',
         'test_neither_flag',
         'tilde_flag',
         'to_htB_flag',
         'two_flag',
         'u_Lucas_flag',
         'u_Lucas_wk_flag',
         'unip_factor_flag',
         'unip_flag',
         'Uparamhash_flag',
         'up_local_test_flag',
         'use_test_coroots_flag',
         'u_wk_flag',
         'u_wk_Lucas_flag',
         'wiggle_flag',
         'write_x_flag',
         'xl_alert_flag']

   log.write("some variables:\n")
   for var in vars:
      atlas_cmd="prints(\"" + var + ": \"," +  var + ")" + "\n"
      #log.write("atlas_cmd: " + atlas_cmd)
      proc.stdin.write(format_cmd(atlas_cmd))
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
      log.write(line)
   starttime=time.time()
   #some initialization
   log.write("\nsome initialization\n")
   atlas_cmd="set G=" + group + "\n"
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   log.write("defined group: " + line + "\n")
   proc.stdout.flush()
   reporting_data=[]
   atlas_cmd="prints(big_unitary_hash.uhash_sizes())\n"
   proc.stdin.write(format_cmd(atlas_cmd))
   log.write("atlas_cmd: " + atlas_cmd + "\n")
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii')
   log.write("starting size of unitary hashes (job number " + str(job_number) + "): " + line)
   log.write("xl_pairs_queue.qsize:" + str(xl_pairs_queue.qsize()) + "\n")
   while xl_pairs_queue.qsize()>0:
      log.write("==================================================================\n")
      log.write("\nJob number: " + str(job_number) + "\n")
      log.write("Job pid: " +  str(pid) + "\n")
      memory_rss=psutil.Process(pid).memory_info().rss / 1024**2 
      memory_vms=psutil.Process(pid).memory_info().vms / 1024**2
      log.write("current memory usage rss:" + str(memory_rss) + "\n")
      log.write("current memory usage vms:" + str(memory_vms) + "\n")
      log.write("elapsed time: " + elapsed_time(starttime))
      try:
         log.write("current queue size: " + str(xl_pairs_queue.qsize()) + "\n")
         log.write("try to get new (x,lambda) pair\n")
         x_lambda_number=xl_pairs_queue.get(block=True, timeout=.001)
         log.write("x_lambda_number: " + str(x_lambda_number) + "\n")
      except Exception as e:
         log.write("Exception\n")
      else:
         log.write("number of new (x,lambda) pair : " + str(x_lambda_number) + "\n")
         atlas_cmd="print_xl_pair(" + group  + " ," + str(x_lambda_number)  + ")\n"
         log.write("atlas_cmd: " + atlas_cmd)
         proc.stdin.write(format_cmd(atlas_cmd))
         proc.stdin.flush()
         line = proc.stdout.readline().decode('ascii').strip()
         vals=line.split(':')
         x_number=vals[0]
         lambda_=vals[1]
         log.write("(x,lambda)=(" + x_number + "," + lambda_ + ")\n")
         data_file=output_dir + "/" + str(job_number) + ".at"
         log.write("data file: " + data_file + "\n")
         x_lambda_start_time=time.time()
         atlas_cmd="prints(big_unitary_hash.uhash_sizes())\n"
         log.write("atlas_cmd: " + atlas_cmd + "\n")
         proc.stdin.write(format_cmd(atlas_cmd))
         proc.stdin.flush()
         line=proc.stdout.readline().decode('ascii')
         log.write("line: " + line + "\n")         
         atlas_cmd="set is_finished=is_finished(" + group + "," + str(x_lambda_number) + ")\n"
         log.write("atlas_cmd: " + atlas_cmd + "\n")
         proc.stdin.write(format_cmd(atlas_cmd))
         proc.stdin.flush()
         line=proc.stdout.readline().decode('ascii')
         log.write("line: " + line + "\n")
         atlas_cmd="prints(is_finished)" + "\n"
         log.write("atlas_cmd: " + atlas_cmd + "\n")
         proc.stdin.write(format_cmd(atlas_cmd))
         proc.stdin.flush()
         line = proc.stdout.readline().decode('ascii')
         log.write("line (value of is_finished): " + line)
         if "true" in line:
            log.write("(x,lambda) is already done\n Go to next pairs")
         elif "false" in line:
            log.write("(x,lambda) not done\n")
            atlas_cmd="set xl_pair=xl_pair(" + group + "," + str(x_lambda_number) + ")\nprints(new_line,\"got_pair\")" + "\n"
            execute_atlas_command(atlas_cmd,"got_pair",log,proc)

            atlas_cmd="print_xl_pair(" + group  + " ," + str(x_lambda_number)  + ")\n"
            log.write("atlas_cmd: " + atlas_cmd + "\n")
            proc.stdin.write(format_cmd(atlas_cmd))
            proc.stdin.flush()
            line = proc.stdout.readline().decode('ascii').strip()
            vals=line.split(':')
            x_number=vals[0]
            lambda_=vals[1]
            log.write("x: " + x_number + " lambda_: " + lambda_ + "\n")
            atlas_cmd="set list=FPP_unitary_hash_bottom_layer(xl_pair);prints(\"done_fpp\")"  + "\n"
            log.write("atlas_cmd: " + atlas_cmd + "\n")
            execute_atlas_command(atlas_cmd,"done_fpp",log,proc)
            log.write("Finished (" + str(x_lambda_number) + "," + x_number + "," + lambda_ + ")\n")
            #get rid of line: Variable list: void...
            line=proc.stdout.readline().decode('ascii').strip()
            #log.write("line: " + line)
         else:
            log.write("Got neither true nor false\n")
            log.write("GOT: " + line + "n")
         #log.write("finished true test\n")
         x_lambda_end_time=time.time()
         x_lambda_total_time=x_lambda_end_time-x_lambda_start_time
         atlas_cmd=">>\"" + data_file  + "\" write_one_pair(" + group + "," + str(x_lambda_number) + ")\n"
         log.write("atlas_cmd: " + atlas_cmd)
         proc.stdin.write(format_cmd(atlas_cmd))
         proc.stdin.flush()

         #log.write("Finished pair x_lambda_number:" + str(x_lambda_number) + "\n")
         log.write("elapsed time: " + nice_time(x_lambda_total_time) + "\ntime from start: " + nice_time(x_lambda_end_time-starttime) + "\n")
         reporting_data.append((job_number,round,x_lambda_number,x_number,lambda_,x_lambda_total_time))
         log.write("added to reporting data: " + str(x_lambda_number))
         log.write("end of loop\n")
         log.write("Get another pair from queue\n")
         proc.stdout.flush()
         #MEMORY TEST:
         #if this process is using more than max_memory: kill this process
         log.write("Memory Test")
         log.write("\nJob number: " + str(job_number) + "\n")
         log.write("Job pid: " +  str(pid) + "\n")
         memory_rss=psutil.Process(pid).memory_info().rss / 1024**2 
         memory_vms=psutil.Process(pid).memory_info().vms / 1024**2
         log.write("checking current memory usage rss:" + str(memory_rss) + "\n")
         log.write("checking current memory usage vms:" + str(memory_vms) + "\n")
         if memory_vms>max_memory:
            log.write("memory_vms>" + str(max_memory)+ ": killing process ")
            log.write("proc.pid: " + str(proc.pid) + "\n")
            log.write("Killing process (sending quit to atlas) " + str(proc.pid) + "\n")
            atlas_cmd="quit" + "\n"
            log.write("atlas_cmd: " + atlas_cmd)
            proc.stdin.write(format_cmd(atlas_cmd))
            proc.stdin.flush()
            log.write("killed the process\n")
            log.write("exiting\n")
            exit()
            log.write("didn't exit, still running...?\n")
         else:
            log.write("memory_vms<" + str(max_memory)+ ": process not killed\n")
            #go back to top of while x_lambdas_todo.qsize()>0 loop
            #lambda_queue.qsize=0: exit while x_lambdas_todo.qsize()>0 loop
   log.write("No more (x,lambda) pairs to do\ntime: " + str(time.ctime()) + "\n")
   log.write("report on times:\n")
   log.write("|job_number|round|pair number|x|lambda:time\n")
   for (job_number,round,x_lambda_number,x_number,lambda_,x_lambda_total_time) in reporting_data:
      log.write("|" + str(job_number) +"|" + round + "|" + str(x_lambda_number) + "|" + str(lambda_) + "|" + nice_time(x_lambda_total_time) +"\n")
#      log.write("xl_pair#:" + str(pair_number) + " time:" + str(nice_time(time)) + "\n" )
#      log.write(str(a) + " " + str(b) + " " + str(nice_time(c)) + "\n", flush=True)
#


   atlas_cmd="big_unitary_hash.uhash_size(" + group + ")\n"
   log.write("atlas_cmd: " + atlas_cmd)
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   line = proc.stdout.readline().decode('ascii').strip()

   if line.find("void")>=0:
      line = proc.stdout.readline().decode('ascii').strip()
      log.write("line: " + line)
   proc.stdout.flush()
   log.write("size of unitary_hash after all pairs: (job number " + str(job_number) +  "): " +  line + "\n")
   stoptime=time.time()
   elapsed = nice_time(stoptime-starttime)

   #         for (x,t) in reporting_data:
   #            log.write("x:" + str(x) + " " + nice_time(t) + "\n")
   #         log.write("Total time for " + str(x_lambda_count) + " x/lambda pairs: "+ elapsed + "\n")
   log.write("finished x_lambda queue\n")
   log.write("Killing process at " + str(time.ctime()) + "\n")
   log.close()
   data.close()
   print("killing process ",job_number, " at " + str(time.ctime()) + "\n")
   proc.kill()
   return()

#MAIN
def main(argv):
   global output_dir,group,group_name,xl_pairs_queue,round,test_var, executable_dir,data_directory, group_definition_file, init_file
   data_directory="./data"
   n_procs=1
   test_var=""
   run_index=-1
   round=""
   group=""
   logfile=""
   executable_dir="/.ccs/u02/jdada11/atlasSoftware/david_facets/"
   files_to_read=[]
   files_to_read_string=""
   default_file_to_read=""
   test_only=False
   init_file=""
   no_init=False
   queue_size=-1
   opts, args = getopt.getopt(argv, "g:d:n:l:x:i:G:t:m:I:h")
   for opt, arg in opts:
      if opt in ('-h'):
         print("\
         \nUsage: FPP.py -g group <options> ([defaults])\n",\
                "-G: group_name (alternate name suitable for files and directories) [group]\n",\
                "-d data directory [data]\n", \
                "-n number of cores [1]\n", \
                "-l logfile [data/group_name/logs/group_name.log.txt]\n", \
                "-i init_file [group_name._init.at]\n",\
                "-x other_files to load (can be used several times)\n", 
                "-m max queue size [# of (x,lambda) pairs]\n", \
                "-I: Ignore (don't load) any init file\n",\
                "-t: testing only (currently no op)\n",\
                "-h: this help file")
         exit();
      elif opt in ('-l'):
         logfile=arg
         print("log file set to: ", logfile)
      elif opt in ('-m'):
         queue_size=int(arg)
      elif opt in ('-I'):
         no_init=True
         print("ignoring init file")
      elif opt in ('-i'):
         if no_init==False:
            init_file=arg
            print("init file set to: ", init_file)
         else:
            print("-I is set, not using init file")
      elif opt in ('-n'):
         n_procs=int(arg)
         print("set number of processes: ", n_procs)
      elif opt in ('-g'):
         group=arg
         group_name=arg
         print("group= " + group)
      elif opt in ('-G'):
         group_name=arg
         print("set group_name: ", group_name)
      elif opt in ('-d'):
         data_directory=arg
         print("parent directory: " + data_directory)
      elif opt in ('-x'):
         files_to_read.append(arg)
      elif opt in ('-t'):
         test_only=True
         print("testing only")
      elif opt in ('-x'):
         files_to_read.append(arg)
   if group=="":
      print("group must be defined with -g (-h for help)");
      exit()
   existing = []
   pattern = re.compile(rf"^{re.escape(group_name)}.*?(\d+)$")
   #print("pattern: ", pattern)

   for filename in os.listdir(data_directory):
      match = pattern.match(filename)
      if match:
         # Extract the integer and store it with the full filename
         integer_value = int(match.group(1))
         existing.append(integer_value)
   if len(existing)==0:
      run_index=1
   else:
      run_index=max(existing)+1
   test_var=str(run_index)
   round=str(run_index)
   print("round: ", round)
   if len(init_file)==0 and no_init==False:
      init_file=group_name + "_init.at"
      print("init file is now: ", init_file)
   output_dir=data_directory + "/" + group_name + "_" + round
   try:
      os.mkdir(output_dir)
      os.mkdir(output_dir + "/logs")
      os.mkdir(output_dir + "/symlinks")
   except OSError as e:
      print(f"Error creating directory: {e}")
   if len(logfile)==0:
      logfile=output_dir + "/logs/" + group_name + "_log_" + round + ".txt"
   print("logfile: ", logfile)
   files_to_read.insert(0,init_file)
   print("files to read: ", files_to_read)
   for file in files_to_read:
      files_to_read_string +=" -x " + file
   print("files to read string: ", files_to_read_string)
   print("Starting at ", time.ctime())
   #print("command line 2: ", " ".join(sys.argv))
   #log_file=output_dir + "/logs/logfile.txt"
   print("redirecting output to ", logfile)
   sys.stdout = open(logfile, 'w')
   sys.stderr = sys.stdout
   print("Starting at ", time.ctime())
   print("command line: ", " ".join(sys.argv))
   print("Number of cores being used: ", n_procs,flush=True)
   #if no_init==True: don't load the init file, and don't initialize the hash
   #if no_init==False:
      #if init_file exists, load it
      #otherwise initialize big_unitary_hash (bad idea for big groups)
   proc=None
   print("proc None",flush=True)
   print("init_file: ", init_file,flush=True)
   #if init_file exists load it
   if len(init_file)>0 and os.path.exists(init_file):
      print("case 1: ", init_file, flush=True)
      atlas_cmd=executable_dir + "atlas"
      proc=subprocess.Popen([atlas_cmd,"all.at",init_file], stdin=PIPE,stdout=PIPE)
      print("Loaded init file ", init_file,flush=True)
      #proc.kill()
   else:
      print("init file is empty, generating big_unitary_hash", flush=True)
      atlas_cmd=executable_dir + "atlas"
      proc=subprocess.Popen([atlas_cmd,'all.at'], stdin=PIPE,stdout=PIPE)
      #atlas_cmd="big_unitary_hash.set_xl_sizes(" + group + ",[]);prints(\"OK2\")\n"
      atlas_cmd="big_unitary_hash.set_xl_sizes(" + group + ",[])\n"
      proc.stdin.write(format_cmd(atlas_cmd+ "\n"))
      proc.stdin.flush()
      print("set xl_sizes in big_unitary_hash\n", flush=True)
      line=proc.stdout.readline().decode('ascii').strip()
      print("IGNORING: ", line)
      #write the init file
      atlas_cmd=">> " + group_name + "_init.at  big_unitary_hash.write()"
      print("atlas: " + atlas_cmd + "\n")
      proc.stdin.write(format_cmd(atlas_cmd+ "\n"))
      proc.stdin.flush()
   print("done with ", flush=True)
   xl_pairs_queue=mp.Queue()
   #print("proc=",proc,flush=True)
   print("created xl_pairs_queue",flush=True)
   #get number of (x,lambda) pairs to make queue of that size
   #if -m was called use that instead (if not, queue_size=-1)
   if queue_size==-1:
      print("setting queue size using atlas")
      atlas_cmd="prints(big_unitary_hash.xl_sizes_cumulative(" + group + ")~[0])\n"
      print("atlas_cmd: ", atlas_cmd.strip(),flush=True)
      proc.stdin.write(format_cmd(atlas_cmd+ "\n"))
      proc.stdin.flush()
      line=proc.stdout.readline().decode('ascii').strip()
      re.sub(".* ","",line)
      print("line:", line,":", flush=True)
      queue_size=int(line)
   else:
      print("queue size set to " + str(queue_size) +  "(-m option)\n")
   for i in range(queue_size):
      xl_pairs_queue.put(i)
   print("xl_pairs_queue size is now: ", xl_pairs_queue.qsize(),flush=True)

   #if x_lambdas_todo_file is given (probably always the case)
   #read file into atlas, just to find the length of the array x_lambdas_todo 

   #launch free.py (monitoring processes)
   freeproc=subprocess.Popen(["./free.py","-d" + output_dir + "/logs"], stdin=PIPE,stdout=PIPE)
   print("----------------------------")
   cpu_count=mp.cpu_count()
   print("Number of cores available: ", cpu_count)
   for file in files_to_copy:
      if len(file)>0:
         print("Copying" + file + "to logs directory")
         shutil.copy(file,output_dir + "/logs")
   #run git log -n 1 to get last log entry and print it out
   print("git log -n 1:")
   git_cmd=['git','log','-n','1']
   proc=subprocess.Popen(git_cmd, stdin=PIPE,stdout=PIPE)
   proc.stdin.flush()
   gitlog=proc.stdout.read().decode('ascii')
   print(gitlog)

   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   #Write group definition file
   print("running one atlas job to write group definition file")
   atlas_cmd=executable_dir + "atlas" 
   print("atlas_cmd: " + atlas_cmd)
   proc=subprocess.Popen([atlas_cmd,"all.at"], stdin=PIPE,stdout=PIPE)
   atlas_cmd="<FPP.at\n"
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   atlas_cmd="write_real_form_plus(" + group + ",\"G_temp\")" + "\n"  #plus: in FPP.at: includes j line
   print("atlas_cmd: ",  atlas_cmd)
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   group_definition_file=data_directory + "/" + group_name + ".at"
   group_definition=open(group_definition_file,"w",buffering=1)
   group_definition.write("<groups.at\n")
   while True:
      line = proc.stdout.readline().decode('ascii').strip()
      group_definition.write(line + "\n")
      if line.find("rf_number")>=0:
         break
         group_definition.close()

   #launch atlas processes
   print("starting ", n_procs, " atlas processes")
   symlinks_dir=data_directory + "/" + group_name + "_" + round + "/symlinks"
   print("symlinks_dir: ", symlinks_dir, flush=True)
   print("going to make symlinks: ", n_procs,flush=True)
   for i in range(n_procs):
      atlas_cmd=symlinks_dir + "/atlas_" + str(i)
      symlink_cmd="ln -s " + executable_dir + "atlas " + atlas_cmd
      #print("symlink_cmd: ", symlink_cmd, flush=True)
      os.system(symlink_cmd)
      myarg=[atlas_cmd,"all.at"]+files_to_read
      #print("myarg: ",  myarg,flush=True)
      proc=subprocess.Popen(myarg, stdin=PIPE,stdout=PIPE)
      procs.append(proc)
      pid=proc.pid
   print("executing ", n_procs, " atlas functions", flush=True)
   with concurrent.futures.ProcessPoolExecutor(n_procs) as P:
      for i in range(n_procs):
         proc=procs[i]
         Q=P.submit(atlas_compute,i,proc.pid)
         print("submitted atlas job ", i, " ", proc.pid)
         T.append(Q)
         #print("TQ",flush=True)
   shutil.rmtree(symlinks_dir)
   print("done with calculation, removed symlinks\n")
   print("now updating init file")
   fpp_init()
   print("finished fpp_init, exiting")
   freeproc.kill()

def fpp_init():
   atlas_executable=executable_dir + "atlas" 
   arg=[executable_dir,"all.at"]

   #read the directories ./data/G2_s_j j=0,1,2,3 (for example)
   #also read G2_s_init.at if it exists
   dirs=[]
   for entry in os.listdir(data_directory):
      print("entry: ", entry)
      if os.path.isdir(os.path.join(data_directory,entry)) and entry.startswith(group_name):
         print("adding", entry)
         dirs.append(os.path.join(data_directory,entry))
         print("dirs: ", dirs)
         files_to_read=["all.at",group_definition_file]
         init_file=group_name + "_init.at"
         print("look for: " + init_file)
   if os.path.exists(init_file):
      files_to_read.append(init_file)
      print("adding " + init_file + " to list of files\n")
   for dir in dirs:
      print("dir: ", dir)
      files=os.listdir(dir)
      print("files: ", files)
      for file in files:
         if file.endswith("at"):
            files_to_read.append(dir + "/" + file)
            #print("dirs: ", dirs)
            #print("files_to_read: ", files_to_read)
            arg=[atlas_executable] + files_to_read
            print("arg: ",arg)
            proc=subprocess.Popen(arg,stdin=PIPE,stdout=PIPE)

   atlas_cmd="big_unitary_hash.set_xl_sizes(" + group +",[]);prints(\"\")\n"  #need prints to get line of output
   print("atlas_cmd: ", atlas_cmd, flush=True)
   proc.stdin.write(format_cmd(atlas_cmd))
   proc.stdin.flush()
   line=proc.stdout.readline().decode('ascii')
   print("line: ", line)

   print("sending output to ", init_file, "\n")
   #    atlas_cmd="> " + init_file + " big_unitary_hash.writeG("+ group + ")"
   atlas_cmd="> " + init_file + " big_unitary_hash.write()"
   print("atlas: " + atlas_cmd + "\n")
   proc.stdin.write(format_cmd(atlas_cmd+ "\n"))
   proc.stdin.flush()
   print("Exiting")
   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



