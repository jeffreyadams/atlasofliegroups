#!/bin/python3

#python script to compute K_characters in parallel
#output to a valid atlas file, of the form
#K_chars=[...]
#sample usage: Kpols.py -g G2_s -d 2 -f facetsG2.txt -c 10 
#command line arguments:
#-g: group
#-d: dimension (for naming output files)
#-f: facet_file
#-o: output_file_stem (optional: override default naming convention)
#-c: number of cores to use
#default output files: directory "group"/"group_dim" + dim + "_" + j + ".at"
#example: G2_s/G2_s_dim1_3.at

#to create a facet file (example):
#in atlas:
#>facetsE6dim4 facets(E6_s,4)   {4 dim facets}
#>facetsE6dim4 facets(E6_s)     {all facets}
import sys, time, os, getopt, subprocess, queue
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing   #only for cpu count

progress_step_size=100 #how often to report progress

#simple reporting function
def report(q,i,proc):
   return("(i)[size of q]<status_i>: (" + str(i) + ")[:" + str(q.qsize()) +  "]<" + str(proc.poll()) + ">")

#call the atlas process proc, with arguments:
#q: queue of facet data, obtained from facet_file
#procs: array of processes
#i: number of process
def atlas_compute(group,output_file_stem,q,procs,i):
   facets_computed_local=0  #facets computed by this process
   global facets_computed #facets computed by all processes
   proc=procs[i]
   outputfile=group + "/" + output_file_stem + "_" + str(i) + ".at"  #atlas requires quotes if filename includes a directory
   f=open(outputfile,"w")
   firstline=True
   while not q.empty():
      if q.qsize()%1000==0:
         print("size of q: "+ str(q.qsize()),end="\r")
      report(q,i,proc)
      try:
         facet = q.get(False)
      except queue.Empty:
         quit_arg='{}'.format("\n quit").encode('utf-8')
         proc.stdin.write(quit_arg)
         print("quitting", i)
         report(q,i,proc)
      else:
         atlas_arg='{}'.format("\n prints(\"(\"," + facet +  "[0]" + ",\",\",K_data(K_char(" + group + "," + facet + "))" + ",\")\")\n" ).encode('utf-8')
         z=proc.stdin.write(atlas_arg)
         proc.stdin.flush()
         line = proc.stdout.readline().decode('ascii')
         #first time write "set K_chars=["
         #all other times write ",(...)"  (no comma at end)
         if firstline:
            line="set K_chars=[" + line
            firstline=False
         else:
            line=","+line
         f.write(line)

      facets_computed+=1  #total facets computed by all processed
      facets_computed_local+=1
   #end of while loop, last line write "]"
   f.write("]")
   return(i,facets_computed_local)

def main(argv):
   global facets_computed
   facets_computed=0
   cpu_count=multiprocessing.cpu_count()
   print("----------------------------")
   print("Number of cores: ", cpu_count)
   opts, args = getopt.getopt(argv, "g:d:c:f:o:")
   output_file_stem_in=""
   if len(opts)==0:
      print("Usage: \n-g: group\n-d: dimension\n-f: facet_file\n-c: number_of_cores\n-o: output_file_stem (optional)")
      exit()
   for opt,arg in opts:
      print(opt, " ", arg)
   for opt, arg in opts:
       if opt in ('-g'):
          group=arg
       elif opt in ('-d'):
          dim=arg
       elif opt in ('-f'):
          facet_file=arg
       elif opt in ('-o'):
          output_file_stem_in=arg
       elif opt in ('-c'):
          max_cores=int(arg)
   print("group: ", group)
   print("reading facets from: ", facet_file)
   if output_file_stem_in !="":
      output_file_stem=output_file_stem_in
   else:
      output_file_stem=group + "_dim" + str(dim)
   if not os.path.exists(group):
      os.makedirs(group)
   q=queue.Queue()
   file=open(facet_file,"r")
   #read facet_file and put each entry on q
   data=file.read().splitlines()
   for d in data:
      q.put(d)
   global facets_input
   facets_input=q.qsize()
   print("number of facets in input file: ", facets_input)
   procs=[]
   T=[]   #array of results from the atlas processes
   print("starting ", max_cores, "atlas processes")
   for i in range(max_cores):
      proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE,stdout=PIPE)
      procs.append(proc)
   print("number of processes started: ", len(procs))
   P=concurrent.futures.ThreadPoolExecutor()
   for i in range(max_cores):
         proc=procs[i]
         T.append(P.submit(atlas_compute,group,output_file_stem,q,procs,i))
   total_facets=0
   print("list of process numbers/number of facets computed")
   for t in T:
      try:
         data = t.result()
      except Exception as exc:
         print('%r generated an exception: %s' % ("x",exc))
      else:
         (p,n)=data
         total_facets+=n
         print(data)
   print("number of facets in input file: ", facets_input)
   print("number of K_polynomials computed (sum of numbers from each process):", total_facets)
   print("number of K_polynomials computed (global counter): ", facets_computed)

#   print("status of processes: ")
#   for i in range(len(procs)):
#      print("i: ", i, procs[i].poll())
#   for p in procs:
#      print(proc.poll())
#   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



