#!/bin/python3

#python script to compute K_characters in parallel
#sample usage: Kpols.py -g G2_s -f facetsG2.txt -o out -c 10
#command line arguments:
#-g: group
#-f: facet_file
#-o: output_file_root
#-c: max_cores
#this will launch max_cores processes
#process #j write output to: group/output_file_root_j
#example above: G2_s/out_1,..., G2_s/out_10

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
def atlas_compute(group,output_file_root,q,procs,i):
   count=0
   proc=procs[i]
   file_name=  group + "/" + output_file_root + "_" + str(i)
   while not q.empty():
      if q.qsize()%1000==0:
         print("size of q: "+ str(q.qsize()),end="\r")
      count+=1
      report(q,i,proc)
      try:
         facet = q.get(False)
      except queue.Empty:
         quit_arg='{}'.format("\n quit").encode('utf-8')
         proc.stdin.write(quit_arg)
#         print("quit job ", i)
         report(q,i,proc)
      else:
#         atlas_arg='{}'.format("\n  prints(\"(\",(" + str(i) + ")," + facet +  "[0]" + ",K_data(K_char(" + group + "," + facet + ")))").encode('utf-8')
         atlas_arg='{}'.format("\n  prints(\"(\"," + facet +  "[0]" + ",\",\",K_data(K_char(" + group + "," + facet + "))" + ",\")\")").encode('utf-8')

         print(atlas_arg)
#         atlas_arg='{}'.format("\n  prints(\"( + facet +  "[0]," + ",K_data(K_char(" + group + "," + facet + ")))").encode('utf-8')
         proc.stdin.write(atlas_arg)
   return(i,count)

def main(argv):
   cpu_count=multiprocessing.cpu_count()
   print("Number of cores: ", cpu_count)
   opts, args = getopt.getopt(argv, "g:c:f:o:")
   if len(opts)==0:
      print("Usage: \n-g: group\n-f: facet_file\n-o: output_file_root\n-c: number_of_cores\n")
      exit()
   for opt,arg in opts:
      print(opt, " ", arg)
   for opt, arg in opts:
       if opt in ('-g'):
          group=arg
       elif opt in ('-f'):
          facet_file=arg
       elif opt in ('-o'):
          output_file_root=arg
       elif opt in ('-c'):
          max_cores=int(arg)
   print("group: ", group)
   print("facet_file: ", facet_file)
   print("output_file_root: ", output_file_root)
   if not os.path.exists(group):
      os.makedirs(group)
   q=queue.Queue()
   file=open(facet_file,"r")
   #read facet_file and put each entry on q
   data=file.read().splitlines()
   for d in data:
      q.put(d)
   facets_input=q.qsize()
   print("number of facets in input file: ", facets_input)
   #initialize array of max_cores atlas processes
   Q=concurrent.futures.ThreadPoolExecutor()
   procs=[]
   T=[]   #array of results from the atlas processes
   #start max_cores atlas processes, taking input from stdin
   #output -> stdout, which writes to group/output_file_roots_i
   print("starting ", max_cores, "atlas processes")
   for i in range(max_cores):
      filename=  group + "/" + output_file_root + "_" + str(i)
      with open(filename, "w") as outfile:
        proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE, stdout=outfile)
        procs.append(proc)
   print("number of processes started: ", len(procs))
   P=concurrent.futures.ThreadPoolExecutor()
   #use P.submit 
   for i in range(max_cores):
#         print("submitting job i=",i)
         proc=procs[i]
         T.append(P.submit(atlas_compute,group,output_file_root,q,procs,i))
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
   print("number of K_polynomials computed: ", total_facets)

   #   print("status of processes: ")
#   for i in range(len(procs)):
#      print("i: ", i, procs[i].poll())
#   for p in procs:
#      print(proc.poll())
#   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



